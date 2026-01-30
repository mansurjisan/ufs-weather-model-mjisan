#!/bin/bash
#==============================================================================
# Build complete nvhpc stack for GPU-enabled UFS on URSA
# This builds: NetCDF-Fortran + ESMF with nvhpc
# Submit: sbatch build_nvhpc_stack.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_nvhpc
#SBATCH --output=build_nvhpc_%j.out
#SBATCH --error=build_nvhpc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:h100:1
#SBATCH --time=04:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building nvhpc stack for GPU UFS on URSA"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date:   $(date)"
echo "=============================================="

# Base directory
INSTALL_BASE=/scratch5/purged/Mansur.Jisan/nvhpc_stack
mkdir -p ${INSTALL_BASE}
cd ${INSTALL_BASE}

# Load modules
module purge
module load nvhpc/24.11
module load cuda/12.8.1
module load hdf5/1.14.3
module load netcdf-c/4.9.2

# Use NVHPC bundled MPI
NVHPC_ROOT=$(dirname $(dirname $(which nvfortran)))
export NVHPC_MPI_ROOT="${NVHPC_ROOT}/comm_libs/mpi"
export PATH="${NVHPC_MPI_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${NVHPC_MPI_ROOT}/lib:${LD_LIBRARY_PATH}"

# Get NetCDF-C paths
export NETCDF_C_ROOT=$(nc-config --prefix)
export CFLAGS="-I${NETCDF_C_ROOT}/include"
export LDFLAGS="-L${NETCDF_C_ROOT}/lib"
export LD_LIBRARY_PATH="${NETCDF_C_ROOT}/lib:${LD_LIBRARY_PATH}"

echo ""
echo "Compiler: $(nvfortran --version | head -2)"
echo "NetCDF-C: ${NETCDF_C_ROOT}"
echo ""

#==============================================================================
# Step 1: Build NetCDF-Fortran with nvhpc
#==============================================================================
NETCDF_F_VERSION=4.6.1
NETCDF_F_INSTALL=${INSTALL_BASE}/netcdf-fortran-${NETCDF_F_VERSION}

if [ ! -f "${NETCDF_F_INSTALL}/lib/libnetcdff.so" ]; then
    echo "=============================================="
    echo "Building NetCDF-Fortran ${NETCDF_F_VERSION}"
    echo "=============================================="

    cd ${INSTALL_BASE}

    # Download if not present
    if [ ! -d "netcdf-fortran-${NETCDF_F_VERSION}" ]; then
        if [ ! -f "netcdf-fortran-${NETCDF_F_VERSION}.tar.gz" ]; then
            # Try to copy from local if download fails
            echo "Downloading NetCDF-Fortran..."
            curl -L -o netcdf-fortran-${NETCDF_F_VERSION}.tar.gz \
                https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NETCDF_F_VERSION}.tar.gz || \
            wget -O netcdf-fortran-${NETCDF_F_VERSION}.tar.gz \
                https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NETCDF_F_VERSION}.tar.gz
        fi
        tar xzf netcdf-fortran-${NETCDF_F_VERSION}.tar.gz
    fi

    cd netcdf-fortran-${NETCDF_F_VERSION}

    # Clean any previous build
    rm -rf build && mkdir build && cd build

    # Configure with nvhpc
    FC=nvfortran CC=nvc CXX=nvc++ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NETCDF_F_INSTALL} \
        -DCMAKE_PREFIX_PATH=${NETCDF_C_ROOT} \
        -DENABLE_TESTS=OFF \
        -DBUILD_SHARED_LIBS=ON

    # Build and install
    make -j${SLURM_NTASKS}
    make install

    echo "NetCDF-Fortran installed to ${NETCDF_F_INSTALL}"
else
    echo "NetCDF-Fortran already installed at ${NETCDF_F_INSTALL}"
fi

# Update paths
export NETCDF_FORTRAN_ROOT=${NETCDF_F_INSTALL}
export LD_LIBRARY_PATH="${NETCDF_FORTRAN_ROOT}/lib:${LD_LIBRARY_PATH}"
export PATH="${NETCDF_FORTRAN_ROOT}/bin:${PATH}"

#==============================================================================
# Step 2: Build ESMF with nvhpc
#==============================================================================
echo ""
echo "=============================================="
echo "Building ESMF with nvhpc"
echo "=============================================="

cd ${INSTALL_BASE}

if [ ! -d "esmf" ]; then
    echo "ERROR: ESMF not cloned. Please clone on login node first:"
    echo "  cd ${INSTALL_BASE}"
    echo "  git clone --branch v8.6.1 --depth 1 https://github.com/esmf-org/esmf.git"
    exit 1
fi

cd esmf

# Clean previous build
rm -rf obj mod lib

# Set ESMF environment
export ESMF_DIR=$(pwd)
export ESMF_INSTALL_PREFIX=${INSTALL_BASE}/esmf-install
export ESMF_COMPILER=nvhpc
export ESMF_COMM=openmpi
export ESMF_BOPT=O
export ESMF_OPTLEVEL=2
export ESMF_ABI=64
export ESMF_OS=Linux
export ESMF_SITE=default

# Point to our nvhpc-built NetCDF-Fortran
export ESMF_NETCDF=nc-config
export ESMF_NFCONFIG=${NETCDF_FORTRAN_ROOT}/bin/nf-config

# Disable features that cause issues
export ESMF_MOAB=OFF
export ESMF_LAPACK=OFF
export ESMF_PIO=OFF

echo "ESMF Configuration:"
echo "  ESMF_DIR:            ${ESMF_DIR}"
echo "  ESMF_INSTALL_PREFIX: ${ESMF_INSTALL_PREFIX}"
echo "  ESMF_COMPILER:       ${ESMF_COMPILER}"
echo "  ESMF_NFCONFIG:       ${ESMF_NFCONFIG}"
echo ""

# Build
echo "Building ESMF (this takes ~30-60 minutes)..."
echo "Start time: $(date)"

make -j${SLURM_NTASKS} 2>&1 | tee ${INSTALL_BASE}/esmf_build.log

echo "Build completed at $(date)"

# Install
echo "Installing ESMF..."
make install 2>&1 | tee ${INSTALL_BASE}/esmf_install.log

#==============================================================================
# Create module file
#==============================================================================
MODULEFILE_DIR=${INSTALL_BASE}/modulefiles
mkdir -p ${MODULEFILE_DIR}/nvhpc-stack

cat > ${MODULEFILE_DIR}/nvhpc-stack/1.0.lua << MODEOF
help([[
nvhpc stack for GPU-enabled UFS
Includes: NetCDF-Fortran 4.6.1, ESMF 8.6.1 (built with nvhpc)
]])

-- Load base modules
load("nvhpc/24.11")
load("cuda/12.8.1")
load("hdf5/1.14.3")
load("netcdf-c/4.9.2")

-- NVHPC MPI
local nvhpc_root = os.getenv("NVHPC_ROOT") or "/opt/nvidia/hpc_sdk/Linux_x86_64/24.11"
local mpi_root = pathJoin(nvhpc_root, "comm_libs/mpi")
prepend_path("PATH", pathJoin(mpi_root, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(mpi_root, "lib"))

-- NetCDF-Fortran (nvhpc-built)
local nf_root = "${NETCDF_FORTRAN_ROOT}"
prepend_path("PATH", pathJoin(nf_root, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(nf_root, "lib"))
prepend_path("CPATH", pathJoin(nf_root, "include"))
prepend_path("CMAKE_PREFIX_PATH", nf_root)
setenv("NETCDF_FORTRAN_ROOT", nf_root)

-- ESMF (nvhpc-built)
local esmf_root = "${ESMF_INSTALL_PREFIX}"
prepend_path("PATH", pathJoin(esmf_root, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(esmf_root, "lib"))
prepend_path("CPATH", pathJoin(esmf_root, "include"))
prepend_path("CMAKE_PREFIX_PATH", esmf_root)
setenv("ESMF_ROOT", esmf_root)
setenv("ESMFMKFILE", pathJoin(esmf_root, "lib/esmf.mk"))

-- Compiler settings
setenv("CC", "nvc")
setenv("CXX", "nvc++")
setenv("FC", "nvfortran")

whatis("nvhpc stack for GPU-enabled UFS")
MODEOF

echo ""
echo "=============================================="
echo "Build completed!"
echo ""
echo "NetCDF-Fortran: ${NETCDF_FORTRAN_ROOT}"
echo "ESMF:           ${ESMF_INSTALL_PREFIX}"
echo ""
echo "To use:"
echo "  module use ${MODULEFILE_DIR}"
echo "  module load nvhpc-stack/1.0"
echo "=============================================="
