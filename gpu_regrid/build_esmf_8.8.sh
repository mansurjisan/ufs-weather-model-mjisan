#!/bin/bash
#==============================================================================
# Build ESMF 8.8.0 with NVHPC for GPU-enabled UFS on URSA
# UFS requires ESMF >= 8.8.0
# Submit: sbatch build_esmf_8.8.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_esmf88
#SBATCH --output=build_esmf88_%j.out
#SBATCH --error=build_esmf88_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:h100:1
#SBATCH --time=03:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building ESMF 8.8.0 with NVHPC on URSA"
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

# Use nvhpc NetCDF-Fortran
export NETCDF_FORTRAN_ROOT=~/local/nvhpc-netcdf

# Use NVHPC bundled MPI
NVHPC_ROOT=$(dirname $(dirname $(which nvfortran)))
export NVHPC_MPI_ROOT="${NVHPC_ROOT}/comm_libs/mpi"
export PATH="${NVHPC_MPI_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${NVHPC_MPI_ROOT}/lib:${LD_LIBRARY_PATH}"

echo ""
echo "Compiler: $(nvfortran --version | head -1)"
echo "MPI:      $(mpifort --version 2>&1 | head -1)"
echo ""

#==============================================================================
# Clone ESMF 8.8.0 if not present
#==============================================================================
ESMF_VERSION=v8.8.0
ESMF_SRC_DIR=${INSTALL_BASE}/esmf-8.8.0

if [ ! -d "${ESMF_SRC_DIR}" ]; then
    echo "Cloning ESMF ${ESMF_VERSION}..."
    # Clone needs to be done on login node, check if exists
    echo "ERROR: ESMF ${ESMF_VERSION} not found."
    echo "Please clone on login node first:"
    echo "  cd ${INSTALL_BASE}"
    echo "  git clone --branch ${ESMF_VERSION} --depth 1 https://github.com/esmf-org/esmf.git esmf-8.8.0"
    exit 1
fi

cd ${ESMF_SRC_DIR}

# Clean previous build
rm -rf obj mod lib

#==============================================================================
# Set ESMF build environment
#==============================================================================
export ESMF_DIR=$(pwd)
export ESMF_INSTALL_PREFIX=${INSTALL_BASE}/esmf-8.8.0-install
export ESMF_COMPILER=nvhpc
export ESMF_COMM=openmpi
export ESMF_NETCDF=nc-config
export ESMF_NFCONFIG=${NETCDF_FORTRAN_ROOT}/bin/nf-config
export ESMF_BOPT=O
export ESMF_OPTLEVEL=2
export ESMF_ABI=64
export ESMF_OS=Linux
export ESMF_SITE=default

# Disable features that cause issues
export ESMF_MOAB=OFF
export ESMF_LAPACK=OFF
export ESMF_PIO=OFF

# Enable OpenACC for GPU
export ESMF_ACC_SOFTWARE_STACK=nvhpc
export ESMF_ACC_DEVICE=nvidia

echo "ESMF Configuration:"
echo "  ESMF_DIR:            ${ESMF_DIR}"
echo "  ESMF_INSTALL_PREFIX: ${ESMF_INSTALL_PREFIX}"
echo "  ESMF_COMPILER:       ${ESMF_COMPILER}"
echo "  ESMF_COMM:           ${ESMF_COMM}"
echo "  ESMF_NFCONFIG:       ${ESMF_NFCONFIG}"
echo ""

#==============================================================================
# Build ESMF (reduced parallelism to avoid OOM)
#==============================================================================
echo "Building ESMF 8.8.0 (this takes ~60-90 minutes)..."
echo "Start time: $(date)"

# Use -j2 to avoid OOM (nvhpc compiler is memory-hungry)
make -j2 2>&1 | tee ${INSTALL_BASE}/esmf_8.8.0_build.log

echo ""
echo "Build completed at $(date)"

#==============================================================================
# Install ESMF
#==============================================================================
echo "Installing ESMF..."
make install 2>&1 | tee ${INSTALL_BASE}/esmf_8.8.0_install.log

#==============================================================================
# Create module file
#==============================================================================
MODULEFILE_DIR=${INSTALL_BASE}/modulefiles/esmf
mkdir -p ${MODULEFILE_DIR}

cat > ${MODULEFILE_DIR}/8.8.0-nvhpc.lua << 'MODEOF'
help([[
ESMF 8.8.0 built with NVHPC for GPU support
]])

local base = "/scratch5/purged/Mansur.Jisan/nvhpc_stack/esmf-8.8.0-install"

prepend_path("PATH", pathJoin(base, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(base, "lib"))
prepend_path("LIBRARY_PATH", pathJoin(base, "lib"))
prepend_path("CPATH", pathJoin(base, "include"))
prepend_path("CMAKE_PREFIX_PATH", base)

setenv("ESMF_ROOT", base)
setenv("ESMFMKFILE", pathJoin(base, "lib/libO/Linux.nvhpc.64.openmpi.default/esmf.mk"))

whatis("ESMF 8.8.0 with NVHPC/GPU support")
MODEOF

echo ""
echo "=============================================="
echo "ESMF 8.8.0 build completed!"
echo "Install location: ${ESMF_INSTALL_PREFIX}"
echo ""
echo "esmf.mk location:"
find ${ESMF_INSTALL_PREFIX} -name "esmf.mk" 2>/dev/null
echo ""
echo "To use:"
echo "  module use ${INSTALL_BASE}/modulefiles"
echo "  module load esmf/8.8.0-nvhpc"
echo ""
echo "Or set directly:"
echo "  export ESMFMKFILE=\$(find ${ESMF_INSTALL_PREFIX} -name 'esmf.mk')"
echo "=============================================="
