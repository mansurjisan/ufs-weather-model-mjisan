#!/bin/bash
#==============================================================================
# Build ESMF with NVHPC for GPU-enabled UFS on URSA
# Submit: sbatch build_esmf_nvhpc.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_esmf
#SBATCH --output=build_esmf_%j.out
#SBATCH --error=build_esmf_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:h100:1
#SBATCH --time=02:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building ESMF with NVHPC on URSA"
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
module load netcdf-fortran/4.6.1

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
# Clone ESMF if not present
#==============================================================================
ESMF_VERSION=v8.6.1
if [ ! -d "esmf" ]; then
    echo "Cloning ESMF ${ESMF_VERSION}..."
    git clone --branch ${ESMF_VERSION} --depth 1 https://github.com/esmf-org/esmf.git
fi

cd esmf

#==============================================================================
# Set ESMF build environment
#==============================================================================
export ESMF_DIR=$(pwd)
export ESMF_INSTALL_PREFIX=${INSTALL_BASE}/esmf-install
export ESMF_COMPILER=nvhpc
export ESMF_COMM=openmpi
export ESMF_NETCDF=nc-config
export ESMF_NFCONFIG=nf-config
export ESMF_BOPT=O
export ESMF_OPTLEVEL=2
export ESMF_ABI=64
export ESMF_OS=Linux
export ESMF_SITE=default

# Enable OpenACC for GPU
export ESMF_ACC_SOFTWARE_STACK=nvhpc
export ESMF_ACC_DEVICE=nvidia

# Parallel build
export ESMF_MOAB=OFF
export ESMF_LAPACK=OFF

echo "ESMF Configuration:"
echo "  ESMF_DIR:            ${ESMF_DIR}"
echo "  ESMF_INSTALL_PREFIX: ${ESMF_INSTALL_PREFIX}"
echo "  ESMF_COMPILER:       ${ESMF_COMPILER}"
echo "  ESMF_COMM:           ${ESMF_COMM}"
echo ""

#==============================================================================
# Build ESMF
#==============================================================================
echo "Building ESMF (this takes ~30-60 minutes)..."
echo "Start time: $(date)"

make -j${SLURM_NTASKS} 2>&1 | tee build.log

echo ""
echo "Build completed at $(date)"

#==============================================================================
# Install ESMF
#==============================================================================
echo "Installing ESMF..."
make install 2>&1 | tee install.log

#==============================================================================
# Create module file
#==============================================================================
MODULEFILE_DIR=${INSTALL_BASE}/modulefiles/esmf
mkdir -p ${MODULEFILE_DIR}

cat > ${MODULEFILE_DIR}/8.6.1-nvhpc.lua << 'MODEOF'
help([[
ESMF 8.6.1 built with NVHPC for GPU support
]])

local base = "/scratch5/purged/Mansur.Jisan/nvhpc_stack/esmf-install"

prepend_path("PATH", pathJoin(base, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(base, "lib"))
prepend_path("LIBRARY_PATH", pathJoin(base, "lib"))
prepend_path("CPATH", pathJoin(base, "include"))
prepend_path("CMAKE_PREFIX_PATH", base)

setenv("ESMF_ROOT", base)
setenv("ESMFMKFILE", pathJoin(base, "lib/esmf.mk"))

whatis("ESMF 8.6.1 with NVHPC/GPU support")
MODEOF

echo ""
echo "=============================================="
echo "ESMF build completed!"
echo "Install location: ${ESMF_INSTALL_PREFIX}"
echo ""
echo "To use:"
echo "  module use ${INSTALL_BASE}/modulefiles"
echo "  module load esmf/8.6.1-nvhpc"
echo "=============================================="
