#!/bin/bash
#==============================================================================
# Build UFS Coastal App with GCC spack-stack + GPU Regrid Library
#
# This uses the hybrid approach:
# - UFS/CMEPS built with GCC from spack-stack
# - GPU regrid as external library (built with nvhpc)
# - Integration via C interoperability
#
# Submit: sbatch build_ufs_gcc_with_gpu_regrid.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_ufs_gcc
#SBATCH --output=build_ufs_gcc_%j.out
#SBATCH --error=build_ufs_gcc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:h100:1
#SBATCH --time=02:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building UFS Coastal with GCC + GPU Regrid"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date:   $(date)"
echo "=============================================="

# Directories
UFS_DIR=/scratch5/purged/Mansur.Jisan/ufs-weather-model-mjisan
GPU_REGRID_ROOT=/scratch5/purged/Mansur.Jisan/nvhpc_stack/gpu_regrid

#==============================================================================
# Load spack-stack environment (GCC-based) - Official URSA config
#==============================================================================
echo ""
echo "Loading URSA spack-stack environment..."

module purge

# Use official URSA spack-stack path
module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-gcc-12.4.0/install/modulefiles/Core

# Load GCC stack
module load stack-gcc/12.4.0
module load stack-openmpi/4.1.6

# Load cmake
module load cmake/3.27.9

# Load UFS dependencies from spack-stack
module load jasper/2.0.32
module load libpng/1.6.37
module load hdf5/1.14.3
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1
module load parallelio/2.6.2
module load esmf/8.8.0
module load fms/2024.02
module load bacio/2.4.1
module load sp/2.5.0
module load w3emc/2.10.0

# Save GCC MPI wrapper paths BEFORE loading nvhpc
GCC_MPI_ROOT=$(dirname $(which mpifort))
export GCC_MPICC="${GCC_MPI_ROOT}/mpicc"
export GCC_MPICXX="${GCC_MPI_ROOT}/mpic++"
export GCC_MPIFORT="${GCC_MPI_ROOT}/mpifort"

# Load nvhpc/cuda for GPU regrid runtime libraries only
module load nvhpc/24.11
module load cuda/12.8.1

# Set compilers to use GCC MPI wrappers (not nvhpc)
export MPI_CC="${GCC_MPICC}"
export MPI_CXX="${GCC_MPICXX}"
export MPI_FC="${GCC_MPIFORT}"
export FC="${GCC_MPIFORT}"
export CC="${GCC_MPICC}"
export CXX="${GCC_MPICXX}"

# Set cmake platform
export CMAKE_Platform=ursa.gnu

echo ""
echo "Environment:"
echo "  Compiler:        $(mpifort --version 2>/dev/null | head -1 || echo 'N/A')"
echo "  ESMFMKFILE:      ${ESMFMKFILE:-not set}"
echo "  GPU_REGRID_ROOT: ${GPU_REGRID_ROOT}"
echo ""
module list

#==============================================================================
# Verify GPU regrid library
#==============================================================================
if [ ! -f "${GPU_REGRID_ROOT}/lib/libgpu_regrid.so" ]; then
    echo "ERROR: GPU regrid library not found at ${GPU_REGRID_ROOT}"
    echo "Please build it first with: sbatch gpu_regrid/lib/build_gpu_regrid_lib.sh"
    exit 1
fi
echo "GPU regrid library found."

# Setup GPU regrid environment
source ${GPU_REGRID_ROOT}/setup_gpu_regrid.sh

#==============================================================================
# Build UFS
#==============================================================================
cd ${UFS_DIR}
rm -rf build_gcc_gpu
mkdir -p build_gcc_gpu
cd build_gcc_gpu

echo ""
echo "=============================================="
echo "Running CMake..."
echo "=============================================="

# Get NVHPC library paths for linking GPU regrid runtime
NVHPC_ROOT=$(dirname $(dirname $(which nvfortran)))
NVHPC_LIB="${NVHPC_ROOT}/compilers/lib"

# Export GPU_REGRID_ROOT for CMake to detect
export GPU_REGRID_ROOT=${GPU_REGRID_ROOT}

cmake .. \
  -DAPP=CSTLS \
  -DUSE_ATMOS=ON \
  -DNO_PARMETIS=OFF \
  -DOLDIO=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DGPU_REGRID=ON \
  -DGPU_REGRID_CINTEROP=ON \
  -DGPU_REGRID_ROOT=${GPU_REGRID_ROOT} \
  -DCMAKE_C_COMPILER="${GCC_MPICC}" \
  -DCMAKE_CXX_COMPILER="${GCC_MPICXX}" \
  -DCMAKE_Fortran_COMPILER="${GCC_MPIFORT}" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${GPU_REGRID_ROOT}/lib -lgpu_regrid -L${NVHPC_LIB} -lnvf -lacchost -laccdevice -Wl,-rpath,${GPU_REGRID_ROOT}/lib -Wl,-rpath,${NVHPC_LIB}" \
  2>&1 | tee cmake_output.log

echo ""
echo "=============================================="
echo "Building UFS..."
echo "=============================================="

make -j${SLURM_CPUS_PER_TASK} 2>&1 | tee build_output.log

echo ""
echo "=============================================="
echo "Build completed!"
echo "Date: $(date)"
echo "=============================================="

# Check result
if [ -f "ufs_model" ]; then
    echo "SUCCESS: ufs_model executable created"
    ls -la ufs_model

    # Verify GPU regrid is linked
    echo ""
    echo "Checking GPU regrid library linkage..."
    ldd ufs_model | grep -i gpu_regrid || echo "Warning: gpu_regrid not in ldd output"
else
    echo "WARNING: ufs_model not found"
fi

echo ""
echo "Runtime requirements:"
echo "  module load nvhpc/24.11"
echo "  source ${GPU_REGRID_ROOT}/setup_gpu_regrid.sh"
