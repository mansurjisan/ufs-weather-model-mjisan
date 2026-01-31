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
# Load spack-stack environment (GCC-based)
#==============================================================================
echo ""
echo "Loading spack-stack environment..."

# URSA spack-stack module path
module purge
module use /apps/modules/spack-stack/modulefiles/core

# Load spack-stack gcc environment
# Note: Adjust this based on what's available on URSA
# Common pattern: spack-stack/1.x.y or compiler/gcc/x.y.z first
module load stack-intel/2022.1.0  2>/dev/null || \
module load stack-gcc/11.3.0      2>/dev/null || \
module load gcc/11.3.0            2>/dev/null || \
echo "Warning: Could not find standard spack-stack module"

# Load required libraries from spack-stack
module load cmake 2>/dev/null || true
module load netcdf-c 2>/dev/null || true
module load netcdf-fortran 2>/dev/null || true
module load esmf 2>/dev/null || true
module load pio 2>/dev/null || true

# Also load nvhpc for GPU regrid runtime
module load nvhpc/24.11
module load cuda/12.8.1

echo ""
echo "Environment:"
echo "  Compiler:        $(gfortran --version 2>/dev/null | head -1 || echo 'N/A')"
echo "  ESMFMKFILE:      ${ESMFMKFILE:-not set}"
echo "  GPU_REGRID_ROOT: ${GPU_REGRID_ROOT}"
echo ""

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

# Get NVHPC library paths for linking GPU regrid
NVHPC_ROOT=$(dirname $(dirname $(which nvfortran)))
NVHPC_LIB="${NVHPC_ROOT}/compilers/lib"

cmake .. \
  -DAPP=CSTLS \
  -DUSE_ATMOS=ON \
  -DNO_PARMETIS=OFF \
  -DOLDIO=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DGPU_REGRID=ON \
  -DGPU_REGRID_ROOT=${GPU_REGRID_ROOT} \
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
