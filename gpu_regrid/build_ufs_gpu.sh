#!/bin/bash
#==============================================================================
# Build UFS Coastal App with GPU Regridding on URSA
# Submit: sbatch build_ufs_gpu.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_ufs_gpu
#SBATCH --output=build_ufs_gpu_%j.out
#SBATCH --error=build_ufs_gpu_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --time=02:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building UFS Coastal App with GPU Regridding"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date:   $(date)"
echo "=============================================="

# Base directories
UFS_DIR=/scratch5/purged/Mansur.Jisan/ufs-weather-model
NVHPC_STACK=/scratch5/purged/Mansur.Jisan/nvhpc_stack
NVHPC_NETCDF=~/local/nvhpc-netcdf

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

# Set ESMF path (nvhpc-built ESMF 8.8.0)
export ESMFMKFILE=${NVHPC_STACK}/esmf-8.8.0-install/lib/libO/Linux.nvhpc.64.openmpi.default/esmf.mk

# Set NetCDF paths
export NetCDF_ROOT=$(nc-config --prefix)
export NETCDF_FORTRAN_ROOT=${NVHPC_NETCDF}
export LD_LIBRARY_PATH="${NVHPC_NETCDF}/lib:${LD_LIBRARY_PATH}"
export PATH="${NVHPC_NETCDF}/bin:${PATH}"

# Verify environment
echo ""
echo "Environment:"
echo "  Compiler:     $(nvfortran --version | head -1)"
echo "  ESMFMKFILE:   ${ESMFMKFILE}"
echo "  NetCDF-C:     ${NetCDF_ROOT}"
echo "  NetCDF-F:     ${NETCDF_FORTRAN_ROOT}"
echo ""

# Verify ESMF file exists
if [ ! -f "${ESMFMKFILE}" ]; then
    echo "ERROR: ESMFMKFILE not found at ${ESMFMKFILE}"
    exit 1
fi
echo "ESMF configuration file found."

# Verify nvhpc NetCDF-Fortran
if [ ! -f "${NVHPC_NETCDF}/lib/libnetcdff.a" ]; then
    echo "ERROR: nvhpc NetCDF-Fortran not found at ${NVHPC_NETCDF}"
    exit 1
fi
echo "nvhpc NetCDF-Fortran found."

# Create build directory
cd ${UFS_DIR}
rm -rf build_gpu
mkdir -p build_gpu
cd build_gpu

echo ""
echo "=============================================="
echo "Running CMake..."
echo "=============================================="

cmake .. \
  -DAPP=CSTLS \
  -DGPU_REGRID=ON \
  -DUSE_ATMOS=ON \
  -DNO_PARMETIS=OFF \
  -DOLDIO=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_C_COMPILER=nvc \
  -DCMAKE_CXX_COMPILER=nvc++ \
  -DNetCDF_ROOT=${NetCDF_ROOT} \
  -DNetCDF_Fortran_INCLUDE_DIRS=${NVHPC_NETCDF}/include \
  -DNetCDF_Fortran_LIBRARY=${NVHPC_NETCDF}/lib/libnetcdff.a \
  2>&1 | tee cmake_output.log

echo ""
echo "=============================================="
echo "Building UFS (this may take 30-60 minutes)..."
echo "=============================================="

make -j${SLURM_CPUS_PER_TASK} 2>&1 | tee build_output.log

echo ""
echo "=============================================="
echo "Build completed!"
echo "Date: $(date)"
echo "=============================================="

# Check if executable was created
if [ -f "ufs_model" ]; then
    echo "SUCCESS: ufs_model executable created"
    ls -la ufs_model
else
    echo "WARNING: ufs_model not found in build directory"
    echo "Searching for executables..."
    find . -name "*.exe" -o -name "ufs_*" 2>/dev/null | head -10
fi
