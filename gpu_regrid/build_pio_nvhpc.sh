#!/bin/bash
#==============================================================================
# Build PIO (Parallel I/O) with NVHPC for GPU-enabled UFS on URSA
# Required by CMEPS
# Submit: sbatch build_pio_nvhpc.sh
#==============================================================================
#SBATCH --partition=u1-h100
#SBATCH --job-name=build_pio
#SBATCH --output=build_pio_%j.out
#SBATCH --error=build_pio_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:h100:1
#SBATCH --time=01:00:00
#SBATCH --account=gpu-nos-surge
#SBATCH --qos=gpu

set -e

echo "=============================================="
echo "Building PIO with NVHPC on URSA"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Date:   $(date)"
echo "=============================================="

# Base directory
INSTALL_BASE=/scratch5/purged/Mansur.Jisan/nvhpc_stack
PIO_VERSION=2.6.9
PIO_SRC=${INSTALL_BASE}/ParallelIO-pio${PIO_VERSION}
PIO_INSTALL=${INSTALL_BASE}/pio-${PIO_VERSION}-install

mkdir -p ${INSTALL_BASE}
cd ${INSTALL_BASE}

# Load modules
module purge
module load nvhpc/24.11
module load cuda/12.8.1
module load cmake
module load hdf5/1.14.3
module load netcdf-c/4.9.2

# Use nvhpc NetCDF-Fortran
export NETCDF_FORTRAN_ROOT=~/local/nvhpc-netcdf

# Use NVHPC bundled MPI
NVHPC_ROOT=$(dirname $(dirname $(which nvfortran)))
export NVHPC_MPI_ROOT="${NVHPC_ROOT}/comm_libs/mpi"
export PATH="${NVHPC_MPI_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${NVHPC_MPI_ROOT}/lib:${LD_LIBRARY_PATH}"

# Get MPI paths from mpicc/mpifort
MPI_INCLUDE=$(mpicc -show | grep -oE '\-I[^ ]+' | head -1 | sed 's/-I//')
MPI_LIB_DIR=$(mpicc -show | grep -oE '\-L[^ ]+' | head -1 | sed 's/-L//')

# Find MPI Fortran module directory (mpi.mod) - it's in lib/ not include/
# Use the same hpcx version as the MPI we're using (12.6)
MPI_FORT_MOD=$(dirname ${MPI_INCLUDE})/lib
if [ ! -f "${MPI_FORT_MOD}/mpi.mod" ]; then
    # Fallback: search for it
    MPI_MOD_FILE=$(find ${NVHPC_ROOT} -path "*12.6*" -name "mpi.mod" 2>/dev/null | head -1)
    if [ -n "${MPI_MOD_FILE}" ]; then
        MPI_FORT_MOD=$(dirname ${MPI_MOD_FILE})
    fi
fi

# NetCDF paths
export NetCDF_C_PATH=$(nc-config --prefix)
export NetCDF_Fortran_PATH=${NETCDF_FORTRAN_ROOT}
export HDF5_PATH=$(h5cc -showconfig | grep "Installation point" | awk '{print $3}')

# Create pkg-config file for nvhpc netcdf-fortran (PIO 2.6.9 uses pkg-config)
mkdir -p ${NETCDF_FORTRAN_ROOT}/lib/pkgconfig
cat > ${NETCDF_FORTRAN_ROOT}/lib/pkgconfig/netcdf-fortran.pc << PKGEOF
prefix=${NETCDF_FORTRAN_ROOT}
exec_prefix=\${prefix}
libdir=\${exec_prefix}/lib
includedir=\${prefix}/include

Name: NetCDF-Fortran
Description: NetCDF Fortran library
Version: 4.6.1
Libs: -L\${libdir} -lnetcdff
Cflags: -I\${includedir}
Requires: netcdf
PKGEOF

export PKG_CONFIG_PATH="${NETCDF_FORTRAN_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"

echo ""
echo "Environment:"
echo "  Compiler:      $(nvfortran --version | head -1)"
echo "  NetCDF-C:      ${NetCDF_C_PATH}"
echo "  NetCDF-F:      ${NetCDF_Fortran_PATH}"
echo "  HDF5:          ${HDF5_PATH}"
echo ""

#==============================================================================
# Clone PIO if not present
#==============================================================================
if [ ! -d "${PIO_SRC}" ]; then
    echo "ERROR: PIO source not found at ${PIO_SRC}"
    echo "Please clone on login node first:"
    echo "  cd ${INSTALL_BASE}"
    echo "  git clone --branch pio2_6_9 --recursive https://github.com/NCAR/ParallelIO.git ParallelIO-pio${PIO_VERSION}"
    exit 1
fi

cd ${PIO_SRC}

#==============================================================================
# Build PIO
#==============================================================================
rm -rf build
mkdir -p build
cd build

# Check for genf90 (must be pre-cloned on login node)
if [ ! -d "${PIO_SRC}/src/flib/genf90/src/genf90" ]; then
    echo "ERROR: genf90 not found. Please clone on login node first:"
    echo "  cd ${PIO_SRC}"
    echo "  git clone https://github.com/PARALLELIO/genf90.git src/flib/genf90/src/genf90"
    exit 1
fi

echo "Configuring PIO..."
echo "  MPI_INCLUDE: ${MPI_INCLUDE}"
echo "  MPI_FORT_MOD: ${MPI_FORT_MOD}"

FC=nvfortran CC=nvc CXX=nvc++ cmake .. \
    -DCMAKE_INSTALL_PREFIX=${PIO_INSTALL} \
    -DUSER_CMAKE_MODULE_PATH=${PIO_SRC}/cmake/CMake_Fortran_utils \
    -DGENF90_PATH=${PIO_SRC}/src/flib/genf90/src/genf90 \
    -DPIO_ENABLE_FORTRAN=ON \
    -DPIO_ENABLE_TIMING=OFF \
    -DPIO_ENABLE_TESTS=OFF \
    -DPIO_ENABLE_EXAMPLES=OFF \
    -DPIO_ENABLE_DOC=OFF \
    -DWITH_PNETCDF=OFF \
    -DNetCDF_C_PATH=${NetCDF_C_PATH} \
    -DNetCDF_Fortran_PATH=${NetCDF_Fortran_PATH} \
    -DMPI_C_COMPILER=mpicc \
    -DMPI_Fortran_COMPILER=mpifort \
    -DCMAKE_C_FLAGS="-fPIC -I${MPI_INCLUDE}" \
    -DCMAKE_Fortran_FLAGS="-fPIC -I${MPI_INCLUDE} -I${MPI_FORT_MOD}"

echo ""
echo "Building PIO..."
make -j${SLURM_CPUS_PER_TASK} 2>&1 | tee build.log

echo ""
echo "Installing PIO..."
make install 2>&1 | tee install.log

echo ""
echo "=============================================="
echo "PIO ${PIO_VERSION} build completed!"
echo "Install location: ${PIO_INSTALL}"
echo ""
echo "For UFS build, set:"
echo "  export PIO=${PIO_INSTALL}"
echo "=============================================="
