#!/bin/bash
#===============================================================================
# Build script for GPU-accelerated regridding on NOAA URSA cluster
#
# Usage: ./build_gpu_regrid.sh [test|cmeps]
#===============================================================================

set -e

BUILD_TYPE=${1:-test}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
CMEPS_DIR="${ROOT_DIR}/CMEPS-interface/CMEPS/mediator"

# Load URSA GPU modules (known working configuration)
echo "Loading GPU modules..."
module purge 2>/dev/null || true

# Load NVIDIA HPC SDK
echo "  Loading nvhpc/24.11..."
module load nvhpc/24.11

# Load CUDA toolkit
echo "  Loading cuda/12.8.1..."
module load cuda/12.8.1

# Load HDF5 and NetCDF (optional but useful)
echo "  Loading hdf5/1.14.3..."
module load hdf5/1.14.3 2>/dev/null || echo "    HDF5 not found (optional)"

echo "  Loading netcdf-c/4.9.2..."
module load netcdf-c/4.9.2 2>/dev/null || echo "    NetCDF-C not found (optional)"

echo "  Loading netcdf-fortran/4.6.1..."
module load netcdf-fortran/4.6.1 2>/dev/null || echo "    NetCDF-Fortran not found (optional)"

echo ""
echo "Loaded modules:"
module list 2>&1

# Compiler settings
FC=nvfortran
FFLAGS_BASE="-O3 -Minfo=all"
FFLAGS_ACC="-acc -gpu=cc80,cc90 -Minfo=accel"
FFLAGS_CUDA="-cuda -gpu=cc80,cc90 -cudalib=cusparse"

echo ""
echo "================================================"
echo " GPU Regrid Build Configuration"
echo "================================================"
echo " Compiler: $FC"
echo " Build type: $BUILD_TYPE"
echo " GPU flags: $FFLAGS_ACC"
echo ""

cd "$CMEPS_DIR"

if [ "$BUILD_TYPE" == "test" ]; then
    #---------------------------------------------------------------------------
    # Build standalone test
    #---------------------------------------------------------------------------
    echo "Building standalone test..."

    # Create build directory
    mkdir -p build_gpu_test
    cd build_gpu_test

    # Compile modules
    echo "  Compiling med_kind_mod (stub)..."
    cat > med_kind_mod.F90 << 'EOF'
module med_kind_mod
  integer, parameter :: SHR_KIND_R8 = selected_real_kind(15,307)
  integer, parameter :: SHR_KIND_I4 = selected_int_kind(9)
  integer, parameter :: SHR_KIND_CL = 256
end module med_kind_mod
EOF
    $FC $FFLAGS_BASE -c med_kind_mod.F90

    echo "  Compiling gpu_regrid_mod..."
    $FC $FFLAGS_BASE $FFLAGS_ACC -c ../gpu_regrid_mod.F90

    echo "  Compiling test program..."
    $FC $FFLAGS_BASE $FFLAGS_ACC -c ../gpu_regrid_test.F90

    echo "  Linking..."
    $FC $FFLAGS_BASE $FFLAGS_ACC -o test_gpu_regrid \
        med_kind_mod.o gpu_regrid_mod.o gpu_regrid_test.o

    echo ""
    echo "Build successful!"
    echo "Test executable: $CMEPS_DIR/build_gpu_test/test_gpu_regrid"
    echo ""
    echo "To run on GPU node:"
    echo "  salloc -p u1-h100 -N 1 --gres=gpu:1 -t 0:30:00"
    echo "  srun ./test_gpu_regrid"

elif [ "$BUILD_TYPE" == "cmeps" ]; then
    #---------------------------------------------------------------------------
    # Build with CMEPS integration
    #---------------------------------------------------------------------------
    echo "Building with CMEPS integration..."
    echo "NOTE: This requires modifying CMEPS CMakeLists.txt"

    # Check for CMakeLists.txt
    if [ ! -f "../CMakeLists.txt" ]; then
        echo "ERROR: CMEPS CMakeLists.txt not found"
        exit 1
    fi

    # Show integration instructions
    cat << 'EOF'

CMEPS Integration Instructions:
===============================

1. Add to CMEPS/CMakeLists.txt:

   option(GPU_REGRID "Enable GPU regridding" OFF)

   if(GPU_REGRID)
     enable_language(CUDA)
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc80")
     list(APPEND CMEPS_SOURCES
       mediator/gpu_regrid_mod.F90
     )
   endif()

2. Modify mediator/med_map_mod.F90:

   Add at top:
     use gpu_regrid_mod, only : gpu_regrid_available, gpu_regrid_apply

   In med_map_field(), before ESMF_FieldRegrid call:
     if (gpu_regrid_available(n1, n2, mapindex)) then
       call ESMF_FieldGet(field_src, farrayPtr=src_ptr, rc=rc)
       call ESMF_FieldGet(field_dst, farrayPtr=dst_ptr, rc=rc)
       call gpu_regrid_apply(n1, n2, mapindex, src_ptr, dst_ptr, rc)
     else
       call ESMF_FieldRegrid(...)
     end if

3. Rebuild UFS with GPU support:

   cmake .. -DGPU_REGRID=ON
   make -j8

EOF

else
    echo "Unknown build type: $BUILD_TYPE"
    echo "Usage: $0 [test|cmeps]"
    exit 1
fi

echo ""
echo "Done!"
