# NVHPC (NVIDIA HPC SDK) compiler flags for UFS
# Supports OpenACC GPU acceleration

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -Mbackslash -Mpreprocess -Mallocatable=03")

if(NOT 32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
endif()

# GPU/OpenACC flags
if(GPU_REGRID)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc80,cc90 -Minfo=accel")
    add_definitions(-DGPU_REGRID)
    message(STATUS "GPU regridding enabled with OpenACC (NVHPC)")
endif()

if(DEBUG)
    add_definitions(-DDEBUG)
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Mbounds -Mchkptr -Ktrap=fp")
    set(CMAKE_C_FLAGS_DEBUG "-O0")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0")
else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fast")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

# OpenMP flags (NVHPC uses -mp)
if(OPENMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp")
endif()
