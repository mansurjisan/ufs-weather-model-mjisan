# GPU-Accelerated Regridding for UFS Coastal Applications

This directory contains GPU-accelerated regridding code for ESMF/CMEPS
to improve performance of coupled coastal simulations (SCHISM + WW3 + DATM).

## Performance

Tested on NVIDIA H100 (NOAA URSA cluster):

| Scenario | Speedup vs CPU |
|----------|---------------|
| Single field (warmed up) | **2.4x** |
| Single field (GPU-resident data) | **3.4x** |

## Overview

ESMF regridding is fundamentally a sparse matrix-vector multiply (SpMV) operation.
This implementation extracts the weight matrix from ESMF and applies it on GPU
using OpenACC for portable GPU acceleration.

## Files

```
gpu_regrid/
├── README.md                          # This file
├── CMEPS_INTEGRATION.md               # Step-by-step CMEPS integration guide
├── lib/                               # C-interop shared library (for hybrid builds)
│   ├── gpu_regrid_c_interface.F90    # Fortran module with C bindings
│   ├── gpu_regrid.h                  # C header file
│   ├── build_gpu_regrid_lib.sh       # Build script for shared library
│   └── test_c_interface.c            # C test program
├── cmeps/
│   └── gpu_regrid_cmeps.F90          # CMEPS integration wrapper
├── cmeps_mediator/
│   ├── gpu_regrid_mod.F90            # OpenACC-based GPU regrid module
│   ├── gpu_regrid_test.F90           # Standalone test program
│   └── gpu_regrid_cusparse.F90       # cuSPARSE version (experimental)
├── build_esmf_8.8.sh                  # ESMF 8.8.0 build with NVHPC
├── build_pio_nvhpc.sh                 # PIO build with NVHPC (may fail)
├── build_ufs_gpu.sh                   # Full NVHPC UFS build
├── build_ufs_gcc_with_gpu_regrid.sh   # Hybrid GCC + GPU library build
└── (see also)
    ├── ../GPU_REGRID_DESIGN.md       # Detailed design document
    ├── ../CMEPS-interface/CMEPS/mediator/gpu_regrid_mod.F90         # OpenACC version
    └── ../CMEPS-interface/CMEPS/mediator/gpu_regrid_mod_cinterop.F90 # C interop version
```

## Installation

### Option 1: Hybrid Build (Recommended for Production)

Build GPU regrid as a separate library with NVHPC, then link from GCC/Intel-compiled UFS.
This is the recommended approach as PIO has compatibility issues with nvfortran.

```bash
# 1. Build the GPU regrid shared library with NVHPC
cd gpu_regrid/lib
sbatch build_gpu_regrid_lib.sh

# 2. Build UFS with GCC spack-stack + GPU regrid library
export GPU_REGRID_ROOT=/scratch5/purged/Mansur.Jisan/nvhpc_stack/gpu_regrid
cd ..
sbatch build_ufs_gcc_with_gpu_regrid.sh
```

### Option 2: Full NVHPC Build (Experimental)

Build entire UFS with NVHPC. Requires building all dependencies (ESMF, PIO, etc.) with NVHPC.
Note: PIO has known issues with nvfortran's int64 handling.

```bash
# Build all dependencies with NVHPC first
sbatch gpu_regrid/build_esmf_8.8.sh
sbatch gpu_regrid/build_pio_nvhpc.sh  # May fail due to nvfortran issues

# Build UFS
sbatch gpu_regrid/build_ufs_gpu.sh
```

### Option 3: Standalone Test Only

```bash
./scripts/build_gpu_regrid.sh test
```

## Quick Start on URSA

```bash
# 1. Build test
./scripts/build_gpu_regrid.sh test

# 2. Run on H100 GPU node
salloc -p u1-h100 -N 1 --gres=gpu:1 -t 0:30:00
cd CMEPS-interface/CMEPS/mediator/build_gpu_test
srun ./test_gpu_regrid
```

## Integration with CMEPS

See **[CMEPS_INTEGRATION.md](CMEPS_INTEGRATION.md)** for detailed step-by-step instructions.

Key steps:
1. Copy `gpu_regrid_mod.F90` to CMEPS mediator directory
2. Modify CMakeLists.txt to add GPU build option
3. Extract weights from `ESMF_FieldRegridStore()` via `factorList` parameter
4. Store weights using `gpu_regrid_store_weights()`
5. Replace `ESMF_FieldRegrid()` calls with `gpu_regrid_apply()`

## Requirements

### For Hybrid Build (Recommended)
- GCC or Intel compiler (for UFS)
- NVIDIA HPC SDK 24.x (for GPU regrid library)
- CUDA Toolkit 12.x
- NVIDIA GPU (tested on H100, supports cc80/cc90)
- spack-stack environment (provides ESMF, PIO, NetCDF, etc.)

### For Full NVHPC Build
- NVIDIA HPC SDK (nvfortran with OpenACC)
- CUDA Toolkit 12.x
- ESMF 8.8.0+ built with NVHPC
- PIO 2.6.x built with NVHPC (has compatibility issues)
- NetCDF-Fortran built with NVHPC

## API Reference

```fortran
! Initialize GPU module
call gpu_regrid_init(rc)

! Store ESMF weights for GPU use
call gpu_regrid_store_weights(comp_src, comp_dst, mapindex, &
     factorList, factorIndexList, src_size, dst_size, rc)

! Apply GPU regridding (single field)
call gpu_regrid_apply(comp_src, comp_dst, mapindex, &
     src_data, dst_data, rc)

! Apply GPU regridding (multiple fields)
call gpu_regrid_apply_batch(comp_src, comp_dst, mapindex, &
     src_data_2d, dst_data_2d, nfields, rc)

! Check if GPU regridding is available
if (gpu_regrid_available(comp_src, comp_dst, mapindex)) then
  ! Use GPU path
end if

! Cleanup
call gpu_regrid_finalize(rc)
```

## Author

GPU acceleration prototype for UFS coastal coupling (SCHISM + WW3 + DATM).
