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
├── cmeps_mediator/
│   ├── gpu_regrid_mod.F90            # OpenACC-based GPU regrid module
│   ├── gpu_regrid_test.F90           # Standalone test program
│   └── gpu_regrid_cusparse.F90       # cuSPARSE version (experimental)
└── (see also)
    ├── ../GPU_REGRID_DESIGN.md       # Detailed design document
    └── ../scripts/build_gpu_regrid.sh # Build script for URSA
```

## Installation

### Option 1: Copy to CMEPS mediator directory

```bash
cp gpu_regrid/cmeps_mediator/*.F90 CMEPS-interface/CMEPS/mediator/
```

### Option 2: Build standalone test

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

- NVIDIA HPC SDK (nvfortran with OpenACC)
- CUDA Toolkit 12.x
- NVIDIA GPU (tested on H100, supports cc80/cc90)

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
