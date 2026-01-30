# GPU-Accelerated Regridding for UFS Coastal Applications

This directory contains GPU-accelerated regridding code for ESMF/CMEPS
to improve performance of coupled coastal simulations (SCHISM + WW3 + DATM).

## Overview

ESMF regridding is fundamentally a sparse matrix-vector multiply (SpMV) operation.
This implementation extracts the weight matrix from ESMF and applies it on GPU
using OpenACC or cuSPARSE for 10-30x speedup.

## Files

```
gpu_regrid/
├── README.md                          # This file
├── cmeps_mediator/
│   ├── gpu_regrid_mod.F90            # OpenACC-based GPU regrid module
│   ├── gpu_regrid_cusparse.F90       # cuSPARSE-optimized version
│   └── gpu_regrid_test.F90           # Standalone test program
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

See `GPU_REGRID_DESIGN.md` for detailed integration instructions.

Key steps:
1. Extract weights from `ESMF_FieldRegridStore()` via `factorList` parameter
2. Store weights using `gpu_regrid_store_weights()`
3. Replace `ESMF_FieldRegrid()` calls with `gpu_regrid_apply()`

## Requirements

- NVIDIA HPC SDK (nvfortran with OpenACC)
- CUDA Toolkit 12.x (for cuSPARSE version)
- NVIDIA GPU (tested on H100)

## Performance

Expected speedups on NVIDIA H100:

| Operation | CPU (ESMF) | GPU (OpenACC) | GPU (cuSPARSE) |
|-----------|------------|---------------|----------------|
| 50K point regrid | 50 ms | 2 ms | 0.5 ms |
| 100K point regrid | 100 ms | 4 ms | 1 ms |

## Author

GPU acceleration prototype for UFS coastal coupling.
