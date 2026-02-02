# GPU-Accelerated Regridding for UFS Coastal Applications

## Technical Documentation for Team Review

**Date:** January 31, 2026
**Platform:** NOAA URSA Cluster (H100 GPUs)
**Author:** GPU Acceleration Team

---

## Executive Summary

This document describes the implementation of GPU-accelerated regridding for the Unified Forecast System (UFS) coastal applications. The work enables GPU acceleration of ESMF regridding operations in CMEPS (Community Mediator for Earth Prediction Systems) for coupled coastal simulations involving SCHISM ocean model and DATM atmospheric data.

**Key Achievement:** Successfully built UFS coastal application with hybrid GCC/NVHPC approach, enabling GPU acceleration while maintaining compatibility with the standard spack-stack build environment.

---

## 1. Problem Statement

### 1.1 Background
ESMF regridding in coupled Earth system models is computationally intensive. For coastal applications coupling atmospheric data (DATM) with ocean models (SCHISM), regridding operations can become a performance bottleneck, especially at high resolutions.

### 1.2 Objective
Accelerate ESMF regridding operations using NVIDIA GPUs (H100) on the NOAA URSA cluster to improve the performance of UFS coastal simulations.

### 1.3 Challenges Encountered
1. **Compiler Incompatibility:** NVIDIA's PIO (Parallel I/O) library has fundamental incompatibilities with `nvfortran` due to Fortran int64 constant expressions
2. **Module File Incompatibility:** Fortran `.mod` files are compiler-specific (GCC vs NVHPC)
3. **Dependency Chain:** UFS requires ESMF 8.8.0, PIO 2.6.x, and numerous other libraries
4. **spack-stack Integration:** URSA's spack-stack provides GCC-built libraries only

---

## 2. Solution Architecture

### 2.1 Hybrid Build Approach

After extensive testing, we developed a **hybrid build strategy** that:
- Builds UFS and all dependencies with **GCC** from spack-stack
- Builds the GPU regrid module as a **separate shared library** with **NVHPC**
- Links the two via **C interoperability** (ISO_C_BINDING)

```
┌─────────────────────────────────────────────────────────────────┐
│                     UFS Coastal Application                      │
│                    (Built with GCC/spack-stack)                  │
├─────────────────────────────────────────────────────────────────┤
│  SCHISM  │  DATM/CDEPS  │  CMEPS Mediator  │  ESMF  │  PIO     │
│  (GCC)   │    (GCC)     │     (GCC)        │ (GCC)  │  (GCC)   │
├──────────┴──────────────┴──────────────────┴────────┴──────────┤
│                    C Interop Interface                          │
│                   (ISO_C_BINDING layer)                         │
├─────────────────────────────────────────────────────────────────┤
│              libgpu_regrid.so (Built with NVHPC)                │
│                     OpenACC GPU Kernels                         │
├─────────────────────────────────────────────────────────────────┤
│                    NVIDIA H100 GPU                              │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 Why This Approach?

| Approach | Pros | Cons | Status |
|----------|------|------|--------|
| Full NVHPC Build | Native GPU support | PIO incompatible with nvfortran | Failed |
| Full GCC Build | spack-stack compatible | No GPU support | Baseline |
| **Hybrid GCC + NVHPC** | Best of both worlds | Requires C interop layer | **Implemented** |

---

## 3. Technical Implementation

### 3.1 GPU Regrid Algorithm

ESMF regridding is fundamentally a **Sparse Matrix-Vector Multiply (SpMV)** operation:

```
dst_field = W × src_field
```

Where `W` is a sparse weight matrix extracted from ESMF's RouteHandle.

**GPU Implementation:**
1. Extract weights from ESMF in COO (Coordinate) format
2. Convert to CSR (Compressed Sparse Row) format for GPU efficiency
3. Transfer CSR matrix to GPU memory (one-time cost)
4. Execute SpMV kernel on GPU using OpenACC
5. Support batch mode for multiple fields (improved GPU utilization)

### 3.2 File Structure

```
gpu_regrid/
├── lib/                                    # NVHPC-compiled shared library
│   ├── gpu_regrid_c_interface.F90         # Fortran module with C bindings
│   ├── gpu_regrid.h                       # C header for external linkage
│   ├── build_gpu_regrid_lib.sh            # Build script (SLURM)
│   └── test_c_interface.c                 # C test program
│
├── cmeps_src/                             # Source for CMEPS integration
│   └── gpu_regrid_mod_cinterop.F90        # CMEPS module (C interop version)
│
├── build_ufs_gcc_with_gpu_regrid.sh       # Main UFS build script
├── build_esmf_8.8.sh                      # ESMF build (if needed)
└── GPU_REGRID_DOCUMENTATION.md            # This document
```

### 3.3 Key Source Files

#### gpu_regrid_c_interface.F90 (NVHPC Library)
```fortran
! C-callable interface for GPU regrid operations
subroutine gpu_regrid_store_weights_c(n_weights, src_size, dst_size, &
     weights, dst_indices, src_indices, map_id, rc) bind(C)
  ! Converts COO to CSR, transfers to GPU
end subroutine

subroutine gpu_regrid_apply_c(map_id, src_size, dst_size, &
     src_data, dst_data, rc) bind(C)
  ! Executes SpMV on GPU using OpenACC
  !$acc parallel loop
  do i = 1, dst_size
    ! CSR SpMV kernel
  end do
end subroutine
```

#### gpu_regrid_mod_cinterop.F90 (CMEPS Integration)
```fortran
! Called from CMEPS mediator, compiled with GCC
module gpu_regrid_mod
  use iso_c_binding

  interface
    subroutine c_gpu_regrid_apply(...) bind(C, name="gpu_regrid_apply")
    end subroutine
  end interface

  ! Wrapper functions that call C interface
end module
```

---

## 4. Build Instructions

### 4.1 Prerequisites
- NOAA URSA cluster access
- Account with GPU partition access (`u1-h100`)
- Access to spack-stack: `/contrib/spack-stack/spack-stack-1.9.2/`

### 4.2 Step-by-Step Build

```bash
# Clone repository
cd /scratch5/purged/$USER
git clone https://github.com/mansurjisan/ufs-weather-model-mjisan.git
cd ufs-weather-model-mjisan
git checkout feature/gpu_regrid

# Step 1: Build GPU regrid shared library (NVHPC)
cd gpu_regrid/lib
sbatch build_gpu_regrid_lib.sh
# Wait for completion (~5 minutes)
# Verify: ls /scratch5/purged/$USER/nvhpc_stack/gpu_regrid/lib/libgpu_regrid.so

# Step 2: Build UFS with hybrid approach (GCC + GPU library)
cd /scratch5/purged/$USER/ufs-weather-model-mjisan/gpu_regrid
sbatch build_ufs_gcc_with_gpu_regrid.sh
# Wait for completion (~30-60 minutes)
# Verify: ls ../build_gcc_gpu/ufs_model
```

### 4.3 Runtime Requirements
```bash
# Before running ufs_model:
module load nvhpc/24.11
source /scratch5/purged/$USER/nvhpc_stack/gpu_regrid/setup_gpu_regrid.sh
```

---

## 5. Accomplishments

### 5.1 Completed Tasks

| Task | Status | Notes |
|------|--------|-------|
| ESMF 8.8.0 build with NVHPC | Completed | For standalone testing |
| GPU regrid standalone test | Completed | 2.4x speedup measured |
| GPU regrid shared library (C interop) | Completed | libgpu_regrid.so |
| CMEPS integration module | Completed | gpu_regrid_mod_cinterop.F90 |
| UFS build with hybrid approach | Completed | ufs_model links GPU library |
| spack-stack integration | Completed | Uses URSA's official modules |

### 5.2 Build Verification

```
SUCCESS: ufs_model executable created
-rwxr-xr-x 1 Mansur.Jisan coastal 31471736 Jan 31 16:25 ufs_model

Checking GPU regrid library linkage...
libgpu_regrid.so => /scratch5/purged/Mansur.Jisan/nvhpc_stack/gpu_regrid/lib/libgpu_regrid.so
```

### 5.3 Preliminary Performance Results (Standalone Test)

| Scenario | Time | Speedup vs CPU |
|----------|------|----------------|
| CPU baseline (single field) | 1.00x | - |
| GPU single field (cold) | ~1.5x | Including transfer overhead |
| GPU single field (warm) | **2.4x** | Weights pre-loaded on GPU |
| GPU batch mode (10 fields) | **2.5x+** | Better GPU utilization |

*Note: Full application benchmarks pending*

---

## 6. Architecture Decisions

### 6.1 Why C Interoperability?

Fortran module (`.mod`) files are **compiler-specific**. A module compiled with NVHPC cannot be `use`d by code compiled with GCC. C interoperability via `ISO_C_BINDING` provides:

1. **ABI compatibility** - C calling conventions are standardized
2. **No .mod dependency** - Only needs C header file
3. **Clean separation** - GPU code isolated in shared library

### 6.2 Why OpenACC?

| Option | Pros | Cons |
|--------|------|------|
| **OpenACC** | Portable, Fortran-native, mature | NVIDIA-specific (in practice) |
| CUDA Fortran | Maximum control | NVIDIA-only, more complex |
| cuSPARSE | Optimized library | C interface, less flexible |

OpenACC was chosen for:
- Native Fortran support
- Simpler code maintenance
- Automatic data management
- Good performance for SpMV patterns

### 6.3 CSR Matrix Format

Compressed Sparse Row (CSR) format chosen for GPU SpMV:
- Efficient memory access patterns
- Standard format for GPU sparse libraries
- One-time conversion cost from ESMF's COO format

---

## 7. Known Limitations

1. **NVHPC Runtime Required:** Even though UFS is built with GCC, NVHPC runtime libraries must be loaded at execution time

2. **GPU Memory:** Weight matrices are stored on GPU memory. Very large grids may require memory management

3. **Single GPU:** Current implementation uses single GPU. Multi-GPU would require additional work

4. **Integration Not Activated:** The GPU regrid module is compiled but not yet called from CMEPS mediator. Requires modification of `med_map_mod.F90` to use GPU path

---

## 8. Next Steps

### 8.1 Immediate
- [ ] Run `coastal_ike_shinnecock_atm2sch` test case
- [ ] Verify GPU library loads correctly at runtime
- [ ] Benchmark GPU vs CPU regridding in full application

### 8.2 Short-term
- [ ] Activate GPU regrid calls in CMEPS `med_map_mod.F90`
- [ ] Add runtime flag to enable/disable GPU regridding
- [ ] Performance profiling with NVIDIA Nsight

### 8.3 Future
- [ ] Multi-GPU support for large domains
- [ ] cuSPARSE backend option for optimized SpMV
- [ ] Integration with other UFS components (WW3, MOM6)

---

## 9. Repository Information

| Item | Details |
|------|---------|
| Repository | https://github.com/mansurjisan/ufs-weather-model-mjisan |
| Branch | `feature/gpu_regrid` |
| Base Branch | `feature/coastal_app` |
| Platform | NOAA URSA (Rocky Linux 9, H100 GPUs) |
| spack-stack | 1.9.2 (GCC 12.4.0, OpenMPI 4.1.6) |
| NVHPC | 24.11 |
| CUDA | 12.8.1 |
| ESMF | 8.8.0 |
| PIO | 2.6.2 |

---

## 10. Contact & Support

For questions about this implementation, contact the GPU Acceleration Team.

---

## Appendix A: Module Dependencies

```
spack-stack modules (GCC):
  stack-gcc/12.4.0
  stack-openmpi/4.1.6
  cmake/3.27.9
  esmf/8.8.0
  parallelio/2.6.2
  netcdf-c/4.9.2
  netcdf-fortran/4.6.1
  hdf5/1.14.3
  fms/2024.02
  bacio/2.4.1
  sp/2.5.0
  w3emc/2.10.0

Runtime additions:
  nvhpc/24.11
  cuda/12.8.1
```

## Appendix B: Error Resolution Summary

| Error | Cause | Solution |
|-------|-------|----------|
| PIO nvfortran int64 errors | nvfortran doesn't support int64 constant expressions | Hybrid build approach |
| ESMF .mod incompatibility | Compiler-specific module files | C interoperability |
| SLURM path issues | Scripts run from spool directory | Absolute paths in scripts |
| nvhpc overriding GCC | Module load order | Save GCC paths before loading nvhpc |
| Missing gpu_regrid_mod_cinterop.F90 | File in submodule | Copy from main repo in build script |
