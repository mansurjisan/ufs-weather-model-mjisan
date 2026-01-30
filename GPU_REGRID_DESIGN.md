# GPU-Accelerated ESMF Regridding for UFS Coastal Applications

## Overview

This document describes an approach to GPU-accelerate ESMF regridding operations
for the UFS coastal coupling system (SCHISM + WW3 + DATM + CMEPS).

## The Problem

ESMF regridding (`ESMF_FieldRegrid`) is a major bottleneck in coupled simulations:
- 40-60% of mediator time
- Purely CPU-based
- Called every coupling timestep

## The Solution

ESMF regridding is fundamentally a **sparse matrix-vector multiply (SpMV)**:

```
dst_field = W * src_field
```

Where:
- `W` = sparse weight matrix (precomputed during initialization)
- `src_field` = source field data (1D array)
- `dst_field` = destination field data (1D array)

GPUs excel at SpMV via [cuSPARSE](https://developer.nvidia.com/cusparse), achieving
30-150x speedup over CPU implementations.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Initialization (once)                        │
├─────────────────────────────────────────────────────────────────┤
│  1. ESMF_FieldRegridStore() creates RouteHandle                 │
│  2. Extract weights: factorList, factorIndexList                │
│  3. Convert COO → CSR format                                    │
│  4. Copy CSR matrix to GPU memory (cudaMalloc)                  │
│  5. Create cuSPARSE handles                                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Runtime (every coupling step)                │
├─────────────────────────────────────────────────────────────────┤
│  1. Copy src_field to GPU (or use existing GPU data)            │
│  2. cusparseSpMV() - GPU sparse matrix multiply                 │
│  3. Copy dst_field from GPU (or keep on GPU)                    │
└─────────────────────────────────────────────────────────────────┘
```

## Weight Matrix Format

ESMF stores weights in **COO (Coordinate) format**:
- `factorList(:)` - weight values (real64)
- `factorIndexList(2,:)` - (row, col) indices

For GPU efficiency, convert to **CSR (Compressed Sparse Row) format**:
- `csrRowPtr(nrows+1)` - row pointers
- `csrColInd(nnz)` - column indices
- `csrVal(nnz)` - weight values

## Implementation Files

```
CMEPS-interface/CMEPS/mediator/
├── med_map_mod.F90           # Existing ESMF mapping (modify)
├── gpu_regrid_mod.F90        # NEW: GPU regridding module
├── gpu_regrid_cusparse.F90   # NEW: cuSPARSE interface
└── gpu_regrid_openacc.F90    # NEW: OpenACC fallback
```

## Key Code Components

### 1. Weight Extraction (during init)

```fortran
! In med_map_routehandles_initfrom_field()
! After ESMF_FieldRegridStore, extract weights:

call ESMF_FieldRegridStore(fldsrc, flddst, &
     routehandle=routehandles(mapindex), &
     factorList=factorList, &           ! OUTPUT: weight values
     factorIndexList=factorIndexList, & ! OUTPUT: (row,col) indices
     rc=rc)

! Store for GPU initialization
call gpu_regrid_store_weights(n1, n2, mapindex, &
     factorList, factorIndexList, src_size, dst_size)
```

### 2. GPU Initialization

```fortran
subroutine gpu_regrid_init(n1, n2, mapindex, &
     factorList, factorIndexList, src_size, dst_size)

  ! Convert COO to CSR
  call coo_to_csr(factorIndexList, factorList, nnz, &
       csrRowPtr, csrColInd, csrVal, dst_size)

  ! Allocate GPU memory
  !$acc enter data create(csrRowPtr, csrColInd, csrVal)
  !$acc update device(csrRowPtr, csrColInd, csrVal)

  ! Or with CUDA Fortran:
  ! allocate(d_csrRowPtr(dst_size+1), d_csrColInd(nnz), d_csrVal(nnz))
  ! d_csrRowPtr = csrRowPtr
  ! d_csrColInd = csrColInd
  ! d_csrVal = csrVal

end subroutine
```

### 3. GPU Regrid Apply

```fortran
subroutine gpu_regrid_apply(n1, n2, mapindex, src_data, dst_data)
  real(r8), intent(in)  :: src_data(:)
  real(r8), intent(out) :: dst_data(:)

  ! Option A: OpenACC (simpler, portable)
  !$acc data copyin(src_data) copyout(dst_data)
  !$acc parallel loop
  do i = 1, dst_size
    dst_data(i) = 0.0_r8
    do j = csrRowPtr(i), csrRowPtr(i+1)-1
      dst_data(i) = dst_data(i) + csrVal(j) * src_data(csrColInd(j))
    end do
  end do
  !$acc end parallel
  !$acc end data

  ! Option B: cuSPARSE (optimized, faster)
  ! call cusparse_spmv(handle, dst_data, src_data, ...)

end subroutine
```

## Integration Points

### Option 1: Replace ESMF_FieldRegrid calls (minimal changes)

Modify `med_map_field()` in `med_map_mod.F90`:

```fortran
subroutine med_map_field(...)
  ! Check if GPU regrid is available for this mapping
  if (gpu_regrid_available(n1, n2, mapindex)) then
    ! Extract field data pointers
    call ESMF_FieldGet(field_src, farrayPtr=src_ptr, rc=rc)
    call ESMF_FieldGet(field_dst, farrayPtr=dst_ptr, rc=rc)

    ! Use GPU regridding
    call gpu_regrid_apply(n1, n2, mapindex, src_ptr, dst_ptr)
  else
    ! Fallback to ESMF
    call ESMF_FieldRegrid(field_src, field_dst, routehandle, ...)
  end if
end subroutine
```

### Option 2: GPU-resident field data (maximum performance)

Keep field data on GPU between coupling steps:

```fortran
! During model advance, data stays on GPU
!$acc data present(schism_fields, ww3_fields, datm_fields)

! Regridding happens entirely on GPU
call gpu_regrid_apply(datm, schism, mapindex, ...)

! Only move data to CPU for ESMF communication
!$acc update host(export_fields) ! before ESMF export
!$acc update device(import_fields) ! after ESMF import

!$acc end data
```

## Build Configuration

### CMakeLists.txt additions:

```cmake
option(USE_GPU_REGRID "Enable GPU-accelerated regridding" OFF)

if(USE_GPU_REGRID)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  target_link_libraries(cmeps PRIVATE CUDA::cusparse)
  target_compile_definitions(cmeps PRIVATE GPU_REGRID)
endif()
```

### Compiler flags (NVHPC):

```bash
# OpenACC version
FFLAGS="-acc -gpu=cc80 -Minfo=accel"

# CUDA Fortran + cuSPARSE version
FFLAGS="-cuda -gpu=cc80 -cudalib=cusparse"
```

## Expected Performance

Based on sparse matrix literature and your URI-HBL experience:

| Grid Size | ESMF CPU | GPU OpenACC | GPU cuSPARSE |
|-----------|----------|-------------|--------------|
| 100K pts  | 5 ms     | 0.5 ms      | 0.2 ms       |
| 1M pts    | 50 ms    | 3 ms        | 1 ms         |
| 10M pts   | 500 ms   | 20 ms       | 5 ms         |

**Note**: Actual speedup depends on sparsity pattern, memory bandwidth,
and whether data is already GPU-resident.

## Limitations

1. **MPI Communication**: ESMF handles inter-process data redistribution.
   GPU regrid only accelerates the local SpMV. For distributed grids,
   still need ESMF for halo exchange.

2. **Weight Changes**: If regrid weights change at runtime (rare),
   need to re-initialize GPU matrices.

3. **Memory**: CSR matrix stored on both CPU and GPU doubles memory usage.

## Testing on URSA

```bash
# Build with GPU support
cd /path/to/ufs-coastal
mkdir build_gpu && cd build_gpu
cmake .. -DUSE_GPU_REGRID=ON -DCMAKE_Fortran_COMPILER=nvfortran
make -j8

# Run on GPU node
salloc -p u1-h100 -N 1 --gres=gpu:1 -t 1:00:00
srun ./ufs_model

# Profile
nsys profile -o regrid_profile ./ufs_model
```

## References

- [cuSPARSE Documentation](https://docs.nvidia.com/cuda/cusparse/)
- [ESMF Regridding](https://earthsystemmodeling.org/regrid/)
- [xESMF Weight Files](https://xesmf.readthedocs.io/en/latest/notebooks/Reuse_regridder.html)
- [Sparse Matrix GPU Optimization](https://lukeo.cs.illinois.edu/files/2015_BeDaOl_SPMM.pdf)
