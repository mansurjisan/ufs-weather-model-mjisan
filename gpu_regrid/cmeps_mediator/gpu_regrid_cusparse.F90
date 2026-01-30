!===============================================================================
! cuSPARSE Interface for GPU-Accelerated Regridding
!
! This module provides optimized sparse matrix-vector multiply using
! NVIDIA cuSPARSE library for maximum performance.
!
! Requires: NVIDIA HPC SDK with CUDA Fortran and cuSPARSE
! Compile with: nvfortran -cuda -cudalib=cusparse
!===============================================================================

#ifdef USE_CUSPARSE

module gpu_regrid_cusparse_mod

  use cudafor
  use cusparse

  use med_kind_mod, only : R8=>SHR_KIND_R8, I4=>SHR_KIND_I4

  implicit none
  private

  public :: cusparse_regrid_init
  public :: cusparse_regrid_setup_matrix
  public :: cusparse_regrid_apply
  public :: cusparse_regrid_finalize

  ! cuSPARSE handles
  type(cusparseHandle) :: cusparse_handle
  logical :: handle_initialized = .false.

  ! Maximum cached matrices
  integer, parameter :: MAX_MATRICES = 100

  ! GPU matrix storage
  type :: gpu_csr_matrix_type
    logical :: initialized = .false.
    integer :: nrows, ncols, nnz
    integer :: comp_src, comp_dst, mapindex

    ! Device arrays
    integer(I4), device, allocatable :: d_rowPtr(:)
    integer(I4), device, allocatable :: d_colInd(:)
    real(R8), device, allocatable :: d_values(:)

    ! cuSPARSE descriptors
    type(cusparseSpMatDescr) :: matDescr
    type(cusparseDnVecDescr) :: vecDescr_src
    type(cusparseDnVecDescr) :: vecDescr_dst

    ! Buffer for SpMV
    integer(c_size_t) :: bufferSize
    type(c_devptr) :: d_buffer
  end type gpu_csr_matrix_type

  type(gpu_csr_matrix_type), target :: gpu_matrices(MAX_MATRICES)
  integer :: num_matrices = 0

contains

  !===============================================================================
  subroutine cusparse_regrid_init(rc)
    integer, intent(out) :: rc
    integer :: istat

    rc = 0

    if (handle_initialized) return

    ! Create cuSPARSE handle
    istat = cusparseCreate(cusparse_handle)
    if (istat /= CUSPARSE_STATUS_SUCCESS) then
      rc = -1
      write(*,'(A,I0)') 'CUSPARSE_REGRID: ERROR creating handle, status=', istat
      return
    end if

    handle_initialized = .true.
    write(*,'(A)') 'CUSPARSE_REGRID: Initialized cuSPARSE handle'

  end subroutine cusparse_regrid_init

  !===============================================================================
  subroutine cusparse_regrid_setup_matrix(comp_src, comp_dst, mapindex, &
       h_rowPtr, h_colInd, h_values, nrows, ncols, nnz, rc)
    ! Setup CSR matrix on GPU with cuSPARSE descriptors

    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer(I4), intent(in) :: h_rowPtr(:)  ! Host row pointers
    integer(I4), intent(in) :: h_colInd(:)  ! Host column indices
    real(R8), intent(in) :: h_values(:)     ! Host values
    integer, intent(in) :: nrows, ncols, nnz
    integer, intent(out) :: rc

    integer :: idx, istat
    type(gpu_csr_matrix_type), pointer :: mat
    real(R8) :: alpha, beta

    rc = 0
    alpha = 1.0_R8
    beta = 0.0_R8

    ! Find or create slot
    idx = find_or_create_matrix_slot(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      return
    end if

    mat => gpu_matrices(idx)
    mat%nrows = nrows
    mat%ncols = ncols
    mat%nnz = nnz
    mat%comp_src = comp_src
    mat%comp_dst = comp_dst
    mat%mapindex = mapindex

    ! Allocate device arrays
    allocate(mat%d_rowPtr(nrows + 1))
    allocate(mat%d_colInd(nnz))
    allocate(mat%d_values(nnz))

    ! Copy to device
    mat%d_rowPtr = h_rowPtr
    mat%d_colInd = h_colInd
    mat%d_values = h_values

    ! Create sparse matrix descriptor (CSR format)
    istat = cusparseCreateCsr(mat%matDescr, &
         int(nrows, c_int64_t), &
         int(ncols, c_int64_t), &
         int(nnz, c_int64_t), &
         c_devloc(mat%d_rowPtr), &
         c_devloc(mat%d_colInd), &
         c_devloc(mat%d_values), &
         CUSPARSE_INDEX_32I, &
         CUSPARSE_INDEX_32I, &
         CUSPARSE_INDEX_BASE_ONE, &  ! Fortran 1-based indexing
         CUDA_R_64F)

    if (istat /= CUSPARSE_STATUS_SUCCESS) then
      rc = -2
      write(*,'(A,I0)') 'CUSPARSE_REGRID: ERROR creating matrix descriptor, status=', istat
      return
    end if

    mat%initialized = .true.
    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
         'CUSPARSE_REGRID: Setup matrix ', comp_src, '->', comp_dst, &
         ' size=', nrows, 'x', ncols, ' nnz=', nnz

  end subroutine cusparse_regrid_setup_matrix

  !===============================================================================
  subroutine cusparse_regrid_apply(comp_src, comp_dst, mapindex, &
       d_src, d_dst, src_size, dst_size, rc)
    ! Apply SpMV using cuSPARSE (data already on GPU)

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), device, intent(in) :: d_src(:)    ! Source vector on GPU
    real(R8), device, intent(out) :: d_dst(:)   ! Dest vector on GPU
    integer, intent(in) :: src_size, dst_size
    integer, intent(out) :: rc

    integer :: idx, istat
    type(gpu_csr_matrix_type), pointer :: mat
    type(cusparseDnVecDescr) :: vecDescr_src, vecDescr_dst
    real(R8) :: alpha, beta
    integer(c_size_t) :: bufferSize
    type(c_devptr) :: d_buffer

    rc = 0
    alpha = 1.0_R8
    beta = 0.0_R8

    ! Find cached matrix
    idx = find_matrix(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      return
    end if

    mat => gpu_matrices(idx)

    ! Create dense vector descriptors
    istat = cusparseCreateDnVec(vecDescr_src, int(src_size, c_int64_t), &
         c_devloc(d_src), CUDA_R_64F)
    if (istat /= CUSPARSE_STATUS_SUCCESS) then
      rc = -2
      return
    end if

    istat = cusparseCreateDnVec(vecDescr_dst, int(dst_size, c_int64_t), &
         c_devloc(d_dst), CUDA_R_64F)
    if (istat /= CUSPARSE_STATUS_SUCCESS) then
      rc = -3
      return
    end if

    ! Get buffer size
    istat = cusparseSpMV_bufferSize(cusparse_handle, &
         CUSPARSE_OPERATION_NON_TRANSPOSE, &
         alpha, mat%matDescr, vecDescr_src, &
         beta, vecDescr_dst, &
         CUDA_R_64F, &
         CUSPARSE_SPMV_ALG_DEFAULT, &
         bufferSize)

    ! Allocate buffer if needed
    if (bufferSize > 0) then
      istat = cudaMalloc(d_buffer, bufferSize)
    end if

    ! Execute SpMV: dst = alpha * A * src + beta * dst
    istat = cusparseSpMV(cusparse_handle, &
         CUSPARSE_OPERATION_NON_TRANSPOSE, &
         alpha, mat%matDescr, vecDescr_src, &
         beta, vecDescr_dst, &
         CUDA_R_64F, &
         CUSPARSE_SPMV_ALG_DEFAULT, &
         d_buffer)

    if (istat /= CUSPARSE_STATUS_SUCCESS) then
      rc = -4
      write(*,'(A,I0)') 'CUSPARSE_REGRID: ERROR in SpMV, status=', istat
    end if

    ! Cleanup
    istat = cusparseDestroyDnVec(vecDescr_src)
    istat = cusparseDestroyDnVec(vecDescr_dst)
    if (bufferSize > 0) then
      istat = cudaFree(d_buffer)
    end if

  end subroutine cusparse_regrid_apply

  !===============================================================================
  subroutine cusparse_regrid_finalize(rc)
    integer, intent(out) :: rc
    integer :: i, istat

    rc = 0

    ! Destroy matrix descriptors and free device memory
    do i = 1, num_matrices
      if (gpu_matrices(i)%initialized) then
        istat = cusparseDestroySpMat(gpu_matrices(i)%matDescr)
        if (allocated(gpu_matrices(i)%d_rowPtr)) deallocate(gpu_matrices(i)%d_rowPtr)
        if (allocated(gpu_matrices(i)%d_colInd)) deallocate(gpu_matrices(i)%d_colInd)
        if (allocated(gpu_matrices(i)%d_values)) deallocate(gpu_matrices(i)%d_values)
        gpu_matrices(i)%initialized = .false.
      end if
    end do

    ! Destroy cuSPARSE handle
    if (handle_initialized) then
      istat = cusparseDestroy(cusparse_handle)
      handle_initialized = .false.
    end if

    num_matrices = 0

  end subroutine cusparse_regrid_finalize

  !===============================================================================
  ! Private helper functions
  !===============================================================================

  integer function find_or_create_matrix_slot(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    ! Check existing
    do i = 1, num_matrices
      if (gpu_matrices(i)%comp_src == comp_src .and. &
          gpu_matrices(i)%comp_dst == comp_dst .and. &
          gpu_matrices(i)%mapindex == mapindex) then
        find_or_create_matrix_slot = i
        return
      end if
    end do

    ! Create new
    if (num_matrices < MAX_MATRICES) then
      num_matrices = num_matrices + 1
      find_or_create_matrix_slot = num_matrices
    else
      find_or_create_matrix_slot = -1
    end if

  end function find_or_create_matrix_slot

  !===============================================================================
  integer function find_matrix(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    find_matrix = -1
    do i = 1, num_matrices
      if (gpu_matrices(i)%comp_src == comp_src .and. &
          gpu_matrices(i)%comp_dst == comp_dst .and. &
          gpu_matrices(i)%mapindex == mapindex .and. &
          gpu_matrices(i)%initialized) then
        find_matrix = i
        return
      end if
    end do

  end function find_matrix

end module gpu_regrid_cusparse_mod

#endif
