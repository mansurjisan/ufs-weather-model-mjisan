!===============================================================================
! GPU-Accelerated Regridding Module for CMEPS
!
! This module provides GPU-accelerated sparse matrix multiply operations
! to replace CPU-based ESMF regridding for improved performance.
!
! Author: GPU Acceleration Prototype
! Target: NVIDIA GPUs via OpenACC or cuSPARSE
!===============================================================================

module gpu_regrid_mod

  use med_kind_mod, only : R8=>SHR_KIND_R8, I4=>SHR_KIND_I4
  use med_kind_mod, only : CL=>SHR_KIND_CL

  implicit none
  private

  ! Public interfaces
  public :: gpu_regrid_init
  public :: gpu_regrid_store_weights
  public :: gpu_regrid_apply
  public :: gpu_regrid_apply_batch   ! New: batch multiple fields
  public :: gpu_regrid_available
  public :: gpu_regrid_finalize

  ! Maximum number of regrid mappings to cache
  integer, parameter :: MAX_MAPPINGS = 100

  ! CSR matrix storage for GPU regridding
  type :: csr_matrix_type
    logical :: initialized = .false.
    integer(I4) :: nrows = 0        ! Number of destination points
    integer(I4) :: ncols = 0        ! Number of source points
    integer(I4) :: nnz = 0          ! Number of non-zeros
    integer(I4), allocatable :: rowPtr(:)   ! CSR row pointers (nrows+1)
    integer(I4), allocatable :: colInd(:)   ! CSR column indices (nnz)
    real(R8), allocatable    :: values(:)   ! CSR values (nnz)
    ! Mapping identifiers
    integer :: comp_src = 0
    integer :: comp_dst = 0
    integer :: mapindex = 0
  end type csr_matrix_type

  ! Cache of CSR matrices for different mappings
  type(csr_matrix_type), target :: csr_cache(MAX_MAPPINGS)
  integer :: num_cached = 0

  ! Module variables
  logical :: gpu_available = .false.
  logical :: module_initialized = .false.

  character(*), parameter :: u_FILE_u = __FILE__

contains

  !===============================================================================
  subroutine gpu_regrid_init(rc)
    ! Initialize the GPU regridding module
    ! Check for GPU availability and initialize runtime
#ifdef _OPENACC
    use openacc
#endif
    integer, intent(out) :: rc
#ifdef _OPENACC
    integer :: num_devices
#endif

    rc = 0

    if (module_initialized) return

#ifdef _OPENACC
    ! Check for available GPUs
    num_devices = acc_get_num_devices(acc_device_nvidia)
    if (num_devices > 0) then
      gpu_available = .true.
      call acc_init(acc_device_nvidia)
      write(*,'(A,I0,A)') 'GPU_REGRID: Found ', num_devices, ' NVIDIA GPU(s)'
    else
      gpu_available = .false.
      write(*,'(A)') 'GPU_REGRID: No GPU available, will use CPU fallback'
    end if
#else
    gpu_available = .false.
    write(*,'(A)') 'GPU_REGRID: Compiled without OpenACC, using CPU fallback'
#endif

    module_initialized = .true.

  end subroutine gpu_regrid_init

  !===============================================================================
  subroutine gpu_regrid_store_weights(comp_src, comp_dst, mapindex, &
       factorList, factorIndexList, src_size, dst_size, rc)
    ! Store ESMF regridding weights for GPU use
    ! Converts COO format to CSR format

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: factorList(:)          ! Weight values
    integer(I4), intent(in) :: factorIndexList(:,:) ! (2, nnz): row, col indices
    integer, intent(in) :: src_size, dst_size
    integer, intent(out) :: rc

    integer :: idx, nnz, i, row, col
    integer, allocatable :: row_counts(:)
    type(csr_matrix_type), pointer :: csr

    rc = 0

    ! Find available slot or existing entry
    idx = find_or_create_slot(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No available cache slots'
      return
    end if

    csr => csr_cache(idx)

    ! Get number of non-zeros
    nnz = size(factorList)

    ! Allocate CSR arrays
    if (allocated(csr%rowPtr)) deallocate(csr%rowPtr)
    if (allocated(csr%colInd)) deallocate(csr%colInd)
    if (allocated(csr%values)) deallocate(csr%values)

    allocate(csr%rowPtr(dst_size + 1))
    allocate(csr%colInd(nnz))
    allocate(csr%values(nnz))

    csr%nrows = dst_size
    csr%ncols = src_size
    csr%nnz = nnz
    csr%comp_src = comp_src
    csr%comp_dst = comp_dst
    csr%mapindex = mapindex

    ! Convert COO to CSR format
    ! Step 1: Count entries per row
    allocate(row_counts(dst_size))
    row_counts = 0
    do i = 1, nnz
      row = factorIndexList(1, i)  ! ESMF uses 1-based indexing
      if (row >= 1 .and. row <= dst_size) then
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    ! Step 2: Build row pointers
    csr%rowPtr(1) = 1
    do i = 1, dst_size
      csr%rowPtr(i+1) = csr%rowPtr(i) + row_counts(i)
    end do

    ! Step 3: Fill column indices and values
    row_counts = 0  ! Reuse as insertion counters
    do i = 1, nnz
      row = factorIndexList(1, i)
      col = factorIndexList(2, i)
      if (row >= 1 .and. row <= dst_size) then
        idx = csr%rowPtr(row) + row_counts(row)
        csr%colInd(idx) = col
        csr%values(idx) = factorList(i)
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    deallocate(row_counts)

    ! Copy CSR data to GPU
#ifdef _OPENACC
    if (gpu_available) then
      !$acc enter data copyin(csr%rowPtr, csr%colInd, csr%values)
      write(*,'(A,I0,A,I0,A,I0)') 'GPU_REGRID: Stored CSR matrix on GPU (', &
           csr%nrows, ' x ', csr%ncols, ', nnz=', csr%nnz
    end if
#endif

    csr%initialized = .true.

    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
         'GPU_REGRID: Created regrid mapping ', comp_src, '->', comp_dst, &
         ' mapindex=', mapindex, ' nnz=', nnz

  end subroutine gpu_regrid_store_weights

  !===============================================================================
  subroutine gpu_regrid_apply(comp_src, comp_dst, mapindex, &
       src_data, dst_data, rc)
    ! Apply GPU regridding (sparse matrix-vector multiply)

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: src_data(:)
    real(R8), intent(out) :: dst_data(:)
    integer, intent(out) :: rc

    integer :: idx, i, j, jstart, jend
    type(csr_matrix_type), pointer :: csr
    real(R8) :: sum

    rc = 0

    ! Find cached CSR matrix
    idx = find_cached_matrix(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No cached matrix found'
      return
    end if

    csr => csr_cache(idx)

#ifdef _OPENACC
    if (gpu_available) then
      ! GPU path: OpenACC SpMV
      ! Use present_or_* to avoid redundant transfers if data already on GPU
      !$acc data present_or_copyin(src_data) present_or_copyout(dst_data) &
      !$acc      present(csr%rowPtr, csr%colInd, csr%values)

      !$acc parallel loop private(sum, jstart, jend)
      do i = 1, csr%nrows
        sum = 0.0_R8
        jstart = csr%rowPtr(i)
        jend = csr%rowPtr(i+1) - 1
        !$acc loop reduction(+:sum)
        do j = jstart, jend
          sum = sum + csr%values(j) * src_data(csr%colInd(j))
        end do
        dst_data(i) = sum
      end do
      !$acc end parallel loop

      !$acc end data
    else
#endif
      ! CPU fallback path
      do i = 1, csr%nrows
        sum = 0.0_R8
        do j = csr%rowPtr(i), csr%rowPtr(i+1) - 1
          sum = sum + csr%values(j) * src_data(csr%colInd(j))
        end do
        dst_data(i) = sum
      end do
#ifdef _OPENACC
    end if
#endif

  end subroutine gpu_regrid_apply

  !===============================================================================
  subroutine gpu_regrid_apply_batch(comp_src, comp_dst, mapindex, &
       src_data, dst_data, nfields, rc)
    ! Apply GPU regridding to multiple fields at once (batch mode)
    ! More efficient than calling gpu_regrid_apply multiple times

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: src_data(:,:)    ! (src_size, nfields)
    real(R8), intent(out) :: dst_data(:,:)   ! (dst_size, nfields)
    integer, intent(in) :: nfields
    integer, intent(out) :: rc

    integer :: idx, i, j, k, jstart, jend
    type(csr_matrix_type), pointer :: csr
    real(R8) :: sum

    rc = 0

    ! Find cached CSR matrix
    idx = find_cached_matrix(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No cached matrix found for batch'
      return
    end if

    csr => csr_cache(idx)

#ifdef _OPENACC
    if (gpu_available) then
      ! GPU path: Process all fields in parallel
      !$acc data present_or_copyin(src_data) present_or_copyout(dst_data) &
      !$acc      present(csr%rowPtr, csr%colInd, csr%values)

      !$acc parallel loop collapse(2) private(sum, jstart, jend)
      do k = 1, nfields
        do i = 1, csr%nrows
          sum = 0.0_R8
          jstart = csr%rowPtr(i)
          jend = csr%rowPtr(i+1) - 1
          !$acc loop reduction(+:sum)
          do j = jstart, jend
            sum = sum + csr%values(j) * src_data(csr%colInd(j), k)
          end do
          dst_data(i, k) = sum
        end do
      end do
      !$acc end parallel loop

      !$acc end data
    else
#endif
      ! CPU fallback path
      do k = 1, nfields
        do i = 1, csr%nrows
          sum = 0.0_R8
          do j = csr%rowPtr(i), csr%rowPtr(i+1) - 1
            sum = sum + csr%values(j) * src_data(csr%colInd(j), k)
          end do
          dst_data(i, k) = sum
        end do
      end do
#ifdef _OPENACC
    end if
#endif

  end subroutine gpu_regrid_apply_batch

  !===============================================================================
  logical function gpu_regrid_available(comp_src, comp_dst, mapindex)
    ! Check if GPU regridding is available for a given mapping

    integer, intent(in) :: comp_src, comp_dst, mapindex

    integer :: idx

    gpu_regrid_available = .false.

    if (.not. module_initialized) return
    if (.not. gpu_available) return

    idx = find_cached_matrix(comp_src, comp_dst, mapindex)
    if (idx > 0) then
      gpu_regrid_available = csr_cache(idx)%initialized
    end if

  end function gpu_regrid_available

  !===============================================================================
  subroutine gpu_regrid_finalize(rc)
    ! Cleanup GPU resources

    integer, intent(out) :: rc
    integer :: i

    rc = 0

    do i = 1, num_cached
      if (csr_cache(i)%initialized) then
#ifdef _OPENACC
        if (gpu_available) then
          !$acc exit data delete(csr_cache(i)%rowPtr, &
          !$acc                  csr_cache(i)%colInd, &
          !$acc                  csr_cache(i)%values)
        end if
#endif
        if (allocated(csr_cache(i)%rowPtr)) deallocate(csr_cache(i)%rowPtr)
        if (allocated(csr_cache(i)%colInd)) deallocate(csr_cache(i)%colInd)
        if (allocated(csr_cache(i)%values)) deallocate(csr_cache(i)%values)
        csr_cache(i)%initialized = .false.
      end if
    end do

    num_cached = 0
    module_initialized = .false.

  end subroutine gpu_regrid_finalize

  !===============================================================================
  ! Private helper functions
  !===============================================================================

  integer function find_or_create_slot(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    ! First check if mapping already exists
    do i = 1, num_cached
      if (csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex) then
        find_or_create_slot = i
        return
      end if
    end do

    ! Create new slot
    if (num_cached < MAX_MAPPINGS) then
      num_cached = num_cached + 1
      find_or_create_slot = num_cached
    else
      find_or_create_slot = -1
    end if

  end function find_or_create_slot

  !===============================================================================
  integer function find_cached_matrix(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    find_cached_matrix = -1
    do i = 1, num_cached
      if (csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex .and. &
          csr_cache(i)%initialized) then
        find_cached_matrix = i
        return
      end if
    end do

  end function find_cached_matrix

end module gpu_regrid_mod
