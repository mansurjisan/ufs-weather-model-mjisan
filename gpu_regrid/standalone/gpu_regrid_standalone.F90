!===============================================================================
! Standalone GPU-Accelerated Regridding Module
! No CMEPS/ESMF dependencies - for testing GPU performance
!===============================================================================

module gpu_regrid_standalone

  implicit none
  private

  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: I4 = selected_int_kind(9)

  ! Public interfaces
  public :: gpu_regrid_init
  public :: gpu_regrid_store_weights
  public :: gpu_regrid_apply
  public :: gpu_regrid_apply_batch
  public :: gpu_regrid_available
  public :: gpu_regrid_finalize

  ! Maximum number of regrid mappings to cache
  integer, parameter :: MAX_MAPPINGS = 100

  ! CSR matrix storage for GPU regridding
  type :: csr_matrix_type
    logical :: initialized = .false.
    integer(I4) :: nrows = 0
    integer(I4) :: ncols = 0
    integer(I4) :: nnz = 0
    integer(I4), allocatable :: rowPtr(:)
    integer(I4), allocatable :: colInd(:)
    real(R8), allocatable    :: values(:)
    integer :: comp_src = 0
    integer :: comp_dst = 0
    integer :: mapindex = 0
  end type csr_matrix_type

  ! Cache of CSR matrices
  type(csr_matrix_type), target :: csr_cache(MAX_MAPPINGS)
  integer :: num_cached = 0

  logical :: gpu_available = .false.
  logical :: module_initialized = .false.

contains

  !===============================================================================
  subroutine gpu_regrid_init(rc)
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
  function find_or_create_slot(comp_src, comp_dst, mapindex) result(idx)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: idx, i

    ! Look for existing entry
    do i = 1, num_cached
      if (csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex) then
        idx = i
        return
      end if
    end do

    ! Create new slot
    if (num_cached < MAX_MAPPINGS) then
      num_cached = num_cached + 1
      idx = num_cached
    else
      idx = -1
    end if
  end function find_or_create_slot

  !===============================================================================
  subroutine gpu_regrid_store_weights(comp_src, comp_dst, mapindex, &
       factorList, factorIndexList, src_size, dst_size, rc)

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: factorList(:)
    integer(I4), intent(in) :: factorIndexList(:,:)
    integer, intent(in) :: src_size, dst_size
    integer, intent(out) :: rc

    integer :: idx, nnz, i, row, col, pos
    integer, allocatable :: row_counts(:)
    type(csr_matrix_type), pointer :: csr

    rc = 0

    idx = find_or_create_slot(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No available cache slots'
      return
    end if

    csr => csr_cache(idx)
    nnz = size(factorList)

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

    ! Convert COO to CSR
    allocate(row_counts(dst_size))
    row_counts = 0
    do i = 1, nnz
      row = factorIndexList(1, i)
      if (row >= 1 .and. row <= dst_size) then
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    csr%rowPtr(1) = 1
    do i = 1, dst_size
      csr%rowPtr(i+1) = csr%rowPtr(i) + row_counts(i)
    end do

    row_counts = 0
    do i = 1, nnz
      row = factorIndexList(1, i)
      col = factorIndexList(2, i)
      if (row >= 1 .and. row <= dst_size) then
        pos = csr%rowPtr(row) + row_counts(row)
        csr%colInd(pos) = col
        csr%values(pos) = factorList(i)
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    deallocate(row_counts)

#ifdef _OPENACC
    if (gpu_available) then
      !$acc enter data copyin(csr%rowPtr, csr%colInd, csr%values)
    end if
#endif

    csr%initialized = .true.
    write(*,'(A,I0,A,I0,A,I0)') 'GPU_REGRID: Stored weights, nnz=', nnz, &
         ', src=', src_size, ', dst=', dst_size

  end subroutine gpu_regrid_store_weights

  !===============================================================================
  function gpu_regrid_available(comp_src, comp_dst, mapindex) result(available)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    logical :: available
    integer :: i

    available = .false.
    if (.not. module_initialized) return

    do i = 1, num_cached
      if (csr_cache(i)%initialized .and. &
          csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex) then
        available = .true.
        return
      end if
    end do
  end function gpu_regrid_available

  !===============================================================================
  subroutine gpu_regrid_apply(comp_src, comp_dst, mapindex, src_data, dst_data, rc)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: src_data(:)
    real(R8), intent(inout) :: dst_data(:)
    integer, intent(out) :: rc

    type(csr_matrix_type), pointer :: csr
    integer :: i, j, jstart, jend, idx
    real(R8) :: sum

    rc = 0

    ! Find the cached CSR matrix
    idx = 0
    do i = 1, num_cached
      if (csr_cache(i)%initialized .and. &
          csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex) then
        idx = i
        exit
      end if
    end do

    if (idx == 0) then
      rc = -1
      return
    end if

    csr => csr_cache(idx)

#ifdef _OPENACC
    if (gpu_available) then
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
      ! CPU fallback
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
       src_batch, dst_batch, nfields, rc)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: src_batch(:,:)
    real(R8), intent(inout) :: dst_batch(:,:)
    integer, intent(in) :: nfields
    integer, intent(out) :: rc

    type(csr_matrix_type), pointer :: csr
    integer :: i, j, f, jstart, jend, idx
    real(R8) :: sum

    rc = 0

    idx = 0
    do i = 1, num_cached
      if (csr_cache(i)%initialized .and. &
          csr_cache(i)%comp_src == comp_src .and. &
          csr_cache(i)%comp_dst == comp_dst .and. &
          csr_cache(i)%mapindex == mapindex) then
        idx = i
        exit
      end if
    end do

    if (idx == 0) then
      rc = -1
      return
    end if

    csr => csr_cache(idx)

#ifdef _OPENACC
    if (gpu_available) then
      !$acc data present_or_copyin(src_batch) present_or_copyout(dst_batch) &
      !$acc      present(csr%rowPtr, csr%colInd, csr%values)
      !$acc parallel loop collapse(2) private(sum, jstart, jend)
      do f = 1, nfields
        do i = 1, csr%nrows
          sum = 0.0_R8
          jstart = csr%rowPtr(i)
          jend = csr%rowPtr(i+1) - 1
          !$acc loop reduction(+:sum)
          do j = jstart, jend
            sum = sum + csr%values(j) * src_batch(csr%colInd(j), f)
          end do
          dst_batch(i, f) = sum
        end do
      end do
      !$acc end parallel loop
      !$acc end data
    else
#endif
      do f = 1, nfields
        do i = 1, csr%nrows
          sum = 0.0_R8
          do j = csr%rowPtr(i), csr%rowPtr(i+1) - 1
            sum = sum + csr%values(j) * src_batch(csr%colInd(j), f)
          end do
          dst_batch(i, f) = sum
        end do
      end do
#ifdef _OPENACC
    end if
#endif

  end subroutine gpu_regrid_apply_batch

  !===============================================================================
  subroutine gpu_regrid_finalize(rc)
    integer, intent(out) :: rc
    integer :: i
    type(csr_matrix_type), pointer :: csr

    rc = 0

    do i = 1, num_cached
      csr => csr_cache(i)
      if (csr%initialized) then
#ifdef _OPENACC
        if (gpu_available) then
          !$acc exit data delete(csr%rowPtr, csr%colInd, csr%values)
        end if
#endif
        if (allocated(csr%rowPtr)) deallocate(csr%rowPtr)
        if (allocated(csr%colInd)) deallocate(csr%colInd)
        if (allocated(csr%values)) deallocate(csr%values)
        csr%initialized = .false.
      end if
    end do

    num_cached = 0
    module_initialized = .false.
    write(*,'(A)') 'GPU_REGRID: Module finalized'

  end subroutine gpu_regrid_finalize

end module gpu_regrid_standalone
