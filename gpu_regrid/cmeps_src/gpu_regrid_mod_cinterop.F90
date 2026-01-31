!===============================================================================
! GPU-Accelerated Regridding Module for CMEPS (C Interop Version)
!
! This version uses C interoperability to call an external GPU library
! (libgpu_regrid.so) built with NVHPC. This allows the rest of UFS to be
! compiled with GCC or Intel while still using GPU acceleration.
!
! Build: Compile with GCC/Intel, link -lgpu_regrid -lnvf -lacchost -laccdevice
!===============================================================================

module gpu_regrid_mod

  use med_kind_mod, only : R8=>SHR_KIND_R8, I4=>SHR_KIND_I4
  use med_kind_mod, only : CL=>SHR_KIND_CL
  use iso_c_binding

  implicit none
  private

  ! Public interfaces
  public :: gpu_regrid_init
  public :: gpu_regrid_store_weights
  public :: gpu_regrid_apply
  public :: gpu_regrid_apply_batch
  public :: gpu_regrid_available
  public :: gpu_regrid_finalize

  ! Maximum number of regrid mappings to cache
  integer, parameter :: MAX_MAPPINGS = 100

  ! C interface to libgpu_regrid.so (built with nvhpc)
  interface
    subroutine c_gpu_regrid_init(rc) bind(C, name="gpu_regrid_init")
      import :: c_int
      integer(c_int), intent(out) :: rc
    end subroutine

    subroutine c_gpu_regrid_finalize(rc) bind(C, name="gpu_regrid_finalize")
      import :: c_int
      integer(c_int), intent(out) :: rc
    end subroutine

    subroutine c_gpu_regrid_store_weights(n_weights, src_size, dst_size, &
         weights, dst_indices, src_indices, map_id, rc) &
         bind(C, name="gpu_regrid_store_weights")
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: n_weights
      integer(c_int), value, intent(in) :: src_size
      integer(c_int), value, intent(in) :: dst_size
      real(c_double), intent(in) :: weights(*)
      integer(c_int), intent(in) :: dst_indices(*)
      integer(c_int), intent(in) :: src_indices(*)
      integer(c_int), intent(out) :: map_id
      integer(c_int), intent(out) :: rc
    end subroutine

    subroutine c_gpu_regrid_apply(map_id, src_size, dst_size, &
         src_data, dst_data, rc) bind(C, name="gpu_regrid_apply")
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: map_id
      integer(c_int), value, intent(in) :: src_size
      integer(c_int), value, intent(in) :: dst_size
      real(c_double), intent(in) :: src_data(*)
      real(c_double), intent(out) :: dst_data(*)
      integer(c_int), intent(out) :: rc
    end subroutine

    subroutine c_gpu_regrid_apply_batch(map_id, src_size, dst_size, n_fields, &
         src_data, dst_data, rc) bind(C, name="gpu_regrid_apply_batch")
      import :: c_int, c_double
      integer(c_int), value, intent(in) :: map_id
      integer(c_int), value, intent(in) :: src_size
      integer(c_int), value, intent(in) :: dst_size
      integer(c_int), value, intent(in) :: n_fields
      real(c_double), intent(in) :: src_data(*)
      real(c_double), intent(out) :: dst_data(*)
      integer(c_int), intent(out) :: rc
    end subroutine

    subroutine c_gpu_regrid_check_gpu(available) bind(C, name="gpu_regrid_check_gpu")
      import :: c_int
      integer(c_int), intent(out) :: available
    end subroutine
  end interface

  ! Mapping info storage (to translate comp_src/comp_dst/mapindex to map_id)
  type :: mapping_info_type
    integer :: comp_src = 0
    integer :: comp_dst = 0
    integer :: mapindex = 0
    integer :: map_id = -1
    integer :: src_size = 0
    integer :: dst_size = 0
    logical :: initialized = .false.
  end type mapping_info_type

  type(mapping_info_type), save :: mappings(MAX_MAPPINGS)
  integer, save :: num_mappings = 0

  ! Module state
  logical, save :: gpu_available = .false.
  logical, save :: module_initialized = .false.

  character(*), parameter :: u_FILE_u = __FILE__

contains

  !===============================================================================
  subroutine gpu_regrid_init(rc)
    ! Initialize the GPU regridding module via C interface
    integer, intent(out) :: rc

    integer(c_int) :: c_available, c_rc

    rc = 0

    if (module_initialized) return

    ! Check GPU availability via C library
    call c_gpu_regrid_check_gpu(c_available)
    gpu_available = (c_available == 1)

    if (gpu_available) then
      call c_gpu_regrid_init(c_rc)
      if (c_rc /= 0) then
        write(*,'(A)') 'GPU_REGRID: WARNING - GPU init failed, using CPU fallback'
        gpu_available = .false.
      else
        write(*,'(A)') 'GPU_REGRID: GPU acceleration enabled via C interop library'
      end if
    else
      write(*,'(A)') 'GPU_REGRID: No GPU available, using CPU fallback'
    end if

    module_initialized = .true.

  end subroutine gpu_regrid_init

  !===============================================================================
  subroutine gpu_regrid_store_weights(comp_src, comp_dst, mapindex, &
       factorList, factorIndexList, src_size, dst_size, rc)
    ! Store ESMF regridding weights for GPU use

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: factorList(:)
    integer(I4), intent(in) :: factorIndexList(:,:)
    integer, intent(in) :: src_size, dst_size
    integer, intent(out) :: rc

    integer :: idx, nnz, i
    integer(c_int) :: c_n_weights, c_src_size, c_dst_size
    integer(c_int) :: c_map_id, c_rc
    real(c_double), allocatable :: c_weights(:)
    integer(c_int), allocatable :: c_dst_indices(:), c_src_indices(:)

    rc = 0

    ! Find or create slot for this mapping
    idx = find_or_create_slot(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No available mapping slots'
      return
    end if

    nnz = size(factorList)

    ! Store GPU weights if available
    if (gpu_available) then
      ! Convert to C arrays
      c_n_weights = int(nnz, c_int)
      c_src_size = int(src_size, c_int)
      c_dst_size = int(dst_size, c_int)

      allocate(c_weights(nnz))
      allocate(c_dst_indices(nnz))
      allocate(c_src_indices(nnz))

      ! Copy data - factorIndexList is (2, nnz) where (1,:) is row and (2,:) is col
      do i = 1, nnz
        c_weights(i) = real(factorList(i), c_double)
        c_dst_indices(i) = int(factorIndexList(1, i), c_int)
        c_src_indices(i) = int(factorIndexList(2, i), c_int)
      end do

      ! Store on GPU via C interface
      call c_gpu_regrid_store_weights(c_n_weights, c_src_size, c_dst_size, &
           c_weights, c_dst_indices, c_src_indices, c_map_id, c_rc)

      deallocate(c_weights, c_dst_indices, c_src_indices)

      if (c_rc /= 0) then
        write(*,'(A)') 'GPU_REGRID: WARNING - GPU weight storage failed'
        gpu_available = .false.
      else
        mappings(idx)%map_id = int(c_map_id)
      end if
    end if

    ! Store mapping info
    mappings(idx)%comp_src = comp_src
    mappings(idx)%comp_dst = comp_dst
    mappings(idx)%mapindex = mapindex
    mappings(idx)%src_size = src_size
    mappings(idx)%dst_size = dst_size
    mappings(idx)%initialized = .true.

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

    integer :: idx, i
    integer(c_int) :: c_map_id, c_src_size, c_dst_size, c_rc
    real(c_double), allocatable :: c_src(:), c_dst(:)

    rc = 0

    ! Find cached mapping
    idx = find_cached_mapping(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No cached mapping found'
      return
    end if

    if (gpu_available .and. mappings(idx)%map_id > 0) then
      ! GPU path via C interface
      c_map_id = int(mappings(idx)%map_id, c_int)
      c_src_size = int(size(src_data), c_int)
      c_dst_size = int(size(dst_data), c_int)

      allocate(c_src(size(src_data)), c_dst(size(dst_data)))
      do i = 1, size(src_data)
        c_src(i) = real(src_data(i), c_double)
      end do

      call c_gpu_regrid_apply(c_map_id, c_src_size, c_dst_size, c_src, c_dst, c_rc)

      if (c_rc == 0) then
        do i = 1, size(dst_data)
          dst_data(i) = real(c_dst(i), R8)
        end do
      else
        rc = -2
      end if

      deallocate(c_src, c_dst)
    else
      ! CPU fallback - should not happen if weights stored correctly
      ! This is a safety net
      rc = -3
      write(*,'(A)') 'GPU_REGRID: WARNING - No GPU mapping, CPU fallback needed'
    end if

  end subroutine gpu_regrid_apply

  !===============================================================================
  subroutine gpu_regrid_apply_batch(comp_src, comp_dst, mapindex, &
       src_data, dst_data, nfields, rc)
    ! Apply GPU regridding to multiple fields at once (batch mode)

    integer, intent(in) :: comp_src, comp_dst, mapindex
    real(R8), intent(in) :: src_data(:,:)
    real(R8), intent(out) :: dst_data(:,:)
    integer, intent(in) :: nfields
    integer, intent(out) :: rc

    integer :: idx, i, k
    integer :: src_size, dst_size
    integer(c_int) :: c_map_id, c_src_size, c_dst_size, c_n_fields, c_rc
    real(c_double), allocatable :: c_src(:,:), c_dst(:,:)

    rc = 0

    ! Find cached mapping
    idx = find_cached_mapping(comp_src, comp_dst, mapindex)
    if (idx < 0) then
      rc = -1
      write(*,'(A)') 'GPU_REGRID: ERROR - No cached mapping found for batch'
      return
    end if

    src_size = size(src_data, 1)
    dst_size = size(dst_data, 1)

    if (gpu_available .and. mappings(idx)%map_id > 0) then
      ! GPU path via C interface
      c_map_id = int(mappings(idx)%map_id, c_int)
      c_src_size = int(src_size, c_int)
      c_dst_size = int(dst_size, c_int)
      c_n_fields = int(nfields, c_int)

      allocate(c_src(src_size, nfields), c_dst(dst_size, nfields))

      ! Copy source data
      do k = 1, nfields
        do i = 1, src_size
          c_src(i, k) = real(src_data(i, k), c_double)
        end do
      end do

      call c_gpu_regrid_apply_batch(c_map_id, c_src_size, c_dst_size, &
           c_n_fields, c_src, c_dst, c_rc)

      if (c_rc == 0) then
        do k = 1, nfields
          do i = 1, dst_size
            dst_data(i, k) = real(c_dst(i, k), R8)
          end do
        end do
      else
        rc = -2
      end if

      deallocate(c_src, c_dst)
    else
      rc = -3
      write(*,'(A)') 'GPU_REGRID: WARNING - No GPU mapping for batch, CPU fallback needed'
    end if

  end subroutine gpu_regrid_apply_batch

  !===============================================================================
  logical function gpu_regrid_available(comp_src, comp_dst, mapindex)
    ! Check if GPU regridding is available for a given mapping

    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: idx

    gpu_regrid_available = .false.

    if (.not. module_initialized) return
    if (.not. gpu_available) return

    idx = find_cached_mapping(comp_src, comp_dst, mapindex)
    if (idx > 0) then
      gpu_regrid_available = (mappings(idx)%initialized .and. &
                              mappings(idx)%map_id > 0)
    end if

  end function gpu_regrid_available

  !===============================================================================
  subroutine gpu_regrid_finalize(rc)
    ! Cleanup GPU resources

    integer, intent(out) :: rc
    integer(c_int) :: c_rc
    integer :: i

    rc = 0

    if (gpu_available) then
      call c_gpu_regrid_finalize(c_rc)
    end if

    ! Clear mapping cache
    do i = 1, num_mappings
      mappings(i)%initialized = .false.
      mappings(i)%map_id = -1
    end do
    num_mappings = 0
    module_initialized = .false.

  end subroutine gpu_regrid_finalize

  !===============================================================================
  ! Private helper functions
  !===============================================================================

  integer function find_or_create_slot(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    ! Check if mapping already exists
    do i = 1, num_mappings
      if (mappings(i)%comp_src == comp_src .and. &
          mappings(i)%comp_dst == comp_dst .and. &
          mappings(i)%mapindex == mapindex) then
        find_or_create_slot = i
        return
      end if
    end do

    ! Create new slot
    if (num_mappings < MAX_MAPPINGS) then
      num_mappings = num_mappings + 1
      find_or_create_slot = num_mappings
    else
      find_or_create_slot = -1
    end if

  end function find_or_create_slot

  !===============================================================================
  integer function find_cached_mapping(comp_src, comp_dst, mapindex)
    integer, intent(in) :: comp_src, comp_dst, mapindex
    integer :: i

    find_cached_mapping = -1
    do i = 1, num_mappings
      if (mappings(i)%comp_src == comp_src .and. &
          mappings(i)%comp_dst == comp_dst .and. &
          mappings(i)%mapindex == mapindex .and. &
          mappings(i)%initialized) then
        find_cached_mapping = i
        return
      end if
    end do

  end function find_cached_mapping

end module gpu_regrid_mod
