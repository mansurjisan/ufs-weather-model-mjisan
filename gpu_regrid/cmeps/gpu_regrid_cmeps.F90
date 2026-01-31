!===============================================================================
! GPU Regrid CMEPS Integration Module
!
! This module provides GPU-accelerated regridding for CMEPS by calling
! the nvhpc-compiled libgpu_regrid.so via C interoperability.
!
! Usage in CMEPS:
!   1. Call gpu_regrid_cmeps_init() during CMEPS initialization
!   2. Call gpu_regrid_cmeps_store_weights() after ESMF route handle created
!   3. Call gpu_regrid_cmeps_apply() or gpu_regrid_cmeps_apply_batch() for regrid
!   4. Call gpu_regrid_cmeps_finalize() during CMEPS finalization
!
! Build: Compile with GCC/Intel (same as UFS), link -lgpu_regrid
!===============================================================================
module gpu_regrid_cmeps

  use iso_c_binding
  use ESMF

  implicit none
  private

  ! Public procedures
  public :: gpu_regrid_cmeps_init
  public :: gpu_regrid_cmeps_finalize
  public :: gpu_regrid_cmeps_store_from_routehandle
  public :: gpu_regrid_cmeps_apply
  public :: gpu_regrid_cmeps_apply_batch
  public :: gpu_regrid_cmeps_is_available

  ! C interface to libgpu_regrid.so
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

  ! Module state
  logical, save :: module_initialized = .false.
  logical, save :: gpu_available = .false.

contains

  !-----------------------------------------------------------------------------
  ! Check if GPU regridding is available
  !-----------------------------------------------------------------------------
  logical function gpu_regrid_cmeps_is_available()
    integer(c_int) :: available

    if (.not. module_initialized) then
      call c_gpu_regrid_check_gpu(available)
      gpu_available = (available == 1)
      module_initialized = .true.
    end if

    gpu_regrid_cmeps_is_available = gpu_available
  end function

  !-----------------------------------------------------------------------------
  ! Initialize GPU regrid module
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_cmeps_init(rc)
    integer, intent(out) :: rc

    integer(c_int) :: c_rc

    rc = ESMF_SUCCESS

    if (.not. gpu_regrid_cmeps_is_available()) then
      rc = ESMF_RC_NOT_IMPL
      return
    end if

    call c_gpu_regrid_init(c_rc)
    if (c_rc /= 0) then
      rc = ESMF_FAILURE
      return
    end if

  end subroutine gpu_regrid_cmeps_init

  !-----------------------------------------------------------------------------
  ! Finalize GPU regrid module
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_cmeps_finalize(rc)
    integer, intent(out) :: rc

    integer(c_int) :: c_rc

    rc = ESMF_SUCCESS

    call c_gpu_regrid_finalize(c_rc)
    if (c_rc /= 0) then
      rc = ESMF_FAILURE
    end if

  end subroutine gpu_regrid_cmeps_finalize

  !-----------------------------------------------------------------------------
  ! Extract weights from ESMF RouteHandle and store on GPU
  ! Returns map_id for use in apply functions
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_cmeps_store_from_routehandle(routehandle, srcField, dstField, &
       map_id, rc)
    type(ESMF_RouteHandle), intent(in) :: routehandle
    type(ESMF_Field), intent(in) :: srcField
    type(ESMF_Field), intent(in) :: dstField
    integer, intent(out) :: map_id
    integer, intent(out) :: rc

    ! ESMF weight extraction variables
    integer :: factorCount, srcTermCount, dstTermCount
    real(ESMF_KIND_R8), allocatable :: factorList(:)
    integer, allocatable :: factorIndexList(:,:)

    ! Grid size variables
    type(ESMF_Grid) :: srcGrid, dstGrid
    integer :: srcSize, dstSize
    integer :: localDeCount

    ! C interface variables
    integer(c_int) :: c_n_weights, c_src_size, c_dst_size
    integer(c_int) :: c_map_id, c_rc
    real(c_double), allocatable :: c_weights(:)
    integer(c_int), allocatable :: c_dst_indices(:), c_src_indices(:)
    integer :: i

    rc = ESMF_SUCCESS
    map_id = -1

    ! Get weight information from route handle
    call ESMF_FieldRegridGetArea(srcField, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    ! Get the sparse matrix (SMM) factors from the route handle
    ! This requires ESMF 8.1.0+
    call ESMF_RouteHandleGet(routehandle, &
         sparseMatMulInfo_factorCount=factorCount, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    if (factorCount == 0) then
      rc = ESMF_RC_NOT_FOUND
      return
    end if

    allocate(factorList(factorCount))
    allocate(factorIndexList(2, factorCount))

    call ESMF_RouteHandleGet(routehandle, &
         sparseMatMulInfo_factorList=factorList, &
         sparseMatMulInfo_factorIndexList=factorIndexList, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      deallocate(factorList, factorIndexList)
      return
    end if

    ! Get grid sizes
    call ESMF_FieldGet(srcField, grid=srcGrid, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call ESMF_FieldGet(dstField, grid=dstGrid, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    ! Get total element counts
    call ESMF_GridGet(srcGrid, localDeCount=localDeCount, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    ! Note: For production, need to get actual grid sizes
    ! This is simplified - actual implementation needs proper size queries
    srcSize = 0
    dstSize = 0
    ! TODO: Get actual sizes from distgrid

    ! For now, infer from index list
    srcSize = maxval(factorIndexList(1,:))
    dstSize = maxval(factorIndexList(2,:))

    ! Convert to C arrays
    c_n_weights = int(factorCount, c_int)
    c_src_size = int(srcSize, c_int)
    c_dst_size = int(dstSize, c_int)

    allocate(c_weights(factorCount))
    allocate(c_src_indices(factorCount))
    allocate(c_dst_indices(factorCount))

    do i = 1, factorCount
      c_weights(i) = real(factorList(i), c_double)
      c_src_indices(i) = int(factorIndexList(1, i), c_int)
      c_dst_indices(i) = int(factorIndexList(2, i), c_int)
    end do

    ! Store on GPU
    call c_gpu_regrid_store_weights(c_n_weights, c_src_size, c_dst_size, &
         c_weights, c_dst_indices, c_src_indices, c_map_id, c_rc)

    map_id = int(c_map_id)

    if (c_rc /= 0) then
      rc = ESMF_FAILURE
    end if

    ! Cleanup
    deallocate(factorList, factorIndexList)
    deallocate(c_weights, c_src_indices, c_dst_indices)

  end subroutine gpu_regrid_cmeps_store_from_routehandle

  !-----------------------------------------------------------------------------
  ! Apply GPU regridding to a single field
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_cmeps_apply(map_id, srcField, dstField, rc)
    integer, intent(in) :: map_id
    type(ESMF_Field), intent(in) :: srcField
    type(ESMF_Field), intent(inout) :: dstField
    integer, intent(out) :: rc

    real(ESMF_KIND_R8), pointer :: srcPtr(:), dstPtr(:)
    integer(c_int) :: c_map_id, c_src_size, c_dst_size, c_rc
    real(c_double), allocatable :: c_src(:), c_dst(:)
    integer :: srcSize, dstSize, i

    rc = ESMF_SUCCESS

    ! Get field data pointers
    call ESMF_FieldGet(srcField, farrayPtr=srcPtr, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    call ESMF_FieldGet(dstField, farrayPtr=dstPtr, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    srcSize = size(srcPtr)
    dstSize = size(dstPtr)

    ! Copy to C arrays
    allocate(c_src(srcSize), c_dst(dstSize))
    do i = 1, srcSize
      c_src(i) = real(srcPtr(i), c_double)
    end do

    ! Call GPU regrid
    c_map_id = int(map_id, c_int)
    c_src_size = int(srcSize, c_int)
    c_dst_size = int(dstSize, c_int)

    call c_gpu_regrid_apply(c_map_id, c_src_size, c_dst_size, c_src, c_dst, c_rc)

    if (c_rc /= 0) then
      rc = ESMF_FAILURE
      deallocate(c_src, c_dst)
      return
    end if

    ! Copy back to Fortran array
    do i = 1, dstSize
      dstPtr(i) = real(c_dst(i), ESMF_KIND_R8)
    end do

    deallocate(c_src, c_dst)

  end subroutine gpu_regrid_cmeps_apply

  !-----------------------------------------------------------------------------
  ! Apply GPU regridding to multiple fields (batch mode)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_cmeps_apply_batch(map_id, srcFields, dstFields, n_fields, rc)
    integer, intent(in) :: map_id
    type(ESMF_Field), intent(in) :: srcFields(:)
    type(ESMF_Field), intent(inout) :: dstFields(:)
    integer, intent(in) :: n_fields
    integer, intent(out) :: rc

    real(ESMF_KIND_R8), pointer :: srcPtr(:), dstPtr(:)
    integer(c_int) :: c_map_id, c_src_size, c_dst_size, c_n_fields, c_rc
    real(c_double), allocatable :: c_src(:,:), c_dst(:,:)
    integer :: srcSize, dstSize, i, k

    rc = ESMF_SUCCESS

    if (n_fields < 1) return

    ! Get sizes from first field
    call ESMF_FieldGet(srcFields(1), farrayPtr=srcPtr, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    srcSize = size(srcPtr)

    call ESMF_FieldGet(dstFields(1), farrayPtr=dstPtr, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    dstSize = size(dstPtr)

    ! Allocate batch arrays (column-major for Fortran)
    allocate(c_src(srcSize, n_fields), c_dst(dstSize, n_fields))

    ! Copy all source fields
    do k = 1, n_fields
      call ESMF_FieldGet(srcFields(k), farrayPtr=srcPtr, rc=rc)
      if (rc /= ESMF_SUCCESS) then
        deallocate(c_src, c_dst)
        return
      end if
      do i = 1, srcSize
        c_src(i, k) = real(srcPtr(i), c_double)
      end do
    end do

    ! Call GPU batch regrid
    c_map_id = int(map_id, c_int)
    c_src_size = int(srcSize, c_int)
    c_dst_size = int(dstSize, c_int)
    c_n_fields = int(n_fields, c_int)

    call c_gpu_regrid_apply_batch(c_map_id, c_src_size, c_dst_size, c_n_fields, &
         c_src, c_dst, c_rc)

    if (c_rc /= 0) then
      rc = ESMF_FAILURE
      deallocate(c_src, c_dst)
      return
    end if

    ! Copy back to destination fields
    do k = 1, n_fields
      call ESMF_FieldGet(dstFields(k), farrayPtr=dstPtr, rc=rc)
      if (rc /= ESMF_SUCCESS) then
        deallocate(c_src, c_dst)
        return
      end if
      do i = 1, dstSize
        dstPtr(i) = real(c_dst(i, k), ESMF_KIND_R8)
      end do
    end do

    deallocate(c_src, c_dst)

  end subroutine gpu_regrid_cmeps_apply_batch

end module gpu_regrid_cmeps
