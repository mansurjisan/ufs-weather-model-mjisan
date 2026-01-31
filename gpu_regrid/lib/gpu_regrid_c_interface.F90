!===============================================================================
! GPU Regrid C Interface - Allows calling from GCC-compiled code
! Build with nvhpc, link dynamically at runtime
!===============================================================================
module gpu_regrid_c_interface
  use iso_c_binding
  implicit none

  integer, parameter :: R8 = c_double
  integer, parameter :: I4 = c_int
  integer, parameter :: I8 = c_long

  ! Maximum number of weight matrices to store
  integer, parameter :: MAX_MAPS = 100

  ! Stored weight data in CSR format (for GPU)
  type :: csr_matrix
    integer(I4) :: n_rows = 0
    integer(I4) :: n_cols = 0
    integer(I4) :: nnz = 0
    real(R8), allocatable :: values(:)
    integer(I4), allocatable :: col_idx(:)
    integer(I4), allocatable :: row_ptr(:)
    logical :: on_device = .false.
  end type csr_matrix

  type(csr_matrix), save, target :: weight_matrices(MAX_MAPS)
  integer, save :: n_matrices = 0
  logical, save :: initialized = .false.

  !$acc declare create(weight_matrices)

contains

  !-----------------------------------------------------------------------------
  ! Initialize GPU regrid module (C-callable)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_init_c(rc) bind(C, name="gpu_regrid_init")
    integer(c_int), intent(out) :: rc

    if (initialized) then
      rc = 0
      return
    end if

    n_matrices = 0
    initialized = .true.

    !$acc init device_type(nvidia)

    rc = 0
  end subroutine gpu_regrid_init_c

  !-----------------------------------------------------------------------------
  ! Finalize GPU regrid module (C-callable)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_finalize_c(rc) bind(C, name="gpu_regrid_finalize")
    integer(c_int), intent(out) :: rc
    integer :: i

    do i = 1, n_matrices
      if (weight_matrices(i)%on_device) then
        !$acc exit data delete(weight_matrices(i)%values, &
        !$acc                   weight_matrices(i)%col_idx, &
        !$acc                   weight_matrices(i)%row_ptr)
      end if
      if (allocated(weight_matrices(i)%values)) deallocate(weight_matrices(i)%values)
      if (allocated(weight_matrices(i)%col_idx)) deallocate(weight_matrices(i)%col_idx)
      if (allocated(weight_matrices(i)%row_ptr)) deallocate(weight_matrices(i)%row_ptr)
      weight_matrices(i)%nnz = 0
      weight_matrices(i)%on_device = .false.
    end do

    n_matrices = 0
    initialized = .false.
    rc = 0
  end subroutine gpu_regrid_finalize_c

  !-----------------------------------------------------------------------------
  ! Store weights from COO format (C-callable)
  ! Returns map_id for later use
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_store_weights_c(n_weights, src_size, dst_size, &
       weights, dst_indices, src_indices, map_id, rc) &
       bind(C, name="gpu_regrid_store_weights")

    integer(c_int), value, intent(in) :: n_weights
    integer(c_int), value, intent(in) :: src_size
    integer(c_int), value, intent(in) :: dst_size
    real(c_double), intent(in) :: weights(n_weights)
    integer(c_int), intent(in) :: dst_indices(n_weights)
    integer(c_int), intent(in) :: src_indices(n_weights)
    integer(c_int), intent(out) :: map_id
    integer(c_int), intent(out) :: rc

    integer :: i, row, nnz_per_row
    integer, allocatable :: row_counts(:)
    type(csr_matrix), pointer :: mat

    if (.not. initialized) then
      rc = -1
      return
    end if

    if (n_matrices >= MAX_MAPS) then
      rc = -2
      return
    end if

    n_matrices = n_matrices + 1
    map_id = n_matrices
    mat => weight_matrices(map_id)

    mat%n_rows = dst_size
    mat%n_cols = src_size
    mat%nnz = n_weights

    ! Convert COO to CSR
    allocate(mat%values(n_weights))
    allocate(mat%col_idx(n_weights))
    allocate(mat%row_ptr(dst_size + 1))
    allocate(row_counts(dst_size))

    ! Count entries per row
    row_counts = 0
    do i = 1, n_weights
      row = dst_indices(i)
      if (row >= 1 .and. row <= dst_size) then
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    ! Build row_ptr
    mat%row_ptr(1) = 1
    do i = 1, dst_size
      mat%row_ptr(i+1) = mat%row_ptr(i) + row_counts(i)
    end do

    ! Reset counts for filling
    row_counts = 0

    ! Fill CSR arrays
    do i = 1, n_weights
      row = dst_indices(i)
      if (row >= 1 .and. row <= dst_size) then
        nnz_per_row = mat%row_ptr(row) + row_counts(row)
        mat%values(nnz_per_row) = weights(i)
        mat%col_idx(nnz_per_row) = src_indices(i)
        row_counts(row) = row_counts(row) + 1
      end if
    end do

    deallocate(row_counts)

    ! Copy to GPU
    !$acc enter data copyin(mat%values, mat%col_idx, mat%row_ptr)
    mat%on_device = .true.

    rc = 0
  end subroutine gpu_regrid_store_weights_c

  !-----------------------------------------------------------------------------
  ! Apply regridding (C-callable)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_apply_c(map_id, src_size, dst_size, src_data, dst_data, rc) &
       bind(C, name="gpu_regrid_apply")

    integer(c_int), value, intent(in) :: map_id
    integer(c_int), value, intent(in) :: src_size
    integer(c_int), value, intent(in) :: dst_size
    real(c_double), intent(in) :: src_data(src_size)
    real(c_double), intent(out) :: dst_data(dst_size)
    integer(c_int), intent(out) :: rc

    type(csr_matrix), pointer :: mat
    integer :: i, j, idx
    real(R8) :: sum_val

    if (.not. initialized .or. map_id < 1 .or. map_id > n_matrices) then
      rc = -1
      return
    end if

    mat => weight_matrices(map_id)

    if (.not. mat%on_device) then
      rc = -2
      return
    end if

    ! GPU SpMV using CSR format
    !$acc data copyin(src_data) copyout(dst_data)
    !$acc parallel loop gang vector private(sum_val, j, idx)
    do i = 1, mat%n_rows
      sum_val = 0.0_R8
      do j = mat%row_ptr(i), mat%row_ptr(i+1) - 1
        idx = mat%col_idx(j)
        if (idx >= 1 .and. idx <= src_size) then
          sum_val = sum_val + mat%values(j) * src_data(idx)
        end if
      end do
      dst_data(i) = sum_val
    end do
    !$acc end parallel loop
    !$acc end data

    rc = 0
  end subroutine gpu_regrid_apply_c

  !-----------------------------------------------------------------------------
  ! Apply batch regridding for multiple fields (C-callable)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_apply_batch_c(map_id, src_size, dst_size, n_fields, &
       src_data, dst_data, rc) bind(C, name="gpu_regrid_apply_batch")

    integer(c_int), value, intent(in) :: map_id
    integer(c_int), value, intent(in) :: src_size
    integer(c_int), value, intent(in) :: dst_size
    integer(c_int), value, intent(in) :: n_fields
    real(c_double), intent(in) :: src_data(src_size, n_fields)
    real(c_double), intent(out) :: dst_data(dst_size, n_fields)
    integer(c_int), intent(out) :: rc

    type(csr_matrix), pointer :: mat
    integer :: i, j, k, idx
    real(R8) :: sum_val

    if (.not. initialized .or. map_id < 1 .or. map_id > n_matrices) then
      rc = -1
      return
    end if

    mat => weight_matrices(map_id)

    if (.not. mat%on_device) then
      rc = -2
      return
    end if

    ! Batch GPU SpMV - process all fields together
    !$acc data copyin(src_data) copyout(dst_data)
    !$acc parallel loop gang collapse(2) private(sum_val, j, idx)
    do k = 1, n_fields
      do i = 1, mat%n_rows
        sum_val = 0.0_R8
        !$acc loop seq
        do j = mat%row_ptr(i), mat%row_ptr(i+1) - 1
          idx = mat%col_idx(j)
          if (idx >= 1 .and. idx <= src_size) then
            sum_val = sum_val + mat%values(j) * src_data(idx, k)
          end if
        end do
        dst_data(i, k) = sum_val
      end do
    end do
    !$acc end parallel loop
    !$acc end data

    rc = 0
  end subroutine gpu_regrid_apply_batch_c

  !-----------------------------------------------------------------------------
  ! Check if GPU is available (C-callable)
  !-----------------------------------------------------------------------------
  subroutine gpu_regrid_check_gpu_c(available) bind(C, name="gpu_regrid_check_gpu")
    integer(c_int), intent(out) :: available

    ! OpenACC check - will be 1 if GPU available
    !$acc init device_type(nvidia)
    available = 1
  end subroutine gpu_regrid_check_gpu_c

end module gpu_regrid_c_interface
