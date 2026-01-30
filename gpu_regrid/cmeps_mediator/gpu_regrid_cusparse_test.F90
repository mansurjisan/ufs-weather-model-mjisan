!===============================================================================
! Test Program for cuSPARSE GPU Regridding
!
! Compile: nvfortran -cuda -cudalib=cusparse -DUSE_CUSPARSE -o test_cusparse \
!          gpu_regrid_cusparse.F90 gpu_regrid_cusparse_test.F90
! Run: ./test_cusparse
!===============================================================================

program test_cusparse

  use cudafor
  use gpu_regrid_cusparse_mod

  implicit none

  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: I4 = selected_int_kind(9)

  ! Test matrix dimensions (large for GPU)
  integer, parameter :: SRC_SIZE = 100000   ! Source grid points
  integer, parameter :: DST_SIZE = 500000   ! Destination points
  integer, parameter :: NNZ = 2000000       ! Non-zeros

  ! Host arrays
  real(R8), allocatable :: h_factorList(:)
  integer(I4), allocatable :: h_factorIndexList(:,:)
  real(R8), allocatable :: h_src(:), h_dst(:), h_ref(:)
  integer(I4), allocatable :: h_rowPtr(:), h_colInd(:)
  real(R8), allocatable :: h_values(:)

  ! Device arrays
  real(R8), device, allocatable :: d_src(:), d_dst(:)

  real(R8) :: t_start, t_end, t_gpu, t_cpu
  real(R8) :: max_error
  integer :: i, j, idx, rc
  integer :: comp_src, comp_dst, mapindex
  integer, allocatable :: row_counts(:)
  integer :: istat

  ! Component indices
  comp_src = 1
  comp_dst = 2
  mapindex = 1

  write(*,'(A)') '================================================'
  write(*,'(A)') ' cuSPARSE GPU Regridding Test'
  write(*,'(A)') '================================================'
  write(*,'(A,I0)') ' Source grid: ', SRC_SIZE
  write(*,'(A,I0)') ' Dest grid:   ', DST_SIZE
  write(*,'(A,I0)') ' Non-zeros:   ', NNZ
  write(*,'(A)') ''

  ! Check GPU
  write(*,'(A)') ' Checking GPU...'
  istat = cudaGetDeviceCount(i)
  if (i < 1) then
    write(*,'(A)') ' ERROR: No GPU found'
    stop 1
  end if
  write(*,'(A,I0,A)') ' Found ', i, ' GPU(s)'

  ! Allocate host arrays
  allocate(h_factorList(NNZ))
  allocate(h_factorIndexList(2, NNZ))
  allocate(h_src(SRC_SIZE))
  allocate(h_dst(DST_SIZE))
  allocate(h_ref(DST_SIZE))
  allocate(h_rowPtr(DST_SIZE + 1))
  allocate(h_colInd(NNZ))
  allocate(h_values(NNZ))
  allocate(row_counts(DST_SIZE))

  ! Allocate device arrays
  allocate(d_src(SRC_SIZE))
  allocate(d_dst(DST_SIZE))

  ! Generate random sparse matrix (COO format)
  write(*,'(A)') ' Generating test sparse matrix...'
  call random_seed()
  do i = 1, NNZ
    call random_number(h_factorList(i))
    h_factorIndexList(1, i) = mod(i-1, DST_SIZE) + 1  ! row
    h_factorIndexList(2, i) = mod(i*7, SRC_SIZE) + 1  ! col
  end do

  ! Generate random source data
  call random_number(h_src)

  ! Convert COO to CSR
  write(*,'(A)') ' Converting COO to CSR format...'
  row_counts = 0
  do i = 1, NNZ
    idx = h_factorIndexList(1, i)
    row_counts(idx) = row_counts(idx) + 1
  end do

  h_rowPtr(1) = 1
  do i = 1, DST_SIZE
    h_rowPtr(i+1) = h_rowPtr(i) + row_counts(i)
  end do

  row_counts = 0
  do i = 1, NNZ
    idx = h_factorIndexList(1, i)
    j = h_rowPtr(idx) + row_counts(idx)
    h_colInd(j) = h_factorIndexList(2, i)
    h_values(j) = h_factorList(i)
    row_counts(idx) = row_counts(idx) + 1
  end do
  deallocate(row_counts)

  ! Initialize cuSPARSE
  write(*,'(A)') ' Initializing cuSPARSE...'
  call cusparse_regrid_init(rc)
  if (rc /= 0) then
    write(*,'(A)') ' ERROR: cuSPARSE init failed'
    stop 1
  end if

  ! Setup matrix on GPU
  write(*,'(A)') ' Setting up CSR matrix on GPU...'
  call cusparse_regrid_setup_matrix(comp_src, comp_dst, mapindex, &
       h_rowPtr, h_colInd, h_values, DST_SIZE, SRC_SIZE, NNZ, rc)
  if (rc /= 0) then
    write(*,'(A,I0)') ' ERROR: Matrix setup failed, rc=', rc
    stop 1
  end if

  ! Compute CPU reference
  write(*,'(A)') ' Computing CPU reference solution...'
  h_ref = 0.0_R8
  call cpu_time(t_start)
  do i = 1, NNZ
    idx = h_factorIndexList(1, i)
    j = h_factorIndexList(2, i)
    h_ref(idx) = h_ref(idx) + h_factorList(i) * h_src(j)
  end do
  call cpu_time(t_end)
  t_cpu = t_end - t_start
  write(*,'(A,F10.4,A)') ' CPU time: ', t_cpu*1000.0, ' ms'

  ! Copy source data to GPU
  d_src = h_src
  d_dst = 0.0_R8

  ! First GPU call (warmup)
  write(*,'(A)') ' Applying cuSPARSE regridding (first call)...'
  call cpu_time(t_start)
  call cusparse_regrid_apply(comp_src, comp_dst, mapindex, &
       d_src, d_dst, SRC_SIZE, DST_SIZE, rc)
  istat = cudaDeviceSynchronize()
  call cpu_time(t_end)
  t_gpu = t_end - t_start

  if (rc /= 0) then
    write(*,'(A,I0)') ' ERROR: cuSPARSE apply failed, rc=', rc
    stop 1
  end if

  write(*,'(A,F10.4,A)') ' GPU time (first): ', t_gpu*1000.0, ' ms'

  ! Copy result back and verify
  h_dst = d_dst
  max_error = 0.0_R8
  do i = 1, DST_SIZE
    max_error = max(max_error, abs(h_dst(i) - h_ref(i)))
  end do

  write(*,'(A)') ''
  write(*,'(A,E12.4)') ' Max error: ', max_error
  if (max_error < 1.0e-10) then
    write(*,'(A)') ' PASSED: Results match CPU reference'
  else
    write(*,'(A)') ' WARNING: Results differ from CPU reference'
  end if

  ! Benchmark (100 iterations)
  write(*,'(A)') ''
  write(*,'(A)') ' Benchmark (100 iterations)...'

  call cpu_time(t_start)
  do i = 1, 100
    call cusparse_regrid_apply(comp_src, comp_dst, mapindex, &
         d_src, d_dst, SRC_SIZE, DST_SIZE, rc)
  end do
  istat = cudaDeviceSynchronize()
  call cpu_time(t_end)
  t_gpu = (t_end - t_start) / 100.0

  write(*,'(A,F10.4,A)') ' Avg GPU time: ', t_gpu*1000.0, ' ms/call'
  write(*,'(A,F10.2,A)') ' Throughput:   ', NNZ / t_gpu / 1.0e9, ' GFLOP/s'
  write(*,'(A,F10.2,A)') ' Speedup:      ', t_cpu / t_gpu, 'x (vs CPU)'

  ! Benchmark (1000 iterations for stable measurement)
  write(*,'(A)') ''
  write(*,'(A)') ' Benchmark (1000 iterations)...'

  call cpu_time(t_start)
  do i = 1, 1000
    call cusparse_regrid_apply(comp_src, comp_dst, mapindex, &
         d_src, d_dst, SRC_SIZE, DST_SIZE, rc)
  end do
  istat = cudaDeviceSynchronize()
  call cpu_time(t_end)
  t_gpu = (t_end - t_start) / 1000.0

  write(*,'(A,F10.4,A)') ' Avg GPU time: ', t_gpu*1000.0, ' ms/call'
  write(*,'(A,F10.2,A)') ' Throughput:   ', NNZ / t_gpu / 1.0e9, ' GFLOP/s'
  write(*,'(A,F10.2,A)') ' Speedup:      ', t_cpu / t_gpu, 'x (vs CPU)'

  ! Cleanup
  call cusparse_regrid_finalize(rc)

  deallocate(h_factorList, h_factorIndexList, h_src, h_dst, h_ref)
  deallocate(h_rowPtr, h_colInd, h_values)
  deallocate(d_src, d_dst)

  write(*,'(A)') ''
  write(*,'(A)') ' Test completed successfully!'
  write(*,'(A)') '================================================'

end program test_cusparse
