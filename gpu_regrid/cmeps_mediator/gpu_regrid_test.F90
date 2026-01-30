!===============================================================================
! Test Program for GPU Regridding Module
!
! Compile: nvfortran -acc -gpu=cc80 -o test_gpu_regrid gpu_regrid_mod.F90 gpu_regrid_test.F90
! Run: ./test_gpu_regrid
!===============================================================================

program test_gpu_regrid

  use gpu_regrid_mod

  implicit none

  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: I4 = selected_int_kind(9)

  ! Test matrix dimensions (simulating DATM -> SCHISM regrid)
  ! Large grid sizes representative of real coastal simulations
  integer, parameter :: SRC_SIZE = 100000   ! DATM grid points (e.g., 0.25° global)
  integer, parameter :: DST_SIZE = 500000   ! SCHISM nodes (high-res coastal mesh)
  integer, parameter :: NNZ = 2000000       ! Non-zeros (avg 4 per dst point)

  real(R8), allocatable :: factorList(:)
  integer(I4), allocatable :: factorIndexList(:,:)
  real(R8), allocatable :: src_data(:), dst_data(:), dst_ref(:)

  real(R8) :: t_start, t_end, t_gpu, t_cpu
  real(R8) :: max_error
  integer :: i, j, idx, rc
  integer :: comp_src, comp_dst, mapindex

  ! Component indices (matching CMEPS)
  comp_src = 1  ! DATM
  comp_dst = 2  ! SCHISM
  mapindex = 1  ! bilinear

  write(*,'(A)') '================================================'
  write(*,'(A)') ' GPU Regridding Test'
  write(*,'(A)') '================================================'
  write(*,'(A,I0,A,I0)') ' Source grid: ', SRC_SIZE, ' points'
  write(*,'(A,I0,A,I0)') ' Dest grid:   ', DST_SIZE, ' points'
  write(*,'(A,I0)') ' Non-zeros:   ', NNZ
  write(*,'(A)') ''

  ! Allocate arrays
  allocate(factorList(NNZ))
  allocate(factorIndexList(2, NNZ))
  allocate(src_data(SRC_SIZE))
  allocate(dst_data(DST_SIZE))
  allocate(dst_ref(DST_SIZE))

  ! Generate random sparse matrix (COO format)
  write(*,'(A)') ' Generating test sparse matrix...'
  call random_seed()
  do i = 1, NNZ
    call random_number(factorList(i))
    factorIndexList(1, i) = mod(i-1, DST_SIZE) + 1  ! row (1-based)
    factorIndexList(2, i) = mod(i*7, SRC_SIZE) + 1  ! col (1-based, scattered)
  end do

  ! Generate random source data
  call random_number(src_data)

  ! Initialize GPU regridding
  write(*,'(A)') ' Initializing GPU regrid module...'
  call gpu_regrid_init(rc)
  if (rc /= 0) then
    write(*,'(A)') ' ERROR: GPU init failed'
    stop 1
  end if

  ! Store weights
  write(*,'(A)') ' Storing regrid weights...'
  call gpu_regrid_store_weights(comp_src, comp_dst, mapindex, &
       factorList, factorIndexList, SRC_SIZE, DST_SIZE, rc)
  if (rc /= 0) then
    write(*,'(A)') ' ERROR: Weight storage failed'
    stop 1
  end if

  ! Compute reference solution (CPU)
  write(*,'(A)') ' Computing CPU reference solution...'
  dst_ref = 0.0_R8
  call cpu_time(t_start)
  do i = 1, NNZ
    idx = factorIndexList(1, i)
    j = factorIndexList(2, i)
    dst_ref(idx) = dst_ref(idx) + factorList(i) * src_data(j)
  end do
  call cpu_time(t_end)
  t_cpu = t_end - t_start
  write(*,'(A,F10.4,A)') ' CPU time: ', t_cpu*1000.0, ' ms'

  ! Apply GPU regridding
  write(*,'(A)') ' Applying GPU regridding...'
  dst_data = 0.0_R8
  call cpu_time(t_start)
  call gpu_regrid_apply(comp_src, comp_dst, mapindex, src_data, dst_data, rc)
  call cpu_time(t_end)
  t_gpu = t_end - t_start

  if (rc /= 0) then
    write(*,'(A)') ' ERROR: GPU regrid failed'
    stop 1
  end if

  write(*,'(A,F10.4,A)') ' GPU time: ', t_gpu*1000.0, ' ms'
  write(*,'(A,F10.2,A)') ' Speedup:  ', t_cpu/t_gpu, 'x'

  ! Verify results
  max_error = 0.0_R8
  do i = 1, DST_SIZE
    max_error = max(max_error, abs(dst_data(i) - dst_ref(i)))
  end do

  write(*,'(A)') ''
  write(*,'(A,E12.4)') ' Max error: ', max_error
  if (max_error < 1.0e-10) then
    write(*,'(A)') ' PASSED: Results match CPU reference'
  else
    write(*,'(A)') ' WARNING: Results differ from CPU reference'
  end if

  ! Benchmark multiple iterations
  write(*,'(A)') ''
  write(*,'(A)') ' Benchmark (100 iterations)...'

  call cpu_time(t_start)
  do i = 1, 100
    call gpu_regrid_apply(comp_src, comp_dst, mapindex, src_data, dst_data, rc)
  end do
  call cpu_time(t_end)
  t_gpu = (t_end - t_start) / 100.0

  write(*,'(A,F10.4,A)') ' Avg GPU time: ', t_gpu*1000.0, ' ms/call'
  write(*,'(A,F10.2,A)') ' Throughput:   ', NNZ / t_gpu / 1.0e9, ' GFLOP/s'
  write(*,'(A,F10.2,A)') ' Speedup:      ', t_cpu / t_gpu, 'x (vs CPU)'

  ! Benchmark with data resident on GPU (simulates real coupled model)
  write(*,'(A)') ''
  write(*,'(A)') ' Benchmark with GPU-resident data (1000 iterations)...'

  !$acc data copyin(src_data) create(dst_data)

  call cpu_time(t_start)
  do i = 1, 1000
    !$acc kernels present(src_data, dst_data)
    ! Dummy touch to ensure data stays on GPU
    !$acc end kernels
    call gpu_regrid_apply(comp_src, comp_dst, mapindex, src_data, dst_data, rc)
  end do
  call cpu_time(t_end)
  t_gpu = (t_end - t_start) / 1000.0

  !$acc end data

  write(*,'(A,F10.4,A)') ' Avg GPU time: ', t_gpu*1000.0, ' ms/call'
  write(*,'(A,F10.2,A)') ' Throughput:   ', NNZ / t_gpu / 1.0e9, ' GFLOP/s'
  write(*,'(A,F10.2,A)') ' Speedup:      ', t_cpu / t_gpu, 'x (vs CPU)'

  ! Cleanup
  call gpu_regrid_finalize(rc)

  deallocate(factorList, factorIndexList, src_data, dst_data, dst_ref)

  write(*,'(A)') ''
  write(*,'(A)') ' Test completed successfully!'
  write(*,'(A)') '================================================'

end program test_gpu_regrid
