program getAtomicPositions
  implicit none
  integer :: nAt = 124, i, j
  double precision :: random, L = 18d0
  double precision, allocatable :: positions(:,:)

  allocate(positions(nAt,3))

  open(1, file='AtomicPositions124SL=18.txt', status='new')
  write(1,*) nAt
  write(1,*) L

  do i = 1, nAt
    do j = 1, 3
      call random_number(random)
      positions(i,j) = random*L
    end do
    write(1,*) positions(i,:)
  end do

  deallocate(positions)
  close(1)

end program getAtomicPositions
