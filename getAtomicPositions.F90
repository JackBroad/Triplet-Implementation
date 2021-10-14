program getAtomicPositions
  implicit none
  integer :: nAt = 800, i, j
  double precision :: random
  double precision, allocatable :: positions(:,:)

  allocate(positions(nAt,3))

  open(1, file='AtomicPositions800.txt', status='new')
  write(1,*) nAt

  do i = 1, nAt
    do j = 1, 3
      call random_number(random)
      positions(i,j) = random*100d0
    end do
    write(1,*) positions(i,:)
  end do

  deallocate(positions)
  close(1)

end program getAtomicPositions
