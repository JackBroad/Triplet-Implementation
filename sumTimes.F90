! sums over columns in slurm files to get results for benchmarking of speed
program main
  implicit none
  integer :: i, j, i_max=150, j_max=6
  double precision :: dataMatrix(150,6), timeVec(6)
  Character(len=300) :: dataFile = 'slurm-14542977.out'

  ! Read in data file
  dataFile = trim(dataFile)
  open(1, file=dataFile, status='old')
  do i = 1, i_max
    read(1,*) (dataMatrix(i,j), j=1,j_max)
  end do
  close(1)

  ! Sum over columns
  do j = 1, j_max
    timeVec(j) = sum(dataMatrix(:,j))
  end do

  print *, timeVec
  
end program main

