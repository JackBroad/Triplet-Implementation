PROGRAM makeDeltaGamma
implicit none

double precision, allocatable :: X_dg(:,:)
integer :: i, j, k, N_a

! Declare number of atoms
N_a = 400

! Allocate array
allocate(X_dg(N_a,N_a))

! Fill array
do i = 1, N_a

  do j = 1, N_a

    if (i .ne. j) then

      X_dg(i,j) = 0.6

    else

      X_dg(i,j) = 0.0

    end if 

  end do

end do

! Write array to file
open(1, file='deltaGamma400x400.txt', status='new')
write(1,*) N_a
do k = 1, N_a
  write(1,*) X_dg(k,:)
end do
close(1)

! De-allocate array
deallocate(X_dg)

END PROGRAM
