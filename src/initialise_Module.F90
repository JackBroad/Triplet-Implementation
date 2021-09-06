module initialise_Module
  use GP_variables
  implicit none


contains


! Subroutine to set up non-additive GP
subroutine initialise_GP_NonAdd(hyperParametersFile, alphaFile, trainingSetFile)
  implicit none
  Character(len=300) :: hyperParametersFile , alphaFile, trainingSetFile
  integer :: i, j, k, l

  ! Read in hyperparameters
  hyperParametersFile = trim( hyperParametersFile )
  open(1, file= hyperParametersFile , status='old')
  do i = 1, 3
    read(1,*) hyperParams(i)
  end do
  close(1)

  ! Get the no of TPs and alpha values for each
  alphaFile = trim( alphaFile)
  open(2, file=alphaFile, status='old')
  read(2,*) N_tp
  if (allocated(alpha)) then
    deallocate(alpha)
  end if
  allocate(alpha(N_tp))
  read(2,*) (alpha(j), j=1,N_tp)
  close(2)

  ! Read in nArgs and the distances for each TP
  trainingSetFile = trim( trainingSetFile )
  open(3, file=trainingSetFile, status='old')
  read(3,*) nArgs
  if (allocated(trainData)) then
    deallocate(trainData)
  end if
  allocate(trainData(N_tp,nArgs))
  do k = 1, N_tp
    read(3,*) (trainData(k,l), l=1,nArgs)
  end do
  close(3)

return
end subroutine initialise_GP_NonAdd


subroutine initialise_Positions(fileName, posAt,N_a)
  implicit none
  character (len=40), intent(in) :: fileName
  double precision, allocatable, intent(out) :: posAt(:,:)
  integer, intent(out) :: N_a
  integer :: m, n

  ! Read in the number of atoms
  open(1, file=fileName, status='old')
  read(1,*) N_a

  ! Get the atomic positions
  allocate(posAt(N_a,nArgs))
  do m = 1, N_a
    read(1,*) (posAt(m,n), n=1,nArgs)
  end do
  close(1)

return
end subroutine initialise_Positions


subroutine initialise_Variables(N_a, N_tri,udSize)
  implicit none
  integer, intent(in) :: N_a
  integer, intent(out) :: N_tri, udSize

  ! Calculate the number of triplets
  N_tri = N_a**3
  N_tri = N_tri - 3*N_a**2
  N_tri = N_tri + 2*N_a
  N_tri = N_tri / 6 ! No. of triplets

  ! Determine size of UD of X_dg
  udSize = ((N_a * N_a) - N_a) / 2

return
end subroutine initialise_Variables


! Moves an atom at random and returns the new atomic positions and the
! index of the moved atom.
! Takes Cartesian atomic positions (pos), the no. of atoms (num)
! and the max. distance to move in any direction (dMax) as arguments
subroutine initialise_Move(pos,num,dMax,  newPos,mover)
  implicit none
  integer, intent(in) :: num
  double precision, intent(in) :: pos(num,3), dMax
  integer, intent(out) :: mover
  double precision, intent(out) :: newPos(num,3)
  integer :: icol, irow
  double precision :: randNo

  ! Pick an atom at random to move (can't generate random int in fortran, hence
  ! the work-around)
  call random_number(randNo)
  mover = 1 + FLOOR(num*randNo)

  ! Change the position of the atom in newPos
  do irow = 1, num
    do icol = 1, 3
      ! Set each element of newPos to be the same as the equiv. in pos
      newPos(irow,icol) = pos(irow,icol)
      ! If in row corresponding to moving atom, change coords of this atom in
      ! newPos
      if (irow .eq. mover) then
        call random_number(randNo)
        newPos(irow,icol) = newPos(irow,icol) + (2.0*randNo-1) * dMax
      end if
    end do
  end do

return
end subroutine initialise_Move


end module initialise_Module
