module initialise_Module
  use GP_variables
  use positionData_Module, only: positionData
  use mpi_variables
  implicit none
  !include 'mpif.h'


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

  Perm = setupPermutationMatrix()

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
subroutine initialise_Move(currentPos,dMax,addSeed,  newPos,mover)
  implicit none
  logical, intent(in) :: addSeed
  double precision, intent(in) :: dMax
  type (positionData), intent(in) :: currentPos
  integer, intent(out) :: mover
  type (positionData), intent(out) :: newPos
  integer :: icol, irow, seed(8) ! Min. size for seed array
  double precision :: randNo

  newPos = currentPos
  if (processRank .eq. root) then
    if (addSeed .eqv. .true.) then
      seed = 168389234
      call random_seed(put=seed)
    end if
  call random_number(randNo)
  mover = 1 + FLOOR(newPos%N_a*randNo)
!  print *, 'moving atom', mover

  ! Change the position of the atom in newPos
  do irow = 1, newPos%N_a
    if (irow .eq. mover) then
      do icol = 1, 3
        call random_number(randNo)
        newPos%posArray(irow,icol) = newPos%posArray(irow,icol) + (2.0*randNo-1) &
                                     * dMax
        !print *, 'moved', (2.0*randNo-1)*dMax, 'in direction', icol
      end do
    end if
  end do
  end if
  
return
end subroutine initialise_Move


function setupPermutationMatrix() result(Perm)
  integer :: Perm(6,3)

  Perm(1,:) = (/1, 2, 3/)
  Perm(2,:) = (/1, 3, 2/)
  Perm(3,:) = (/2, 1, 3/)
  Perm(4,:) = (/2, 3, 1/)
  Perm(5,:) = (/3, 1, 2/)
  Perm(6,:) = (/3, 2, 1/)

end function setupPermutationMatrix


function initialise_positionDataStruct(fileName) result(currentPosition)
  type (positionData) :: currentPosition
  character(len=40) :: fileName
  call initialise_Positions(fileName, currentPosition%posArray,currentPosition%N_a)
  call initialise_Variables(currentPosition%N_a, currentPosition%N_tri, &
                            currentPosition%N_distances)
end function


end module initialise_Module
