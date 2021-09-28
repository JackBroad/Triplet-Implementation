program main
  use mpi_variables
  use expShare_variables
  use GP_Variables
  use triplet_mod
  use tmpi_calcFullSimBoxEnergy_mod, only: tmpi_calcFullSimBoxEnergy
  use tmpi_calcAtomMoveEnergy_mod, only: tmpi_calcAtomMoveEnergy
  use toyMove_Module, only: toyMove
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use updateData
  use assert_module
  use initialise_Module
  implicit none
  include 'mpif.h'


  integer :: move, i
  logical :: seed=.false., acceptMove
  double precision :: dist
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  Character(len=300) :: positionFile = 'AtomicPositions400.txt'
  type (energiesData) :: currentEnergies, proposedEnergies
  type (positionData) :: currentPosition, proposedPosition

  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls for full calc
  call initialise_GP_NonAdd(hyperParametersFile, alphaFile, trainingSetFile)
  currentPosition = initialise_positionDataStruct(positionFile)  


  ! Calculate energy for full sim box
  currentEnergies = tmpi_calcFullSimBoxEnergy(currentPosition)
  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  ! Atom move
  dist = 1.5d0
  acceptMove = .true.
  do i = 1, 150
    call initialise_Move(currentPosition,currentEnergies,dist,seed, &
                         proposedPosition,proposedEnergies,move)
    call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(proposedPosition%posArray(move,:), 3, MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierror)

    !proposedEnergies = tmpi_calcAtomMoveEnergy(move,proposedPosition,currentEnergies)
    proposedEnergies = toyMove(move,proposedPosition,currentEnergies)

    if (acceptMove .eqv. .true.) then
      call updateDataAfterMove(proposedEnergies,proposedPosition, &
                               currentEnergies,currentPosition)
    end if
  end do

  deallocate(alpha)
  deallocate(trainData)
  deallocate(expUpdate)
  deallocate(expUpdateInd)
  deallocate(expUpdateNoRepeat)
  deallocate(expUpdateIndNoRepeat)
  deallocate(changeExpData)


  call MPI_FINALIZE(ierror)
  

end program main

