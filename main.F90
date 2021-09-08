program main
  use mpi_variables
  use GP_Variables
  use triplet_mod
  use tmpi_calcFullSimBoxEnergy_mod, only: tmpi_calcFullSimBoxEnergy
  use tmpi_calcAtomMoveEnergy_mod, only: tmpi_calcAtomMoveEnergy
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use assert_module
  use initialise_Module
  implicit none
  include 'mpif.h'


  integer :: move, i
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
  do i = 1, 3
  call initialise_Move(currentPosition,dist, proposedPosition,move)
  call MPI_Bcast(proposedPosition%posArray, 3*proposedPosition%N_a, &
                 MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)


  proposedEnergies = tmpi_calcAtomMoveEnergy(1,move,proposedPosition,currentEnergies)
  end do

  deallocate(alpha)
  deallocate(trainData)


  call MPI_FINALIZE(ierror)
  

end program main

