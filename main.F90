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


  integer :: move
  double precision :: dist
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  Character(len=300) :: positionFile = 'AtomicPositions5.txt'
  type (energiesData) :: currentEnergies, proposedEnergies
  type (positionData) :: currentPosition, proposedPosition

  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls for full calc
  call initialise_GP_NonAdd(hyperParametersFile, alphaFile, trainingSetFile)
  call initialise_Positions(positionFile, currentPosition%posArray, &
                            currentPosition%N_a)
  call initialise_Variables(currentPosition%N_a, currentPosition%N_tri, &
                            currentPosition%N_distances)
  

  currentEnergies = tmpi_calcFullSimBoxEnergy(currentPosition)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  ! Set-up calls for atom move
  proposedPosition = currentPosition
  dist = 1.5d0
  call initialise_Move(currentPosition%posArray,currentPosition%N_a,dist, &
                       proposedPosition%posArray,move)


  call tmpi_calcAtomMoveEnergy(20,move,proposedPosition,currentEnergies, &
                               proposedEnergies)


  deallocate(alpha)
  deallocate(trainData)


  call MPI_FINALIZE(ierror)
  

end program main

