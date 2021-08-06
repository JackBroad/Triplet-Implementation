program main
  use mpi_variables
  use GP_Variables
  use triplet_mod
  use tmpi_calcFullSimBoxEnergy_mod, only: tmpi_calcFullSimBoxEnergy
  use tmpi_calcAtomMoveEnergy_mod, only: tmpi_calcAtomMoveEnergy
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use assert_module
  implicit none
  include 'mpif.h'


  !integer :: N_a, N_tri, udSize
  !double precision, allocatable :: posArray(:,:)
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  type (energiesData) :: currentEnergies, proposedEnergies
  type (positionData) :: currentPosition, proposedPosition

  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls
  call initialise_GP(hyperParametersFile, alphaFile, trainingSetFile)
  call initialise_Positions('AtomicPositions5.txt', currentPosition%posArray, &
                             currentPosition%N_a)
  call initialise_Variables(currentPosition%N_a, currentPosition%N_tri, &
                            currentPosition%N_distances)
  

  currentEnergies = tmpi_calcFullSimBoxEnergy(currentPosition%N_a, &
                                              currentPosition%N_tri, &
                                              currentPosition%N_distances, &
                                              currentPosition%posArray)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call tmpi_calcAtomMoveEnergy(20,1.5d0,currentPosition%N_a,currentPosition%N_distances, &
                               currentPosition%N_tri,currentEnergies, &
                               currentPosition%posArray,proposedEnergies)


  deallocate(alpha)
  deallocate(trainData)


  call MPI_FINALIZE(ierror)
  

end program main

