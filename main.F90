program main
  use mpi_variables
  use GP_Variables
  use tmpi_calcFullSimBoxEnergy_mod
  use energiesData_Module, only: energiesData
  use assert_module
  implicit none


  integer :: N_a, N_tri, udSize
  double precision, allocatable :: posArray(:,:)!, X_dg(:,:)
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  type( energiesData ) :: currentEnergies, proposedEnergies

  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls
  call initialise_GP(hyperParametersFile, alphaFile, trainingSetFile)
  !  call initialise_Positions('AtomicPositions400.txt', posArray,N_a)
  call initialise_Positions('AtomicPositions5.txt', posArray,N_a)
  call initialise_Variables(N_a, N_tri,udSize)
  

  currentEnergies = tmpi_calcFullSimBoxEnergy(N_a,N_tri,udSize,posArray)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call tmpi_calcAtomMoveEnergy(20,1.5d0,N_a,udSize,N_tri,currentEnergies, &
                               posArray,proposedEnergies)


  deallocate(alpha)
  deallocate(trainData)


  call MPI_FINALIZE(ierror)
  

end program main

