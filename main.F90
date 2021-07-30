program main
  !use mpi
  use mpi_variables
  use GP_Variables
  use triplet_mpi_mod
  use energiesData_Module, only: energiesData
  implicit none

  integer :: N_a, N_tri, udSize
  double precision, allocatable :: posArray(:,:)!, X_dg(:,:)
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  type( energiesData ) :: currentEnergies
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls
  call initialise_GP(hyperParametersFile, alphaFile, trainingSetFile)
  call initialise_Positions('AtomicPositions5.txt', posArray,N_a)
  call initialise_Variables(N_a, N_tri,udSize)
  

!  call triplet_mpi_fullNonAdd(N_a,N_tri,udSize,posArray, X_dg, &
!       disIntMat,expMat,uFinal,uVecFinal)

  call triplet_mpi_fullNonAdd(N_a,N_tri,udSize,posArray, currentEnergies)!%interatomicDistances, &
!       currentEnergies%distancesIntMat,currentEnergies%expMatrix, &
 !      currentEnergies%Utotal, currentEnergies%tripletEnergies)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call triplet_mpi_moveNonAdd(20,1.5d0,N_a,N_tri,udSize,posArray,  &
       currentEnergies%interatomicDistances, &
       currentEnergies%distancesIntMat,currentEnergies%expMatrix, &
       currentEnergies%Utotal, currentEnergies%tripletEnergies)


  deallocate(alpha)
  deallocate(trainData)

  call MPI_FINALIZE(ierror)
  

end program main

