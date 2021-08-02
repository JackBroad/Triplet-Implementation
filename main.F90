program main
  !use mpi
  use mpi_variables
  use GP_Variables
  use triplet_mpi_mod
  use energiesData_Module, only: energiesData
  use assert_module
  implicit none

  integer :: N_a, N_tri, udSize, i
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
  call initialise_Positions('AtomicPositions5.txt', posArray,N_a)
  call initialise_Variables(N_a, N_tri,udSize)
  

  currentEnergies = triplet_mpi_fullNonAdd(N_a,N_tri,udSize,posArray)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call triplet_mpi_moveNonAdd(1,1.5d0,N_a,N_tri,udSize,currentEnergies,posArray, &
                              proposedEnergies)


  if (processRank .eq. 0) then
  print *, currentEnergies%Utotal
  print *, proposedEnergies%Utotal
  print *, ' '
  do i = 1, N_a
    print *, currentEnergies%interatomicDistances(i,:)
  end do
  print *, ' '
  do i = 1, N_a
    print *, proposedEnergies%interatomicDistances(i,:)
  end do
  !do i = 1, udSize
  !  print *, currentEnergies%expMatrix(:,:,i)
  !  print *, proposedEnergies%expMatrix(:,:,i)
  !end do
  end if


  deallocate(alpha)
  deallocate(trainData)

  call MPI_FINALIZE(ierror)
  

end program main

