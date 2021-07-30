program main
  !use mpi
  use mpi_variables
  use GP_Variables
  use triplet_mpi_mod
  implicit none

  integer :: N_a, N_tri, udSize
  double precision :: uFinal, uChange
  integer, allocatable :: disIntMat(:,:)
  double precision, allocatable :: posArray(:,:), X_dg(:,:)
  double precision, allocatable :: expMat(:,:,:), uVecFinal(:)
  double precision, allocatable :: uVecChange(:)
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set-up calls
  call initialise_GP()
  call initialise_Positions('AtomicPositions5.txt', posArray,N_a)
  call initialise_Variables(N_a, N_tri,udSize)
  

  call triplet_mpi_fullNonAdd(N_a,N_tri,udSize,posArray, X_dg, &
                              disIntMat,expMat,uFinal,uVecFinal)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call triplet_mpi_moveNonAdd(20,1.5d0,N_a,N_tri,udSize,posArray,X_dg, &
                              disIntMat,expMat,uChange,uVecChange)


  deallocate(alpha)
  deallocate(trainData)
  deallocate(X_dg)
  deallocate(expMat)
  deallocate(disIntMat)
  deallocate(posArray)
  !deallocate(uVecFinal)
  !deallocate(uVecChange)


  call MPI_FINALIZE(ierror)
  

end program main

