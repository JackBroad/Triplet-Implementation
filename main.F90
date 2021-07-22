program main
  !use mpi
  use mpi_variables
  use triplet_mpi_mod
  implicit none

  integer :: N_at, nArg, N_t, N_trip, nUD, Permu(6,3), N_per
  double precision :: hyperParam(3), uFinal, uChange
  integer, allocatable :: disIndMat(:,:)
  double precision, allocatable :: posiArray(:,:), Xdg(:,:), trainDat(:,:)
  double precision, allocatable :: alph(:), expMat(:,:,:), uVecFinal(:)
  double precision, allocatable :: uVecChange(:)
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)
  

  call triplet_mpi_fullNonAdd('AtomicPositions5.txt', N_at,nArg,N_trip, &
                              N_per,nUD,posiArray,Xdg,disIndMat,expMat,uFinal,uVecFinal)


  call MPI_BARRIER(MPI_COMM_WORLD, barError)


  call triplet_mpi_moveNonAdd(5,1.5d0,N_at,nArg,N_trip,N_per,nUD, &
                              posiArray,Xdg,disIndMat,expMat,uChange,uVecChange)


  deallocate(alpha)
  deallocate(trainData)
  deallocate(Xdg)
  deallocate(expMat)
  deallocate(disIndMat)
  deallocate(posiArray)


  call MPI_FINALIZE(ierror)
  

end program main

