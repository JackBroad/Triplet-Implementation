program tryMPI
  !use mpi
  use mpi_variables
  use triplet_mpi_mod
  implicit none
  
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)
  

  call triplet_mpi(fileName='AtomicPositions5.txt')


  call MPI_FINALIZE(ierror)
  
END PROGRAM tryMPI

