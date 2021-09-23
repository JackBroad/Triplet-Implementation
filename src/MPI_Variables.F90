module mpi_variables
  implicit none

  integer :: processRank, clusterSize, ierror
  integer :: root=0, barError

end module mpi_variables


module expShare_variables
  implicit none

  integer, allocatable :: expUpdateInd(:,:),expUpdateIndNoRepeat(:,:)
  double precision, allocatable :: expUpdate(:), expUpdateNoRepeat(:)
  double precision, allocatable :: changeExpData(:,:,:)

end module expShare_variables
