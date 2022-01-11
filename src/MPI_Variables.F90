module mpi_variables
  implicit none

  integer :: processRank, clusterSize, ierror
  integer :: root=0, barError

end module mpi_variables


module expShare_variables
  implicit none

  integer, allocatable :: expUpdateInd(:,:),expUpdateIndNoRepeat(:,:)
  integer, allocatable :: changedTriInd(:)
  double precision, allocatable :: expUpdate(:), expUpdateNoRepeat(:)
  double precision, allocatable :: changeExpData(:,:,:), expArray(:,:,:)
  double precision, allocatable :: oldExpData(:,:,:)

end module expShare_variables

module dataStructure_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none

  type (energiesData) :: currentEnergyData, proposedEnergyData
  type (positionData) :: currentPositionData, proposedPositionData

end module dataStructure_variables
