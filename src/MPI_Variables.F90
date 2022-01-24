module mpi_variables
  use, intrinsic :: ISO_C_BINDING, only : C_PTR
  use mpi
  implicit none

  ! MPI variables
  integer :: processRank, clusterSize, ierror
  integer :: root=0, barError

  ! Shared mem variables
  integer :: hostComm, sharedSize, hostRank
  integer :: win, disp_unit, shapeArray(3)
  type(C_PTR) :: baseptr
  integer, pointer :: dummy(:,:,:)
  double precision, pointer :: expArray(:,:,:)
  integer(KIND=MPI_ADDRESS_KIND) :: windowsize

end module mpi_variables


module expShare_variables
  implicit none

  integer, allocatable :: expUpdateInd(:,:),expUpdateIndNoRepeat(:,:)
  integer, allocatable :: changedTriInd(:)
  double precision, allocatable :: expUpdate(:), expUpdateNoRepeat(:)
  double precision, allocatable :: changeExpData(:,:,:)
  double precision, allocatable :: oldExpData(:,:,:)

end module expShare_variables

module dataStructure_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none

  type (energiesData) :: currentEnergyData, proposedEnergyData
  type (positionData) :: currentPositionData, proposedPositionData

end module dataStructure_variables
