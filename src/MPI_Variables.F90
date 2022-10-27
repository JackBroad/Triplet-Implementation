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


module pbcAndMic_variables
  implicit none

  integer :: nExplicit
  logical :: calculate
  double precision :: sideLength, Rcut
  double precision :: density
  integer, parameter :: dp = SELECTED_REAL_KIND(15)
  double precision, allocatable :: typeTwoTable(:,:)
  double precision, allocatable :: typeThreeTable(:,:)

end module pbcAndMic_variables


module expShare_variables
  implicit none

  integer :: N_changed_exp_per_host, N_exp_per_host
  integer, allocatable :: expUpdateInd(:,:),expUpdateIndNoRepeat(:,:)
  integer, allocatable :: changedTriInd(:), hostIndices(:,:)
  integer, allocatable :: fullHostInds(:,:)
  double precision, allocatable :: expUpdate(:), expUpdateNoRepeat(:)
  double precision, allocatable :: changeExpData(:,:,:), hostDists(:)
  double precision, allocatable :: oldExpData(:,:,:), fullHostDists(:)

end module expShare_variables


module dataStructure_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none

  type (energiesData) :: currentEnergyData, proposedEnergyData
  type (positionData) :: currentPositionData, proposedPositionData

end module dataStructure_variables


module time_variables
  implicit none

  double precision :: moveTime, expTime, sumTime, setTime
  double precision :: gatherTime, xTime, tripTime, partialSumTime
  double precision :: extractTime, tripSumTime, rootSumTime

end module time_variables
