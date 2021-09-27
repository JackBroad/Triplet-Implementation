module updateData
  use mpi_variables
  use expShare_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none
  include 'mpif.h'

  private
  public updateDataAfterMove

contains

  subroutine broadcastEnergyData(proposedEnergyData,proposedPosition)
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergyData

    call MPI_Bcast(proposedEnergyData%tripletEnergies, proposedPosition%N_tri, &
                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(proposedEnergyData%Utotal, 1, MPI_DOUBLE_PRECISION,  root, &
                   MPI_COMM_WORLD, ierror)

  end subroutine broadcastEnergyData

  subroutine shareChangedExponentials(proposedEnergyData)
    implicit none
    integer :: length, lengthVec(clusterSize), sumLength
    integer :: maxLength, reLength, j, disp(clusterSize)
    type (energiesData) :: proposedEnergyData
    integer, allocatable :: changeExpInd(:,:)
    integer, allocatable :: expUpdateNoRepeatTrans(:,:)
    double precision, allocatable :: changeExpMat(:,:,:)

    ! Only need to do anything if >1 process present
    if (clusterSize .gt. 1) then

      ! Find number of changed exps across all processors
      length = size(expUpdateNoRepeat)
      call MPI_gather(length, 1, MPI_INT, lengthVec, 1, MPI_INT, root, &
                      MPI_COMM_WORLD, ierror)
      call MPI_Bcast(lengthVec, clusterSize, MPI_INT, root, MPI_COMM_WORLD, &
                     ierror)
      call MPI_BARRIER(MPI_COMM_WORLD, barError)

      ! Find displacement of each process when gathering exps
      disp(1) = 0
      do j = 2, clusterSize
        disp(j) = sum(lengthVec(1:j-1))
      end do
      sumLength = sum(lengthVec)

      ! Gather changed exponentials to root
      allocate(changeExpMat(nArgs,N_tp,sumLength))
      call MPI_gatherv(changeExpData, N_tp*nArgs*length, MPI_DOUBLE_PRECISION, &
                       changeExpMat, N_tp*nArgs*lengthVec, N_tp*nArgs*disp, &
                       MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

      ! Gather indices if changed exponentials to root
      allocate(changeExpInd(2,sumLength))
      allocate(expUpdateNoRepeatTrans(2,length))
      expUpdateNoRepeatTrans = transpose(expUpdateIndNoRepeat)
      call MPI_gatherv(expUpdateNoRepeatTrans, 2*length, MPI_INT, changeExpInd, &
                       2*lengthVec, 2*disp, MPI_INT, 1, MPI_COMM_WORLD, &
                       ierror)

      ! Broadcast updated exps and indices
      !=======Would trimming the repeats prior to Bcasting be good here?=======
      call MPI_Bcast(changeExpMat, N_tp*nArgs*sumLength, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(changeExpInd, 2*sumLength, MPI_INT, 1, MPI_COMM_WORLD, &
                     ierror)

      ! Update exp matrix on all processes
      changeExpInd = transpose(changeExpInd)
      call updateExpMatrix(proposedEnergyData,changeExpMat,changeExpInd, &
                           sumLength)
      deallocate(changeExpMat,changeExpInd)

    end if

  return
  end subroutine shareChangedExponentials


  subroutine updateDataAfterMove(proposedEnergyData,proposedPositionData, &
                                 currentEnergyData,currentPositionData)
    implicit none
    type (energiesData) :: proposedEnergyData, currentEnergyData
    type (positionData) :: proposedPositionData, currentPositionData

    call shareChangedExponentials(proposedEnergyData)
    call broadcastEnergyData(proposedEnergyData,proposedPositionData)

    currentEnergyData = proposedEnergyData
    currentPositionData = proposedPositionData

  return
  end subroutine updateDataAfterMove
  
end module updateData
