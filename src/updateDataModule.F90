module updateData
  use mpi_variables
  use expShare_variables
  use dataStructure_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  implicit none
  include 'mpif.h'


private
public updateCurrentDataStructures, resetProposedDataStructures


contains


  subroutine updateCurrentDataStructures(move)
    implicit none
    integer :: move
    double precision :: energyTime, xdgTime, tripletTime
    double precision :: expTime, posTime

    energyTime = MPI_Wtime()
    currentEnergyData%Utotal = currentEnergyData%Utotal + &
                               proposedEnergyData%Utotal
    energyTime = MPI_Wtime() - energyTime
    xdgTime = MPI_Wtime()
    call updateDistanceData(currentEnergyData,proposedEnergyData,expUpdateIndNoRepeat, &
                            currentPositionData%N_a-1)
    xdgTime = MPI_Wtime() - xdgTime
    tripletTime = MPI_Wtime()
    call updateTripletEnergies(currentEnergyData,proposedEnergyData)
    tripletTime = MPI_Wtime() - tripletTime
    expTime = MPI_Wtime()
    !call updateExpMatrix(changeExpData,expUpdateIndNoRepeat,currentPositionData%N_a-1)
    expTime = MPI_Wtime() - expTime
    posTime = MPI_Wtime()
    currentPositionData%posArray(move,:) = proposedPositionData%posArray(move,:)
    posTime = MPI_Wtime() - posTime

    !if (processRank .eq. root) then
    !  print *, energyTime, xdgTime, tripletTime, expTime, posTime, 0d0, 0d0, 0d0
    !end if

  return
  end subroutine updateCurrentDataStructures


  subroutine resetProposedDataStructures(move)
    implicit none
    integer :: move
    double precision :: xdgTime, tripletTime
    double precision :: expTime, posTime

    xdgTime = MPI_Wtime()
    call updateDistanceData(proposedEnergyData,currentEnergyData,expUpdateIndNoRepeat, &
                            proposedPositionData%N_a-1)
    xdgTime = MPI_Wtime() - xdgTime
    tripletTime = MPI_Wtime()
    call updateTripletEnergies(proposedEnergyData,currentEnergyData)
    tripletTime = MPI_Wtime() - tripletTime
    expTime = MPI_Wtime()
    call resetExpMatrix(proposedEnergyData,expUpdateIndNoRepeat,proposedPositionData%N_a-1)
    expTime = MPI_Wtime() - expTime
    posTime = MPI_Wtime()
    proposedPositionData%posArray(move,:) = currentPositionData%posArray(move,:)
    posTime = MPI_Wtime() - posTime

    !if (processRank .eq. root) then
    !  print *, 0d0, xdgTime, tripletTime, expTime, posTime, 0d0, 0d0, 0d0
    !end if

  return
  end subroutine resetProposedDataStructures


  subroutine resetExpMatrix(proposedEnergyData,expInd,length)
    integer :: indj, length, expInd(length,2), j
    type (energiesData) :: proposedEnergyData

    do j = 1, length
      indj = proposedEnergyData%distancesIntMat(expInd(j,1), expInd(j,2))
      expArray(1:N_tp,1:nArgs,indj) = oldExpData(1:N_tp,1:nArgs,j)
    end do

  return
  end subroutine resetExpMatrix

  subroutine updateDistanceData(changedEnData,unchangedEnData,expInd,length)
    integer :: length, expInd(length,2), j
    type (energiesData) :: changedEnData, unchangedEnData

    do j = 1, length
      changedEnData%interatomicDistances(expInd(j,1),expInd(j,2)) = &
      unchangedEnData%interatomicDistances(expInd(j,1),expInd(j,2))
    end do

  return
  end subroutine updateDistanceData

  subroutine updateTripletEnergies(changedEnData,unchangedEnData)
    integer :: j, nUpdates, enUpdate
    type (energiesData) :: changedEnData, unchangedEnData

    nUpdates = size(changedTriInd)

    do j = 1, nUpdates
      enUpdate = changedTriInd(j)
      changedEnData%tripletEnergies(enUpdate) = &
      unchangedEnData%tripletEnergies(enUpdate)
    end do

  return
  end subroutine updateTripletEnergies

end module updateData
