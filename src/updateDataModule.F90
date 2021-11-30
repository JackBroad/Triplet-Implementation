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
  public updateDataAfterMove

contains


  subroutine updateTripletEnergies()
    implicit none
    integer :: allInds(proposedPositionData%N_changed_triplets), length
    integer :: i, disp(clusterSize), lengthVec(clusterSize)
    double precision :: allEnergies(proposedPositionData%N_changed_triplets)

    ! Set up arrays and constants needed to gather energies
    length = size(proposedEnergyData%changedTriU)
    call MPI_gather(length, 1, MPI_INT, lengthVec, 1, MPI_INT, root, &
                    MPI_COMM_WORLD, ierror)
    call MPI_Bcast(lengthVec, clusterSize, MPI_INT, root, MPI_COMM_WORLD, &
                   ierror)
    disp(1) = 0
    do i = 2, clusterSize
       disp(i) = sum(lengthVec(1:i-1))
    end do

    ! Gather all new trip energies and their indices on root
    call MPI_gatherv(proposedEnergyData%changedTriInd, length, MPI_INT, &
                     allInds, lengthVec, disp, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_gatherv(proposedEnergyData%changedTriU, length, MPI_DOUBLE_PRECISION, &
                     allEnergies, lengthVec, disp, MPI_DOUBLE_PRECISION, root, &
                     MPI_COMM_WORLD, ierror)

    ! Broadcast the complete vectors required for the update (do this rather
    ! than sum on root and broadcast whole tripletEnergy vector as this is
    ! faster for > 5 atoms)
    call MPI_Bcast(allInds, proposedPositionData%N_changed_triplets, MPI_INT, root, &
                   MPI_COMM_WORLD, ierror)
    call MPI_Bcast(allEnergies, proposedPositionData%N_changed_triplets, MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierror)

    ! Update triplet energies on each process
    do i = 1, proposedPositionData%N_changed_triplets
      proposedEnergyData%tripletEnergies(allInds(i)) = allEnergies(i)
    end do

  return
  end subroutine updateTripletEnergies


  subroutine updateCurrentDataStructures()
    implicit none

    currentEnergyData%Utotal = currentEnergyData%Utotal + &
                               proposedEnergyData%Utotal
    currentEnergyData%interatomicDistances = proposedEnergyData%&
                                             interatomicDistances
    currentEnergyData%tripletEnergies = proposedEnergyData%&
                                        tripletEnergies
    currentEnergyData%expMatrix = proposedEnergyData%expMatrix
    currentPositionData = proposedPositionData

  return
  end subroutine


  subroutine updateDataAfterMove()
    implicit none

    call updateTripletEnergies()
    call updateCurrentDataStructures()

  return
  end subroutine updateDataAfterMove
  
end module updateData
