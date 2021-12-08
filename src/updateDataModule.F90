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

    call updateCurrentDataStructures()

  return
  end subroutine updateDataAfterMove
  
end module updateData
