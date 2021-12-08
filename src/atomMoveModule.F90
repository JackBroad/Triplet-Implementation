module atomMoveModule
  use GP_variables
  use mpi_variables
  use dataStructure_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use triplet_mod
  implicit none
  !include 'mpif.h'


contains


  subroutine getTripletEnergiesAtomMove(move,nDists,nTrips, uVec)
    implicit none
    integer :: i, j, counter, al, be, triplet(3)
    double precision :: U(1)
    integer, intent(in) :: move, nDists, nTrips
    double precision, intent(inout) :: uVec(nTrips)

    counter = 0
    do i = 1, nDists
      al = proposedEnergyData%alphaBetaPairs(i,1)
      be = proposedEnergyData%alphaBetaPairs(i,2)
      if (be .lt. move) then
        counter = counter + (move - be)
        triplet = (/ al, be, move /)
        call tripletEnergiesNonAdd(triplet,proposedEnergyData%distancesIntMat, &
                                   1,N_tp,proposedPositionData%N_a,N_p,nArgs, &
                                   Perm,proposedPositionData%N_distances, &
                                   proposedEnergyData%expMatrix,alpha, &
                                   hyperParams(2), U)
        uVec(counter) = U(1)
        counter = counter + (proposedPositionData%N_a - move)
      else if ( (al .eq. move) .or. (be .eq. move) ) then
        do j = be+1, proposedPositionData%N_a
          counter = counter + 1
          triplet = (/ al, be, j /)
          call tripletEnergiesNonAdd(triplet,proposedEnergyData%distancesIntMat, &
                                     1,N_tp,proposedPositionData%N_a,N_p,nArgs, &
                                     Perm,proposedPositionData%N_distances, &
                                     proposedEnergyData%expMatrix,alpha, &
                                     hyperParams(2), U)
          uVec(counter) = U(1)
        end do
      else
        counter = counter + (proposedPositionData%N_a - be) ! Skip all trips
                                                            ! involving current 
                                                            ! alpha and beta
      end if
    end do

  return
  end subroutine getTripletEnergiesAtomMove


end module atomMoveModule
