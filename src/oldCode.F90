module oldCode

contains

  subroutine updateChangedTripletEnergies()

    proposedEnergyData%tripletEnergies = currentEnergyData%tripletEnergies
    do j = 1, proposedPositionData%N_changed_triplets

      proposedEnergyData%tripletEnergies(tripIndex(j)) = newUfull(j)

    end do

  end subroutine updateChangedTripletEnergies

  subroutine getAffectedTripletDistances(move,proposedEnergy)
    implicit none
    integer, intent(in) :: move
    type (energiesData) :: proposedEnergy

    counter = 0
    do j = 1, triPerProc
      if (scatterTrip(1,j) .eq. move) then

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(move,scatterTrip(2,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(2,j)
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(move,scatterTrip(3,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(3,j)

      else if (scatterTrip(2,j) .eq. move) then

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(scatterTrip(1,j),move)
        expUpdateInd(counter,1) = scatterTrip(1,j)
        expUpdateInd(counter,2) = move
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(move,scatterTrip(3,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(3,j)

      else

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(scatterTrip(1,j),move)
        expUpdateInd(counter,1) = scatterTrip(1,j)
        expUpdateInd(counter,2) = move
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergy%interatomicDistances(scatterTrip(2,j),move)
        expUpdateInd(counter,1) = scatterTrip(2,j)
        expUpdateInd(counter,2) = move

      end if
    end do

    ! Remove repeat distances
    allocate(mask(2*triPerProc))
    mask = .TRUE.
    do j = 2*triPerProc,2,-1
      mask(j) = .NOT.(ANY(expUpdate(:j-1)==expUpdate(j)))
    end do
    allocate(indexVector, source=PACK([(j,j=1,2*triPerProc)], mask))
    if (allocated(expUpdateNoRepeat)) then
      deallocate(expUpdateNoRepeat)
    end if
    allocate(expUpdateNoRepeat, source=expUpdate(indexVector))
    deallocate(indexVector)

    ! Do the same for the distance indices
    mask = .TRUE.
    do j = 2*triPerProc,2,-1
      mask(j) = .NOT.(ANY(expUpdateInd(:j-1,1)==expUpdateInd(j,1) .AND. &
                expUpdateInd(:j-1,2)==expUpdateInd(j,2)))
    end do
    allocate(indexVector, source=PACK([(j,j=1,2*triPerProc)], mask))
    if (allocated(expUpdateIndNoRepeat)) then
      deallocate(expUpdateIndNoRepeat)
    end if
    allocate(expUpdateIndNoRepeat, source=expUpdateInd(indexVector,:))
    deallocate(mask,indexVector)

  return
  end subroutine getAffectedTripletDistances

end module oldCode
