module atomMoveModule
  use GP_variables
  use mpi_variables
  use expShare_variables
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
    integer :: tmpChangedTri(nTrips), sizeCounter
    double precision :: U(1)
    integer, intent(in) :: move, nDists, nTrips
    double precision, intent(inout) :: uVec(nTrips)

    counter = 0
    sizeCounter = 0
    do i = 1, nDists
      al = proposedEnergyData%alphaBetaPairs(i,1)
      be = proposedEnergyData%alphaBetaPairs(i,2)
      if (be .lt. move) then
        sizeCounter = sizeCounter + 1
        counter = counter + (move - be)
        triplet = (/ al, be, move /)
        call tripletEnergiesNonAdd(triplet,proposedEnergyData%distancesIntMat, &
                                   1,N_tp,proposedPositionData%N_a,N_p,nArgs, &
                                   Perm,proposedPositionData%N_distances, &
                                   expArray,alpha,hyperParams(2), U)
        uVec(counter) = U(1)
        tmpChangedTri(sizeCounter) = counter
        counter = counter + (proposedPositionData%N_a - move)
      else if ( (al .eq. move) .or. (be .eq. move) ) then
        do j = be+1, proposedPositionData%N_a
          counter = counter + 1
          sizeCounter = sizeCounter + 1
          triplet = (/ al, be, j /)
          call tripletEnergiesNonAdd(triplet,proposedEnergyData%distancesIntMat, &
                                     1,N_tp,proposedPositionData%N_a,N_p,nArgs, &
                                     Perm,proposedPositionData%N_distances, &
                                     expArray,alpha,hyperParams(2), U)
          uVec(counter) = U(1)
          tmpChangedTri(sizeCounter) = counter
        end do
      else
        counter = counter + (proposedPositionData%N_a - be) ! Skip all trips
                                                            ! involving current 
                                                            ! alpha and beta
      end if
    end do

    if (allocated(changedTriInd)) then
      deallocate(changedTriInd)
    end if
    allocate(changedTriInd(sizeCounter))
    changedTriInd = tmpChangedTri(1:sizeCounter)

  return
  end subroutine getTripletEnergiesAtomMove


  function getTriPerAtom(nAt) result(nPer)
    implicit none
    integer, intent(in) :: nAt
    integer :: nPer
    integer :: a, b

    ! Determine the number of triplets each atom is involved in
    nPer = 0
    do a = 2, nAt-1
      do b = a+1, nAt
        nPer = nPer + 1
      end do
    end do

  return
  end function getTriPerAtom


  ! Uses the number of atoms (num) and index of moved atom (atom) to determine
  ! which exponentials need re-calculating after a move and return a matrix of
  ! their postions in X_dg
  subroutine extractChangedExps(num,atInd,Xdg, change,dists)
    implicit none
    integer, intent(in) :: num, atInd
    double precision, intent(in) :: Xdg(num,num)
    integer, intent(out) :: change(num-1,2)
    double precision, intent(out) :: dists(num-1)
    integer :: i, j

    do i = 1, num
      if (i .lt. atInd) then

        change(i,1) = i
        change(i,2) = atInd

      else if (i .gt. atInd) then

        change(i-1,1) = atInd
        change(i-1,2) = i

      end if
    end do

    do j = 1, num-1

      dists(j) = Xdg(change(j,1),change(j,2))

    end do

  return
  end subroutine extractChangedExps


  subroutine updateXdg(move,N_a,positions, X)
    implicit none
    ! Arguments
    integer, intent(in) :: move, N_a
    double precision, intent(in) :: positions(N_a,3)
    double precision, intent(inout) :: X(N_a,N_a)
    ! Local variables
    integer :: i
    double precision :: changedPosition(3)

    ! Identify the position of the moved atom
    changedPosition = positions(move,:)

    ! Find the distances between this atom and all others
    do i = 1, N_a
      if (i .lt. move) then

        X(i,move) = (positions(i,1) - changedPosition(1))**2 + &
                    (positions(i,2) - changedPosition(2))**2 + &
                    (positions(i,3) - changedPosition(3))**2
        X(i,move) = (X(i,move))**0.5
        X(i,move) = 1 / X(i,move)
        X(move,i) = X(i,move)

      else if (i .gt. move) then

        X(move,i) = (positions(i,1) - changedPosition(1))**2 + &
                    (positions(i,2) - changedPosition(2))**2 + &
                    (positions(i,3) - changedPosition(3))**2
        X(move,i) = (X(move,i))**0.5
        X(move,i) = 1 / X(move,i)
        X(i,move) = X(move,i)

      end if
    end do

    return
  end subroutine updateXdg


  subroutine getChangedTriplets(atom, changedTriplets,tripIndex)
    implicit none
    integer, intent(in) :: atom
    integer, intent(out) :: changedTriplets(3,proposedPositionData%N_changed_triplets)
    integer, intent(out) :: tripIndex(proposedPositionData%N_changed_triplets)
    integer :: al, be, ga, counter, indCounter

    ! Fill changed triplet array and vector of indices of changed triplets
    counter = 0
    indCounter = 0
    do while (counter .lt. proposedPositionData%N_changed_triplets)
    do al = 1, proposedPositionData%N_a-2
      do be = al+1, proposedPositionData%N_a-1
        do ga = be+1, proposedPositionData%N_a
          indCounter = indCounter + 1
          if (al .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga
            tripIndex(counter) = indCounter

          else if (be .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga
            tripIndex(counter) = indCounter

          else if (ga .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga
            tripIndex(counter) = indCounter

          end if
        end do
      end do
    end do
    end do

  return
  end subroutine getChangedTriplets


  subroutine saveOldExponentials()
    implicit none
    integer :: i, j

    do i = 1,proposedPositionData%N_a-1
      j = proposedEnergyData%distancesIntMat(expUpdateIndNoRepeat(i,1), &
                                             expUpdateIndNoRepeat(i,2))
      oldExpData(1:N_tp,1:nArgs,i) = expArray(1:N_tp,1:nArgs,j)
    end do

  return
  end subroutine saveOldExponentials


end module atomMoveModule
