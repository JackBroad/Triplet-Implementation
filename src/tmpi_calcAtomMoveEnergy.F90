module tmpi_calcAtomMoveEnergy_mod
  use mpi_variables
  use expShare_variables
  use dataStructure_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use global_Flags, only: textOutput
  use assert_module
  implicit none
  include 'mpif.h'

  private
  public tmpi_calcAtomMoveEnergy, extractChangedExps, getTriPerAtom, &
         updateXdg, getNPerProcNonAdd, getVarrays, getTripletScatterData, &
         getChangedTriplets


  integer :: N_changed_triplets, nTriMax, nTriRe,triPerProc
  integer :: j, counter
  double precision :: moveTime, expTime, sumTime, setTime
  double precision :: gatherAndSumTime, xTime, tripTime
  double precision :: extractTime, tripSumTime, newTripU
  double precision :: oldTripU, partialDeltaU
  integer, allocatable :: changedTriplets(:,:), tripIndex(:)
  integer, allocatable :: scounts(:), displs(:), newExpInt(:,:)
  integer, allocatable :: scatterTrip(:,:), indexVector(:)
  integer, allocatable :: scatterTripInd(:)
  logical, allocatable :: mask(:)
  double precision, allocatable :: newDists(:), newUvec(:)
  double precision, allocatable :: newUfull(:)


contains


  double precision function tmpi_calcAtomMoveEnergy(atomToMove)
    ! Inputs
    integer, intent(in) :: atomToMove


    ! Set up calculation
    moveTime = MPI_Wtime()
    setTime = MPI_Wtime()
    call setUpCalculation(atomToMove)
    setTime = MPI_Wtime() - setTime


    ! Re-calculate changed exponentials
    expTime = MPI_Wtime()
    call changedExponentialCalculation(atomToMove)
    expTime = MPI_Wtime() - expTime


    ! Find the new energy after the move by summing over triplet energies
    sumTime = MPI_Wtime()
    call changedTripletEnergyCalc(atomToMove)
    tmpi_calcAtomMoveEnergy = proposedEnergyData%Utotal
    sumTime = MPI_Wtime() - sumTime


    ! Deallocate all arrays
    call deallocateArrays()
    moveTime = MPI_Wtime() - moveTime


    ! Finalise MPI and print times taken for each step of calculation
    if (processRank .eq. root) then
       call finalTextOutput()
    end if
    call finalAsserts(proposedPositionData%N_a)


  return
  end function tmpi_calcAtomMoveEnergy


  subroutine setUpCalculation(mover)
    implicit none
    integer, intent(in) :: mover

    ! Initial set-up and asserts
    call initialAsserts(proposedPositionData%N_a,proposedPositionData%N_tri, &
                        proposedPositionData%N_distances)
    N_changed_triplets = getTriPerAtom(proposedPositionData%N_a)
    call allocateArrays(proposedPositionData,proposedEnergyData)

    ! Re-calculate interatomicDistances for the new atomic positions
    xTime = MPI_Wtime()
    call updateXdg(mover,proposedPositionData%N_a,proposedPositionData%posArray, &
                   proposedEnergyData%interatomicDistances)
    xTime = MPI_Wtime() - xTime

  return
  end subroutine setUpCalculation


  subroutine changedExponentialCalculation(mover)
    implicit none
    integer, intent(in) :: mover

    ! Find changed exps
    extractTime = MPI_Wtime()
    call extractChangedExps(proposedPositionData%N_a,mover,proposedEnergyData%interatomicDistances, &
                            newExpInt,newDists)
    expUpdateNoRepeat = newDists
    expUpdateIndNoRepeat = newExpInt
    extractTime = MPI_Wtime() - extractTime

    ! Calculate new values of changed exps
    if (allocated(changeExpData)) then
      deallocate(changeExpData)
    end if
    allocate(changeExpData(nArgs,N_tp,proposedPositionData%N_a-1))
    call calculateExponentialsNonAdd(proposedPositionData%N_a-1,N_tp,nArgs,trainData, &
                                     hyperParams(1),newDists, changeExpData)

    ! Update exp array
    call updateExpMatrix(proposedEnergyData,changeExpData, newExpInt, &
                         proposedPositionData%N_a-1)

  return
  end subroutine changedExponentialCalculation


  subroutine changedTripletEnergyCalc(mover)
    implicit none
    integer, intent(in) :: mover

    ! Determine which triplets have changed
    tripTime = MPI_Wtime()
    call getChangedTriplets(mover, changedTriplets,tripIndex)

    ! Prepare trip. data for scattering
    call getTripletScatterData()
    call allocateTripletScatterArrays()

    ! Scatter triplets across processes
    scatterTrip = changedTriplets(1:3,1+displs(processRank+1):displs(processRank+1)+&
                                  scounts(processRank+1))
    scatterTripInd = tripIndex(1+displs(processRank+1):displs(processRank+1)+scounts(processRank+1))
    tripTime = MPI_Wtime() - tripTime

    ! Calculate the non-additive energies for the changed triplets
    tripSumTime = MPI_Wtime()
    call tripletEnergiesNonAdd(scatterTrip,proposedEnergyData%distancesIntMat,triPerProc,N_tp, &
                               proposedPositionData%N_a,N_p,nArgs,Perm,proposedPositionData%N_distances, &
                               proposedEnergyData%expMatrix,alpha,hyperParams(2), newUvec)
    newTripU = sum(newUvec)
    oldTripU = findOldTripU()
    partialDeltaU = newTripU - oldTripU
    tripSumTime = MPI_Wtime() - tripSumTime

    ! Gather the change in triplet energies on the root process for summation
    gatherAndSumTime = MPI_Wtime()
    call MPI_gather(partialDeltaU, 1, MPI_DOUBLE_PRECISION, newUfull, 1, &
                    MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

    ! Find new total non-additive energy after moving an atom
    proposedEnergyData%Utotal = 0d0
    if (processRank .eq. root) then
      proposedEnergyData%Utotal = sum(newUfull)
    end if
    gatherAndSumTime = MPI_Wtime() - gatherAndSumTime

  return
  end subroutine changedTripletEnergyCalc


  double precision function findOldTripU()
    implicit none
    double precision :: oldUvec(triPerProc)

    counter = 1
    do j = 1, triPerProc
      oldUvec(j) = currentEnergyData%tripletEnergies(scatterTripInd(j))
    end do

    findOldTripU = sum(oldUvec)

  return
  end function findOldTripU


  subroutine initialAsserts(N_a,N_tri,N_distances)
    integer, intent(in) :: N_a, N_tri, N_distances
    integer :: N_tri_ex

    if ( processRank == root ) then

       N_tri_ex = N_a**3
       N_tri_ex = N_tri_ex - 3*N_a**2
       N_tri_ex = N_tri_ex + 2*N_a
       N_tri_ex = N_tri_ex / 6

       call assertTrue(N_a > 0, 'tmpi_calcAtomMoveEnergy: should have N_a > 0')
       call assertTrue(N_tri .eq. N_tri_ex, &
       'tmpi_calcAtomMoveEnergy: should have N_tri = (N_a^3 - 3N_a^2 + 2N_a) / 6')
       call assertTrue(N_distances .eq. ((N_a * N_a) - N_a) / 2, &
       'tmpi_calcAtomMoveEnergy: should have N_distances = (N_a^2 - N_a) / 2')

    end if

  end subroutine initialAsserts


  subroutine finalAsserts(N_a)
    integer, intent(in) :: N_a

    if ( processRank == root ) then

       call assertTrue( N_a>0 , 'tmpi_calcAtomMoveEnergy: should have N_a > 0')

    end if

  end subroutine finalAsserts


  subroutine allocateArrays(proposedPosition,proposedEnergy)
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergy
 
    allocate(newExpInt(proposedPosition%N_a-1,2))
    allocate(newDists(proposedPosition%N_a-1))
    allocate(changedTriplets(3,N_changed_triplets))
    allocate(newUfull(clusterSize))
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    if (allocated(proposedEnergy%interatomicDistances)) then
      deallocate(proposedEnergy%interatomicDistances)
    end if
    allocate(proposedEnergy%interatomicDistances(proposedPosition%N_a, &
             proposedPosition%N_a))
    allocate(tripIndex(N_changed_triplets))

  end subroutine allocateArrays


  subroutine getTripletScatterData()

    call getNPerProcNonAdd(N_changed_triplets,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)

  end subroutine getTripletScatterData


  subroutine allocateTripletScatterArrays()

    allocate(scatterTripInd(triPerProc))
    allocate(scatterTrip(3,triPerProc))
    allocate(newUvec(triPerProc))
    if (allocated(expUpdate)) then
      deallocate(expUpdate)
    end if
    !allocate(expUpdate(2*triPerProc))
    allocate(expUpdate(proposedPositionData%N_a-1))
    if (allocated(expUpdateInd)) then
      deallocate(expUpdateInd)
    end if
    !allocate(expUpdateInd(2*triPerProc,2))
    allocate(expUpdateInd(proposedPositionData%N_a-1,2))

  end subroutine allocateTripletScatterArrays


  subroutine updateChangedTripletEnergies()

    proposedEnergyData%tripletEnergies = currentEnergyData%tripletEnergies
    do j = 1, N_changed_triplets

      proposedEnergyData%tripletEnergies(tripIndex(j)) = newUfull(j)

    end do

  end subroutine updateChangedTripletEnergies


  subroutine deallocateArrays()

    deallocate(scatterTripInd)
    deallocate(scatterTrip)
    deallocate(newUvec)
    deallocate(newExpInt)
    deallocate(newDists)
    deallocate(changedTriplets)
    deallocate(newUfull)
    deallocate(scounts)
    deallocate(displs)
    deallocate(tripIndex)

  end subroutine deallocateArrays


  subroutine finalTextOutput()

    if (textOutput) then
      !print *, moveTime, setTime, expTime, sumTime
      print *, sumTime, tripTime, tripSumTime, gatherAndSumTime
    end if

  end subroutine finalTextOutput


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


  subroutine getChangedTriplets(atom, changedTriplets,tripIndex)
    implicit none
    integer, intent(in) :: atom
    integer, intent(out) :: changedTriplets(3,N_changed_triplets)
    integer, intent(out) :: tripIndex(N_changed_triplets)
    integer :: al, be, ga, counter, indCounter

    ! Fill changed triplet array and vector of indices of changed triplets
    counter = 0
    indCounter = 0
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

  return
  end subroutine getChangedTriplets


  subroutine findChangedTriIndex(nPerAt,nAt,atom, triIndex)
    implicit none
    integer, intent(in) :: nPerAt, nAt, atom
    integer, intent(out) :: triIndex(nPerAt)
    integer :: al, be, ga, i, counter

    counter = 0
    i = 0
    do al = 1, nAt-2
      do be = al+1, nAt-1
        do ga = be+1, nAt

          i = i + 1

          if (al .eq. atom) then

            counter = counter + 1
            triIndex(counter) = i

          else if (be .eq. atom) then

            counter = counter + 1
            triIndex(counter) = i

          else if (ga .eq. atom) then

            counter = counter + 1
            triIndex(counter) = i

          end if
        end do
      end do
    end do

  return
  end subroutine findChangedTriIndex


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


end module tmpi_calcAtomMoveEnergy_mod
