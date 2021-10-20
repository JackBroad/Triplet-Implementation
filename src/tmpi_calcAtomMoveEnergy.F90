module tmpi_calcAtomMoveEnergy_mod
  use mpi_variables
  use expShare_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use global_Flags, only: textOutput
  use assert_module
  implicit none
  include 'mpif.h'

  private
  public tmpi_calcAtomMoveEnergy,getChangedTripletData,extractChangedExps,getTriPerAtom, &
         findChangedTriIndex, getChangedTriplets, updateXdg, getNPerProcNonAdd, getVarrays, &
         getTripletScatterData


  integer :: N_changed_triplets, nTriMax, nTriRe,triPerProc, j, counter
  double precision :: totTime, moveTime, expTime, sumTime, setTime
  double precision :: gatherTripTime, xTime, shareExpTime, tripTime
  double precision :: extractTime
  integer, allocatable :: changedTriplets(:,:), tripIndex(:)
  integer, allocatable :: scounts(:), displs(:), newExpInt(:,:)
  integer, allocatable :: scatterTrip(:,:), indexVector(:)
  logical, allocatable :: mask(:)
  double precision, allocatable :: newDists(:), newUvec(:)
  double precision, allocatable :: newUfull(:)


contains


  function tmpi_calcAtomMoveEnergy(atomToMove,proposedPosition,currentEnergyData) result(proposedEnergyData)
    ! Inputs
    integer, intent(in) :: atomToMove!, N_move
    type (energiesData), intent(in) :: currentEnergyData
    type (positionData), intent(in) :: proposedPosition

    ! Output
    type (energiesData) :: proposedEnergyData


    call initialAsserts(proposedPosition%N_a,proposedPosition%N_tri, &
                        proposedPosition%N_distances)
    textOutput = .false.
    call firstTextOutput()
    moveTime = MPI_Wtime()
    N_changed_triplets = getTriPerAtom(proposedPosition%N_a)
    call allocateArrays(proposedPosition,proposedEnergyData)
    proposedEnergyData = currentEnergyData


    ! Set up
    setTime = MPI_Wtime()
    if (processRank .eq. root) then

       call moveTextOutput(atomToMove)

    end if


    ! Re-calculate interatomicDistances for the new atomic positions
    xTime = MPI_Wtime()
    call updateXdg(atomToMove,proposedPosition%N_a,proposedPosition%posArray, proposedEnergyData%interatomicDistances)
    xTime = MPI_Wtime() - xTime

    ! Find the indices of the affected exponentials
    extractTime = MPI_Wtime()
    call extractChangedExps(proposedPosition%N_a,atomToMove,proposedEnergyData%interatomicDistances, &
                            newExpInt,newDists)
    extractTime = MPI_Wtime() - extractTime

    ! Determine which triplets have undergone a change
    tripTime = MPI_Wtime()
    call getChangedTripletInfo(atomToMove,proposedPosition)
    !call getChangedTripletData(atomToMove,proposedPosition%N_a,proposedPosition%N_tri,triPerAt, &
    !                           proposedEnergyData%triMat, tripIndex,changedTriplets)
    tripTime = MPI_Wtime() - tripTime
    setTime = MPI_Wtime() - setTime
    expTime = MPI_Wtime()

    ! Prepare trip. data for scattering
    call getTripletScatterData()
    call allocateTripletScatterArrays()

    ! Scatter triplets across processes
    scatterTrip = changedTriplets(1:3,1+displs(processRank+1):displs(processRank+1)+&
                  scounts(processRank+1))

    ! Find all distances required on each process for the triplet energy calc. for
    ! which the exponentials have changed due to the move
    call getAffectedTripletDistances(atomToMove,proposedEnergyData)

    ! Calculate exponentials for the changed, un-repeated distances
    if (allocated(changeExpData)) then
      deallocate(changeExpData)
    end if
    allocate(changeExpData(nArgs,N_tp,size(expUpdateNoRepeat)))
    call calculateExponentialsNonAdd(size(expUpdateNoRepeat),N_tp,nArgs,trainData, &
                                     hyperParams(1),expUpdateNoRepeat, changeExpData)

    ! Update the exponential matrix
    call updateExpMatrix(proposedEnergyData,changeExpData,expUpdateIndNoRepeat, & 
                         size(expUpdateNoRepeat))
    expTime = MPI_Wtime() - expTime

    ! Calculate the non-additive energies for the changed triplets
    sumTime = MPI_Wtime()
    call tripletEnergiesNonAdd(scatterTrip,proposedEnergyData%distancesIntMat,triPerProc,N_tp, &
                               proposedPosition%N_a,N_p,nArgs,Perm,proposedPosition%N_distances, &
                               proposedEnergyData%expMatrix,alpha,hyperParams(2), newUvec)

    ! Gather the changed triplet energies on the root process for summation
    gatherTripTime = MPI_Wtime()
    call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, scounts, &
                     displs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    gatherTripTime = MPI_Wtime() - gatherTripTime

    ! Find total change in non-additive energy from moving an atom
    if (processRank .eq. root) then

       ! Update energies of changed triplets
       call updateChangedTripletEnergies(currentEnergyData, proposedEnergyData)

       ! Evaluate total non-add energy after changes and print it to screen
       call totalEnergyNonAdd(proposedEnergyData%tripletEnergies,proposedPosition%N_tri, &
                              proposedEnergyData%Utotal)
       call energyTextOutput(proposedEnergyData)

    end if
    sumTime = MPI_Wtime() - sumTime


    ! Deallocate all arrays
    call deallocateArrays()


    moveTime = MPI_Wtime() - moveTime


    ! Finalise MPI and print times taken for each step of calculation
    if (processRank .eq. root) then

       call finalTextOutput()
       print *, moveTime, setTime, expTime, sumTime

    end if
    call finalAsserts(proposedPosition%N_a)


  return
  end function tmpi_calcAtomMoveEnergy


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


  subroutine allocateArrays(proposedPosition,proposedEnergyData)
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergyData
 
    allocate(newExpInt(2,proposedPosition%N_a-1))
    allocate(newDists(proposedPosition%N_a-1))
    allocate(changedTriplets(3,N_changed_triplets))
    Allocate(newUfull(N_changed_triplets))
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    if (allocated(proposedEnergyData%interatomicDistances)) then
      deallocate(proposedEnergyData%interatomicDistances)
    end if
    allocate(proposedEnergyData%interatomicDistances(proposedPosition%N_a, &
             proposedPosition%N_a))
    allocate(tripIndex(N_changed_triplets))

  end subroutine allocateArrays


  subroutine firstTextOutput()

    if (processRank .eq. root) then
       if (textOutput) then

         !print *, ' '
         !print *, ' '
         !print *, '========================'
         !print *, 'Beginning non-additive calculation for atom move'
         !print *, ' '

      end if
    end if

  end subroutine firstTextOutput


  subroutine moveTextOutput(move)
    integer :: move

    if (textOutput) then

      !print *, '------------------------'
      !print *, "Moving atom", move
      !print *, "                 "

    end if

  end subroutine moveTextOutput


  subroutine getTripletScatterData()

    call getNPerProcNonAdd(N_changed_triplets,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)

  end subroutine getTripletScatterData


  subroutine allocateTripletScatterArrays()

    allocate(scatterTrip(3,triPerProc))
    allocate(newUvec(triPerProc))
    if (allocated(expUpdate)) then
      deallocate(expUpdate)
    end if
    allocate(expUpdate(2*triPerProc))
    if (allocated(expUpdateInd)) then
      deallocate(expUpdateInd)
    end if
    allocate(expUpdateInd(2*triPerProc,2))

  end subroutine allocateTripletScatterArrays


  subroutine updateChangedTripletEnergies(currentEnergyData, proposedEnergyData)
    type (energiesData) :: currentEnergyData, proposedEnergyData

    proposedEnergyData%tripletEnergies = currentEnergyData%tripletEnergies
    do j = 1, N_changed_triplets

      proposedEnergyData%tripletEnergies(tripIndex(j)) = newUfull(j)

    end do

  end subroutine updateChangedTripletEnergies


  subroutine energyTextOutput(proposedEnergyData)
    type (energiesData) :: proposedEnergyData

    if (textOutput) then

      !print *, "The non-additive energy after the move is", &
      !         proposedEnergyData%Utotal
      !print *, '------------------------'
      !print *, ' '

    end if

  end subroutine energyTextOutput


  subroutine deallocateArrays()

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

      !print *, "The time taken to do all moves was", moveTime, &
      !         "seconds"
      !print *, "The total time for the program to run was", totTime, &
      !         "seconds"
      !print *, "The time for the exp calc was", expTime
      !print *, "The time for the sum was", sumTime
      !print *, "The time to gather the triplet energies was", gatherTripTime
      !print *, "The time for the set-up was", setTime
      !print *, "The time for the Xdg set-up was", xTime
      !print *, "The time for the exp extraction was", extractTime
      !print *, "The time for the triplet set up was ", tripTime
      !print *, ' '
      !print *, 'Non-additive calculation for atom move complete'
      print *, moveTime,expTime,sumTime,shareExpTime,gatherTripTime,setTime
      !print *, '========================'
      !print *, ' '
      print *, ' '

    end if

  end subroutine finalTextOutput


  subroutine getAffectedTripletDistances(move,proposedEnergyData)
    implicit none
    integer, intent(in) :: move
    type (energiesData) :: proposedEnergyData

    counter = 0
    do j = 1, triPerProc
      if (scatterTrip(1,j) .eq. move) then

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(move,scatterTrip(2,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(2,j)
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(move,scatterTrip(3,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(3,j)

      else if (scatterTrip(2,j) .eq. move) then

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(scatterTrip(1,j),move)
        expUpdateInd(counter,1) = scatterTrip(1,j)
        expUpdateInd(counter,2) = move
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(move,scatterTrip(3,j))
        expUpdateInd(counter,1) = move
        expUpdateInd(counter,2) = scatterTrip(3,j)

      else

        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(scatterTrip(1,j),move)
        expUpdateInd(counter,1) = scatterTrip(1,j)
        expUpdateInd(counter,2) = move
        counter = counter + 1
        expUpdate(counter) = &
        proposedEnergyData%interatomicDistances(scatterTrip(2,j),move)
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
    integer, intent(out) :: change(2,num-1)
    double precision, intent(out) :: dists(num-1)
    integer :: i, j

    do i = 1, num
      if (i .lt. atInd) then

        change(1,i) = i
        change(2,i) = atInd

      else if (i .gt. atInd) then

        change(1,i-1) = atInd
        change(2,i-1) = i

      end if
    end do

    do j = 1, num-1

      dists(j) = Xdg(change(1,j),change(2,j))

    end do

  return
  end subroutine extractChangedExps


  subroutine getChangedTripletData(atom,nAt,nTri,nPerAt,tripMat, tripletIndex,changedTriplets)
    implicit none
    integer, intent(in) :: atom, nAt, nTri, nPerAt, tripMat(3,nTri)
    integer, intent(out) :: changedTriplets(3,nPerAt), tripletIndex(nPerAt)
    integer :: i, counter
    logical :: mask(nTri)

    counter = 0
    mask = any(tripMat .eq. atom, 1) ! Vector that's 'true' for any col of tripMat containing moved atom
    do i = 1, nTri
      if (mask(i) .eqv. .true.) then
        counter = counter + 1
        tripletIndex(counter) = i
        changedTriplets(:,counter) = tripMat(:,i)
      end if
    end do

  return
  end subroutine getChangedTripletData


  subroutine getChangedTripletInfo(move,proposedPosition)
    integer :: move
    type (positionData) :: proposedPosition

    call getChangedTriplets(move,proposedPosition%N_a,N_changed_triplets, changedTriplets)
    call findChangedTriIndex(N_changed_triplets,proposedPosition%N_a,move, tripIndex)

  end subroutine getChangedTripletInfo


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


  subroutine getChangedTriplets(atom,nAt,nPerAt, changedTriplets)
    implicit none
    integer, intent(in) :: atom, nAt, nPerAt
    integer, intent(out) :: changedTriplets(3,nPerAt)
    integer :: al, be, ga, counter

    ! Fill changed triplet array
    counter = 0
    if (atom .eq. 1) then
      al = atom
      do be = al+1, nAt-1
        do ga = be+1, nAt

          counter = counter + 1
          changedTriplets(1,counter) = al
          changedTriplets(2,counter) = be
          changedTriplets(3,counter) = ga

        end do
      end do
    else if (atom .eq. nAt) then
      ga = atom
      do al = 1, nAt-2
        do be = al+1, nAt-1

          counter = counter + 1
          changedTriplets(1,counter) = al
          changedTriplets(2,counter) = be
          changedTriplets(3,counter) = ga

        end do
      end do
    else
      do al = 1, nAt-2
        do be = al+1, nAt-1
          do ga = be+1, nAt
            if (al .eq. atom) then

              counter = counter + 1
              changedTriplets(1,counter) = al
              changedTriplets(2,counter) = be
              changedTriplets(3,counter) = ga

            else if (be .eq. atom) then

              counter = counter + 1
              changedTriplets(1,counter) = al
              changedTriplets(2,counter) = be
              changedTriplets(3,counter) = ga

            else if (ga .eq. atom) then

              counter = counter + 1
              changedTriplets(1,counter) = al
              changedTriplets(2,counter) = be
              changedTriplets(3,counter) = ga

            end if
          end do
        end do
      end do
    end if

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
    
  
end module tmpi_calcAtomMoveEnergy_mod
