module tmpi_calcAtomMoveEnergy_mod
  use mpi_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use global_Flags, only: textOutput
  use assert_module
  implicit none
  include 'mpif.h'

  private
  public tmpi_calcAtomMoveEnergy,getChangedTriplets,extractChangedExps,getTriPerAtom, &
         findChangedTriIndex


  integer :: triPerAt, nTriMax, nTriRe,triPerProc, j, counter
  double precision :: totTime, moveTime, expTime, sumTime, setTime
  double precision :: gatherTripTime, xTime, shareExpTime, tripTime
  double precision :: extractTime
  integer, allocatable :: changedTriplets(:,:), tripIndex(:)
  integer, allocatable :: scounts(:), displs(:), newExpInt(:,:)
  integer, allocatable :: scatterTrip(:,:), indexVector(:)
  integer, allocatable :: expUpdateInd(:,:),expUpdateIndNoRepeat(:,:)
  logical, allocatable :: mask(:)
  double precision, allocatable :: newDists(:), changeExpData(:,:,:)
  double precision, allocatable :: newUfull(:), newUvec(:)
  double precision, allocatable :: expUpdate(:), expUpdateNoRepeat(:)


contains


  function tmpi_calcAtomMoveEnergy(move,proposedPosition,currentEnergyData) result(proposedEnergyData)
    ! Inputs
    integer, intent(in) :: move!, N_move
    type (energiesData), intent(in) :: currentEnergyData
    type (positionData), intent(in) :: proposedPosition

    ! Output
    type (energiesData) :: proposedEnergyData


    call initialAsserts(proposedPosition%N_a,proposedPosition%N_tri, &
                        proposedPosition%N_distances)
    root = 0
    call firstTextOutput()
    totTime = MPI_Wtime()
    call getTriPerAtom(proposedPosition%N_a, triPerAt)
    call allocateArrays(proposedPosition,proposedEnergyData)
    proposedEnergyData = currentEnergyData


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()

       ! Set up
       setTime = MPI_Wtime()
       if (processRank .eq. root) then

          call moveTextOutput(move)

       end if

       ! Re-calculate interatomicDistances for the new atomic positions
       xTime = MPI_Wtime()
       call updateXdg(move,proposedPosition%N_a,proposedPosition%posArray, proposedEnergyData%interatomicDistances)
       xTime = MPI_Wtime() - xTime

       ! Find the indices of the affected exponentials
       extractTime = MPI_Wtime()
       call extractChangedExps(proposedPosition%N_a,move,proposedEnergyData%interatomicDistances, &
                               newExpInt,newDists)
       extractTime = MPI_Wtime() - extractTime

       ! Determine which triplets have undergone a change
       tripTime = MPI_Wtime()
       call changedTripletInfo(move,proposedPosition)
       tripTime = MPI_Wtime() - tripTime
       setTime = MPI_Wtime() - setTime
       expTime = MPI_Wtime()

       ! Prepare trip. data for scattering
       call getTripletScatterData()
       call allocateTripletScatterArrays()

       ! Scatter triplets across processes
       scatterTrip = changedTriplets(1:3,1+displs(processRank+1):displs(processRank+1)+&
                     scounts(processRank+1))
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Find all distances required on each process for the triplet energy calc. for
       ! which the exponentials have changed due to the move
       call getAffectedTripletDistances(move,proposedEnergyData)

       ! Calculate exponentials for the changed, un-repeated distances
       allocate(changeExpData(nArgs,N_tp,size(expUpdateNoRepeat)))
       call calculateExponentialsNonAdd(size(expUpdateNoRepeat),N_tp,nArgs,trainData, &
                                        hyperParams(1),expUpdateNoRepeat, changeExpData)
       expTime = MPI_Wtime() - expTime

       ! Update the exponential matrix
       sumTime = MPI_Wtime()
       call updateExpMatrix(proposedEnergyData,changeExpData,expUpdateIndNoRepeat, & 
                            size(expUpdateNoRepeat))

       ! Calculate the non-additive energies for the changed triplets
       call tripletEnergiesNonAdd(scatterTrip,proposedEnergyData%distancesIntMat,triPerProc,N_tp, &
                                  proposedPosition%N_a,N_p,nArgs,Perm,proposedPosition%N_distances, &
                                  proposedEnergyData%expMatrix,alpha,hyperParams(2), newUvec)
       sumTime = MPI_Wtime() - sumTime

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

       shareExpTime = MPI_Wtime()
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call shareChangedExponentials(proposedEnergyData)
       shareExpTime = MPI_Wtime() - shareExpTime

    moveTime = MPI_Wtime() - moveTime


    ! Deallocate all arrays
    call deallocateArrays()


    ! Finalise MPI and print times taken for each step of calculation
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       call finalTextOutput()

    end if
    call finalAsserts(proposedPosition%N_a)


    ! Broadcast the aspects of proposedEnergyData that are unique to root to all other processes
    call broadcastEnergyData(proposedEnergyData,proposedPosition)


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
    allocate(changedTriplets(3,triPerAt))
    allocate(newUfull(triPerAt))
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    if (allocated(proposedEnergyData%interatomicDistances)) then
      deallocate(proposedEnergyData%interatomicDistances)
    end if
    allocate(proposedEnergyData%interatomicDistances(proposedPosition%N_a, &
             proposedPosition%N_a))
    allocate(tripIndex(triPerAt))

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


  subroutine changedTripletInfo(move,proposedPosition)
    integer :: move
    type (positionData) :: proposedPosition

    call getChangedTriplets(move,proposedPosition%N_a,triPerAt, changedTriplets)
    call findChangedTriIndex(triPerAt,proposedPosition%N_a,move, tripIndex)

  end subroutine changedTripletInfo


  subroutine broadcastEnergyData(proposedEnergyData,proposedPosition)
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergyData

!    call MPI_Bcast(proposedEnergyData%tripletEnergies, proposedPosition%N_tri, &
!                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(proposedEnergyData%Utotal, 1, MPI_DOUBLE_PRECISION,  root, &
                   MPI_COMM_WORLD, ierror)

  end subroutine broadcastEnergyData


  subroutine updateExpMatrix(proposedEnergyData,changeData,expInd,length)
    integer :: indj, length, expInd(length,2)
    double precision :: changeData(nArgs,N_tp,length)
    type (energiesData) :: proposedEnergyData

    do j = 1, length

      indj = proposedEnergyData%distancesIntMat(expInd(j,1), expInd(j,2))
      proposedEnergyData%expMatrix(1:nArgs,1:N_tp,indj) = changeData(1:nArgs,1:N_tp,j)

    end do

  end subroutine updateExpMatrix


  subroutine getTripletScatterData()

    call getNPerProcNonAdd(triPerAt,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)

  end subroutine getTripletScatterData


  subroutine allocateTripletScatterArrays()

    allocate(scatterTrip(3,triPerProc))
    allocate(newUvec(triPerProc))
    allocate(expUpdate(2*triPerProc))
    allocate(expUpdateInd(2*triPerProc,2))

  end subroutine allocateTripletScatterArrays


  subroutine updateChangedTripletEnergies(currentEnergyData, proposedEnergyData)
    type (energiesData) :: currentEnergyData, proposedEnergyData

    proposedEnergyData%tripletEnergies = currentEnergyData%tripletEnergies
    do j = 1, triPerAt

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
    deallocate(changeExpData)
    deallocate(newExpInt)
    deallocate(newDists)
    deallocate(changedTriplets)
    deallocate(newUfull)
    deallocate(scounts)
    deallocate(displs)
    deallocate(expUpdate)
    deallocate(expUpdateInd)
    deallocate(expUpdateNoRepeat)
    deallocate(expUpdateIndNoRepeat)
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
    allocate(expUpdateNoRepeat, source=expUpdate(indexVector))
    deallocate(indexVector)

    ! Do the same for the distance indices
    mask = .TRUE.
    do j = 2*triPerProc,2,-1
      mask(j) = .NOT.(ANY(expUpdateInd(:j-1,1)==expUpdateInd(j,1) .AND. &
                expUpdateInd(:j-1,2)==expUpdateInd(j,2)))
    end do
    allocate(indexVector, source=PACK([(j,j=1,2*triPerProc)], mask))
    allocate(expUpdateIndNoRepeat, source=expUpdateInd(indexVector,:))
    deallocate(mask,indexVector)

  return
  end subroutine getAffectedTripletDistances


  subroutine shareChangedExponentials(proposedEnergyData)
    implicit none
    integer :: length, lengthVec(clusterSize), sumLength
    integer :: maxLength, reLength
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
      displs(1) = 0
      do j = 2, clusterSize
        displs(j) = sum(lengthVec(1:j-1))
      end do
      sumLength = sum(lengthVec)

      ! Gather changed exponentials to root
      allocate(changeExpMat(nArgs,N_tp,sumLength))
      call MPI_gatherv(changeExpData, N_tp*nArgs*length, MPI_DOUBLE_PRECISION, &
                       changeExpMat, N_tp*nArgs*lengthVec, N_tp*nArgs*displs, &
                       MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

      ! Gather indices if changed exponentials to root
      allocate(changeExpInd(2,sumLength))
      allocate(expUpdateNoRepeatTrans(2,length))
      expUpdateNoRepeatTrans = transpose(expUpdateIndNoRepeat)
      call MPI_gatherv(expUpdateNoRepeatTrans, 2*length, MPI_INT, changeExpInd, &
                       2*lengthVec, 2*displs, MPI_INT, 1, MPI_COMM_WORLD, ierror)
    
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


  subroutine getTriPerAtom(nAt, nPer)
    implicit none
    integer, intent(in) :: nAt
    integer, intent(out) :: nPer
    integer :: a, b

    ! Determine the number of triplets each atom is involved in
    nPer = 0
    do a = 2, nAt-1
      do b = a+1, nAt
        nPer = nPer + 1
      end do
    end do

  return
  end subroutine getTriPerAtom


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
