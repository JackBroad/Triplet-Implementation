module toyMove_Module
  use mpi
  use mpi_variables
  use expShare_variables
  use time_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use tmpi_calcFullSimBoxEnergy_mod
  use tmpi_calcAtomMoveEnergy_mod
  use global_Flags, only: textOutput
  use fullBoxModule
  use atomMoveModule
  use assert_module
  implicit none

  integer :: triPerAt, nTriMax, nTriRe,triPerProc, j, counter
  double precision :: gatherTripTime, shareExpTime
  type (energiesData) :: toyEnergyData
  type (positionData) :: movePositionData
  integer, allocatable :: changedTriplets(:,:), tripIndex(:)
  integer, allocatable :: scounts(:), displs(:), newExpInt(:,:)
  integer, allocatable :: scatterTrip(:,:), indexVector(:)
  logical, allocatable :: mask(:)
  double precision, allocatable :: newDists(:), newUvec(:)
  double precision, allocatable :: newUfull(:)


contains


  function toyMoveMinimalScatter(oldPositionData) result(moveEnergyData)
    implicit none
    type (positionData) :: oldPositionData
    type (energiesData) :: moveEnergyData
    double precision :: uTotPerProc, U, randomNo
    double precision, allocatable :: uVector(:)
    integer :: nChangedDists, i
    integer :: atomToMove=1

    ! Set up
    call allocateToyDataStructArrays(oldPositionData,moveEnergyData)
    call setUpEnergyData(oldPositionData,moveEnergyData)
    call makeDisIntMatNonAdd(oldPositionData%N_a,moveEnergyData%distancesIntMat)
    do i = 1, oldPositionData%N_distances
      call random_number(randomNo)
      randomNo = 10*randomNo
      expArray(:,:,i) = randomNo
    end do
    moveEnergyData%triMat = makeTripletMatrix(oldPositionData%N_a,oldPositionData%N_tri)
    moveTime = MPI_Wtime()
    triPerAt = getTriPerAtom(oldPositionData%N_a)
    nChangedDists = oldPositionData%N_a-1

    ! Do exp calculation
    if (allocated(newDists)) then
      deallocate(newDists)
    end if
    if (allocated(newExpInt)) then
      deallocate(newExpInt)
    end if
    allocate(newDists(nChangedDists),newExpInt(nChangedDists,2))
    call extractChangedExps(oldPositionData%N_a,atomToMove,moveEnergyData%interatomicDistances, &
                            newExpInt,newDists)
    if (allocated(changeExpData)) then
      deallocate(changeExpData)
    end if
    allocate(changeExpData(nArgs,N_tp,nChangedDists))
    call calculateExponentialsNonAdd(nChangedDists,N_tp,nArgs,trainData, &
                                     hyperParams(1),newDists, changeExpData)
    call updateExpMatrix(changeExpData,newExpInt,nChangedDists)

    ! Find no. of triplets per proc
    setTime = MPI_Wtime()
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    call getNPerProcNonAdd(triPerAt,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)

    ! Scatter triplets
    allocate(scatterTrip(3,triPerProc))
    do i = 1, triPerProc
       scatterTrip(:,i) = moveEnergyData%triMat(:,i)
    end do
    allocate(newUvec(triPerProc))
    allocate(newUfull(oldPositionData%N_tri))
    setTime = MPI_Wtime() - setTime

    ! Sum over the triplets
    sumTime = MPI_Wtime()
    call toyTripletSum(scatterTrip,triPerProc,expArray, &
                       oldPositionData%N_distances,moveEnergyData%distancesIntMat, &
                       oldPositionData%N_a, newUvec)
    uTotPerProc = sum(newUvec)
    sumTime = MPI_Wtime() - sumTime
    gatherTripTime = MPI_Wtime()
    allocate(uVector(clusterSize))
    call MPI_gather(uTotPerProc, 1, MPI_DOUBLE_PRECISION, uVector, &
                    1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                    ierror)
    if (processRank .eq. root) then
      U = sum(uVector)
    end if
    gatherTripTime = MPI_Wtime() - gatherTripTime
    moveTime = MPI_Wtime() - moveTime
    if (processRank .eq. root) then
      print *, moveTime, setTime, sumTime, gatherTripTime
    end if

    call deallocateAllArrays(moveEnergyData)

  return
  end function toyMoveMinimalScatter


  function toyMoveDistScatter(oldPositionData) result(moveEnergyData)
    implicit none
    type (positionData) :: oldPositionData
    type (energiesData) :: moveEnergyData
    double precision :: randomNo
    double precision, allocatable :: changedDists(:), scatterDists(:)
    double precision, allocatable :: changeExpData(:,:,:), changeExpMat(:,:,:)
    integer :: nChangedDists, distsPerProc, i

    ! Set up
    call allocateToyDataStructArrays(oldPositionData,moveEnergyData)
    call setUpEnergyData(oldPositionData,moveEnergyData)
    call makeDisIntMatNonAdd(oldPositionData%N_a,moveEnergyData%distancesIntMat)
    do i = 1, oldPositionData%N_distances
      call random_number(randomNo)
      randomNo = 10*randomNo
      expArray(:,:,i) = randomNo
    end do
    moveTime = MPI_Wtime()
    !expArray = randomNo
    triPerAt = getTriPerAtom(oldPositionData%N_a)
    nChangedDists = oldPositionData%N_a-1

    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    allocate(changedDists(nChangedDists))
    call getNPerProcNonAdd(nChangedDists,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    changedDists = 1d0
    distsPerProc = scounts(processRank+1)
    allocate(scatterDists(distsPerProc))

    ! Do the scattering
    scatterDists = changedDists(1+displs(processRank+1): &
                   displs(processRank+1)+scounts(processRank+1))

    ! Dummy exp update
    allocate(changeExpData(nArgs,N_tp,distsPerProc))
    allocate(changeExpMat(nArgs,N_tp,nChangedDists))
    do i = 1, distsPerProc
      changeExpData(:,:,i) = i*2d0 + i
    end do
    do i = 1, nChangedDists
      changeExpMat(:,:,i) = i*2d0 + i
    end do

    ! Gather in changed exps
    call MPI_gatherv(changeExpData, N_tp*nArgs*distsPerProc, MPI_DOUBLE_PRECISION, &
                     changeExpMat, N_tp*nArgs*scounts, N_tp*nArgs*displs, &
                     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(changeExpMat, N_tp*nArgs*nChangedDists, &
                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

    ! Update expMatrix
    expArray(nArgs,N_tp,1:nChangedDists) = &
    changeExpMat(nArgs,N_tp,:)

    call getNPerProcNonAdd(triPerAt,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)
    allocate(scatterTrip(3,triPerProc)) ! must allocate after above call
    allocate(newUvec(triPerProc))
    allocate(newUfull(oldPositionData%N_tri))

    ! Do the scattering
    scatterTrip = moveEnergyData%triMat(1:3,1+displs(processRank+1):&
                  displs(processRank+1)+scounts(processRank+1))

    ! Sum over the triplets
    call toyTripletSum(scatterTrip,triPerProc,expArray, &
                       oldPositionData%N_distances,moveEnergyData%distancesIntMat, &
                       oldPositionData%N_a, newUvec)
    call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, &
                     scounts, displs, MPI_DOUBLE_PRECISION, root, &
                     MPI_COMM_WORLD, ierror)

    ! Change triPerAt triplets in each triplet vector on root
    if (processRank .eq. root) then
      moveEnergyData%tripletEnergies(1:triPerAt) = 3d0
      moveEnergyData%Utotal = sum(moveEnergyData%tripletEnergies)
    end if
    moveTime = MPI_Wtime() - moveTime
    if (processRank .eq. root) then
      print *, moveTime, 0, 0, 0
    end if

    call deallocateAllArrays(moveEnergyData)

  return
  end function toyMoveDistScatter


  function toyMove(oldPositionData) result(moveEnergyData)
    implicit none
    double precision :: randomNo
    type (positionData) :: oldPositionData
    type (energiesData) :: moveEnergyData
    integer :: nChangedDists, i

    ! Set up
    call allocateToyDataStructArrays(oldPositionData,moveEnergyData)
    call setUpEnergyData(oldPositionData,moveEnergyData)
    call makeDisIntMatNonAdd(oldPositionData%N_a,moveEnergyData%distancesIntMat)
    call random_number(randomNo)
    randomNo = 10*randomNo

    ! Each move changes triPerAt triplets and N_a-1 distances,
    ! so we only need to share the first triPerAt cols of the
    ! triplet matrix and alter N_a-1 sets of exponentials in
    ! expMat
    moveTime = MPI_Wtime()
    expArray = randomNo
    triPerAt = getTriPerAtom(oldPositionData%N_a)
    nChangedDists = oldPositionData%N_a-1

    do i = 1, nChangedDists
      expArray(:,:,i) = i*2d0 + i
    end do

    ! Set up arrays for scatter
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    call getNPerProcNonAdd(triPerAt,clusterSize, nTriMax,nTriRe)
    call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
    triPerProc = scounts(processRank+1)
    allocate(scatterTrip(3,triPerProc)) ! must allocate after above call
    allocate(newUvec(triPerProc))
    allocate(newUfull(oldPositionData%N_tri))

    ! Do the scattering
    scatterTrip = moveEnergyData%triMat(1:3,1+displs(processRank+1):&
                  displs(processRank+1)+scounts(processRank+1))

    ! Sum over the triplets
    call toyTripletSum(scatterTrip,triPerProc,expArray, &
                       oldPositionData%N_distances,moveEnergyData%distancesIntMat, &
                       oldPositionData%N_a, newUvec)
    call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, &
                     scounts, displs, MPI_DOUBLE_PRECISION, root, &
                     MPI_COMM_WORLD, ierror)

    ! Change triPerAt triplets in each triplet vector on root
    if (processRank .eq. root) then
      moveEnergyData%tripletEnergies(1:triPerAt) = 3d0
      moveEnergyData%Utotal = sum(moveEnergyData%tripletEnergies)
    end if
    moveTime = MPI_Wtime() - moveTime
    if (processRank .eq. root) then
      print *, moveTime, 0, 0, 0
    end if

    call deallocateAllArrays(moveEnergyData)

  return
  end function toyMove


  subroutine deallocateAllArrays(moveEnergyData)
    implicit none
    type (energiesData) :: moveEnergyData

    deallocate(scatterTrip,newUvec,newUfull)
    deallocate(scounts,displs)
    deallocate(toyEnergyData%distancesIntMat,moveEnergyData%distancesIntMat)
    deallocate(toyEnergyData%triMat,moveEnergyData%triMat)
    deallocate(toyEnergyData%interatomicDistances,moveEnergyData%interatomicDistances)
    deallocate(toyEnergyData%tripletEnergies,moveEnergyData%tripletEnergies)
    deallocate(expArray)

  end subroutine deallocateAllArrays


  subroutine toyTripletSum(tripletData,nTrip,expMat,nDists,intMat,nAt, uVec)
    implicit none
    integer, intent(in) :: nTrip, nAt, tripletData(3,nTrip), nDists, intMat(nAt,nAt)
    double precision, intent(in) :: expMat(nArgs,N_tp,nDists)
    double precision, intent(out) :: uVec(nTrip)
    integer :: i, k, al, be, ga, alDis, beDis, gaDis, counter
    double precision :: firstSum, secndSum, expProd
    character (len=10) :: clusterString, frmt

    frmt = '(I5.5)'
    write(clusterString,frmt) clusterSize
    N_p = 6

    do i = 1, nTrip
 
      counter = counter + 1

      al = tripletData(1,i)
      be = tripletData(2,i)
      ga = tripletData(3,i)

      alDis = intMat(al,be)
      beDis = intMat(al,ga)
      gaDis = intMat(be,ga)

      firstSum = 0d0

      do j = 1, N_tp

        counter = counter + 1

        secndSum = 0d0

        do k = 1, N_p

          counter = counter + 1

          expProd = expMat(Perm(k,1),j,alDis) * expMat(Perm(k,2),j,beDis) * &
                    expMat(Perm(k,3),j,gaDis)

          secndSum = secndSum + expProd

        end do

        firstSum = firstSum + (secndSum*0.5)

      end do

      uVec(i) = 0.5 * firstSum

    end do

  end subroutine toyTripletSum


  subroutine setUpEnergyData(oldPositionData,moveEnergyData)
    implicit none
    type (positionData), intent(in) :: oldPositionData
    type (energiesData) :: moveEnergyData

    movePositionData = oldPositionData
    toyEnergyData%Utotal = 1d0
    toyEnergyData%triMat = makeTripletMatrix(oldPositionData%N_a,oldPositionData%N_tri)
    toyEnergyData%interatomicDistances = 1d0
    toyEnergyData%tripletEnergies = 1d0
    expArray = 1d0
    moveEnergyData = toyEnergyData

  return
  end subroutine setUpEnergyData


  subroutine allocateToyDataStructArrays(oldPositionData,moveEnergyData)
    implicit none
    type (positionData), intent(in) :: oldPositionData
    type (energiesData) :: moveEnergyData

    allocate(toyEnergyData%distancesIntMat(oldPositionData%N_a,oldPositionData%N_a))
    !allocate(moveEnergyData%distancesIntMat(oldPositionData%N_a,oldPositionData%N_a))
    allocate(toyEnergyData%triMat(3,oldPositionData%N_tri))
    allocate(moveEnergyData%triMat(3,oldPositionData%N_tri))
    allocate(toyEnergyData%interatomicDistances(oldPositionData%N_a,oldPositionData%N_a))
    allocate(moveEnergyData%interatomicDistances(oldPositionData%N_a,oldPositionData%N_a))
    allocate(toyEnergyData%tripletEnergies(oldPositionData%N_tri))
    allocate(moveEnergyData%tripletEnergies(oldPositionData%N_tri))
    allocate(expArray(nArgs,N_tp,oldPositionData%N_distances))

  return
  end subroutine


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
    do j = 1, triPerAt

      proposedEnergyData%tripletEnergies(tripIndex(j)) = newUfull(j)

    end do

  end subroutine updateChangedTripletEnergies


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


end module toyMove_Module
