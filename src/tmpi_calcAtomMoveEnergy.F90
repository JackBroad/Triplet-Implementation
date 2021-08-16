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


  integer :: triPerAt, nPerProc, nExpMax, nExpRe, nTriMax, nTriRe
  integer :: triPerProc, j
  double precision :: totTime, moveTime, expTime, sumTime, setTime
  double precision :: gatherExpTime, gatherTripTime, xTime
  double precision :: extractTime, tripTime
  integer, allocatable :: changedTriplets(:,:), tripIndex(:)
  integer, allocatable :: scounts(:), displs(:), newExpInt(:,:)
  integer, allocatable :: scatterTrip(:,:)
  double precision, allocatable :: newDists(:), changeExpMat(:,:,:)
  double precision, allocatable :: newUfull(:), changeExpData(:,:,:)
  double precision, allocatable :: newUvec(:), scatterDists(:)


contains


  function tmpi_calcAtomMoveEnergy(N_move,move,proposedPosition,currentEnergyData) result(proposedEnergies)
    ! Inputs
    integer, intent(in) :: N_move, move
    type (energiesData), intent(in) :: currentEnergyData
    type (positionData), intent(in) :: proposedPosition

    ! Output
    type (energiesData) :: proposedEnergies

    ! Local variables
    integer :: i


    call initialAsserts(proposedPosition%N_a,proposedPosition%N_tri, &
                        proposedPosition%N_distances)
    root = 0
    call firstTextOutput()
    totTime = MPI_Wtime()
    call getTriPerAtom(proposedPosition%N_a, triPerAt)
    call allocateArrays(proposedPosition,proposedEnergies)
    proposedEnergies = currentEnergyData


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()
    do i = 1, N_move

       ! Set up on root
       setTime = MPI_Wtime()
       if (processRank .eq. root) then

          call moveTextOutput(move)

       end if

       xTime = MPI_Wtime()
       ! Re-calculate interatomicDistances for the new atomic positions
       call makeXdgNonAdd(proposedPosition%N_a,proposedPosition%posArray, proposedEnergies%interatomicDistances)
       !call updateXdg(move,proposedPosition%N_a,proposedPosition%posArray, proposedEnergies%interatomicDistances)
       xTime = MPI_Wtime() - xTime

       extractTime = MPI_Wtime()
       ! Find the indices of the affected exponentials
       call extractChangedExps(proposedPosition%N_a,move,proposedEnergies%interatomicDistances, &
                               newExpInt,newDists)
       extractTime = MPI_Wtime() - extractTime

       tripTime = MPI_Wtime()
       ! Determine which triplets have undergone a change
       call changedTripletInfo(move,proposedPosition)
       tripTime = MPI_Wtime() - tripTime
       setTime = MPI_Wtime() - setTime

       ! Prepare dist. data for scattering
       expTime = MPI_Wtime()
       call getDistScatterData(proposedPosition)

       ! Allocate arrays for dist. scattering on first loop
       if (i .eq. 1) then

         call allocateDistScatterArrays()

       end if

       ! Read distances for each process to take for exp calculation
       scatterDists = newDists(1+displs(processRank+1):displs(processRank+1)+&
                      scounts(processRank+1))

       ! Update exponentials of distances affected by the move
       call calculateExponentialsNonAdd(nPerProc,N_tp,nArgs,trainData, &
                                        hyperParams(1),scatterDists, &
                                        changeExpData)
       expTime = MPI_Wtime() - expTime
       gatherExpTime = MPI_Wtime()

       ! Gather in all updated exps and broadcast the resultant matrix to all procs
       call MPI_gatherv(changeExpData, N_tp*nArgs*nPerProc, MPI_DOUBLE_PRECISION, &
                        changeExpMat, N_tp*nArgs*scounts, N_tp*nArgs*displs, &
                        MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(changeExpMat, N_tp*nArgs*(proposedPosition%N_a-1), &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Update the exponential matrix
       call updateExpMatrix(proposedPosition, proposedEnergies)
       gatherExpTime = MPI_Wtime() - gatherExpTime

       ! Prepare trip. data from scattering
       sumTime = MPI_Wtime()
       call getTripletScatterData()

       ! Allocate requisite arrays for trip. sactter in first loop
       if (i .eq. 1) then

         call allocateTripletScatterArrays()

       end if

       ! Scatter triplets across processes
       scatterTrip = changedTriplets(1:3,1+displs(processRank+1):displs(processRank+1)+&
                     scounts(processRank+1))

       ! Calculate the non-additive energies for the changed triplets and gather on root
       call tripletEnergiesNonAdd(scatterTrip,proposedEnergies%distancesIntMat,triPerProc,N_tp, &
                                  proposedPosition%N_a,N_p,nArgs,Perm,proposedPosition%N_distances, &
                                  proposedEnergies%expMatrix,alpha,hyperParams(2), newUvec)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       sumTime = MPI_Wtime() - sumTime
       gatherTripTime = MPI_Wtime()
       call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, scounts, &
                        displs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Find total change in non-add energy from moving atom
       if (processRank .eq. root) then

          ! Update energies of changed triplets
          call updateChangedTripletEnergies(currentEnergyData, proposedEnergies)

          ! Evaluate total non-add energy after changes and print it to screen
          call totalEnergyNonAdd(proposedEnergies%tripletEnergies,proposedPosition%N_tri, &
                                 proposedEnergies%Utotal)
          call energyTextOutput(proposedEnergies)

       end if
       gatherTripTime = MPI_Wtime() - gatherTripTime

    end do
    moveTime = MPI_Wtime() - moveTime


    ! Deallocate all arrays
    call deallocateArrays()


    ! Finalise MPI and print times taken for each step of calculation
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       call finalTextOutput(N_move)

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


  subroutine allocateArrays(proposedPosition,proposedEnergies)
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergies
 
    allocate(newExpInt(2,proposedPosition%N_a-1))
    allocate(newDists(proposedPosition%N_a-1))
    allocate(changeExpMat(nArgs,N_tp,proposedPosition%N_a-1))
    allocate(changedTriplets(3,triPerAt))
    allocate(newUfull(triPerAt))
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    allocate(proposedEnergies%interatomicDistances(proposedPosition%N_a, &
             proposedPosition%N_a))
    allocate(tripIndex(triPerAt))

  end subroutine allocateArrays


  subroutine firstTextOutput()

    if (processRank .eq. root) then
       if (textOutput) then

         print *, ' '
         print *, ' '
         print *, '========================'
         print *, 'Beginning non-additive calculation for atom move'
         print *, ' '

      end if
    end if

  end subroutine firstTextOutput


  subroutine moveTextOutput(move)
    integer :: move

    if (textOutput) then

      print *, '------------------------'
      print *, "Moving atom", move
      print *, "                 "

    end if

  end subroutine moveTextOutput


  subroutine changedTripletInfo(move,proposedPosition)
    integer :: move
    type (positionData) :: proposedPosition

    call getChangedTriplets(move,proposedPosition%N_a,triPerAt, changedTriplets)
    call findChangedTriIndex(triPerAt,proposedPosition%N_a,move, tripIndex)

  end subroutine changedTripletInfo


!  subroutine broadcastRootData(proposedPosition,proposedEnergies,move)
!    type (positionData) :: proposedPosition
!    type (energiesData) :: proposedEnergies
!    integer :: move

!    call MPI_Bcast(newExpInt, 2*(proposedPosition%N_a-1), MPI_INT, root, &
!                   MPI_COMM_WORLD, ierror)
!    call MPI_Bcast(proposedEnergies%interatomicDistances, &
!                   proposedPosition%N_a*proposedPosition%N_a, &
!                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
!    call MPI_Bcast(proposedPosition%posArray, 3*proposedPosition%N_a, &
!                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
!    call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)

!  end subroutine broadcastRootData


  subroutine getDistScatterData(proposedPosition)
    type (positionData) :: proposedPosition

     call getNPerProcNonAdd(proposedPosition%N_a-1,clusterSize, nExpMax,nExpRe)
     call getVarrays(clusterSize,nExpMax,nExpRe, scounts,displs)
     nPerProc = scounts(processRank+1)

  end subroutine getDistScatterData


  subroutine allocateDistScatterArrays()

    allocate(scatterDists(nPerProc))
    allocate(changeExpData(nArgs,N_tp,nPerProc))

  end subroutine allocateDistScatterArrays


  subroutine updateExpMatrix(proposedPosition, proposedEnergies)
    integer :: indj
    type (positionData) :: proposedPosition
    type (energiesData) :: proposedEnergies

    do j = 1, proposedPosition%N_a-1

      indj = proposedEnergies%distancesIntMat(newExpInt(1,j),newExpInt(2,j))
      proposedEnergies%expMatrix(1:nArgs,1:N_tp,indj) = changeExpMat(1:nArgs,1:N_tp,j)

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

  end subroutine allocateTripletScatterArrays


  subroutine updateChangedTripletEnergies(currentEnergyData, proposedEnergies)
    type (energiesData) :: currentEnergyData, proposedEnergies

    proposedEnergies%tripletEnergies = currentEnergyData%tripletEnergies
    do j = 1, triPerAt

      proposedEnergies%tripletEnergies(tripIndex(j)) = newUfull(j)

    end do

  end subroutine updateChangedTripletEnergies


  subroutine energyTextOutput(proposedEnergies)
    type (energiesData) :: proposedEnergies

    if (textOutput) then

      print *, "The non-additive energy after the move is", &
               proposedEnergies%Utotal
      print *, '------------------------'
      print *, ' '

    end if

  end subroutine energyTextOutput


  subroutine deallocateArrays()

    deallocate(scatterDists)
    deallocate(scatterTrip)
    deallocate(newUvec)
    deallocate(changeExpData)
    deallocate(changeExpMat)
    deallocate(newExpInt)
    deallocate(newDists)
    deallocate(changedTriplets)

  end subroutine deallocateArrays


  subroutine finalTextOutput(N_move)
    integer, intent(in) :: N_move

    if (textOutput) then

      print *, "The time taken to do", N_move, "moves was", moveTime, &
               "seconds"
      print *, "The total time for the program to run was", totTime, &
               "seconds"
      print *, "The time for the exp calc was", expTime
      print *, "The time for the exp scatter/gather was", gatherExpTime
      print *, "The time for the sum was", sumTime
      print *, "The time for the sum scatter/gather was", gatherTripTime
      print *, "The time for the set-up was", setTime
      print *, "The time for the Xdg set-up was", xTime
      print *, "The time for the exp extraction was", extractTime
      print *, "The time for the triplet set up was ", tripTime
      print *, ' '
      print *, 'Non-additive calculation for atom move complete'
      print *, '========================'
      print *, ' '
      print *, ' '

    end if

  end subroutine finalTextOutput


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
        X(move,i) = X(i,move)

      else if (i .gt. move) then

        X(move,i) = (positions(i,1) - changedPosition(1))**2 + &
                    (positions(i,2) - changedPosition(2))**2 + &
                    (positions(i,3) - changedPosition(3))**2
        X(move,i) = (X(move,i))**0.5
        X(i,move) = X(move,i)

      end if
    end do

    return
  end subroutine updateXdg
    
  
end module tmpi_calcAtomMoveEnergy_mod
