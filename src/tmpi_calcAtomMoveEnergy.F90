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
  public tmpi_calcAtomMoveEnergy, getChangedTriplets, extractChangedExps, getTriPerAtom, &
  findChangedTriIndex


contains


  function tmpi_calcAtomMoveEnergy(N_move,move,proposedPosition,currentEnergyData) result(proposedEnergies)
    ! Inputs
    integer, intent(in) :: N_move, move
    type (energiesData), intent(in) :: currentEnergyData
    type (positionData), intent(in) :: proposedPosition

    ! Output
    type (energiesData) :: proposedEnergies

    ! Local variables
    integer :: triPerProc, i, j, indj, nPerProc, triPerAt, nExpMax, nExpRe
    integer :: nTriMax, nTriRe
    double precision :: totTime, moveTime
    integer, allocatable :: newExpInt(:,:), changedTriplets(:,:), scounts(:)
    integer, allocatable :: scatterTrip(:,:), displs(:), tripIndex(:)
    double precision, allocatable :: newDists(:), newUvec(:), scatterDists(:)
    double precision, allocatable :: changeExpData(:,:,:), changeExpMat(:,:,:)
    double precision, allocatable :: newUfull(:)


    root = 0
    if (processRank .eq. root) then

       if (textOutput) then
         print *, ' '
         print *, ' '
         print *, '========================'
         print *, 'Beginning non-additive calculation for atom move'
         print *, ' '
      end if

    end if


    totTime = MPI_Wtime()
    call getTriPerAtom(proposedPosition%N_a, triPerAt)
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


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()
    do i = 1, N_move

       proposedEnergies%expMatrix = currentEnergyData%expMatrix

       ! Set up on root
       if (processRank .eq. root) then

          if (textOutput) then
            print *, '------------------------'
            print *, "Moving atom", move
            print *, "                 "
          end if

          ! Re-calculate interatomicDistances for the new atomic positions
          call makeXdgNonAdd(proposedPosition%N_a,proposedPosition%posArray, proposedEnergies%interatomicDistances)

          ! Find the indices of the affected exponentials
          call extractChangedExps(proposedPosition%N_a,move,proposedEnergies%interatomicDistances, &
                                  newExpInt,newDists)

          ! Determine which triplets have undergone a change
          call getChangedTriplets(move,proposedPosition%N_a,triPerAt, changedTriplets)
          call findChangedTriIndex(triPerAt,proposedPosition%N_a,move, tripIndex)

       end if

       ! Scatter all requisite data from move set-up on root to all procs
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_Bcast(newExpInt, 2*(proposedPosition%N_a-1), MPI_INT, root, MPI_COMM_WORLD, &
                      ierror)
       call MPI_Bcast(proposedEnergies%interatomicDistances, proposedPosition%N_a*proposedPosition%N_a, &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(proposedPosition%posArray, 3*proposedPosition%N_a, MPI_DOUBLE_PRECISION, root, &
                      MPI_COMM_WORLD, ierror)
       call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)

       ! Determine no. of distances to scatter to each process for exp re-calc
       call getNPerProcNonAdd(proposedPosition%N_a-1,clusterSize, nExpMax,nExpRe)
       call getVarrays(clusterSize,nExpMax,nExpRe, scounts,displs)
       nPerProc = scounts(processRank+1)

       ! Allocate arrays on first loop
       if (i .eq. 1) then

         allocate(scatterDists(nPerProc))
         allocate(changeExpData(nArgs,N_tp,nPerProc))

       end if

       ! Scatter distances across processes
       call MPI_scatterv(newDists, scounts, displs, MPI_DOUBLE_PRECISION, &
                         scatterDists, nPerProc, MPI_DOUBLE_PRECISION, &
                         root, MPI_COMM_WORLD, ierror)

       ! Update exponentials in changed triplets
       call calculateExponentialsNonAdd(nPerProc,N_tp,nArgs,trainData, &
                                        hyperParams(1),scatterDists, &
                                        changeExpData)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Gather in all updated exps and broadcast the resultant matrix to all procs
       call MPI_gatherv(changeExpData, N_tp*nArgs*nPerProc, MPI_DOUBLE_PRECISION, &
                        changeExpMat, N_tp*nArgs*scounts, N_tp*nArgs*displs, &
                        MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(changeExpMat, N_tp*nArgs*(proposedPosition%N_a-1), &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Update the exp matrix
       do j = 1, proposedPosition%N_a-1

          indj = currentEnergyData%distancesIntMat(newExpInt(1,j),newExpInt(2,j))
          proposedEnergies%expMatrix(1:nArgs,1:N_tp,indj) = changeExpMat(1:nArgs,1:N_tp,j)

       end do

       ! Determine how many triplets to send to each proc
       call getNPerProcNonAdd(triPerAt,clusterSize, nTriMax,nTriRe)
       call getVarrays(clusterSize,nTriMax,nTriRe, scounts,displs)
       triPerProc = scounts(processRank+1)

       ! Allocate requisite arrasy in first loop
       if (i .eq. 1) then

         allocate(scatterTrip(3,triPerProc))
         allocate(newUvec(triPerProc))

       end if

       ! Scatter triplets across processes
       call MPI_scatterv(changedTriplets, scounts*3, displs*3, MPI_INT, scatterTrip, &
                         triPerProc*3, MPI_INT, root, MPI_COMM_WORLD, ierror)

       ! Calculate the non-additive energies for the changed triplets and gather
       call tripletEnergiesNonAdd(scatterTrip,currentEnergyData%distancesIntMat,triPerProc,N_tp, &
                                  proposedPosition%N_a,N_p,nArgs,Perm,proposedPosition%N_distances, &
                                  proposedEnergies%expMatrix,alpha,hyperParams(2), newUvec)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, scounts, &
                        displs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Find total change in non-add energy from moving atom
       if (processRank .eq. root) then

          proposedEnergies%tripletEnergies = currentEnergyData%tripletEnergies
          do j = 1, triPerAt

            proposedEnergies%tripletEnergies(tripIndex(j)) = newUfull(j)

          end do

          call totalEnergyNonAdd(proposedEnergies%tripletEnergies,proposedPosition%N_tri, &
                                 proposedEnergies%Utotal)

          if (textOutput) then          
            print *, "The non-additive energy after the move is", proposedEnergies%Utotal
            print *, '------------------------'
            print *, ' '
          end if

       end if

       ! Update data that hasn't yet been
       proposedEnergies%distancesIntMat = currentEnergyData%distancesIntMat

    end do
    moveTime = MPI_Wtime() - moveTime


    ! Deallocate all arrays
    deallocate(scatterDists)
    deallocate(scatterTrip)
    deallocate(newUvec)
    deallocate(changeExpData)
    deallocate(changeExpMat)
    deallocate(newExpInt)
    deallocate(newDists)
    deallocate(changedTriplets)


    ! Finalise MPI and print times taken for each step of calculation
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       if (textOutput) then
         print *, "The time taken to do", N_move, "moves was", moveTime, "seconds"
         print *, "The total time for the program to run was", totTime, "seconds"
         print *, ' '
         print *, 'Non-additive calculation for atom move complete'
         print *, '========================'
         print *, ' '
         print *, ' '
       end if

    end if

    return
  end function tmpi_calcAtomMoveEnergy


  subroutine initialAsserts(N_a)
    ! Input variables
    integer, intent(in) :: N_a


    if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    endif

  end subroutine initialAsserts


  subroutine finalAsserts(N_a)
    ! Input variables
    integer, intent(in) :: N_a

     if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    endif

  end subroutine finalAsserts


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
    
  
end module tmpi_calcAtomMoveEnergy_mod
