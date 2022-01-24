module tmpi_calcAtomMoveEnergy_mod
  use mpi
  use mpi_variables
  use expShare_variables
  use dataStructure_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use global_Flags, only: textOutput
  use assert_module
  use atomMoveModule
  implicit none


  private
  public tmpi_calcAtomMoveEnergy


  double precision :: moveTime, expTime, sumTime, setTime
  double precision :: gatherTime, xTime, tripTime, partialSumTime
  double precision :: extractTime, tripSumTime, newTripU
  double precision :: oldTripU, partialDeltaU, rootSumTime
  integer, allocatable :: newExpInt(:,:)
  double precision, allocatable :: newDists(:), newUfull(:)


contains


  double precision function tmpi_calcAtomMoveEnergy(atomToMove)
    ! Inputs
    integer, intent(in) :: atomToMove


    ! Set up calculation
    moveTime = MPI_Wtime()
    setTime = MPI_Wtime()
    call setUpCalculation(atomToMove)
    setTime = MPI_Wtime() - setTime ! b


    ! Re-calculate changed exponentials
    expTime = MPI_Wtime()
    call changedExponentialCalculation(atomToMove)
    expTime = MPI_Wtime() - expTime ! c


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
    !print *, processRank, moveTime, size(proposedEnergyData%tripletEnergies), &
    !         size(changedTriInd)
    call finalAsserts(proposedPositionData%N_a)


  return
  end function tmpi_calcAtomMoveEnergy


  subroutine setUpCalculation(mover)
    implicit none
    integer, intent(in) :: mover

    ! Initial set-up and asserts
    call initialAsserts(proposedPositionData%N_a,proposedPositionData%N_tri, &
                        proposedPositionData%N_distances)
    proposedPositionData%N_changed_triplets = getTriPerAtom(proposedPositionData%N_a)
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
    allocate(changeExpData(N_tp,nArgs,proposedPositionData%N_a-1))
    if (allocated(oldExpData)) then
      deallocate(oldExpData)
    end if
    allocate(oldExpData(N_tp,nArgs,proposedPositionData%N_a-1))
    call calculateExponentialsNonAdd(proposedPositionData%N_a-1,N_tp,nArgs,trainData, &
                                     hyperParams(1),newDists, changeExpData)

    ! Update exp array
    call saveOldExponentials()
    call updateExpMatrix(changeExpData,newExpInt,proposedPositionData%N_a-1)

  return
  end subroutine changedExponentialCalculation


  subroutine changedTripletEnergyCalc(mover)
    implicit none
    integer :: N_dists_per_proc, N_tri_per_proc
    integer, intent(in) :: mover

    ! Determine which triplets have changed
    tripTime = MPI_Wtime()
    N_dists_per_proc = size(proposedEnergyData%processDists)
    N_tri_per_proc = size(proposedEnergyData%tripletEnergies)
    call getTripletEnergiesAtomMove(mover,N_dists_per_proc,N_tri_per_proc, &
                                    proposedEnergyData%tripletEnergies)
    tripTime = MPI_Wtime() - tripTime ! d

    ! Calculate the non-additive energies for the changed triplets
    tripSumTime = MPI_Wtime()
    tripSumTime = MPI_Wtime() - tripSumTime ! e

    ! Find the change in U on each process
    partialSumTime = MPI_Wtime()
    newTripU = sum(proposedEnergyData%tripletEnergies)
    oldTripU = sum(currentEnergyData%tripletEnergies)
    partialDeltaU = newTripU - oldTripU
    partialSumTime = MPI_Wtime() - partialSumTime

    ! Gather the change in triplet energies on the root process for summation
    gatherTime = MPI_Wtime()
    call MPI_gather(partialDeltaU, 1, MPI_DOUBLE_PRECISION, newUfull, 1, &
                    MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    gatherTime = MPI_Wtime() - gatherTime ! f = gatherTime+partialSumTime

    ! Find new total non-additive energy after moving an atom
    rootSumTime = MPI_Wtime()
    proposedEnergyData%Utotal = 0d0
    if (processRank .eq. root) then
      proposedEnergyData%Utotal = sum(newUfull)
    end if
    rootSumTime = MPI_Wtime() - rootSumTime ! g

  return
  end subroutine changedTripletEnergyCalc


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
 
    allocate(newUfull(clusterSize))
    allocate(newDists(proposedPosition%N_a-1))
    allocate(newExpInt(proposedPosition%N_a-1,2))
    if (allocated(proposedEnergy%interatomicDistances)) then
      deallocate(proposedEnergy%interatomicDistances)
    end if
    allocate(proposedEnergy%interatomicDistances(proposedPosition%N_a, &
             proposedPosition%N_a))

  end subroutine allocateArrays


  subroutine deallocateArrays()

    deallocate(newUfull)
    deallocate(newDists)
    deallocate(newExpInt)

  end subroutine deallocateArrays


  subroutine finalTextOutput()

    if (textOutput) then
      print *, moveTime, setTime, expTime, tripTime, tripSumTime, &
               partialSumTime, gatherTime, rootSumTime
    end if

  end subroutine finalTextOutput


end module tmpi_calcAtomMoveEnergy_mod
