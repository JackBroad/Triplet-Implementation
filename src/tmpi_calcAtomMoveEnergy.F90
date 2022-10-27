module tmpi_calcAtomMoveEnergy_mod
  use, intrinsic :: ISO_C_BINDING, only : C_F_POINTER
  use mpi
  use mpi_variables
  use expShare_variables
  use pbcAndMic_variables
  use dataStructure_variables
  use time_variables
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


  double precision :: partialDeltaU
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


    ! Final asserts
    !call finalAsserts(proposedPositionData%N_a)

  return
  end function tmpi_calcAtomMoveEnergy


  subroutine setUpCalculation(mover)
    implicit none
    integer, intent(in) :: mover

    ! Initial set-up and asserts
    call initialAsserts(proposedPositionData%N_a,proposedPositionData%N_tri, &
                        proposedPositionData%N_distances)
    proposedPositionData%N_changed_triplets = getTriPerAtom(proposedPositionData%N_a)
    call allocateArrays()

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
    call distributeDistsSharedMem()
    extractTime = MPI_Wtime() - extractTime

    ! Calculate new values of changed exps
    if (allocated(changeExpData)) then
      deallocate(changeExpData)
    end if
    allocate(changeExpData(N_tp,nArgs,N_changed_exp_per_host))
    if (allocated(oldExpData)) then
      deallocate(oldExpData)
    end if
    allocate(oldExpData(N_tp,nArgs,N_changed_exp_per_host))
    call calcChangedExposNonAdd(N_changed_exp_per_host,N_tp,nArgs,trainData, &
                                hyperParams(1),hostDists, changeExpData)
    call MPI_WIN_FENCE(0,win,ierror)

    ! Update exp array separately on each host using the exps they calculated
    call saveOldExponentials()
    call MPI_WIN_FENCE(0,win,ierror)
    call updateExpMatrix(changeExpData,hostIndices,N_changed_exp_per_host)
    call MPI_WIN_FENCE(0,win,ierror)

  return
  end subroutine changedExponentialCalculation


  subroutine changedTripletEnergyCalc(mover)
    implicit none
    integer :: N_dists_per_proc, N_tri_per_proc
    integer, intent(in) :: mover

    ! Find the energies of the changed trips and the change in U over the trips
    ! on each proc
    tripTime = MPI_Wtime()
    N_dists_per_proc = size(proposedEnergyData%processDists)
    N_tri_per_proc = size(proposedEnergyData%tripletEnergies)
    call getTripletEnergiesAtomMove(mover,N_dists_per_proc,N_tri_per_proc, &
                                    proposedEnergyData%tripletEnergies,&
                                    partialDeltaU)
    tripTime = MPI_Wtime() - tripTime - partialSumTime ! d (partial sum now done
                                                       ! in triplet energy calc)

    ! Gather the change in triplet energies on the root process for summation
    gatherTime = MPI_Wtime()
    call MPI_gather(partialDeltaU, 1, MPI_DOUBLE_PRECISION, newUfull, 1, &
                    MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    gatherTime = MPI_Wtime() - gatherTime ! f = gatherTime

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


  subroutine allocateArrays()
 
    allocate(newUfull(clusterSize))
    allocate(newDists(proposedPositionData%N_a-1))
    allocate(newExpInt(proposedPositionData%N_a-1,2))

  return
  end subroutine allocateArrays


  subroutine deallocateArrays()

    deallocate(newUfull)
    deallocate(newDists)
    deallocate(newExpInt)

  return
  end subroutine deallocateArrays


subroutine distributeDistsSharedMem()
  implicit none
  integer :: N_dists, N_even, N_spare, start, stopp
  integer :: countVec(sharedSize), dispVec(sharedSize)
  integer :: disp

  ! Get the number of changed distances and build an array to hold the number
  ! assigned to each host
  N_dists = proposedPositionData%N_a - 1
  call getNPerProcNonAdd(N_dists,sharedSize, N_even,N_spare)
  call getVarrays(sharedSize,N_even,N_spare, countVec,dispVec)

  ! Have each host read the number of dists it is assigned and allocate vec to
  ! hold them
  N_changed_exp_per_host = countVec(hostRank+1)
  disp = dispVec(hostRank+1)
  allocate(hostDists(N_changed_exp_per_host))
  allocate(hostIndices(N_changed_exp_per_host,2))

  ! Fill the host arrays from the full arrays of interatomic distances and
  ! indices
  start = disp+1
  stopp = start+N_changed_exp_per_host-1
  hostDists = expUpdateNoRepeat(start:stopp)
  hostIndices = expUpdateIndNoRepeat(start:stopp,:)

  call MPI_BARRIER(MPI_COMM_WORLD, barError)

return
end subroutine distributeDistsSharedMem


end module tmpi_calcAtomMoveEnergy_mod
