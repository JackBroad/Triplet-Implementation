module tmpi_calcFullSimBoxEnergy_mod
  use, intrinsic :: ISO_C_BINDING, only : C_F_POINTER
  use mpi
  use mpi_variables
  use dataStructure_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use global_Flags, only: textOutput
  use assert_module
  use fullBoxModule
  implicit none


  private
  public  tmpi_calcFullSimBoxEnergy,makeDisIntMatNonAdd,makeUDdgNonAdd, &
          makeTripletMatrix, makeXdg


  integer :: N_dists_per_proc, N_tri_per_proc
  double precision :: totTime, setUpTime, expTime
  double precision :: shareTime, rootTime, dataTime
  double precision :: sumTime, partTime, partialU
  integer, allocatable :: scounts(:), displs(:)
  double precision, allocatable :: UD_dg(:)
  double precision, allocatable :: allPartEnergies(:)


contains


  double precision function tmpi_calcFullSimBoxEnergy()
    implicit none

    call initialAsserts(currentPositionData%N_a,currentPositionData%N_tri, &
                        currentPositionData%N_distances)

    ! Set up
    totTime = MPI_Wtime()
    setUpTime = MPI_Wtime()
    call setUpFullBoxCalculation()
    setUpTime = MPI_Wtime() - setUpTime

    ! Calculate exponentials and assign distances to each proc
    expTime = MPI_Wtime()
    call setUpExpCalculation()
    call calculateExponentialsNonAdd(currentPositionData%N_distances,N_tp,nArgs,trainData,&
                                     hyperParams(1),UD_dg, expArray)
    expTime = MPI_Wtime() - expTime

    call findTripletEnergies()

    call getPartialSums()

    call getTotalBoxEnergy()
    tmpi_calcFullSimBoxEnergy = currentEnergyData%Utotal

    call deallocateLocalArrays()

    dataTime = MPI_Wtime()
    call instantiateProposedDataStructs()
    dataTime = MPI_Wtime() - dataTime

    totTime = MPI_Wtime() - totTime

    if (processRank .eq. root) then
       call finalTextOutput() ! Prints times taken by each segment of code
    end if

    call finalAsserts(currentPositionData%N_a)
    
    return
  end function tmpi_calcFullSimBoxEnergy


  subroutine initialAsserts(N_a,N_tri,N_distances)
    integer, intent(in) :: N_a, N_tri, N_distances
    integer :: N_tri_ex

    if ( processRank == root ) then

       N_tri_ex = N_a**3
       N_tri_ex = N_tri_ex - 3*N_a**2
       N_tri_ex = N_tri_ex + 2*N_a
       N_tri_ex = N_tri_ex / 6

       call assertTrue(N_a > 0, 'tmpi_calcFullSimBoxEnergy: should have N_a > 0')
       call assertTrue(N_tri .eq. N_tri_ex, &
       'tmpi_calcFullSimBoxEnergy: should have N_tri = (N_a^3 - 3N_a^2 + 2N_a) / 6')
       call assertTrue(N_distances .eq. ((N_a * N_a) - N_a) / 2, &
       'tmpi_calcFullSimBoxEnergy: should have N_distances = (N_a^2 - N_a) / 2')

    end if

  end subroutine initialAsserts


  subroutine finalAsserts(N_a)
    integer, intent(in) :: N_a

    if ( processRank == root ) then

       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy: should have N_a > 0')

    end if

  end subroutine finalAsserts


  subroutine setUpExpCalculation()
    implicit none

    N_dists_per_proc = getNdistsPerProcFullBox() ! Get no. of dists on each proc
    allocate(currentEnergyData%processDists(N_dists_per_proc))
    !if (allocated(expArray)) then
    !  deallocate(expArray)
    !end if
    !allocate(expArray(N_tp,nArgs,currentPositionData%N_distances))
    currentEnergyData%processDists = distributeDistances(N_dists_per_proc,UD_dg)

    ! Set up shared memory window for exp calc
    shapeArray = (/ N_tp,nArgs,currentPositionData%N_distances /) 
    if (hostrank .eq. 0) then
      windowsize = int(1*10**9,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
    else
      windowsize = 0_MPI_ADDRESS_KIND
    end if
    disp_unit = 1
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, hostComm, &
                                 baseptr, win, ierror)
    if (hostrank .ne. 0) then
       CALL MPI_WIN_SHARED_QUERY(win, 0, windowsize, disp_unit, baseptr, ierror)
    end if
    CALL C_F_POINTER(baseptr, expArray, shapeArray)

    return
  end subroutine setUpExpCalculation


  subroutine setUpFullBoxCalculation()
    implicit none

    root = 0
    N_p = 6
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    currentEnergyData = setupCurrentEnergyData(currentPositionData)
    call broadcastRootData()

    return
  end subroutine setUpFullBoxCalculation


  function setupCurrentEnergyData(currentPosition) result(currentEnergy)
    ! Input variables
    type (positionData), intent(in) :: currentPosition

    ! Output variables
    type (energiesData) :: currentEnergy

    !allocate(currentEnergy%tripletEnergies(currentPosition%N_tri))
    allocate(currentEnergy%interatomicDistances(currentPosition%N_a, &
             currentPosition%N_a))

    ! Set up the arrays required for the non-additive calculation
    call makeXdg(currentPosition%N_a,currentPosition%posArray, &
                 currentEnergy%interatomicDistances)
    call makeDisIntMatNonAdd(currentPosition%N_a, currentEnergy%distancesIntMat)
    call makeUDdgNonAdd(currentPosition%N_a,currentPosition%N_distances, &
                        currentEnergy%interatomicDistances, UD_dg)

    ! Set up array of all possible triplets
    allocate(currentEnergy%triMat(3,currentPosition%N_tri))
    currentEnergy%triMat = makeTripletMatrix(currentPosition%N_a,currentPosition%N_tri)

  end function setupCurrentEnergyData


  subroutine broadcastRootData()

    call MPI_Bcast(hyperParams, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                   ierror)
    call MPI_Bcast(N_tp, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nArgs, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(Perm, 18, MPI_INT, root, MPI_COMM_WORLD, ierror)

  end subroutine broadcastRootData


  subroutine deallocateLocalArrays()

    deallocate(scounts)
    deallocate(displs)
    deallocate(UD_dg)
    deallocate(allPartEnergies)

  end subroutine deallocateLocalArrays


  subroutine finalTextOutput()

    if (textOutput) then

       !print *, "The time taken for the exponentials was", expTime, "seconds"
       !print *, "The time taken for the sum was", sumTime, "seconds"
       !print *, "The time taken to set up was", setUpTime, "seconds"
       !print *, "The total time for the program to run was", totTime, "seconds"
       !print *, ' '
       !print *, 'Non-additive calculation for full sim box complete'
       !print *, '========================'
       !print *, ' '
       !print *, ' '
       print *, totTime, setUpTime, expTime, sumTime, partTime, shareTime, rootTime, dataTime

    end if

  end subroutine finalTextOutput


  subroutine getPartialSums()
    implicit none

    partTime = MPI_Wtime()
    partialU = sum(currentEnergyData%tripletEnergies)
    partTime = MPI_Wtime() - partTime

  return
  end subroutine getPartialSums


  subroutine getTotalBoxEnergy()
    implicit none

    shareTime = MPI_Wtime()
    allocate(allPartEnergies(clusterSize))
    call MPI_gather(partialU, 1, MPI_DOUBLE_PRECISION, allPartEnergies, 1, &
                    MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    shareTime = MPI_Wtime() - shareTime

    rootTime = MPI_Wtime()
    if (processRank .eq. root) then
      currentEnergyData%Utotal = sum(allPartEnergies)
    end if
    rootTime = MPI_Wtime() - rootTime

  return
  end subroutine getTotalBoxEnergy


  subroutine findTripletEnergies()
    implicit none

    ! Determine no. of triplets to send to each process for triplet calc.
    ! then find triplet energies
    sumTime = MPI_Wtime()
    allocate(currentEnergyData%alphaBetaPairs(N_dists_per_proc,2))
    currentEnergyData%alphaBetaPairs = getAlphaBetaPairs(N_dists_per_proc, &
                                                         currentEnergyData%interatomicDistances)
    N_tri_per_proc = getNtripsPerProcFullBox(N_dists_per_proc)
    allocate(currentEnergyData%tripletEnergies(N_tri_per_proc))
    currentEnergyData%tripletEnergies = getTripletEnergiesFullBox(N_dists_per_proc, &
                                                                  N_tri_per_proc)
    sumTime = MPI_Wtime() - sumTime

  return
  end subroutine findTripletEnergies

  
end module tmpi_calcFullSimBoxEnergy_mod
