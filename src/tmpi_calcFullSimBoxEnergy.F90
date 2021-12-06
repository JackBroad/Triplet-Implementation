module tmpi_calcFullSimBoxEnergy_mod
  use mpi_variables
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
  public  tmpi_calcFullSimBoxEnergy,makeDisIntMatNonAdd,makeUDdgNonAdd, &
          makeTripletMatrix, makeXdg


  integer :: nSum, maxTriPerProc, remTriPerProc
  integer :: N_dists_per_proc
  double precision :: totTime, setUpTime, expTime
  double precision :: shareTime, rootTime, dataTime
  double precision :: sumTime
  integer, allocatable :: scounts(:), displs(:)
  integer, allocatable :: triScatter(:,:)
  double precision, allocatable :: UD_dg(:), uVec(:)


contains


  double precision function tmpi_calcFullSimBoxEnergy()
    implicit none


    call initialAsserts(currentPositionData%N_a,currentPositionData%N_tri, &
                        currentPositionData%N_distances)
    call declareConstantsAndRowsOfPermutationMatrix()


    ! Set up
    totTime = MPI_Wtime()
    setUpTime = MPI_Wtime()
    call setUpFullBoxCalculation()
    setUpTime = MPI_Wtime() - setUpTime


    ! Calculate exponentials and assign distances to each proc
    expTime = MPI_Wtime()
    call setUpExpCalculation()
    call calculateExponentialsNonAdd(currentPositionData%N_distances,N_tp,nArgs,trainData,&
                                     hyperParams(1),UD_dg, currentEnergyData%expMatrix)
    expTime = MPI_Wtime() - expTime

    ! Determine no. of triplets to send to each process for triplet calc.
    ! then find triplet energies
    sumTime = MPI_Wtime()
    call setUpTripletSum()
    call tripletEnergiesNonAdd(triScatter,currentEnergyData%distancesIntMat,nSum,N_tp,currentPositionData%N_a, &
                               N_p,nArgs,Perm,currentPositionData%N_distances,currentEnergyData%expMatrix,alpha, &
                               hyperParams(2), uVec)
    sumTime = MPI_Wtime() - sumTime


    ! Gather in the triplet energies and broadcast triplet energies
    shareTime = MPI_Wtime()
    call MPI_gatherv(uVec, nSum, MPI_DOUBLE_PRECISION, currentEnergyData%tripletEnergies, scounts, displs, &
                     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(currentEnergyData%tripletEnergies,currentPositionData%N_tri,MPI_DOUBLE_PRECISION, &
                   root, MPI_COMM_WORLD, ierror)
    shareTime = MPI_Wtime() - shareTime


    ! Find the total non-additive energy for the system by summing triplet
    ! energies
    rootTime = MPI_Wtime()
    call sumEnergiesNonAdd()
    tmpi_calcFullSimBoxEnergy = currentEnergyData%Utotal
    rootTime = MPI_Wtime() - rootTime
    call deallocateLocalArrays()


    dataTime = MPI_Wtime()
    call instantiateProposedDataStructs()
    dataTime = MPI_Wtime() - dataTime


    ! Print times taken for each part of subroutine to run
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then
       !call finalTextOutput()
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


  subroutine declareConstantsAndRowsOfPermutationMatrix()

    totTime = MPI_Wtime()
    root = 0
    N_p = 6
    setUpTime = MPI_Wtime()
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))

  end subroutine declareConstantsAndRowsOfPermutationMatrix
  

  subroutine finalAsserts(N_a)
    integer, intent(in) :: N_a

    if ( processRank == root ) then

       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy: should have N_a > 0')

    end if

  end subroutine finalAsserts


  subroutine setUpExpCalculation()
    implicit none

    N_dists_per_proc = getNdistsPerProc() ! Get no. of dists on each proc
    allocate(currentEnergyData%processDists(N_dists_per_proc))
    allocate(currentEnergyData%expMatrix(N_tp,nArgs,currentPositionData%N_distances))
    currentEnergyData%processDists = distributeDistances(N_dists_per_proc,UD_dg)

    return
  end subroutine setUpExpCalculation


  subroutine setUpFullBoxCalculation()
    implicit none

    !call declareConstantsAndRowsOfPermutationMatrix()
    currentEnergyData = setupCurrentEnergyDataAndArrays(currentPositionData)
    call broadcastRootData()

    return
  end subroutine setUpFullBoxCalculation


  subroutine setUpTripletSum()
    implicit none

    call getNPerProcNonAdd(currentPositionData%N_tri,clusterSize, maxTriPerProc,remTriPerProc)
    call getVarrays(clusterSize,maxTriPerProc,remTriPerProc, scounts,displs)
    call getnSumAndTripletArrays()


    ! Scatter the triplet matrix
    triScatter = currentEnergyData%triMat(1:3,1+displs(processRank+1):displs(processRank+1)+&
                                          scounts(processRank+1))

    return
  end subroutine setUpTripletSum


  subroutine sumEnergiesNonAdd()
    implicit none

    currentEnergyData%Utotal = 0d0
    if (processRank .eq. root) then
      call totalEnergyNonAdd(currentEnergyData%tripletEnergies,currentPositionData%N_tri, &
                             currentEnergyData%Utotal)
    end if

    return
  end subroutine sumEnergiesNonAdd


  subroutine instantiateProposedDataStructs()
    implicit none

    proposedPositionData = currentPositionData
    proposedEnergyData = currentEnergyData

  end subroutine instantiateProposedDataStructs


  function setupCurrentEnergyDataAndArrays(currentPosition) result(currentEnergy)
    ! Input variables
    type (positionData), intent(in) :: currentPosition
       
    ! Output variables
    type (energiesData) :: currentEnergy

    ! Read in all necessary info from files
    allocate(currentEnergy%tripletEnergies(currentPosition%N_tri))
    allocate(currentEnergy%interatomicDistances(currentPosition%N_a, &
             currentPosition%N_a))

    ! Set up the arrays required for the non-additive calculation
    call makeXdg(currentPosition%N_a,currentPosition%posArray, &
                 currentEnergy%interatomicDistances)
    call makeDisIntMatNonAdd(currentPosition%N_a, currentEnergy%distancesIntMat)
    call makeUDdgNonAdd(currentPosition%N_a,currentPosition%N_distances, &
                        currentEnergy%interatomicDistances, UD_dg)

    ! Set up array of a all possible triplets
    allocate(currentEnergy%triMat(3,currentPosition%N_tri))
    currentEnergy%triMat = makeTripletMatrix(currentPosition%N_a,currentPosition%N_tri)

  end function setupCurrentEnergyDataAndArrays


  subroutine broadcastRootData()

    call MPI_Bcast(hyperParams, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                   ierror)
    call MPI_Bcast(N_tp, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nArgs, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(Perm, 18, MPI_INT, root, MPI_COMM_WORLD, ierror)

  end subroutine broadcastRootData



  subroutine getnSumAndTripletArrays()

    nSum = scounts(processRank+1)
    allocate(triScatter(3,nSum))
    allocate(uVec(nSum))

  end subroutine getnSumAndTripletArrays


  subroutine deallocateLocalArrays()

    deallocate(uVec)
    deallocate(triScatter)
    deallocate(scounts)
    deallocate(displs)
    deallocate(UD_dg)

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
       print *, totTime, setUpTime, expTime, sumTime, shareTime, rootTime, dataTime, 0d0

    end if

  end subroutine finalTextOutput


  ! Makes index equivalent to X_dg
  subroutine makeDisIntMatNonAdd(nAt, disIntMat)
    implicit none
    integer, intent(in) :: nAt
    integer, allocatable, intent(out) :: disIntMat(:,:)
    integer :: ind, indi, indj

    ! Do the same for the integer equivalent of X_dg
    allocate(disIntMat(nAt,nAt))
    ind = 0
    do indi = 1, nAt
      do indj = 1, nAt
        if (indi .eq. indj) then
                  disIntMat(indi,indj) = 0
        else if (indi .lt. indj) then
          ind = ind + 1
          disIntMat(indi,indj) = ind
        else
          disIntMat(indi,indj) = disIntMat(indj,indi)
        end if
      end do
    end do

  return
  end subroutine makeDisIntMatNonAdd


  ! Extracts UD of X_dg
  subroutine makeUDdgNonAdd(nAt,udSize,X_dg, UD_dg)
    implicit none
    integer, intent(in) :: nAt, udSize
    double precision, intent(in) :: X_dg(nAt,nAt)
    double precision, allocatable, intent(out) :: UD_dg(:)
    integer :: nMinus, ele, m, n

    allocate(UD_dg(udSize))
    nMinus = nAt - 1
    ele = 0
    do m = 1, nMinus
      do n = m+1, nAt
        ele = ele + 1
        UD_dg(ele) = X_dg(m,n)
      end do
    end do

  return
  end subroutine makeUDdgNonAdd


  ! Uses the number of atoms (nAt) and number of triplets (nTri) to build a matrix
  ! of all possible triplets
  function makeTripletMatrix(nAt,nTri) result(tripletMatrix)
    implicit none
    integer, intent(in) :: nAt, nTri
    integer :: tripletMatrix(3,nTri)
    integer :: al, be, ga, counter

    counter = 0

    do al = 1, nAt-2
      do be = al+1, nAt-1
        do ga = be+1, nAt

          counter = counter+1
          tripletMatrix(1,counter) = al
          tripletMatrix(2,counter) = be
          tripletMatrix(3,counter) = ga

        end do
      end do
    end do

  return
  end function makeTripletMatrix


  subroutine makeXdg(nAt,posArray, X_dg)
    implicit none
    integer, intent(in) :: nAt
    double precision, intent(in) :: posArray(nAt,3)
    double precision, intent(out) :: X_dg(nAt,nAt)
    integer :: i, j

    ! Find X_dg for the atomic positions in posArray
    do i = 1, nAt
      do j = 1, nAt
        if (i .eq. j) then

          X_dg(i,j) = 0

        else if (i .lt. j) then

          X_dg(i,j) = (posArray(i,1)-posArray(j,1))**2 + &
                      (posArray(i,2)-posArray(j,2))**2 + &
                      (posArray(i,3)-posArray(j,3))**2
          X_dg(i,j) = (X_dg(i,j))**0.5
          X_dg(i,j) = 1 / X_dg(i,j) ! Convert to inverse distance

        else

          X_dg(i,j) = X_dg(j,i)

        end if
      end do
    end do

  return
  end subroutine makeXdg
    
  
end module tmpi_calcFullSimBoxEnergy_mod
