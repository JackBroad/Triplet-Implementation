module tmpi_calcFullSimBoxEnergy_mod
  use mpi_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use global_Flags, only: textOutput
  use assert_module
  implicit none
  include 'mpif.h'

  private
  public  tmpi_calcFullSimBoxEnergy

  integer :: dataSize
  double precision :: totTime, setUpTime
  integer, allocatable :: scounts(:), displs(:)
  integer, allocatable :: triMat(:,:)
  double precision, allocatable ::  UD_dg(:), scatterData(:), expData(:,:,:)

contains


  function tmpi_calcFullSimBoxEnergy(N_a,N_tri,N_distances,posArray) result(currentEnergyData)
    ! Input variables
    integer, intent(in) :: N_a, N_tri, N_distances
    double precision, intent(in) :: posArray(N_a,3)
   
    ! Output variables
    type (energiesData) :: currentEnergyData
    
    ! Local variables
    integer :: nSum, maxnSum, maxDataSize
    integer :: reNsum, reDataSize
    double precision :: expTime, sumTime
    integer, allocatable :: triScatter(:,:) 
    double precision, allocatable :: uVec(:)

    call initialAsserts(N_a)
    call declareConstantsAndRowsOfPermutationMatrix()


    ! Set up on root
    if (processRank .eq. root) then

       call initialTextOutput()
       currentEnergyData = setupCurrentEnergyDataAndArrays(N_a,N_tri,N_distances,posArray) 
       Perm = setupPermutationMatrix()

    end if


    ! Hold all processes here until root process has finished setting up
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Broadcasts
    call MPI_Bcast(N_a, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(N_tri, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(N_distances, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call broadcastRootData()
    call broadcastCurrentEnergyData(currentEnergyData,N_a)


    ! Determine max no. of elements of UD_dg to send to each process for exp
    ! calculations
    call getNPerProcNonAdd(N_distances,clusterSize, maxDataSize,reDataSize)
    setUpTime = MPI_Wtime() - setUpTime


    ! Determine actual no. of elements to send to each process for exp calc.
    expTime = MPI_Wtime()
    call getVarrays(clusterSize,maxDataSize,reDataSize, scounts,displs)
    call getDataSizeAndExpArrays()


    ! Scatter the interatomic distances in U_dg to all processes
    call MPI_scatterv(UD_dg, scounts, displs, MPI_DOUBLE_PRECISION, scatterData, dataSize, &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)


    ! Calculate the exponentials for each distance on each process
    call calculateExponentialsNonAdd(dataSize,N_tp,nArgs,trainData,hyperParams(1), &
                                     scatterData, expData)


    ! Allocate an array to hold all exps
    allocate(currentEnergyData%expMatrix(nArgs,N_tp,N_distances))


    ! Gather expData arrays from the other processes and add them to expMatrix on
    ! the root process
    call MPI_gatherv(expData, N_tp*nArgs*dataSize, MPI_DOUBLE_PRECISION, currentEnergyData%expMatrix, &
                     N_tp*nArgs*scounts, N_tp*nArgs*displs, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)
    expTime = MPI_Wtime() - expTime


    ! Broadcast expMatrix to all processes so that sum can be parallelised
    sumTime = MPI_Wtime()
    call MPI_Bcast(currentEnergyData%expMatrix, N_tp*nArgs*N_distances, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                   ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Determine actual no. of elements to send to each process for triplet calc.
    call getNPerProcNonAdd(N_tri,clusterSize, maxnSum,reNsum)
    call getVarrays(clusterSize,maxNsum,reNsum, scounts,displs)
    nSum = scounts(processRank+1)


    ! Scatter the triplet matrix
    allocate(triScatter(3,nSum))
    call MPI_scatterv(triMat, scounts*3, displs*3, MPI_INT, triScatter, nSum*3, MPI_INT, &
                      root, MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Find the energies of the triplets assigned to each process
    allocate(uVec(nSum))
    call tripletEnergiesNonAdd(triScatter,currentEnergyData%distancesIntMat,nSum,N_tp,N_a,N_p,nArgs,Perm, &
                               N_distances,currentEnergyData%expMatrix,alpha,hyperParams(2), uVec)


    ! Gather in the triplet energies and sum them to get total non-add energy
    call MPI_gatherv(uVec, nSum, MPI_DOUBLE_PRECISION, currentEnergyData%tripletEnergies, scounts, displs, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)


    ! Find the total non-additive energy for the system
    currentEnergyData%Utotal = 0d0
    if (processRank .eq. root) then

       call totalEnergyNonAdd(currentEnergyData%tripletEnergies,N_tri, currentEnergyData%Utotal)

       if (textOutput ) then
          print *, "The total non-additive energy is", currentEnergyData%Utotal
          print *, "              "
       end if

    end if
    sumTime = MPI_Wtime() - sumTime

    ! De-allocate arrays not passed to atom-move subroutine
    deallocate(scatterData)
    deallocate(expData)
    deallocate(uVec)
    deallocate(triScatter)
    deallocate(scounts)
    deallocate(displs)
    if (processRank .eq. root) then

       deallocate(UD_dg)
       deallocate(triMat)

    end if


    ! Print times taken for each part of subroutine to run
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       if (textOutput) then
          print *, "The time taken for the exponentials was", expTime, "seconds"
          print *, "The time taken for the sum was", sumTime, "seconds"
          print *, "The time taken to set up was", setUpTime, "seconds"
          print *, "The total time for the program to run was", totTime, "seconds"
          print *, ' '
          print *, 'Non-additive calculation for full sim box complete'
          print *, '========================'
          print *, ' '
          print *, ' '
       end if

    end if

    call finalAsserts(N_a)
    
    return
  end function tmpi_calcFullSimBoxEnergy



  subroutine initialAsserts(N_a)
    ! Input variables
    integer, intent(in) :: N_a


    if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    endif

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
    ! Input variables
    integer, intent(in) :: N_a

    if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    end if

  end subroutine finalAsserts


  subroutine initialTextOutput()
    if (textOutput) then
       print *, ' '
       print *, ' '
       print *, '========================'
       print *, 'Beginning non-additive calculation for whole sim box'
       print *, ' '
    end if
  end subroutine initialTextOutput


  function setupCurrentEnergyDataAndArrays(N_a,N_tri,N_distances,posArray) result(currentEnergyData)
    ! Input variables
    integer, intent(in) :: N_a, N_tri, N_distances
    double precision, intent(in) :: posArray(N_a,3)
       
    ! Output variables
    type( energiesData):: currentEnergyData

    ! Read in all necessary info from files
    allocate(currentEnergyData%tripletEnergies(N_tri))
    allocate(currentEnergyData%interatomicDistances(N_a,N_a))

    ! Set up the arrays required for the non-additive calculation
    call makeXdgNonAdd(N_a,posArray, currentEnergyData%interatomicDistances)
    call makeDisIntMatNonAdd(N_a, currentEnergyData%distancesIntMat)
    call makeUDdgNonAdd(N_a,N_distances,currentEnergyData%interatomicDistances, UD_dg)

    ! Set up array of a all possible triplets
    allocate(triMat(3,N_tri))
    call makeTripletMatrix(N_a,N_tri, triMat)

  end function setupCurrentEnergyDataAndArrays


  function setupPermutationMatrix() result(Perm)
    integer :: Perm(6,3)

    Perm(1,:) = (/1, 2, 3/)
    Perm(2,:) = (/1, 3, 2/)
    Perm(3,:) = (/2, 1, 3/)
    Perm(4,:) = (/2, 3, 1/)
    Perm(5,:) = (/3, 1, 2/)
    Perm(6,:) = (/3, 2, 1/)

  end function setupPermutationMatrix

  subroutine broadcastRootData()

    call MPI_Bcast(hyperParams, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                   ierror)
    call MPI_Bcast(N_tp, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nArgs, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(Perm, 18, MPI_INT, root, MPI_COMM_WORLD, ierror)

  end subroutine broadcastRootData


  subroutine broadcastCurrentEnergyData(currentEnergyData,N_a)
    integer, intent(in) :: N_a
    type (energiesData) :: currentEnergyData

    if (processRank .ne. root) then

       allocate(currentEnergyData%distancesIntMat(N_a,N_a))
       allocate(currentEnergyData%interatomicDistances(N_a,N_a))

    end if

    ! Broadcast new arrays from root
    call MPI_Bcast(currentEnergyData%distancesIntMat, N_a*N_a, MPI_INT, root, &
                   MPI_COMM_WORLD, ierror)
    call MPI_Bcast(currentEnergyData%interatomicDistances, N_a*N_a, &
                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

  end subroutine broadcastCurrentEnergyData


  subroutine getDataSizeAndExpArrays()

    dataSize = scounts(processRank+1)
    allocate(scatterData(dataSize))
    allocate(expData(nArgs,N_tp,dataSize))

  end subroutine getDataSizeAndExpArrays
    
  
end module tmpi_calcFullSimBoxEnergy_mod
