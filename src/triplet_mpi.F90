module triplet_mpi_mod
  use mpi_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use global_Flags, only: textOutput
  implicit none
  include 'mpif.h'


contains


  function tmpi_calcFullSimBoxEnergy(N_a,N_tri,udSize,posArray) result(currentEnergies)
    ! Input variables
    integer, intent(in) :: N_a, N_tri, udSize
    double precision, intent(in) :: posArray(N_a,3)
   
    ! Output variables
    type( energiesData):: currentEnergies
    
    ! Local variables
    integer :: eCols, nSum, maxnSum, dataSize, maxDataSize
    integer :: reNsum, reDataSize
    double precision :: expTime, sumTime, totTime, setUpTime
    integer, allocatable :: triMat(:,:), triScatter(:,:) 
    integer, allocatable :: scounts(:), displs(:)
    double precision, allocatable :: scatterData(:), UD_dg(:)
    double precision, allocatable :: expData(:,:,:), uVec(:)



    ! Declare constants and rows of permutation matrix
    totTime = MPI_Wtime()
    root = 0
    N_p = 6
    setUpTime = MPI_Wtime()
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))


    ! Set up on root
    if (processRank .eq. root) then

       if (textOutput) then
          print *, ' '
          print *, ' '
          print *, '========================'
          print *, 'Beginning non-additive calculation for whole sim box'
          print *, ' '
       end if

       ! Read in all necessary info from files
       allocate(currentEnergies%tripletEnergies(N_tri))
       allocate(currentEnergies%interatomicDistances(N_a,N_a))

       ! Set up the arrays required for the non-additive calculation
       call makeXdgNonAdd(N_a,posArray, currentEnergies%interatomicDistances)
       call makeDisIntMatNonAdd(N_a, currentEnergies%distancesIntMat)
       call makeUDdgNonAdd(N_a,udSize,currentEnergies%interatomicDistances, UD_dg)

       ! Set up array of a all possible triplets
       allocate(triMat(3,N_tri))
       call makeTripletMatrix(N_a,N_tri, triMat)

       ! Declare permutation matrix
       Perm(1,:) = (/1, 2, 3/)
       Perm(2,:) = (/1, 3, 2/)
       Perm(3,:) = (/2, 1, 3/)
       Perm(4,:) = (/2, 3, 1/)
       Perm(5,:) = (/3, 1, 2/)
       Perm(6,:) = (/3, 2, 1/)

    end if


    ! Hold all processes here until root process has finished setting up
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Broadcast all other requisite data that is already allocated everywhere from root to 
    ! all processes
    call MPI_Bcast(hyperParams, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(N_a, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(N_tri, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(udSize, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(N_tp, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(nArgs, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(Perm, 18, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Use info from last broadcast to allocate arrays on other processes
    if (processRank .ne. root) then

       allocate(currentEnergies%distancesIntMat(N_a,N_a))
       allocate(currentEnergies%interatomicDistances(N_a,N_a))

    end if

    ! Broadcast new arrays from root
    call MPI_Bcast(currentEnergies%distancesIntMat, N_a*N_a, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(currentEnergies%interatomicDistances, N_a*N_a, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Determine max no. of elements of UD_dg to send to each process for exp
    ! calculations
    call getNPerProcNonAdd(udSize,clusterSize, maxDataSize,reDataSize)
    setUpTime = MPI_Wtime() - setUpTime
    expTime = MPI_Wtime()


    ! Determine actual no. of elements to send to each process for exp calc.
    call getVarrays(clusterSize,maxDataSize,reDataSize, scounts,displs)
    dataSize = scounts(processRank+1)
    allocate(scatterData(dataSize)) ! Allocate array to send exponentials


    ! Scatter the interatomic distances in U_dg to all processes
    call MPI_scatterv(UD_dg, scounts, displs, MPI_DOUBLE_PRECISION, scatterData, dataSize, &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)


    ! Calculate the exponentials for each distance on each process
    allocate(expData(nArgs,N_tp,dataSize))
    call calculateExponentialsNonAdd(dataSize,N_tp,nArgs,trainData,hyperParams(1), &
                                     scatterData, expData)


    ! Allocate an array to hold all exps
    eCols = nArgs*N_tri
    allocate(currentEnergies%expMatrix(nArgs,N_tp,udSize))


    ! Gather expData arrays from the other processes and add them to expMatrix on
    ! the root process
    call MPI_gatherv(expData, N_tp*nArgs*dataSize, MPI_DOUBLE_PRECISION, currentEnergies%expMatrix, &
                     N_tp*nArgs*scounts, N_tp*nArgs*displs, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)
    expTime = MPI_Wtime() - expTime


    ! Broadcast expMatrix to all processes so that sum can be parallelised
    sumTime = MPI_Wtime()
    call MPI_Bcast(currentEnergies%expMatrix, N_tp*nArgs*udSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
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
    call tripletEnergiesNonAdd(triScatter,currentEnergies%distancesIntMat,nSum,N_tp,N_a,N_p,nArgs,Perm, &
                               udSize,currentEnergies%expMatrix,alpha,hyperParams(2), uVec)


    ! Gather in the triplet energies and sum them to get total non-add energy
    call MPI_gatherv(uVec, nSum, MPI_DOUBLE_PRECISION, currentEnergies%tripletEnergies, scounts, displs, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)


    ! Find the total non-additive energy for the system
    currentEnergies%Utotal = 0d0
    if (processRank .eq. root) then

       call totalEnergyNonAdd(currentEnergies%tripletEnergies,N_tri, currentEnergies%Utotal)

       if (textOutput ) then
          print *, "The total non-additive energy is", currentEnergies%Utotal
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

    return
  end function tmpi_calcFullSimBoxEnergy


  subroutine tmpi_calcAtomMoveEnergy(N_move,dist,N_a,udSize,N_tri,currentEnergies,posArray, &
                                     proposedEnergies)
    ! Input variables
    integer, intent(in) :: N_a, udSize, N_move, N_tri
    double precision, intent(in) :: dist
    type( energiesData ), intent(in) :: currentEnergies

    ! In/out variables
    double precision, intent(inout) :: posArray(N_a,nArgs)

    ! Output variables
    type( energiesData ), intent(out) :: proposedEnergies

    ! Local variables
    integer :: triPerProc, i, j, indj, nPerProc, triPerAt, move, nExpMax, nExpRe
    integer :: nTriMax, nTriRe
    double precision :: newPosAt(N_a,nArgs), totTime, moveTime
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
    call getTriPerAtom(N_a, triPerAt)
    allocate(newExpInt(2,N_a-1))
    allocate(newDists(N_a-1))
    allocate(changeExpMat(nArgs,N_tp,N_a-1))
    allocate(changedTriplets(3,triPerAt))
    allocate(newUfull(triPerAt))
    allocate(scounts(clusterSize))
    allocate(displs(clusterSize))
    allocate(proposedEnergies%interatomicDistances(N_a,N_a))
    allocate(tripIndex(triPerAt))


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()
    do i = 1, N_move

       proposedEnergies%expMatrix = currentEnergies%expMatrix

       ! Set up on root
       if (processRank .eq. root) then

          ! Move an atom
          call moveAt(posArray,N_a,dist, newPosAt,move)
          !if (textOutput) then
          !  print *, '------------------------'
          !  print *, "Moving atom", move
          !  print *, "                 "
          !end if

          ! Re-calculate interatomicDistances for the new atomic positions
          call makeXdgNonAdd(N_a,newPosAt, proposedEnergies%interatomicDistances)

          ! Find the indices of the affected exponentials
          call extractChangedExps(N_a,move,proposedEnergies%interatomicDistances, &
                                  newExpInt,newDists)

          ! Determine which triplets have undergone a change
          call getChangedTriplets(move,N_a,triPerAt, changedTriplets)
          call findChangedTriIndex(triPerAt,N_a,move, tripIndex)

       end if

       ! Scatter all requisite data from move set-up on root to all procs
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_Bcast(newExpInt, 2*(N_a-1), MPI_INT, root, MPI_COMM_WORLD, &
                      ierror)
       call MPI_Bcast(proposedEnergies%interatomicDistances, N_a*N_a, &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(newPosAt, 3*N_a, MPI_DOUBLE_PRECISION, root, &
                      MPI_COMM_WORLD, ierror)
       call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)

       ! Determine no. of distances to scatter to each process for exp re-calc
       call getNPerProcNonAdd(N_a-1,clusterSize, nExpMax,nExpRe)
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
       call MPI_Bcast(changeExpMat, N_tp*nArgs*(N_a-1), &
                      MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Update the exp matrix
       do j = 1, N_a-1

          indj = currentEnergies%distancesIntMat(newExpInt(1,j),newExpInt(2,j))
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
       call tripletEnergiesNonAdd(scatterTrip,currentEnergies%distancesIntMat,triPerProc, &
                                  N_tp,N_a,N_p,nArgs,Perm,udSize,proposedEnergies%expMatrix, &
                                  alpha,hyperParams(2), newUvec)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, scounts, &
                        displs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Find total change in non-add energy from moving atom
       if (processRank .eq. root) then

          proposedEnergies%tripletEnergies = currentEnergies%tripletEnergies
          do j = 1, triPerAt

            proposedEnergies%tripletEnergies(tripIndex(j)) = newUfull(j)

          end do

          call totalEnergyNonAdd(proposedEnergies%tripletEnergies,N_tri, &
                                 proposedEnergies%Utotal)

          !if (textOutput) then          
          !  print *, "The non-additive energy after the move is", proposedEnergies%Utotal
          !  print *, '------------------------'
          !  print *, ' '
          !end if

       end if

       ! Update data that hasn't yet been
       posArray = newPosAt
       proposedEnergies%distancesIntMat = currentEnergies%distancesIntMat

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
  end subroutine tmpi_calcAtomMoveEnergy
end module triplet_mpi_mod
