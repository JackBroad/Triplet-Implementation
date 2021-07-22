module triplet_mpi_mod
  use mpi_variables
  use triplet_mod
  use GP_variables, only: hyperParams, alpha, Perm, trainData, N_tp, nArgs,N_p 
  implicit none
  include 'mpif.h'


contains

  subroutine triplet_mpi_fullNonAdd(N_a,N_tri,udSize,posArray, X_dg, &
                                    disIntMat,expMatrix,U,uFull)
    ! Input variables
    integer, intent(in) :: N_a, N_tri, udSize
    double precision, intent(in) :: posArray(N_a,3)
    !currently Na readfrom inputfile and other calculated but this
    !needs to moved to outside this function
    ! Can be calculated from Na so should be local variables (ie not returned)
   
    ! Output variables 
    double precision, allocatable, intent(out) :: X_dg(:,:), uFull(:)
    double precision, allocatable, intent(out) :: expMatrix(:,:,:)
    integer, allocatable, intent(out) :: disIntMat(:,:)
    double precision, intent(out):: U
    
    ! Local variables
    double precision, allocatable ::  scatterData(:), UD_dg(:)
    double precision, allocatable :: expData(:,:,:), uVec(:)
    double precision ::  expTime, sumTime, totTime, setUpTime
    integer, allocatable :: triMat(:,:), triScatter(:,:)
    integer ::  newSize, eCols, i, kP(3), nSum, dataSize


    ! Declare constants and rows of permutation matrix
    totTime = MPI_Wtime()
    root = 0
    N_p = 6
    setUpTime = MPI_Wtime()


    ! Set up on root
    if (processRank .eq. root) then

       print *, ' '
       print *, ' '
       print *, '========================'
       print *, 'Beginning non-additive calculation for whole sim box'
       print *, ' '

       ! Read in all necessary info from files
       !call initialise(fileName, posArray,trainData,alpha,hyperParams,N_tp,nArgs, &
       !                N_a,N_tri,udSize)
       allocate(uFull(N_tri))
       allocate(X_dg(N_a,N_a))

       ! Set up the arrays required for the non-additive calculation
       call makeXdgNonAdd(N_a,posArray, X_dg)
       call makeDisIntMatNonAdd(N_a, disIntMat)
       call makeUDdgNonAdd(N_a,udSize,X_dg, UD_dg)

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

  !     allocate(alpha(N_tp))
  !     allocate(trainData(N_tp,nArgs))
       allocate(disIntMat(N_a,N_a))
       allocate(X_dg(N_a,N_a))
       !allocate(posArray(N_a,3))

    end if
    print *, 'end if'

    ! Broadcast new arrays from root
    !call MPI_Bcast(alpha, N_tp, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    !call MPI_Bcast(trainData, N_tp*nArgs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
    !               ierror)
    call MPI_Bcast(disIntMat, N_a*N_a, MPI_INT, root, MPI_COMM_WORLD, ierror)
    print *, 'broadcasted disIntMat'
    call MPI_Bcast(X_dg, N_a*N_a, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    print *, 'broadcasted Xdg'
    !call MPI_Bcast(posArray, N_a*3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)
    print *, 'broadcast'


    ! Determine no. of elements of UD_dg to send to each process for exp
    ! calculations and no. of triplets to send to each for energy calculations
    call getDistsAndTripletsPerProcNonAdd(udSize,N_tri,clusterSize, dataSize,nSum)
    allocate(scatterData(dataSize)) ! Allocate array to send exponentials
    print *, dataSize


    ! Scatter the interatomic distances in U_dg to all processes
    call MPI_Scatter(UD_dg, dataSize, MPI_DOUBLE_PRECISION, scatterData, dataSize, &
                     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)


    ! Scatter the triplet matrix
    allocate(triScatter(3,nSum))
    call MPI_Scatter(triMat, nSum*3, MPI_INT, triScatter, nSum*3, MPI_INT, root, &
                     MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Calculate the exponentials for each distance on each process
    setUpTime = MPI_Wtime() - setUpTime
    expTime = MPI_Wtime()
    allocate(expData(nArgs,N_tp,dataSize))
    call calculateExponentialsNonAdd(dataSize,N_tp,nArgs,trainData,hyperParams(1), &
                                     scatterData,N_a, expData)


    ! Allocate an array to hold all exps
    eCols = nArgs*N_tri !udSize
    allocate(expMatrix(nArgs,N_tp,udSize))


    ! Gather expData arrays from the other processes and add them to expMatrix on
    ! the root process
    call MPI_Gather(expData, N_tp*nArgs*dataSize, MPI_DOUBLE_PRECISION, expMatrix, &
                    N_tp*nArgs*dataSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                    ierror)
    expTime = MPI_Wtime() - expTime


    ! Broadcast expMatrix to all processes so that sum can be parallelised
    sumTime = MPI_Wtime()
    call MPI_Bcast(expMatrix, N_tp*nArgs*udSize, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
                   ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Find the energies of the triplets assigned to each process
    allocate(uVec(nSum))
    call tripletEnergiesNonAdd(triScatter,disIntMat,nSum,N_tp,N_a,N_p,nArgs,Perm, &
                               udSize,expMatrix,alpha,hyperParams(2), uVec)


    ! Gather in the triplet energies and sum them to get total non-add energy
    call MPI_Gather(uVec, nSum, MPI_DOUBLE_PRECISION, uFull, nSum, MPI_DOUBLE_PRECISION, &
                    root, MPI_COMM_WORLD, ierror)


    ! Find the total non-additive energy for the system
    U = 0d0
    if (processRank .eq. root) then

       call totalEnergyNonAdd(uFull,N_tri, U)
       print *, uFull
       print *, "The total non-additive energy is", U
       print *, "              "

    end if


    ! De-allocate arrays not passed to atom-move subroutine
    deallocate(scatterData)
    deallocate(expData)
    deallocate(uVec)
    deallocate(triScatter)
    if (processRank .eq. root) then

       deallocate(UD_dg)
       deallocate(triMat)
       deallocate(uFull)

    end if


    ! Print times taken for each part of subroutine to run
    sumTime = MPI_Wtime() - sumTime
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

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

    return
  end subroutine triplet_mpi_fullNonAdd



  subroutine triplet_mpi_moveNonAdd(N_move,dist,N_a,N_tri,udSize,posArray,X_dg,disIntMat,&
                                    expMatrix,deltaU,newUfull)
    integer :: N_a, root=0, newSize,  triPerProc, j, indj, nPerProc
    integer :: dataSize, barError, udSize, kP(3), N_tri
    integer :: triInt, i, N_move, move, triPerAt, triInd, disIntMat(N_a,N_a)
    double precision :: X_dg(N_a,N_a)
    double precision :: expMatrix(nArgs,N_tp,udSize), posArray(N_a,nArgs)
    double precision :: newPosAt(N_a,nArgs), newX_dg(N_a,N_a), deltaU
    double precision :: totTime, moveTime, dist
    integer, allocatable :: newExpInt(:,:), changeTriMat(:,:), tripIndex(:)
    integer, allocatable :: indPerTrip(:,:), scatterTripInd(:), scatterTrip(:,:)
    integer, allocatable :: scatterIndPerTrip(:,:), changedTriplets(:,:)
    integer, allocatable :: scatterExpInd(:,:)
    double precision, allocatable :: newDists(:), newUfull(:), newUvec(:)
    double precision, allocatable :: changedTriDists(:,:), newExpMat(:,:,:)
    double precision, allocatable :: scatterDists(:), changeExpData(:,:,:)
    double precision, allocatable :: changeExpMat(:,:,:)


    if (processRank .eq. root) then

       print *, ' '
       print *, ' '
       print *, '========================'
       print *, 'Beginning non-additive calculation for atom move'
       print *, ' '

    end if


    totTime = MPI_Wtime()
    call getTriPerAtom(N_a, triPerAt)
    triPerProc = triPerAt / clusterSize
    nPerProc = (N_a-1)/clusterSize
    allocate(scatterDists(nPerProc))
    allocate(scatterTrip(3,triPerProc))
    allocate(scatterExpInd(2,nPerProc))
    allocate(newExpInt(2,N_a-1))
    allocate(newDists(N_a-1))
    allocate(indPerTrip(2,triPerAt))
    allocate(changeExpData(nArgs,N_tp,nPerProc))
    allocate(changeExpMat(nArgs,N_tp,N_a-1))
    allocate(newExpMat(nArgs,N_tp,udSize))
    allocate(changedTriplets(3,triPerAt))
    allocate(newUvec(triPerProc))
    allocate(newUfull(triPerAt))
    allocate(changedTriDists(3,triPerAt))
    allocate(tripIndex(triPerAt))


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()
    do i = 1, N_move

       newExpMat = expMatrix

       ! Set up on root
       if (processRank .eq. root) then

          ! Move an atom
          call moveAt(posArray,N_a,dist, newPosAt,move)

          ! Re-calculate X_dg for the new atomic positions
          call makeXdgNonAdd(N_a,newPosAt, newX_dg)

          ! Find the indices of the affected exponentials
          call extractChangedExps(N_a,move,newX_dg, newExpInt,newDists)

          ! Determine which triplets have undergone a change
          call getChangedTriplets(move,N_a,newX_dg,triPerAt, &
                                  changedTriplets,changedTriDists)
          call findChangedTriIndex(triPerAt,N_a,move, tripIndex)
          call findChangedDistsPerTrip(triPerAt,changedTriplets,move, indPerTrip)

       end if

       ! Scatter all requisite data from move set-up on root to all procs
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_Scatter(changedTriplets, triPerProc*3, MPI_INT, scatterTrip, &
            triPerProc*3, MPI_INT, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(newExpInt, 2*(N_a-1), MPI_INT, root, MPI_COMM_WORLD, &
            ierror)
       call MPI_Bcast(newX_dg, N_a*N_a, MPI_DOUBLE_PRECISION, root, &
            MPI_COMM_WORLD, ierror)
       call MPI_Bcast(newPosAt, 3*N_a, MPI_DOUBLE_PRECISION, root, &
            MPI_COMM_WORLD, ierror)
       call MPI_Scatter(newDists, nPerProc, MPI_DOUBLE_PRECISION, &
            scatterDists, nPerProc, MPI_DOUBLE_PRECISION, &
            root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(move, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(tripIndex, triPerAt, MPI_INT, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(indPerTrip, 2*triPerAt, MPI_INT, root, MPI_COMM_WORLD, &
            ierror)
       call MPI_Bcast(changedTriDists, 3*triPerAt, MPI_DOUBLE_PRECISION, root, &
            MPI_COMM_WORLD, ierror)

       ! Update exponentials in changed triplets
       call calculateExponentialsNonAdd(nPerProc,N_tp,nArgs,trainData, &
            hyperParams(1),scatterDists,N_a, &
            changeExpData)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Gather in all updated exps and broadcast the resultant matrix to all procs
       call MPI_Gather(changeExpData, N_tp*nArgs*nPerProc, &
            MPI_DOUBLE_PRECISION, changeExpMat, N_tp*nArgs*nPerProc, &
            MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_Bcast(changeExpMat, N_tp*nArgs*(N_a-1), &
            MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)

       ! Update the exp matrix
       do j = 1, N_a-1

          indj = disIntMat(newExpInt(1,j),newExpInt(2,j))
          newExpMat(1:nArgs,1:N_tp,indj) = changeExpMat(1:nArgs,1:N_tp,j)

       end do

       ! Calculate the non-additive energies for the changed triplets and gather
       call tripletEnergiesNonAdd(scatterTrip,disIntMat,triPerProc,N_tp,N_a,N_p,nArgs,Perm, &
            N_a-1,newExpMat,alpha,hyperParams(2), newUvec)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_Gather(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, triPerProc, &
            MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Find total change in non-add energy from moving atom
       deltaU = 0d0
       if (processRank .eq. root) then

          call totalEnergyNonAdd(newUfull,triPerAt, deltaU)
          print *, "The change in non-additive energy after the move is", deltaU
          print *, '------------------------'
          print *, ' '

       end if

       ! Update data (accept all moves for now so auto-update in each loop)
       posArray = newPosAt
       X_dg = newX_dg
       do j = 1, N_a-1

          indj = disIntMat(newExpInt(1,j),newExpInt(2,j))
          expMatrix(1:nArgs,1:N_tp,indj) = changeExpMat(1:nArgs,1:N_tp,j)

       end do

    end do
    moveTime = MPI_Wtime() - moveTime


    ! Deallocate all arrays on root process and any shared across processes
    deallocate(scatterTrip)
    deallocate(scatterExpInd)
    deallocate(indPerTrip)
    deallocate(changeExpData)
    deallocate(changeExpMat)
    deallocate(tripIndex)
    deallocate(newExpInt)
    deallocate(newDists)
    deallocate(newExpMat)
    deallocate(changedTriplets)
    deallocate(changedTriDists)


    ! Finalise MPI and print times taken for each step of calculation
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       print *, "The time taken to do", N_move, "moves was", moveTime, "seconds"
       print *, "The total time for the program to run was", totTime, "seconds"
       print *, ' '
       print *, 'Non-additive calculation for atom move complete'
       print *, '========================'
       print *, ' '
       print *, ' '

    end if

    return
  end subroutine triplet_mpi_moveNonAdd
end module triplet_mpi_mod
