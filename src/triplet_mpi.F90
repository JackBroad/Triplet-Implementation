module triplet_mpi_mod
  use triplet_mod
  implicit none
  include 'mpif.h'

contains


  subroutine triplet_mpi()
    double precision, allocatable :: X_dg(:,:), scatterData(:), UD_dg(:), alpha(:)
    double precision, allocatable :: trainData(:,:), expData(:,:), expMatrix(:,:)
    double precision, allocatable :: posArray(:,:), uVec(:), uFull(:), newPosAt(:,:)
    double precision, allocatable :: newX_dg(:,:), newDists(:), newExpons(:)
    double precision :: hyperParams(3), expTime, sumTime, totTime, setUpTime, U
    double precision :: moveTime, dist
    integer, allocatable :: triMat(:,:), triScatter(:,:), disIntMat(:,:)
    integer, allocatable :: newExpInt(:,:), changeTriMat(:,:)
    integer :: processRank, clusterSize, ierror, N_a, root, newSize, N_tp, N_p
    integer :: dataSize, barError, udSize, nArgs, Perm(6,3), kP(3), N_tri, nSum
    integer :: triInt, eCols, i, N_move, move, triPerAt


    ! Call functions to initialise MPI
    totTime = MPI_Wtime()
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


    ! Declare constants and rows of permutation matrix
    root = 0
    N_p = 6
    N_move = 10
    dist = 1.5
    setUpTime = MPI_Wtime()


    ! Set up on root
    if (processRank .eq. root) then

       ! Read in all necessary info from files
       call initialise(posArray, trainData, alpha, hyperParams, N_tp, nArgs, N_a, N_tri, udSize)
       allocate(uFull(N_tri))

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

       allocate(alpha(N_tp))
       allocate(trainData(N_tp,nArgs))
       allocate(disIntMat(N_a,N_a))

    end if


    ! Broadcast new arrays from root
    call MPI_Bcast(alpha, N_tp, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(trainData, N_tp*nArgs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
         ierror)
    call MPI_Bcast(disIntMat, N_a*N_a, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Determine no. of elements of UD_dg to send to each process for exp
    ! calculations and no. of triplets to send to each for energy calculations
    call getDistsAndTripletsPerProcNonAdd(udSize,N_tri,clusterSize, dataSize,nSum)
    allocate(scatterData(dataSize)) ! Allocate array to send exponentials


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
    allocate(expData(N_tp,nArgs*dataSize))
    call calculateExponentialsNonAdd(dataSize,N_tp,nArgs,scatterData,trainData, &
         hyperParams(1), expData) ! Only expData is new


    ! Allocate an array to hold all exps
    eCols = nArgs*udSize
    allocate(expMatrix(N_tp,eCols))


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
         eCols,expMatrix,alpha,hyperParams(2), uVec)


    ! Gather in the triplet energies and sum them to get total non-add energy
    call MPI_Gather(uVec, nSum, MPI_DOUBLE_PRECISION, uFull, nSum, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, ierror)


    ! Find the total non-additive energy for the system
    U = 0
    if (processRank .eq. root) then

       call totalEnergyNonAdd(uFull,N_tri, U)
       print *, "The total non-additive energy is", U
       print *, "              "

    end if
    sumTime = MPI_Wtime() - sumTime


    ! Loop over N moves, moving an atom and re-calculating the energy each time
    moveTime = MPI_Wtime()
    allocate(newPosAt(N_a,nArgs))
    allocate(newExpInt(N_a-1,2))
    allocate(newDists(N_a-1))
    do i = 1, N_move
       if (processRank .eq. root) then

          !  print *, 'Moving an atom...'
          ! Move an atom
          call moveAt(posArray,N_a,dist, newPosAt,move)
          !  print *, 'Done'

          !  print *, 'Finding new X_dg...'
          ! Re-calculate X_dg for the new atomic positions
          call makeXdgNonAdd(N_a,newPosAt, newX_dg)
          !  print *, 'Done'

          !  print *, 'Extracting changed distances...'
          ! Find the indices of the affected exponentials
          call extractChangedExps(N_a,move,newX_dg, newExpInt,newDists)
          !  print *, 'Done'

          !  print *, 'Finding the triplets that were affected...'
          ! Determine which triplets have undergone a change
          call getChangedTriplets(move,N_a, changeTriMat,triPerAt)
          !  print *, 'Done'

          ! Update the values of the changed exponentials in expMatrix
          !  call updateExponentialsNonAdd()

          !  print *, 'Updating the triplet energies'
          !  call updateTripletsNonAdd()
          !  print *, 'Done'

       end if
    end do
    !moveTime = MPI_Wtime() - moveTime


    ! Deallocate all arrays on root process and any shared across processes
    deallocate(scatterData)
    deallocate(alpha)
    deallocate(trainData)
    deallocate(expData)
    deallocate(expMatrix)
    deallocate(uVec)
    deallocate(triScatter)
    deallocate(disIntMat)
    !deallocate(newPosAt)
    !deallocate(newExpInt)
    !deallocate(newDists)
    if (processRank .eq. root) then

       deallocate(X_dg)
       deallocate(UD_dg)
       deallocate(triMat)
       deallocate(uFull)
       !  deallocate(newX_dg)

    end if


    ! Finalise MPI and print times taken for each step of calculation
    call MPI_FINALIZE(ierror)
    totTime = MPI_Wtime() - totTime
    if (processRank .eq. root) then

       print *, "The time taken for the exponentials was", expTime, "seconds"
       print *, "The time taken for the sum was", sumTime, "seconds"
       print *, "The time taken to set up was", setUpTime, "seconds"
       ! print *, "The time taken to do", N_move, "moves was", moveTime, "seconds"
       print *, "The total time for the program to run was", totTime, "seconds"

    end if
  end subroutine triplet_mpi
END module triplet_mpi_mod
