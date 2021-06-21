
! Call functions to initialise MPI
totTime = MPI_Wtime()
call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


! Declare constants and rows of permutation matrix
root = 0
N_p = 6


setUpTime = MPI_Wtime()
! Read in all files to root process only
if (processRank .eq. root) then

  ! Hyperparameters
  open(1, file='hyperParam.txt', status='old')
  do i = 1, 3
    read(1,*) hyperParams(i)
  end do
  close(1)
  !print *, "The hyperparameters are", hyperParams
  !print *, "                         "

  ! Delta-Gamma matrix and no. of atoms
  open(2, file='deltaGamma400x400.txt', status='old')
  read(2,*) N_a
  allocate(X_dg(N_a,N_a))
  read(2,*) ((X_dg(j,k), k=1,N_a), j=1,N_a)
  close(2)
  N_tri = N_a**3
  N_tri = N_tri - 3*N_a**2 
  N_tri = N_tri + 2*N_a 
  N_tri = N_tri / 6 ! No. of triplets
  allocate(uFull(N_tri))
  print *, "The number of atoms is", N_a
  !print *, "The array of distances is", X_dg
  print *, "The number of triplets is", N_tri
  !print *, "                         "

  ! Find the size of the UD of X_dg and the number of distances to be sent to
  ! each process
  udSize = ((N_a * N_a) - N_a) / 2
  !print *, "The number of interatomic distances is", udSize
  !print *, "                         "
  allocate(UD_dg(udSize))

  ! Read the UD of X_dg into a separate vector
  nMinus = N_a - 1
  ele = 0
  do m = 1, nMinus
    do n = m+1, N_a
      ele = ele + 1!m + n - 2
      UD_dg(ele) = X_dg(m,n)
    end do
  end do
  !print *, "The upper diagonal of the array of distances is", UD_dg
  !print *, "                         "

  ! Read in the no of TPs and alpha values for each
  open(3, file='alpha.txt', status='old')
  read(3,*) N_tp
  allocate(alpha(N_tp))
  read(3,*) (alpha(o), o=1,N_tp)
  close(3)
  !print *, "The number of TPs is", N_tp
  !print *, "The first value in alpha is", alpha(1)
  !print *, "The final value in alpha is", alpha(N_tp)
  !print *, "                         "

  ! Read in nArgs and the distances for each TP
  open(4, file='trainingSet.txt', status='old')
  read(4,*) nArgs
  allocate(trainData(N_tp,nArgs))
  do p = 1, N_tp
    read(4,*) (trainData(p,q), q=1,nArgs)
  end do
  close(4)
  !print *, "nArgs is", nArgs
  !print *, "The first TP has inv distances", trainData(1,1:nArgs)
  !print *, "The second TP has inv distances", trainData(2,1:nArgs)
  !print *, "The last TP has inv distances", trainData(N_tp,1:nArgs)
  !print *, "                         "

  ! Declare permutation matrix
  Perm(1,:) = (/1, 2, 3/)
  Perm(2,:) = (/1, 3, 2/)
  Perm(3,:) = (/2, 1, 3/)
  Perm(4,:) = (/2, 3, 1/)
  Perm(5,:) = (/3, 1, 2/)
  Perm(6,:) = (/3, 2, 1/)
  !print *, "The permutation matrix is"
  !print *, Perm(1,1:3)
  !print *, Perm(2,1:3)
  !print *, Perm(3,1:3)
  !print *, Perm(4,1:3)
  !print *, Perm(5,1:3)
  !print *, Perm(6,1:3)
  !print *, "                         "

end if


! Hold all processes here until root process has finished reading data
call MPI_BARRIER(MPI_COMM_WORLD, barError)


! Broadcast all requisite data that is already allocated everywhere from root to 
! all processes
call MPI_Bcast(hyperParams, 3, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(N_a, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(N_tri, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(udSize, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(N_tp, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(nArgs, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(Perm, 18, MPI_INT, root, MPI_COMM_WORLD, ierror)
call MPI_BARRIER(MPI_COMM_WORLD, barError)


! Use info that was just broadcasted to allocate arrays on other processes and
! then broadcast equivalent arrays from root
if (processRank .ne. root) then

  allocate(alpha(N_tp))
  allocate(trainData(N_tp,nArgs))

end if

call MPI_Bcast(alpha, N_tp, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
call MPI_Bcast(trainData, N_tp*nArgs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, &
               ierror)
call MPI_BARRIER(MPI_COMM_WORLD, barError)


! Set up matrix of the indices of each IA distance (same rank as X_dg)
allocate(disIntMat(N_a,N_a))
ind = 0
do indi = 1, N_a

  do indj = 1, N_a

    if (indi .eq. indj) then

      disIntMat(indi,indj) = 0

    else if (indi .lt. indj) then

      ind = ind + 1
      disIntMat(indi,indj) = ind

    else

      disIntMat(indi,indj) = disIntMat(indj,indi)

    end if

  end do

  !if (processRank .eq. root) then

   ! print *, disIntMat(indi,1:N_a)

  !end if

end do
!print *, "               "



! Determine no. of elements of UD_dg to send to each process and allocate vector
! to do it
dataSize = udSize / clusterSize
allocate(scatterData(dataSize))
setUpTime = MPI_Wtime() - setUpTime
expTime = MPI_Wtime()


! Scatter the interatomic distances in U_dg to all processes
call MPI_Scatter(UD_dg, dataSize, MPI_DOUBLE_PRECISION, scatterData, dataSize, &
                 MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)
!print *, "Process ", processRank, " received ", scatterData
!print *, "                         "


! Calculate the exponentials for each distance on each process
allocate(expData(N_tp,nArgs*dataSize))
do dsInt = 1, dataSize

  do neRow = 1, N_tp

    do neCol = 1, nArgs

      expon = (scatterData(dsInt) - trainData(neRow,neCol))**2 / 2/hyperParams(1)**2
      expData(neRow,neCol+(dsInt-1)*nArgs) = exp(-1*expon) 

    end do

  end do

end do
call MPI_BARRIER(MPI_COMM_WORLD, barError)
!print *, "The first row of the exp matrix for process", processRank,"is", &
!expData(1,:)
!print *, "                         "


! Allocate an array to hold all exps
allocate(expMatrix(N_tp,nArgs*udSize))

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


! Declare how many triplet sums each process must undertake
nSum = N_tri / clusterSize
allocate(triScatter(3,nSum))
allocate(uVec(nSum))
!print *, "The number of triplets per cluster is", nSum
!print *, "                  "


! On the root loop over alpha, beta and gamma to give a matrix of all triplets
if (processRank .eq. root) then

  allocate(triMat(3,N_tri))
  t = 0

  do al = 1, N_a-2

    do be = al+1, N_a-1

      do ga = be+1, N_a

      t = t+1
      triMat(1,t) = al
      triMat(2,t) = be
      triMat(3,t) = ga

      end do

    end do

  end do

  !print *, "The triplet matrix is"
  !print *, triMat(1,:)
  !print *, triMat(2,:)
  !print *, triMat(3,:)
  !print *, "                  "

end if


! Scatter the triplet matrix
call MPI_Scatter(triMat, nSum*3, MPI_INT, triScatter, nSum*3, MPI_INT, root, &
                 MPI_COMM_WORLD, ierror)
call MPI_BARRIER(MPI_COMM_WORLD, barError)
!print *, "Process", processRank, "received"
!print *, triScatter(:,1)
!print *, triScatter(:,2)
!print *, "                   "


! Find the products of the exponential for each training point and sum over the
! number of permutations to get the 'symmetric kernel exponential'; multiply
! this by the alpha value for the TP, sum the results for each TP and multiply
! by the signal variance to get the energy for each triplet
do triInt = 1, nSum

  ! Get values of alpha, beta and gamma
  alp = triScatter(1,triInt)
  bet = triScatter(2,triInt)
  gam = triScatter(3,triInt)
  !print *, "The triplet comprises atoms", alp, bet, gam

  ! Convert these into the indices of the IA distances that describe the triplet
  alDis = disIntMat(alp,bet)
  beDis = disIntMat(alp,gam)
  gaDis = disIntMat(bet,gam)
  !print *, "This corresponds to interatomic distances",alDis, beDis, gaDis
  !print *, "                    "

  ! Need to add to alpha sum for each TP so set to zero out of TP loop
  alphaSum = 0

  do r = 1, N_tp

    ! The sum of the exp products is different for each TP so reset for
    ! each

    expSum = 0

    do s = 1, N_p

      ! Take sth row of permutation matrix
      kP = Perm(s,1:3)

      ! Find where to look in expMatrix for exponentials
      alInt = kP(1) + (alDis-1)*nArgs
      beInt = kP(2) + (beDis-1)*nArgs
      gaInt = kP(3) + (gaDis-1)*nArgs
      !if (r .eq. 1) then

      !    print *, "The current row of the perm matrix is", kP
      !    print *, "Looking in column", alInt, "for exps for first IA dist"
      !    print *, "Looking in column", beInt, "for exps for second IA dist"
      !    print *, "Looking in column", gaInt, "for exps for third IA dist"
      !    print *, "               "

      !end if

      ! Find the product of the relevant exps under this permutation
      expProd = expMatrix(r,alInt) * expMatrix(r,beInt) * expMatrix(r,gaInt)

      ! Add the product to the sum for this training point
      expSum = expSum + expProd

    end do

    ! Multiply the sum by the relevant alpha value and add the result to all
    ! previous results
    alphaSum = alphaSum + expSum * alpha(r)

  end do

  uVec(triInt) = hyperParams(2) * alphaSum

  !print *, "The energy of triplet", triInt, "on process", processRank, &
  !"is", uVec(triInt)

end do


! Gather in the triplet energies and sum them to get total non-add energy
call MPI_Gather(uVec, nSum, MPI_DOUBLE_PRECISION, uFull, nSum, MPI_DOUBLE_PRECISION, &
                root, MPI_COMM_WORLD, ierror)
U = 0
if (processRank .eq. root) then

  do sumU = 1, N_tri

    U = U + uFull(sumU)

  end do

  print *, "The total non-additive energy is", U
  print *, "              "

end if
sumTime = MPI_Wtime() - sumTime


! Deallocate all arrays on root process and any shared across processes
deallocate(scatterData)
deallocate(alpha)
deallocate(trainData)
deallocate(expData)
deallocate(expMatrix)
deallocate(uVec)
deallocate(triScatter)
deallocate(disIntMat)

if (processRank .eq. root) then

  deallocate(X_dg)
  deallocate(UD_dg)
  deallocate(triMat)
  deallocate(uFull)

end if


! Finalise MPI and print times taken for each step of calculation
call MPI_FINALIZE(ierror)
totTime = MPI_Wtime() - totTime
if (processRank .eq. root) then

  print *, "The time taken for the exponentials was", expTime, "seconds"
  print *, "The time taken for the sum was", sumTime, "seconds"
  print *, "The time taken to set up was", setUpTime, "seconds"
  print *, "The total time for the program to run was", totTime, "seconds"

end if
