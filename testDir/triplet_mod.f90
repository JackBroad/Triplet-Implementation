module triplet_mod
  implicit none


contains


! Subroutine to read in matrix of atomic positions and return all requisite
! constants
subroutine initialise(posAt, trainData, alpha, hyperParams, N_tp, nArgs, N_a, N_tri, udSize)
  implicit none
  double precision, allocatable, intent(out) :: posAt(:,:), alpha(:)
  double precision, allocatable, intent(out) :: trainData(:,:)
  double precision, intent(out) :: hyperParams(3)
  integer :: i, j, k, l, m, n
  integer, intent(out) :: N_tp, nArgs, N_a, N_tri, udSize

  ! Read in hyperparameters
  open(1, file='hyperParam.txt', status='old')
  do i = 1, 3
    read(1,*) hyperParams(i)
  end do
  close(1)

  ! Get the no of TPs and alpha values for each
  open(2, file='alpha.txt', status='old')
  read(2,*) N_tp
  allocate(alpha(N_tp))
  read(2,*) (alpha(j), j=1,N_tp)
  close(2)

  ! Read in nArgs and the distances for each TP
  open(3, file='trainingSet.txt', status='old')
  read(3,*) nArgs
  allocate(trainData(N_tp,nArgs))
  do k = 1, N_tp
    read(3,*) (trainData(k,l), l=1,nArgs)
  end do
  close(3)

  ! Read in the number of atoms
  open(4, file='AtomicPositions.txt', status='old')
  read(4,*) N_a

  ! Get the atomic positions
  allocate(posAt(N_a,nArgs))
  do m = 1, N_a
    read(4,*) (posAt(m,n), n=1,nArgs)
  end do
  close(4)

  ! Calculate the number of triplets
  N_tri = N_a**3
  N_tri = N_tri - 3*N_a**2
  N_tri = N_tri + 2*N_a
  N_tri = N_tri / 6 ! No. of triplets

  ! Determine size of UD of X_dg
  udSize = ((N_a * N_a) - N_a) / 2

return
end subroutine initialise


! Subroutine to make all arrays for the non-additive calculation
subroutine makeArraysNonAdd(nAt,udSize,posArray, X_dg,disIntMat,UD_dg)
  implicit none
  integer, intent(in) :: nAt, udSize
  double precision, intent(in) :: posArray(nAt,nAt)
  double precision, allocatable, intent(out) :: X_dg(:,:), UD_dg(:)
  integer, allocatable, intent(out) :: disIntMat(:,:)
  integer :: i, j, ind, indi, indj, nMinus, ele, m, n

  ! Find X_dg for the atomic positions in posArray
  allocate(X_dg(nAt,nAt))
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

  ! Read the UD of X_dg into a separate vector
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
end subroutine makeArraysNonAdd


! Gets distances per process that need exponential calculations for and number
! of triplets needed per node for non-additve energy calculation
subroutine getDistsAndTripletsPerProcNonAdd(nDist,nTrips,nProc, &
                                            distsPerProc,tripsPerProc)
  implicit none
  integer, intent(in) :: nDist, nProc, nTrips
  integer, intent(out) :: distsPerProc, tripsPerProc

  distsPerProc = nDist / nProc
  tripsPerProc = nTrips / nProc

return
end subroutine getDistsAndTripletsPerProcNonAdd


! Calculates the exponentials required to find the non-additive energy
subroutine calculateExponentialsNonAdd(nJobs,nTP,nArguments,procData,trainingData, &
                                       lengthscale, exponentials)
  implicit none
  integer, intent(in) :: nJobs, nTP, nArguments
  double precision, intent(in) :: procData(nJobs), trainingData(nTP,nArguments)
  double precision, intent(in) :: lengthscale
  double precision, intent(out) :: exponentials(nTP,nJobs*nArguments)
  integer :: i, j, k
  double precision :: expon

  do i = 1, nJobs
    do j = 1, nTP
      do k = 1, nArguments

        expon = (procData(i) - trainingData(j,k))**2 / 2/lengthscale**2
        exponentials(j,k+(i-1)*nArguments) = exp(-1*expon)

      end do
    end do
  end do

return
end subroutine calculateExponentialsNonAdd


subroutine makeTripletMatrix(nAt,nTri, tripletMatrix)
  implicit none
  integer, intent(in) :: nAt, nTri
  integer, intent(out) :: tripletMatrix(3,nTri)
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
end subroutine makeTripletMatrix


! Calculates the non-additive energy for each triplet in the cluster
subroutine tripletEnergiesNonAdd(triData,intMat,nJob,nTP,nAt,nPerm,nArg,permMat, &
                                 expCols,expMat,alphaVec,sigVar, uVector)
  implicit none
  integer, intent(in) :: nJob, nTP, nAt, nPerm, nArg, expCols
  integer, intent(in) :: triData(3,nJob), intMat(nAt,nAt), permMat(nPerm,nArg)
  double precision, intent(in) :: expMat(nTP,expCols), alphaVec(nTP), sigVar
  double precision, intent(out) :: uVector(nJob)
  integer :: triInt, alp, bet, gam, alDis, beDis, gaDis, r, s, kP(nArg)
  integer :: alInt, beInt, gaInt
  double precision :: alphaSum, expSum, expProd

  do triInt = 1, nJob

    ! Get values of alpha, beta and gamma
    alp = triData(1,triInt)
    bet = triData(2,triInt)
    gam = triData(3,triInt)

    ! Convert these into the indices of the IA distances that describe the triplet
    alDis = intMat(alp,bet)
    beDis = intMat(alp,gam)
    gaDis = intMat(bet,gam)

    ! Need to add to alpha sum for each TP so set to zero out of TP loop
    alphaSum = 0

    do r = 1, nTP

      ! The sum of the exp products is different for each TP so reset for
      ! each
      expSum = 0

      do s = 1, nPerm

        ! Take sth row of permutation matrix
        kP = permMat(s,1:3)

        ! Find where to look in expMatrix for exponentials
        alInt = kP(1) + (alDis-1)*nArg
        beInt = kP(2) + (beDis-1)*nArg
        gaInt = kP(3) + (gaDis-1)*nArg

        ! Find the product of the relevant exps under this permutation
        expProd = expMat(r,alInt) * expMat(r,beInt) * expMat(r,gaInt)

        ! Add the product to the sum for this training point
        expSum = expSum + expProd

      end do

      ! Multiply the sum by the relevant alpha value and add the result to all
      ! previous results
      alphaSum = alphaSum + expSum * alphaVec(r)

    end do

    uVector(triInt) = sigVar * alphaSum

  end do

return
end subroutine tripletEnergiesNonAdd


subroutine totalEnergyNonAdd(uVector,nTri, uFinal)
  implicit none
  integer, intent(in) :: nTri
  double precision, intent(in) :: uVector(nTri)
  double precision, intent(out) :: uFinal
  integer :: i

  uFinal = 0.0
  do i = 1, nTri

    uFinal = uFinal + uVector(i)

  end do

return
end subroutine totalEnergyNonAdd


end module triplet_mod
