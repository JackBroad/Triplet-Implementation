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
  open(4, file='AtomicPositions5.txt', status='old')
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


! Subroutine to make X_dg for the non-additive calculation
subroutine makeXdgNonAdd(nAt,posArray, X_dg)
  implicit none
  integer, intent(in) :: nAt
  double precision, intent(in) :: posArray(nAt,nAt)
  double precision, allocatable, intent(out) :: X_dg(:,:)
  integer :: i, j

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

return
end subroutine makeXdgNonAdd


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


! Uses the number of atoms (nAt) and number of triplets (nTri) to build a matrix
! of all possible triplets
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


! Sums the vector of triplet energies (uVector) over the no. of triplets (nTri)
! to deteremine the total non-additive energy
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


! Subroutine to move an atom at random and return the new atomic
! positions and the index of the moved atom.
! Takes Cartesian atomic positions (pos), the no. of atoms (num)
! and the max. distance to move in any direction (dMax) as its arguments
subroutine moveAt(pos,num,dMax,  newPos,mover)
  implicit none
  integer, intent(in) :: num
  double precision, intent(in) :: pos(num,3), dMax
  integer, intent(out) :: mover
  double precision, intent(out) :: newPos(num,3)
  integer :: icol, irow
  double precision :: randNo

  ! Pick an atom at random to move (can't generate random int in fortran, hence
  ! the work-around)
  call random_number(randNo)
  mover = 1 + FLOOR(num*randNo)

  ! Change the position of the atom in newPos
  do irow = 1, num
    do icol = 1, 3
      ! Set each element of newPos to be the same as the equiv. in pos
      newPos(irow,icol) = pos(irow,icol)
      ! If in row corresponding to moving atom, change coords of this atom in
      ! newPos
      if (irow .eq. mover) then
        call random_number(randNo)
        newPos(irow,icol) = newPos(irow,icol) + (2.0*randNo-1) * dMax
      end if
    end do
  end do

  !print *, "Moving atom", mover
  !print *, "                 "

  !print *, "Its old position was:"
  !print *, pos(mover,:)
  !print *, "                 "

  !print *, "Its new position is:"
  !print *, newPos(mover,:)
  !print *, "                 "

return
end subroutine moveAt


! Uses the number of atoms (num) and index of moved atom (atom) to determine
! which exponentials need re-calculating after a move and return a matrix of
! their postions in X_dg
subroutine extractChangedExps(num,atInd,Xdg, change,dists)
  implicit none
  integer, intent(in) :: num, atInd
  double precision, intent(in) :: Xdg(num,num)
  integer, intent(out) :: change(num-1,2)
  double precision, intent(out) :: dists(num-1)
  integer :: i, j

  print *, '============================================'
  print *, 'Moved atom', atInd
  print *, '                     '

  do i = 1, num

    if (i .lt. atInd) then

      change(i,1) = i
      change(i,2) = atInd

    else if (i .gt. atInd) then

      change(i-1,1) = atInd
      change(i-1,2) = i

    end if

  end do

  do j = 1, num-1

    dists(j) = Xdg(change(j,1),change(j,2))

  end do

  print *, 'This changed the following distances:'
  print *, change(1,:)
  print *, change(2,:)
  print *, change(3,:)
  print *, change(4,:)
  print *, '                   '

  !print *, 'These were extracted from the following X_dg:'
  !print *, Xdg(1,:)
  !print *, Xdg(2,:)
  !print *, Xdg(3,:)
  !print *, Xdg(4,:)
  !print *, Xdg(5,:)
  !print *, '                   '

  !print *, 'The distances extracted were:'
  !print *, dists
  !print *, '                   '

return
end subroutine extractChangedExps


subroutine getChangedTriplets(atom,nAt, changedTriplets,nPerAt)
  implicit none
  integer, intent(in) :: atom, nAt
  integer, intent(out) :: nPerAt
  integer, allocatable, intent(out) :: changedTriplets(:,:)
  integer :: a, b, al, be, ga, counter

  ! Determine the number of triplets each atom is involved in
  nPerAt = 0
  do a = 2, nAt-1
    do b = a+1, nAt
      nPerAt = nPerAt + 1
    end do
  end do
  print *, 'The number of triplets each atom is part of is', nPerAt
  print *, '                     '

  ! Allocate array of changed triplets
  allocate(changedTriplets(nPerAt,3))

  ! Fill changed triplet array
  counter = 0
  if (atom .eq. 1) then
    al = atom
    do be = al+1, nAt-1
      do ga = be+1, nAt

        counter = counter + 1
        changedTriplets(counter,1) = al
        changedTriplets(counter,2) = be
        changedTriplets(counter,3) = ga

      end do
    end do
  else if (atom .eq. nAt) then
    ga = atom
    do al = 1, nAt-2
      do be = al+1, nAt-1

        counter = counter + 1
        changedTriplets(counter,1) = al
        changedTriplets(counter,2) = be
        changedTriplets(counter,3) = ga

      end do
    end do
  else
    do al = 1, nAt-2
      do be = al+1, nAt-1
        do ga = be+1, nAt
          if (al .eq. atom) then

            counter = counter + 1
            changedTriplets(counter,1) = al
            changedTriplets(counter,2) = be
            changedTriplets(counter,3) = ga

          else if (be .eq. atom) then

            counter = counter + 1
            changedTriplets(counter,1) = al
            changedTriplets(counter,2) = be
            changedTriplets(counter,3) = ga

          else if (ga .eq. atom) then

            counter = counter + 1
            changedTriplets(counter,1) = al
            changedTriplets(counter,2) = be
            changedTriplets(counter,3) = ga

          end if
        end do
      end do
    end do
  end if

  print *, 'The affected triplets are:'
  do counter = 1, nPerAt
    print *, changedTriplets(counter,:)
  end do
  print *, '============================================'
  print *, '                          '

return
end subroutine getChangedTriplets


subroutine updateExponentialsNonAdd(nJob,nTP,nArgument,nAt,ecol,ind,change,dists, &
                                    expData,trainingData,triplets,nPer,perms, &
                                    updateExp)
  implicit none
  integer, intent(in) :: nJob, nTP, nArgument, nAt, ind(nAt,nAt), change(nJob,2)
  integer, intent(in) :: triplets(nJob,3), nPer, perms(nPer,3), ecol
  double precision, intent(in) :: dists(nJob), expData(nTP,ecol)
  double precision, intent(in) :: trainingData(nTP,nArgument)
  double precision, intent(out) :: updateExp(nTP,ecol)
  integer :: i, j, k, l, a, b, g, aIn, bIn, gIn, rowPerm(3)
  double precision :: expon, lengthscale

  do i = 1, nJob

    ! Find indices of the two atoms that comprise the ith changed pair
    a = change(i,1)
    b = change(i,2)
    !g = triplets(i,3)

    ! Find the index of this pairwise interaction
    aIn = ind(a,b)
    !bIn = ind(a,g)
    !gIn = ind(b,g)

    do j = 1, nTP
      do k = 1, nArgument

        ! Calculate the new exponenitial for this interaction
        expon = (dists(i) - trainingData(j,k))**2 / 2/lengthscale**2
        expon = exp(-1*expon)

        ! Update expMatrix for all triplets where this interaction appears
        do l = 1, nPer

          rowPerm = perms(l,:)

        end do
        !updateExp(j,k+(i-1)*nArguments) = exp(-1*expon)

      end do
    end do
  end do

return
end subroutine updateExponentialsNonAdd


subroutine updateTripletsNonAdd(nJob,nAt,triplets,intMat)
  implicit none
  integer, intent(in) :: nJob, nAt, intMat(nAt,nAt)
  double precision, intent(in) :: triplets(nJob,3)
  integer :: i, a, b, g, aIn, bIn, gIn

  do i = 1, nJob

    a = triplets(i,1)
    b = triplets(i,2)
    g = triplets(i,3)

    aIn = intMat(a,b)
    bIn = intMat(a,g)
    gIn = intMat(b,g)

  end do

return
end subroutine updateTripletsNonAdd


end module triplet_mod
