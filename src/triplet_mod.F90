module triplet_mod
  implicit none


contains


! Subroutine to read in matrix of atomic positions and return all requisite
! constants
subroutine initialise(fileName, posAt,trainData,alpha,hyperParams,N_tp,nArgs,N_a,N_tri,udSize)
  implicit none
  character (len=40), intent(in) :: fileName
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
  open(2, file='smallAlpha.txt', status='old')
  read(2,*) N_tp
  allocate(alpha(N_tp))
  read(2,*) (alpha(j), j=1,N_tp)
  close(2)

  ! Read in nArgs and the distances for each TP
  open(3, file='smallTrainingSet.txt', status='old')
  read(3,*) nArgs
  allocate(trainData(N_tp,nArgs))
  do k = 1, N_tp
    read(3,*) (trainData(k,l), l=1,nArgs)
  end do
  close(3)

  ! Read in the number of atoms
  open(4, file=fileName, status='old')
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
subroutine calculateExponentialsNonAdd(nJobs,nTP,nArguments,trainingData, &
                                       lengthscale,dists,nAt, exponentials)
  implicit none
  integer, intent(in) :: nJobs, nTP, nArguments, nAt
  double precision, intent(in) :: trainingData(nTP,nArguments)
  double precision, intent(in) :: lengthscale, dists(nJobs)
  double precision, intent(out) :: exponentials(nArguments,nJobs,nTP)
  integer :: i, j, k, atOne, atTwo
  double precision :: expon, dist, num, denom

  do i = 1, nJobs
      do j = 1, nTP
        do k = 1, nArguments

          num = dists(i) - trainingData(j,k)
          num = num**2

          denom = lengthscale**2
          denom = 2*denom

          expon = num / denom
          exponentials(k,i,j) = exp(-1*expon)

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
  double precision, intent(in) :: expMat(nArg,expCols,nTP), alphaVec(nTP), sigVar
  double precision, intent(out) :: uVector(nJob)
  integer :: triInt, alp, bet, gam, alDis, beDis, gaDis, r, s, kP(nArg)
  integer :: alInt, beInt, gaInt
  double precision :: alphaSum, expSum, expProd

  do triInt = 1, nJob

    ! Get values of alpha, beta and gamma
    alp = triData(1,triInt)
    bet = triData(2,triInt)
    gam = triData(3,triInt)
    !print *, '=================='
    !print *, alp, bet, gam
    !print*, ' '

    ! Convert these into the indices of the IA distances that describe the triplet
    alDis = intMat(alp,bet)
    beDis = intMat(alp,gam)
    gaDis = intMat(bet,gam)
    !print *, alDis, beDis, gaDis
    !print *, ' '

    ! Need to add to alpha sum for each TP so set to zero out of TP loop
    alphaSum = 0

    do r = 1, nTP

      ! The sum of the exp products is different for each TP so reset for
      ! each
      expSum = 0

      do s = 1, nPerm

        ! Take sth row of permutation matrix
        kP = permMat(s,1:3)

        ! Find the product of the relevant exps under this permutation
        expProd = expMat(kP(1),alDis,r) * expMat(kP(2),beDis,r) * expMat(kP(3),gaDis,r)

        ! Add the product to the sum for this training point
        expSum = expSum + expProd

      end do

      ! Multiply the sum by the relevant alpha value and add the result to all
      ! previous results
      alphaSum = alphaSum + (expSum*alphaVec(r))

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

  uFinal = 0d0
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

  print *, '========================'
  print *, "Moving atom", mover
  print *, "                 "

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
  integer, intent(out) :: change(2,num-1)
  double precision, intent(out) :: dists(num-1)
  integer :: i, j

  do i = 1, num

    if (i .lt. atInd) then

      change(1,i) = i
      change(2,i) = atInd

    else if (i .gt. atInd) then

      change(1,i-1) = atInd
      change(2,i-1) = i

    end if

  end do

  do j = 1, num-1

    dists(j) = Xdg(change(1,j),change(2,j))

  end do

return
end subroutine extractChangedExps


subroutine getTriPerAtom(nAt, nPer)
  implicit none
  integer, intent(in) :: nAt
  integer, intent(out) :: nPer
  integer :: a, b

  ! Determine the number of triplets each atom is involved in
  nPer = 0
  do a = 2, nAt-1
    do b = a+1, nAt
      nPer = nPer + 1
    end do
  end do

return
end subroutine getTriPerAtom


subroutine getChangedTriplets(atom,nAt,Xdg,nPerAt, changedTriplets, &
                              changedTriDists)
  implicit none
  integer, intent(in) :: atom, nAt, nPerAt
  double precision, intent(in) :: Xdg(nAt,nAt)
  integer, intent(out) :: changedTriplets(3,nPerAt)
  double precision, intent(out) :: changedTriDists(3,nPerAt)
  integer :: a, b, al, be, ga, counter, i

  ! Fill changed triplet array
  counter = 0
  if (atom .eq. 1) then
    al = atom
    do be = al+1, nAt-1
      do ga = be+1, nAt

        counter = counter + 1
        changedTriplets(1,counter) = al
        changedTriplets(2,counter) = be
        changedTriplets(3,counter) = ga

      end do
    end do
  else if (atom .eq. nAt) then
    ga = atom
    do al = 1, nAt-2
      do be = al+1, nAt-1

        counter = counter + 1
        changedTriplets(1,counter) = al
        changedTriplets(2,counter) = be
        changedTriplets(3,counter) = ga

      end do
    end do
  else
    do al = 1, nAt-2
      do be = al+1, nAt-1
        do ga = be+1, nAt
          if (al .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga

          else if (be .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga

          else if (ga .eq. atom) then

            counter = counter + 1
            changedTriplets(1,counter) = al
            changedTriplets(2,counter) = be
            changedTriplets(3,counter) = ga

          end if
        end do
      end do
    end do
  end if

  do i = 1, nPerAt

    changedTriDists(1,i) = Xdg(changedTriplets(1,i),changedTriplets(2,i))
    changedTriDists(2,i) = Xdg(changedTriplets(1,i),changedTriplets(3,i))
    changedTriDists(3,i) = Xdg(changedTriplets(2,i),changedTriplets(3,i))

  end do

return
end subroutine getChangedTriplets


subroutine findChangedTriIndex(nPerAt,nAt,atom, triIndex)
  implicit none
  integer, intent(in) :: nPerAt, nAt, atom
  integer, intent(out) :: triIndex(nPerAt)
  integer :: al, be, ga, i, counter

  counter = 0
  i = 0
  do al = 1, nAt-2
    do be = al+1, nAt-1
      do ga = be+1, nAt

        i = i + 1

        if (al .eq. atom) then

          counter = counter + 1
          triIndex(counter) = i

        else if (be .eq. atom) then

          counter = counter + 1
          triIndex(counter) = i

        else if (ga .eq. atom) then

          counter = counter + 1
          triIndex(counter) = i

        end if
      end do
    end do
  end do

return
end subroutine findChangedTriIndex


subroutine findChangedDistsPerTrip(nPerAt,changedTri,atom, indPerTrip)
  implicit none
  integer, intent(in) :: nPerAt, changedTri(3,nPerAt), atom
  integer, intent(out) :: indPerTrip(2,nPerAt)
  integer :: i, al, be, ga, j

  do i = 1, nPerAt

    al = changedTri(1,i)
    be = changedTri(2,i)
    ga = changedTri(3,i)

    if (al .eq. atom) then

      indPerTrip(1,i) = 1
      indPerTrip(2,i) = 2

    else if (be .eq. atom) then

      indPerTrip(1,i) = 1
      indPerTrip(2,i) = 3

    else !if (ga .eq. atom) then

      indPerTrip(1,i) = 2
      indPerTrip(2,i) = 3

    end if
  end do

return
end subroutine findChangedDistsPerTrip


subroutine updateExponentialsNonAdd(nJob,nTP,nArgument,dists,ind,triIndVec, &
                                    trainingData,lengthscale,triPerAt,nDat, &
                                    pR, updateExp)
  implicit none
  integer, intent(in) :: nJob, nTP, nArgument, triPerAt, nDat, pR
  integer, intent(in) :: ind(2,nJob), triIndVec(triPerAt)
  double precision, intent(in) :: dists(3,nJob), lengthscale
  double precision, intent(in) :: trainingData(nTP,nArgument)
  double precision, intent(out) :: updateExp(nTP,nDat*nArgument)
  double precision :: expon
  integer :: i, j, k, l, triInd, m, n, minInd, maxInd, countDone, countJobs, o
  integer :: first

  minInd = pR*nDat + 1
  maxInd = (pR+1)*nDat

  countDone = 0
  do m = 1, triPerAt
    if (triIndVec(m) .lt. minInd) then

      countDone = countDone + 1

    end if
  end do
  first = triIndVec(countDone+1)

  countJobs = 0
  do n = 1, triPerAt
    if (triIndVec(n) .ge. minInd) then
      if (triIndVec(n) .le. maxInd) then

        countJobs = countJobs + 1

      end if
    end if
  end do

  if (pR .eq. 0) then
  print *, ' '
!  print *, 'The indices of the affected triplets are:'
!  print *, triIndVec
!  print *, 'minInd is', minInd
!  print *, 'maxInd is', maxInd
!  print *, 'Max no. of triplets per node is', triPerAt
!  print *, 'No. of preceding triplets already re-calculated on other procs is', &
!           countDone
!  print *, 'No. of changed triplets to be calced by this proc is', countJobs
!  print *, 'No. of changed triplets left to be calced by procs w/ higher rank is', &
!           triPerAt-(countDone+countJobs) 
!  print *, ' '
  end if

  do i = countDone+1, countDone+countJobs

    ! Find col index of last element of first altered triplet in the full expMatrix
    triInd = triIndVec(i)

    ! Convert to first index of the same triplet in same matrix
    triInd = 3*triInd - 2
    
    ! Convert to index of same triplet in expData for current process
    triInd = triInd - (pR*nArgument*nDat)

    do j = 1, nTP
      do k = 1, 2

        l = ind(k,i)

        if (pR .eq. 0) then
        if (i .eq. first) then
          if (j .eq. 1) then
            if (k .eq. 1) then
            !print *, ' '
            !print *, 'First triplet was', first
            !print *, 'On triplet', triIndVec(i)
            !print *, 'Looking at distances', ind(1,i), 'and', ind(2,i), 'in this triplet'
            !print *, 'These had new values of', dists(ind(1,i),i), 'and', dists(ind(2,i),i)
            !print *, 'The old exponentials (from old dists not those above) were:'
            !print *, updateExp(j,(triInd-1)*nArgument+1), updateExp(j,(triInd-1)*nArgument+2) &
            !       , updateExp(j,(triInd-1)*nArgument+3)
            print *, ' '
            end if
          end if
        end if
        end if

        ! Calculate the new exponenitial and add to updated expMatrix
        expon = (dists(l,i) - trainingData(j,l))**2 / 2/lengthscale**2
        updateExp(j,(triInd-1)*nArgument+l) = exp(-1*expon)

        if (pR .eq. 0) then
        if (i .eq. first) then
          if (j .eq. 1) then
            if (k .eq. 2) then
            print *, ' '
            !print *, 'The new exponentials are:'
            !print *, updateExp(j,(triInd-1)*nArgument+1), updateExp(j,(triInd-1)*nArgument+2) &
            !       , updateExp(j,(triInd-1)*nArgument+3)
            !print *, ' '
            end if
          end if
        end if
        end if

      end do
    end do
  end do

return
end subroutine updateExponentialsNonAdd


end module triplet_mod