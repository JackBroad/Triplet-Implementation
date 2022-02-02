module triplet_mod
  !use mpi
  use, intrinsic :: ISO_C_BINDING, only : C_PTR
  use GP_variables
  use mpi_variables
  use dataStructure_variables
  use expShare_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData  
  implicit none


  public


contains


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


! Get amount of data to send to each process for a scatterv call (or receive in
! a gatehrv call)
subroutine getNPerProcNonAdd(nDist,nProc, maxDat,reDat)
  implicit none
  integer, intent(in) :: nDist, nProc
  integer, intent(out) :: maxDat, reDat
  integer :: counter

  counter = 1

  do while (nDist .ge. counter*nProc)

    maxDat = counter
    counter = counter + 1

  end do

  reDat = nDist - maxDat*(nProc-1)

  if (reDat .le. 0) then

    reDat = maxDat

  end if

return
end subroutine getNPerProcNonAdd


! Fill arrays needed for scatterv/gatherv calls
subroutine getVarrays(nProc,maxSize,reSize, counts,stride)
  implicit none
  integer, intent(in) :: nProc, maxSize, reSize
  integer, intent(out) :: counts(nProc), stride(nProc)
  integer :: i

  do i = 1, nProc
    if (i .lt. nProc) then
      counts(i) = maxSize
    else
      counts(i) = reSize
    end if
    stride(i) = (i-1) * maxSize
  end do

return
end subroutine getVarrays


! Calculates the exponentials required to find the non-additive energy
subroutine calculateExponentialsNonAdd(nJobs,nTP,nArguments,trainingData, &
                                       lengthscale,dists, exponentials)
  implicit none
  integer, intent(in) :: nJobs, nTP, nArguments
  double precision, intent(in) :: trainingData(nTP,nArguments)
  double precision, intent(in) :: lengthscale, dists(nJobs)
  double precision, intent(out) :: exponentials(nArguments,nTP,nJobs)
  integer :: i, j, k
  double precision :: expon, num, denom

  do i = 1, nJobs
      do j = 1, nTP
        do k = 1, nArguments

          num = dists(i) - trainingData(j,k)
          num = num**2

          denom = lengthscale**2
          denom = 2*denom

          expon = num / denom
          exponentials(j,k,i) = exp(-1*expon)

      end do
    end do
  end do

return
end subroutine calculateExponentialsNonAdd


subroutine calcAllExposNonAddSharedMem(nTP,nArguments,trainingData,lengthscale)
  implicit none
  integer, intent(in) :: nTP, nArguments
  double precision, intent(in) :: trainingData(nTP,nArguments)
  double precision, intent(in) :: lengthscale
  integer :: i, j, k, nJobs, atOne, atTwo, distIndex
  double precision :: expon, num, denom, dist

  nJobs = N_exp_per_host 
  do i = 1, nJobs
    dist = fullHostDists(i)
    atOne = fullHostInds(i,1)
    atTwo = fullHostInds(i,2)
    distIndex = currentEnergyData%distancesIntMat(atOne,atTwo)
      do j = 1, nTP
        do k = 1, nArguments
          num = dist - trainingData(j,k)
          num = num**2

          denom = lengthscale**2
          denom = 2*denom

          expon = num / denom
          expArray(j,k,distIndex) = exp(-1*expon)
      end do
    end do
  end do

return
end subroutine calcAllExposNonAddSharedMem


! Calculates the non-additive energy for each triplet in the cluster
subroutine tripletEnergiesNonAdd(triData,intMat,nJob,nTP,nAt,nPerm,nArg,permMat, &
                                 expCols,expMat,alphaVec,sigVar, uVector)
  implicit none
  integer, intent(in) :: nJob, nTP, nAt, nPerm, nArg, expCols
  integer, intent(in) :: triData(3,nJob), intMat(nAt,nAt), permMat(nPerm,nArg)
  double precision, intent(in) :: expMat(nArg,nTP,expCols), alphaVec(nTP), sigVar
  double precision, intent(out) :: uVector(nJob)
  integer :: triInt, alp, bet, gam, alDis, beDis, gaDis, r, s, kP(nArg)
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

        ! Find the product of the relevant exps under this permutation
        expProd = expMat(r,kP(1),alDis) * expMat(r,kP(2),beDis) * expMat(r,kP(3),gaDis)
        !expProd = triInt*kP(1) + r*kP(2) + s*kP(3) + triInt*r*s + kP(1)*kP(2)*kP(3)

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


! Returns indices of changed distances in each triplet
subroutine findChangedIndPerTrip(nPerAt,changedTri,atom, indPerTrip)
  implicit none
  integer, intent(in) :: nPerAt, changedTri(3,nPerAt), atom
  integer, intent(out) :: indPerTrip(2,nPerAt)
  integer :: i, al, be, ga

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
end subroutine findChangedIndPerTrip


! Returns the distances in all triplets after a move
subroutine findTripletDistances(nAt,nPerAt,changedTri,X, indPerTrip)
  implicit none
  integer, intent(in) :: nAt, nPerAt, changedTri(3,nPerAt)
  double precision, intent(in) :: X(nAt,nAt)
  double precision, intent(out) :: indPerTrip(3,nPerAt)
  integer :: i, al, be, ga

  do i = 1, nPerAt

    al = changedTri(1,i)
    be = changedTri(2,i)
    ga = changedTri(3,i)

    indPerTrip(1,i) = X(changedTri(1,i),changedTri(2,i))
    indPerTrip(2,i) = X(changedTri(1,i),changedTri(3,i))
    indPerTrip(3,i) = X(changedTri(2,i),changedTri(3,i))

  end do

return
end subroutine findTripletDistances


function energyCheckCalc(xStar) result(PES_GP)
  implicit none
  double precision :: xStar(3), alpha(3)
  double precision :: PES_GP, xTraining(3,3), xTrainingPerm(3,6,3)
  double precision :: lScale=1.010933217823705432e-01
  double precision :: expVar=5.065189615948498530e-05
  integer :: i, j, k, nTraining, nPerms, nDim, perm(3,6), l, m, n
  double precision :: kSqExpAllPerms, kSqExpJthPerm, kKernTotal

  !xStar = (/ 4.9530743964810102E-002, 2.1838616589923775E-002,
  !3.0229538987517568E-002/)
  perm(:,1) = (/1, 2, 3/)
  perm(:,2) = (/1, 3, 2/)
  perm(:,3) = (/2, 1, 3/)
  perm(:,4) = (/2, 3, 1/)
  perm(:,5) = (/3, 1, 2/)
  perm(:,6) = (/3, 2, 1/)

  ! Get the no of TPs and alpha values for each
  open(2, file='smallAlpha.txt', status='old')
  read(2,*) nTraining
  read(2,*) (alpha(l), l=1,nTraining)
  close(2)

  ! Read in nArgs and the distances for each TP
  open(3, file='smallTrainingSet.txt', status='old')
  read(3,*) nDim
  do m = 1, nTraining
    read(3,*) (xTraining(n,m), n=1,nDim)
  end do
  close(3)

  nPerms = 6
  do i=1,nDim
     do j=1,nPerms
        do k=1,nTraining
           xTrainingPerm(i,j,k)=xTraining(perm(i,j),k)
        end do
     end do
  end do
  do m = 1, nPerms
    !print *, xTrainingPerm(:,m,3)
  end do

  kKernTotal=0
  do i=1,nTraining
     kSqExpAllPerms=0
     do j=1,nPerms
        kSqExpJthPerm=1
        do k=1,nDim
           kSqExpJthPerm  =  kSqExpJthPerm * &
           (exp( - (xStar(k)-xTrainingPerm(k,j,i))**2 /2.0/lScale**2))
        end do !Dimensions (k)
        kSqExpAllPerms = kSqExpAllPerms + kSqExpJthPerm
     end do !Permuations (j)
     kKernTotal = kKernTotal + alpha(i) * kSqExpAllPerms
  end do !Training points (i)

  PES_GP=kKernTotal * expVar
!  print *, 'Non-additive E:', PES_GP
end function energyCheckCalc


  subroutine updateExpMatrix(changeData,expInd,length)
    integer :: indj, length, expInd(length,2), j
    double precision :: changeData(N_tp,nArgs,length)

    do j = 1, length
      indj = proposedEnergyData%distancesIntMat(expInd(j,1), expInd(j,2))
      expArray(1:N_tp,1:nArgs,indj) = changeData(1:N_tp,1:nArgs,j)
    end do

  return
  end subroutine updateExpMatrix


end module triplet_mod
