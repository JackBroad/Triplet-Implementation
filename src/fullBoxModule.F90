module fullBoxModule
  !use mpi
  use GP_variables
  use mpi_variables
  use expShare_variables
  use pbcAndMic_variables
  use dataStructure_variables
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use triplet_mod
  implicit none


contains


  subroutine instantiateProposedDataStructs()
    implicit none

    proposedPositionData = currentPositionData
    proposedEnergyData = currentEnergyData

  end subroutine instantiateProposedDataStructs


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
  subroutine makeUDdgNonAdd(nAt,udSize,X_dg, UD_dg,distInds)
    implicit none
    integer, intent(in) :: nAt, udSize
    double precision, intent(in) :: X_dg(nAt,nAt)
    integer, allocatable, intent(out) :: distInds(:,:)
    double precision, allocatable, intent(out) :: UD_dg(:)
    integer :: nMinus, ele, m, n

    allocate(UD_dg(udSize))
    allocate(distInds(udSize,2))
    nMinus = nAt - 1
    ele = 0
    do m = 1, nMinus
      do n = m+1, nAt
        ele = ele + 1
        UD_dg(ele) = X_dg(m,n)
        distInds(ele,1) = m
        distInds(ele,2) = n
      end do
    end do

  return
  end subroutine makeUDdgNonAdd


  ! Uses the number of atoms (nAt) and number of triplets (nTri) to build a
  ! matrix
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
    integer :: i, j, k
    double precision :: halfLength, delK, delX, delY, delZ

    halfLength = sideLength/2d0

    ! Find X_dg for the atomic positions in posArray
    do i = 1, nAt
      do j = 1, nAt
        if (i .eq. j) then
          X_dg(i,j) = 0
        else if (i .lt. j) then
          do k = 1, 3 ! Loop over x,y,z components
            delK = posArray(i,k) - posArray(j,k)
            if (delK .gt. halfLength) then
              delK = delK - sideLength ! If diff between kth coords is more
                                       ! than +ve L/2 then use image of atom
                                       ! j that is shifted DOWN along kth axis
                                       ! by L
            else if (delK .lt. -1d0*halfLength) then
              delK = delK + sideLength ! If diff between kth coords is less
                                       ! than -ve L/2 then use image of atom
                                       ! j that is shifted UP along kth axis
                                       ! by L
            end if
            if (k .eq. 1) then
              delX = delK
            else if (k .eq. 2) then
              delY = delK
            else if (k .eq. 3) then
              delZ = delK
            end if
          end do

          X_dg(i,j) = delX**2 + delY**2 + delZ**2
          X_dg(i,j) = (X_dg(i,j))**0.5
          X_dg(i,j) = 1 / X_dg(i,j) ! Convert to inverse distance

        else

          X_dg(i,j) = X_dg(j,i)

        end if
      end do
    end do

  return
  end subroutine makeXdg


integer function getNtripsPerProcFullBox(nDists)
  implicit none
  integer :: nDists, counter, i, j
  integer :: al, be

  counter = 0
  do i = 1, nDists
    al = currentEnergyData%alphaBetaPairs(i,1)
    be = currentEnergyData%alphaBetaPairs(i,2)
    do j = 1, currentPositionData%N_a
      if (be .lt. j) then
        counter = counter + 1
      end if
    end do
  end do

  getNtripsPerProcFullBox = counter

return
end function getNtripsPerProcFullBox


function getTripletEnergiesFullBox(nDists,nTrips) result(uVec)
  implicit none
  integer :: nDists, nTrips, triplet(3)
  integer :: i, j, al, be, counter
  double precision :: U(1), uVec(nTrips)

  counter = 0
  do i = 1, nDists
    al = currentEnergyData%alphaBetaPairs(i,1)
    be = currentEnergyData%alphaBetaPairs(i,2)
    do j = 1, currentPositionData%N_a
      if (be .lt. j) then
        counter = counter + 1
        triplet = (/ al, be, j /)
        call tripletEnergiesNonAdd(triplet,1,N_tp,N_p,nArgs,Perm, &
                                   currentPositionData,expArray, &
                                   alpha,hyperParams(2), &
                                   currentEnergyData, U)
        uVec(counter) = U(1)
      end if
    end do
  end do

return
end function getTripletEnergiesFullBox


integer function getNdistsPerProcFullBox()
  implicit none
  integer :: counter, nLoops, i, j
  logical :: countBack

  counter = 0
  i = 1
  nLoops = 1
  countBack = .false.
  do while (i .le. currentPositionData%N_distances)
    do j = 1, clusterSize
      if (processRank .eq. j-1) then
        if (countBack .eqv. .false.) then
          if (((nLoops-1)*clusterSize)+processRank+1 .le. &
              currentPositionData%N_distances) then
            counter = counter + 1
          end if
        else
          if ((nLoops*clusterSize)-processRank .le. &
              currentPositionData%N_distances) then
            counter = counter + 1
          end if
        end if
      end if
      i = i + 1
      if (j .eq. clusterSize) then
        nLoops = nLoops + 1
        if (countBack .eqv. .false.) then
          countBack = .true.
        else
          countBack = .false.
        end if
      end if
    end do
  end do
  getNdistsPerProcFullBox = counter

return
end function getNdistsPerProcFullBox


function distributeDistances(nDists,allDists) result(procDists)
  implicit none
  integer :: counter, nLoops, i, j, nDists
  logical :: countBack
  double precision :: allDists(currentPositionData%N_distances)
  double precision :: procDists(nDists)

  counter = 0
  i = 1
  nLoops = 1
  countBack = .false.

  do while (i .le. currentPositionData%N_distances)
    do j = 1, clusterSize
      if (processRank .eq. j-1) then
        if (countBack .eqv. .false.) then
          if (((nLoops-1)*clusterSize)+processRank+1 .le. &
              currentPositionData%N_distances) then
            counter = counter + 1
            procDists(counter) = allDists(i)
          end if
        else
          if ((nLoops*clusterSize)-processRank .le. &
              currentPositionData%N_distances) then
            counter = counter + 1
            procDists(counter) = allDists((nLoops*clusterSize)-&
                                          processRank)
          end if
        end if
      end if
      i = i + 1
      if (j .eq. clusterSize) then
        nLoops = nLoops + 1
        if (countBack .eqv. .false.) then
          countBack = .true.
        else
          countBack = .false.
        end if
      end if
    end do
  end do

return
end function distributeDistances


function getAlphaBetaPairs(nDists,Xdg) result(alphaBeta)
  implicit none
  integer :: alpha, beta, counter, nDists, alphaBeta(nDists,2)
  double precision :: Xdg(currentPositionData%N_a,currentPositionData%N_a)

  counter = 1
  do while (counter .le. nDists)
    do alpha = 1, currentPositionData%N_a-1
      do beta = alpha+1, currentPositionData%N_a
        if (currentEnergyData%processDists(counter) .eq. Xdg(alpha,beta)) then
          alphaBeta(counter,1) = alpha
          alphaBeta(counter,2) = beta
          counter = counter + 1
        end if
      end do
    end do
  end do

return
end function getAlphaBetaPairs


end module fullBoxModule
