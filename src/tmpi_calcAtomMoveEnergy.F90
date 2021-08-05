module tmpi_calcAtomMoveEnergy_mod
  use mpi_variables
  use triplet_mod
  use GP_variables, only: hyperParams,alpha,Perm,trainData,N_tp,nArgs,N_p
  use energiesData_Module, only: energiesData
  use global_Flags, only: textOutput
  use assert_module
  implicit none
  include 'mpif.h'

  private
  public tmpi_calcAtomMoveEnergy


contains


  subroutine tmpi_calcAtomMoveEnergy(N_move,dist,N_a,N_distances,N_tri,currentEnergyData,posArray, &
                                     proposedEnergies)
    ! Input variables
    integer, intent(in) :: N_a, N_distances, N_move, N_tri
    double precision, intent(in) :: dist
    type( energiesData ), intent(in) :: currentEnergyData

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

       proposedEnergies%expMatrix = currentEnergyData%expMatrix

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

          indj = currentEnergyData%distancesIntMat(newExpInt(1,j),newExpInt(2,j))
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
       call tripletEnergiesNonAdd(scatterTrip,currentEnergyData%distancesIntMat,triPerProc, &
                                  N_tp,N_a,N_p,nArgs,Perm,N_distances,proposedEnergies%expMatrix, &
                                  alpha,hyperParams(2), newUvec)
       call MPI_BARRIER(MPI_COMM_WORLD, barError)
       call MPI_gatherv(newUvec, triPerProc, MPI_DOUBLE_PRECISION, newUfull, scounts, &
                        displs, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

       ! Find total change in non-add energy from moving atom
       if (processRank .eq. root) then

          proposedEnergies%tripletEnergies = currentEnergyData%tripletEnergies
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
       proposedEnergies%distancesIntMat = currentEnergyData%distancesIntMat

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


  subroutine initialAsserts(N_a)
    ! Input variables
    integer, intent(in) :: N_a


    if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    endif

  end subroutine initialAsserts

  subroutine finalAsserts(N_a)
    ! Input variables
    integer, intent(in) :: N_a

     if ( processRank == root ) then
       call assertTrue( N_a>0 , 'tmpi_calcFullSimBoxEnergy argument N_a should be >0')
    endif

  end subroutine finalAsserts
    
  
end module tmpi_calcAtomMoveEnergy_mod
