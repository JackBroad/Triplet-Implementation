program main
  use mpi
  use mpi_variables
  use expShare_variables
  use dataStructure_variables
  use time_variables
  use GP_Variables
  use triplet_mod
  use tmpi_calcFullSimBoxEnergy_mod, only: tmpi_calcFullSimBoxEnergy
  use tmpi_calcAtomMoveEnergy_mod, only: tmpi_calcAtomMoveEnergy
  use toyMove_Module, only: toyMove, toyMoveDistScatter, toyMoveMinimalScatter
  use energiesData_Module, only: energiesData
  use positionData_Module, only: positionData
  use updateData
  use assert_module
  use initialise_Module
  implicit none


  integer :: movedAtom, i
  logical :: setSeed=.false., acceptMove=.true., useToyCode=.false.
  logical :: moveFlag=.true.
  double precision :: dist, time, fullEnergy, moveEnergy, check
  double precision :: atomMoveTime, acceptTime, rejectTime
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  Character(len=300) :: positionFile = 'AtomicPositions500.txt'

  ! Set up MPI 
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, clusterSize, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, processRank, ierror)


  ! Set up shared memory
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
                           0, MPI_INFO_NULL, hostComm, ierror)
  call MPI_COMM_SIZE(hostComm, sharedSize, ierror)
  call MPI_COMM_RANK(hostComm, hostRank, ierror)


  ! Set-up calls for full calc
  call initialise_GP_NonAdd(hyperParametersFile, alphaFile, trainingSetFile)
  currentPositionData = initialise_positionDataStruct(positionFile)


  dist = 1.5d0


  if (useToyCode .eqv. .false.) then

    ! Calculate energy for full sim box
    fullEnergy = tmpi_calcFullSimBoxEnergy()
    call MPI_BARRIER(MPI_COMM_WORLD, barError)


    ! Atom move
    if (moveFlag .eqv. .true.) then
    do i = 1, 150
      call initialise_Move(dist,setSeed, movedAtom)
      time = MPI_Wtime()
      atomMoveTime = MPI_Wtime()
      call MPI_Bcast(movedAtom, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(proposedPositionData%posArray(movedAtom,:), 3, MPI_DOUBLE_PRECISION, &
                     root, MPI_COMM_WORLD, ierror)

      !***********True move code***********
      moveEnergy = tmpi_calcAtomMoveEnergy(movedAtom)
      atomMoveTime = MPI_Wtime() - atomMoveTime
      !if (processRank .eq. root) then
      !  check = -79272925.072343603d0
      !  call assertEqual_double(fullEnergy, check, 1d0*10d0**(-5d0), 'incorrect energy found')
      !end if

      if (acceptMove .eqv. .true.) then
        acceptTime = MPI_Wtime()
        call updateCurrentDataStructures(movedAtom)
        fullEnergy = fullEnergy + moveEnergy
        acceptTime = MPI_Wtime() - acceptTime
        rejectTime = 0d0
        acceptMove = .false.
      else if (acceptMove .eqv. .false.) then
        rejectTime = MPI_Wtime()
        call resetProposedDataStructures(movedAtom)
        rejectTime = MPI_Wtime() - rejectTime
        acceptTime = 0d0
        acceptMove = .true.
      end if
      time = MPI_Wtime() - time
      if (processRank .eq. root) then
        print *, time, atomMoveTime, setTime, expTime, tripTime, partialSumTime, acceptTime, rejectTime
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, barError)
      deallocate(changeExpData)
      deallocate(oldExpData)
      deallocate(hostDists)
      deallocate(hostIndices)

    end do
    end if

    deallocate(alpha)
    deallocate(trainData)

  else

    if (processRank .eq. root) then
      print *, 1, 1, 1, 1, 1, 1, 1, 1! Keeps output in format the b/marking script expects
    end if

    do i = 1, 500

      !***********Toy move code***********
      !proposedEnergyData = toyMove(currentPositionData)
      !proposedEnergyData = toyMoveDistScatter(currentPositionData)
      proposedEnergyData = toyMoveMinimalScatter(currentPositionData)

    end do

    deallocate(alpha)
    deallocate(trainData)

  end if

  call MPI_WIN_FREE(win,ierror)
  call MPI_FINALIZE(ierror)

  contains

  double precision function checkSimBoxEnergy()
    implicit none
    integer :: i
    double precision, allocatable :: tripDist(:,:), tripEnergies_Ex(:)

    allocate(tripDist(3,currentPositionData%N_tri))
    allocate(tripEnergies_Ex(currentPositionData%N_tri))

    call findTripletDistances(currentPositionData%N_a,currentPositionData%N_tri, &
                              proposedEnergyData%triMat,proposedEnergyData%interatomicDistances, &
                              tripDist)

    do i = 1, currentPositionData%N_tri
      tripEnergies_Ex(i) = energyCheckCalc(tripDist(:,i)) ! Energies of each triplet
    end do

    checkSimBoxEnergy = sum(tripEnergies_Ex)

  return
  end function checkSimBoxEnergy

end program main

