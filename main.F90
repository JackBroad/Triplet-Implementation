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
  logical :: setSeed=.false., acceptMove=.true., moveFlag=.true.
  double precision :: dist, time, fullEnergy, moveEnergy, check
  double precision :: atomMoveTime, acceptTime, rejectTime
  Character(len=300) :: hyperParametersFile = 'hyperParam.txt'
  Character(len=300) :: alphaFile = 'alpha.txt'
  Character(len=300) :: trainingSetFile = 'trainingSet.txt'
  Character(len=300) :: positionFile = 'AtomicPositions500SL=18.txt'


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

  ! Set step-size for atom-move calc
  dist = 2d0

  ! Calculate energy for full sim box
  fullEnergy = tmpi_calcFullSimBoxEnergy()
  call MPI_BARRIER(MPI_COMM_WORLD, barError)

  ! Atom move
  if (moveFlag .eqv. .true.) then
  do i = 1, 150

    call initialise_Move(dist,setSeed, movedAtom)
    time = MPI_Wtime()
    atomMoveTime = MPI_Wtime()
    call broadcastMoveData()

    moveEnergy = tmpi_calcAtomMoveEnergy(movedAtom)
    atomMoveTime = MPI_Wtime() - atomMoveTime

    if (acceptMove .eqv. .true.) then
      acceptTime = MPI_Wtime()
      call acceptMoveRoutine()
      acceptTime = MPI_Wtime() - acceptTime
      rejectTime = 0d0
    else if (acceptMove .eqv. .false.) then
      rejectTime = MPI_Wtime()
      call rejectMoveRoutine()
      rejectTime = MPI_Wtime() - rejectTime
      acceptTime = 0d0
    end if
    time = MPI_Wtime() - time

    call MPI_BARRIER(MPI_COMM_WORLD, barError)
    if (processRank .eq. root) then
      print *, time, atomMoveTime, setTime, expTime, tripTime, &
               partialSumTime, acceptTime, rejectTime
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, barError)

    call deallocateMoveArrays()

  end do
  end if

  call deallocateArraysGP()

  ! Close shared memory window and MPI
  call MPI_WIN_FREE(win,ierror)
  call MPI_FINALIZE(ierror)


  contains


  subroutine broadcastMoveData()
    implicit none

    call MPI_Bcast(movedAtom, 1, MPI_INT, root, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(proposedPositionData%posArray(movedAtom,:), 3, &
                   MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierror)

  return
  end subroutine broadcastMoveData


  subroutine acceptMoveRoutine()
    implicit none

    call updateCurrentDataStructures(movedAtom)
    acceptMove = .false.
    fullEnergy = fullEnergy + moveEnergy

  return
  end subroutine acceptMoveRoutine


  subroutine rejectMoveRoutine()
    implicit none

    call resetProposedDataStructures(movedAtom)
    acceptMove = .true.

  return
  end subroutine rejectMoveRoutine


  subroutine deallocateMoveArrays()
    implicit none

    deallocate(changeExpData)
    deallocate(oldExpData)
    deallocate(hostDists)
    deallocate(hostIndices)

  return
  end subroutine deallocateMoveArrays


  subroutine deallocateArraysGP()
    implicit none

    deallocate(alpha)
    deallocate(trainData)

  return
  end subroutine deallocateArraysGP


end program main

