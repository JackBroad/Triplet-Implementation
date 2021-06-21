module test
  use triplet_mod
  implicit none
  double precision:: tolerance=1e-15

contains


!@test
subroutine test_init()
  implicit none
  integer :: atom, trip, args, tp, ud
  integer :: trueAtom, trueArgs, trueTP
  double precision, allocatable :: pos(:,:), alp(:), train(:,:)
  double precision :: hypers(3)

  ! Assign known values to true variables
  trueAtom = 400
  trueArgs = 3
  trueTP = 337

  ! Call initialisation subroutine to get values to test
  call initialise(pos, train, alp, hypers, trip, args, atom, &
                  trip, ud)

  ! Check that all values returned by the call to initialise 
  ! match the true values
  call assertEqual(trueAtom, atom, tolerance,'test for &
                   number of atoms')
  if (anyExceptions()) then
    print *, 'The number of atoms was read in incorrectly'
  end if
  call assertEqual(trueArgs, args, tolerance,'test for &
                   number of arguments')
  if (anyExceptions()) then
    print *, 'The number of arguments was read in incorrectly'
  end if
  call assertEqual(trueTP, tp, tolerance,'test for &
                   number of training points')
  if (anyExceptions()) then
    print *, 'The number of TPs was read in incorrectly'
  end if

return
end subroutine test_init


end module test
