module Assert_Module
  implicit none
  public assertEqual_double, assertLessThanOrEqual_double

contains
  
  subroutine assertEqual_double( actual, expected, tolerance, message)
    double precision, intent(in):: actual, expected, tolerance
    Character(len = *), intent(in):: message

#ifdef ASSERTS
    if(abs(actual - expected) > tolerance) then
       print *,'Failed assertion'
       print *, message
       print *,'Actual  : ',actual
       print *,'Expected: ',expected
       STOP 100
    endif
#endif
  end subroutine assertEqual_double


  subroutine assertLessThanOrEqual_double( lower, higher, tolerance, message)
    double precision, intent(in):: lower, higher, tolerance
    Character(len = *), intent(in):: message

#ifdef ASSERTS
    if( lower-higher > tolerance ) then
       print *,'Failed assertion'
       print *, message
       print *,'Lower : ',lower
       print *,'Higher: ',higher
       STOP 100
    endif
#endif
  end subroutine assertLessThanOrEqual_double

  
  subroutine assertEqual_int( actual, expected, message)
    integer, intent(in):: actual, expected
    Character(len = *), intent(in):: message

#ifdef ASSERTS
    if( actual .NE. expected ) then
       print *,'Failed assertion'
       print *, message
       print *,'Actual  : ',actual
       print *,'Expected: ',expected
       STOP 100
    endif
#endif
  end subroutine assertEqual_int

  
  subroutine assertTrue( expression,  message)
    logical, intent(in):: expression
    Character(len = *), intent(in):: message

#ifdef ASSERTS
    if( .not.expression ) then
       print *,'Failed assertTrue'
       print *, message
       STOP 100
    endif
#endif
  end subroutine assertTrue


  subroutine assertFalse( expression,  message)
    logical, intent(in):: expression
    Character(len = *), intent(in):: message

#ifdef ASSERTS
    if( expression ) then
       print *,'Failed assertFalse'
       print *, message
       STOP 100
    endif
#endif
  end subroutine assertFalse

  
end module Assert_Module
