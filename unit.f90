module unit
    implicit none

contains

  subroutine assertTrue(assertion,msg)
    implicit none
    logical,intent(in) :: assertion
    character (len=*),intent(in) :: msg

    if(assertion) then
      write(0,*) "Unit test passed"
    else
      write(0,'"Unit test failed: ",(a)') msg
      stop 31
    endif
  end subroutine assertTrue


  subroutine assertFalse(assertion,msg)
    implicit none
    logical,intent(in) :: assertion
    character (len=*),intent(in) :: msg

    if(.not. assertion) then
      write(0,*) "Unit test passed"
    else
      write(0,'"Unit test failed: ",(a)') msg
      stop 31
    endif
  end subroutine assertFalse

  subroutine assertEqualInt(int1,int2,msg)
    implicit none
    integer,intent(in) :: int1,int2
    character (len=*),intent(in) :: msg

    if(int1 == int2) then
      write(0,*) "Unit test passed"
    else
      write(0,'"Unit test failed: ",(a)') msg
      stop 31
    endif
  end subroutine assertEqualInt

  subroutine assertEqualReal(num1,num2,msg)
    implict none
    real,intent(in) :: num1,num2
    character (len=*),intent(in) :: msg

    if(num1 == num2) then
      write(0,*) "Unit test passed"
    else
      write(0,'"Unit test failed: ",(a)') msg
      stop 31
    endif
  end subroutine assertEqualReal

  subroutine assertWithinTol(num1,num2,tolerance,msg)
    implicit none
    real,intent(in) :: num1,num2,tolerance
    character (len=*),intent(in) :: msg

    if(abs(num1-num2) <= tolerance) then
      write(0,*) "Unit test passed"
    else
      write(0,'"Unit test failed: ",(a)') msg
      stop 31
    endif
  end subroutine assertWithinTol

end module unit
