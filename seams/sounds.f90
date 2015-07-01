
module Sounds

contains

  subroutine Woof()
    print *,'Woof'
  end subroutine

  subroutine Meouw() 
    print *,'Meouw'
  end subroutine

  subroutine MakeSoundTenTimes(soundFunc)
    integer :: i
    interface
      subroutine soundFunc()
      end subroutine
    end interface
    do i = 1, 10 
      call soundFunc()
    enddo
  end subroutine

end module
