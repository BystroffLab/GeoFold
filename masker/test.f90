program test
  character(len=1000) :: aline
  character(len=10),dimension(3) :: alist=(/"ONE       ","TWO       ","THREE     "/)
  write(*,*) junk_var()
contains
  integer function junk_var()
    junk_var=9
  end function junk_var
end program test
