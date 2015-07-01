program makecollars
use masker
character(len=200) :: cfile
!!
call getenv("COLLARS",cfile)
if (cfile==" ") then
  write(*,*) 'Please set the environment variable COLLARS'
  stop
else
  write(*,*) "Creating collars file for masker: ",trim(cfile)
endif
call masker_initcollars()
call masker_initmasks()
call masker_deallo()
end program makecollars

