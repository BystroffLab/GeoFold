program hb2cij
  use contact_obj

!generate a contact map based only on backbone hydrogen bonds
implicit none
character(len=:),allocatable :: contactfile,outputfile
character(len=1000) :: tmp
type(contact_list) :: contacts
integer :: i

if(iargc() /= 2) then
  write (0,*) "Usage: xhb2cij contactfile outputfile"
  write(0,*) "Number of arguments given",iargc()
  do i = 1, iargc()
    call getarg(i,tmp)
    write(0,*) trim(adjustl(tmp))
  enddo
  stop "hb2cij.f90"
endif

call getarg(1,tmp)
contactfile = trim(adjustl(tmp))
write(*,*) contactfile
call getarg(2,tmp)
outputfile = trim(adjustl(tmp))
write(*,*) outputfile

call read_from_file(contacts,contactfile)
!test
! write(0,*) "TESTING"
! call print(contacts%head_)
! call print(contacts%tail_)
call write_output(contacts,outputfile)
call contacts%cl_dtor()
CONTAINS

!------------------------------------------------------------------------------!
!-------------------- Read Hbonds file into a contact list --------------------!
subroutine read_from_file(contacts,contactfile)
  type(contact_list),intent(inout) :: contacts
  character(len=:),allocatable,intent(in) :: contactfile
  character(len=:),allocatable :: format_str
  character(len=2000) :: aline
  integer :: dunit,ios,res1,res2,tflag
  character(len=3) :: atom1,atom2,flag
  ! type(contact_node),target :: tmp

  contacts = contact_list()
  format_str = '(i7,x,a3,x,i7,x,a3,x,a3)'
  open(newunit=dunit,file=contactfile,status="old",form="formatted",iostat=ios)
  if(ios /= 0) stop "Error opening hydrogen bond contact file"
  do
    read(dunit,'(a)',iostat=ios) aline
    if(ios /= 0) exit
    if(aline(1:1) == "!") cycle
    read(aline,format_str,iostat=ios) res1, atom1, res2, atom2, flag
    if(ios /= 0) stop 'Bad line in Hydrogen bond contact file'
    ! tmp = contact_node()
    tflag = flag2tflag(flag,atom1,atom2)
    ! call tmp%cn_ctor(res1,res2,tflag)
    call contacts%append(res1,res2,tflag)
  enddo
  close(dunit)
endsubroutine read_from_file
!---------- Determine the contact_node flag type from the hbond file ----------!
integer function flag2tflag(flag,atom1,atom2) result(tflag)
  character(len=3),intent(in) :: flag,atom1,atom2
  character(len=:),allocatable :: tmp

  tmp = trim(adjustl(flag))
  if(tmp == 'S') then
    tflag = 1
  elseif(tmp == 'H') then
    tmp = trim(adjustl(atom1))
    if(tmp /= "N" .and. tmp /= "O" .and. tmp /= "OXT") then
      tflag = 3
    else
      tmp = trim(adjustl(atom2))
      if(tmp /= "N" .and. tmp /= "O" .and. tmp /= "OXT") then
        tflag = 3
      else
        tflag = 2
      endif
    endif
  endif
end function flag2tflag
!----------------------------- Write out contacts -----------------------------!
subroutine write_output(contacts,outputfile)
  type(contact_list),intent(in) :: contacts
  character(len=*),intent(in) :: outputfile
  integer :: dunit,ios,count,i
  integer :: res1,res2,tflag

  count = 0
  open(newunit=dunit,file=outputfile,form="formatted",status="replace",action="write",iostat=ios)
  if(ios /= 0) then
    write(0,*) "IOError: Could not open file for writing: ",outputfile
    stop 'hb2cij::write_output'
  endif
  if(contacts%size() == 0) then
    write(0,*) "Warning: Empty contact list"
  else
    do i = 1, contacts%size()
      call contacts%get_indx_values(i,res1,res2,tflag)
      if (tflag /= 2) cycle
      write(dunit,'(2i6,f5.2)') res1,res2,1.0
    enddo
    close(dunit)
  endif
endsubroutine write_output

end program
