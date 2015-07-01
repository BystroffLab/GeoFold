program maskerserver
use masker
implicit none
character(len=200) :: servedpdb="submitted.pdb",surfacefile="junk.pdb"
character(len=200) :: aline,webpage="surf.html"
integer :: i,j,k,ios=0,jarg,iargc,nat
logical :: isthere,isopen,drawingatall=.false.
real,dimension(:,:),allocatable :: xyz
real :: ses
integer :: dtvalues(8),dtmin,dtsec,dtms
real :: before,after

rendering = .false.
lbox = 0.
jarg = iargc()
if (jarg < 3) then
  write(*,*) 'Usage: xmaskerserver servedpdbfile webpage yes/no (whether to draw surface)'
  write(*,*) 'This program starts a deamon to run MASKER and output a web page.'
  stop 'maskerserver.f90 v.Fri Jun  7 08:27:22 EDT 2002'
endif
call getarg(1,servedpdb)
call getarg(2,webpage)
call getarg(3,aline)
if (aline(1:1) == "y") then
  drawingatall = .true.
else
  drawingatall = .false.
endif
drawing = drawingatall

call masker_initmasks    !! this assumes that the proper environment variables have been
                 !! set and the masker files are present. See masker.f90 for details.

do 
  call sleep(1)
  inquire(file=servedpdb,exist=isthere,opened=isopen)
  if (isthere.and..not.isopen) then
    !! write(*,*) 'File exists'
    call getsurface(servedpdb,ses)
    !! close(12,status='delete')
    call writewebpage(webpage,ses,sasa,ssasa,bsasa)
    close(14)
  endif
enddo

CONTAINS
!!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine getsurface(pfile,ses)
implicit none
real,intent(out) :: ses
character(len=*),intent(in) :: pfile
character(len=200) :: aline
real,dimension(:,:),allocatable :: xyz
real,dimension(3) :: vec
integer :: i,j,k,ios,jj ,punit=12
!! integer, global :: nat
real :: x
open(punit,file=pfile,status='old',form='formatted',iostat=ios)
if (ios/=0) then
  write(*,*) 'File is missing'
  ses = -1    !! flag for error
  return
endif
!! Read coordinates
i = 0
ses = 0
do
  read(punit,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)=='PROBE') then
    read(aline(6:),*,iostat=ios) rw
    if (ios/=0) then
      write(*,'("Bad PROBE line:",a)') trim(aline)
      rw = 1.4
    else
      write(*,*) "Probe radius=",rw
    endif
  elseif (aline(1:6)=='PDBOUT') then
    surfacefile = trim(aline(7:))
    if (surfacefile==servedpdb) then
      write(*,*) 'Error: cant recycle pdb.'
      surfacefile = "junk.pdb"
    endif
  elseif (aline(1:7)=='WEBPAGE') then
    webpage = trim(aline(8:))
    if (webpage==servedpdb) then
      write(*,*) 'Error: cant recycle pages.'
      webpage = "surf.html"
    endif
  endif
  if (aline(1:5)=='ATOM ') i = i + 1
  if (aline(1:6)=='HETATM') i = i + 1
enddo
nat = i
!! diagnostic
!! write(*,*) 'Number of atoms:', nat
rewind(punit)
allocate(xyz(3,nat),stat=ios)
if (ios/=0) then
  write(*,*) 'Error allocating xyz'
  ses = -1    !! flag for error
  return
endif
call masker_allocatoms(nat)
vec = 0.
i = 0
do
  read(punit,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) then
      write(*,*) 'Error reading coordinates'
      ses = -1    !! flag for error
      return
    endif
    jj = 0
    do j=1,nattype
      if (aline(14:17)==atomlib(j)%name(1:4)) then
        jj = j
        exit
      endif
    enddo
    if (jj==0) then
      write(*,'("Unknown atom type: ",a4)') aline(14:17)
      atype(i) = 1
    else
      atype(i) = jj
    endif
    !! diagnostic
    !! write(*,*) 'Atom ',i,jj,xyz(1:3,i)
    vec = vec + xyz(1:3,i)
  endif
enddo
close(12,status='delete')
vec = vec/real(nat)
!! move coordinates to the origin
do i=1,nat
  xyz(1:3,i) = xyz(1:3,i) - vec(1:3)
enddo
!! Use this code to get VDW energy
  ! do i=1,nat-1
  !   do j=i+1,nat
  !     d = distmod(xyznew(1:3,i),xyznew(1:3,j),lbox)
  !     r = atomlib(atype(i))%r + atomlib(atype(j))%r
  !     if (d < tooclose) then
  !       write(*,'("Atoms too close: ",i5,3f8.3,i5,3f8.3)') i,xyznew(1:3,i),j,xyznew(1:3,j)
  !     endif
  !     !! diagnostic
  !     !!  write(*,*) 'r=',r,  i, atype(i), atomlib(atype(i))%r,  j,  atype(j), atomlib(atype(j))%r
  !     v = v + vdwnrg(d,r)
  !   enddo
  ! enddo
  !! get molecular surface area
  !! if (mod(icycle,1)==0) then
  !! Set this to true to draw the surface
if (drawingatall) then
  drawing = .true.
  open(13,file=surfacefile,status='replace',form='formatted')
  !! must use unit=13 here.
  if (lbox > 0.) then
    x = lbox*BOXSCALE
    write(13,'("CRYST1",3f9.3,3f7.2," P 1")') x,x,x,90.,90.,90.
  endif
  do i=1,nat
    vec = xyz(1:3,i)
    if (lbox > 0.) then
      vec = vec*BOXSCALE
    endif
    write(13,'("ATOM  ",i5,2x,a1,3x,a3,1x,a1,i4,4x,3f8.3,2f6.2)') &
              i,atomlib(atype(i))%name(1:1), &
              atomlib(atype(i))%name," ",i,vec,0., 0.
  enddo
endif
!! get the SES
CALL DATE_AND_TIME(VALUES=dtvalues)
before = dtvalues(6)*60 + dtvalues(7) + dtvalues(8)/1000.
call masker_getms(xyz,atype,nat,x)
CALL DATE_AND_TIME(VALUES=dtvalues)
after = dtvalues(6)*60 + dtvalues(7) + dtvalues(8)/1000.
ses = sasa + ssasa + bsasa
if (drawing) then
  close(13)
  if (rendering) close(runit)
  drawing = .false.
endif

deallocate(xyz)

end subroutine getsurface
!----------------------------------------------------!
subroutine writewebpage(webpage,ses,s1,s2,s3)
implicit none
real,intent(in) :: ses,s1,s2,s3
character(len=*),intent(in) :: webpage
integer :: ios=0
logical :: isthere
inquire(file=webpage,exist=isthere)
if (isthere) then
  open(14,file=webpage,status='old',position='append',form='formatted',iostat=ios)
else
  open(14,file=webpage,status='new',form='formatted',iostat=ios)
endif
if (ios/=0) stop 'cant open webpage'
write(14,'("<p>Atoms=",i6," SES=",f9.2," contact=",f9.2," toroidal=",f9.2," reentrant=",f9.2," probe_radius=",f6.2)') &
  nat,ses,s1,s2,s3,rw
write(14,'("<br>time elapsed=",f12.3," sec")') (after - before)
close(14)

end subroutine writewebpage
!----------------------------------------------------!
real function vdwnrg(d,r)
real,intent(in) :: d,r
real,parameter :: eps=1.2264
real :: x,y
x = (r/d)**6
y = eps*(x*x - x)
vdwnrg = y
end function vdwnrg
!----------------------------------------------------!
!!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end program maskerserver
