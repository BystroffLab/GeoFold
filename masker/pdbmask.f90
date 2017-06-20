program pdbmask
use masker
!! INPUT: coordinates in PDB format:
!!  Each atom field must match one of the lines in ATOMLIB
!!  ATOM lines and HETATM lines are accepted.
!! Coordinates are translated to a cubic periodic box of 
!! length lbox.  Symmatry operations are checked
!! when computing the MS.
!! MC moves are random vectors drawn from a Guassian
!! distribution, truncated at 3 sigma.
!! Equation for length of vector: (for now)
!!   (ran()^3)*3*sigma
!! 
implicit none
character(len=200) :: xyzfile, aline,outfile, dockfile,r3dfile,dotfile
real,dimension(:,:),allocatable :: xyz,xyznew
real,dimension(:),allocatable :: bfac
integer,dimension(:),allocatable :: resseq
character(len=3),dimension(:),allocatable :: resname
integer,dimension(:),pointer :: atype_ptr
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,jj,nat=0
real :: tt, sig=0.3,x,y,z,nrg,scl,xmin,xmax,ymin,ymax,zmin,zmax
real,dimension(3) :: vec
real :: dockdist,lbox=0.0
integer :: iargc, time, runit=15, ios=0
logical :: drawingatall = .false., docking=.false., rendering=.false.,sasonly=.true.
xmin = 999.; ymin = 999.; zmin = 999.
xmax = -999.; ymax = -999.; zmax = -999.
masker_saying = .true.
dockdist=6.5
wzs = .false.
!! ---------------------------------------------------------------
!! Read command line
!! ---------------------------------------------------------------
tt = 100.
ncycle = 1
nseed = time()
jarg = iargc()
lbox = 0.00
outfile = " "
dotfile = 'n'
if (jarg < 1) then
  write(*,'("Usage: xpdbmask coordfile [output boxsize proberadius surfdots r3dfile ]")')
  write(*,'(a,f7.2,a)') 'boxsize=',lbox,'   Set to -1 for no box.'
  write(*,'(a,f7.2)') 'surfdots      Draw surface dots to file in PDB format for debugging. n=no dots'
  write(*,'(a,f7.2)') 'proberadius rw=',rw
  call masker_usage()
  stop 'pdbmask.f90 v.16-MAY-11'
endif
call getarg(1,xyzfile)
if (jarg >=2 ) then
  call getarg(2,outfile)
endif
lbox = 0.
if (jarg >=3 ) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) lbox
  if (ios/=0) stop 'Bad value for lbox'
  write(*,'("Periodic boundary: ",f8.2)') lbox
  call masker_setlbox(lbox)
endif
if (jarg >=4 ) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) rw
  if (ios/=0) stop 'Bad value for rw'
  write(*,'("Probe radius: ",f8.2)') rw
endif
if (jarg >=5 ) then
  call getarg(5,aline)
  if (aline(1:1)/='n') then
    write(*,*) "Drawing surface dots"
    drawingatall = .true.
    dotfile = trim(aline)
  endif
endif
if (jarg >=6) then
  call getarg(6,r3dfile)
  rendering = .true.
else
  rendering = .false.
endif
tt = 1.
ncycle = 1
!! ---------------------------------------------------------------
!! Initialize masks
!! ---------------------------------------------------------------
call masker_plot(drawingatall)
call masker_initmasks
!! ---------------------------------------------------------------
!! Read atom coordinates
!! ---------------------------------------------------------------
open(12,file=xyzfile,status='old',form='formatted')
i = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)=='ATOM ') i = i + 1
  if (aline(1:6)=='HETATM') i = i + 1
enddo
nat = i
!! diagnostic
write(*,*) 'Number of atoms:', nat
rewind(12)
allocate(xyz(3,nat),stat=ios)
if (ios/=0) stop 'pdbmask:: Error allocating xyz'
allocate(bfac(nat),stat=ios)
if (ios/=0) stop 'pdbmask:: Error allocating bfac'
allocate(xyznew(3,nat),stat=ios)
if (ios/=0) stop 'pdbmask:: Error allocating xyznew'
allocate(resseq(nat),stat=ios)
if (ios/=0) stop 'pdbmask:: Error allocating resseq'
allocate(resname(nat),stat=ios)
if (ios/=0) stop 'pdbmask:: Error allocating resname'
call masker_allocatoms(nat)
!! ---------------------------------------------------------------
!! Get atoms and atom types
!! ---------------------------------------------------------------
i = 0
vec = 0.
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    jj = 0
    do j=1,nattype
      if (trim(adjustl(atomlib(j)%name(1:4)))==trim(adjustl(aline(13:16)))) then
        jj = j
        exit
      endif
    enddo
    if (jj==0) then
      write(*,'("pdbmask:: Unknown atom type: >>",a4,"<< ignored")') aline(13:16)
      !! diagnostic
      do j=1,nattype;write(*,'(a4,"|",$)') atomlib(j)%name(1:4);enddo;write(*,*)
      !atype(i) = -1   !! dont draw the surface
      cycle
    elseif (aline(14:14)=="H") then
      write(*,'("pdbmask:: Hydrogen ignored >>",a,"<<")') trim(aline)
      cycle
    endif
    !! diagnostic
    ! write(*,'(a,4x,a5,f9.3)') trim(aline),atomlib(jj)%name,atomlib(jj)%r
    i = i + 1
    atype(i) = jj
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'pdbmask:: Format error parsing xyz'
    read(aline(23:26),*,iostat=ios) resseq(i)
    if (ios/=0) stop 'pdbmask:: Format error parsing resseq'
    resname(i) = aline(18:20)
    vec = vec + xyz(1:3,i)
  endif
enddo
nat = i
vec = -vec/real(nat)
!! ---------------------------------------------------------------
!! Enforce periodic boundary conditions, if present.
!! ---------------------------------------------------------------
if (lbox > 0.) then   !! enforce the periodic box, if present
  do i=1,nat
    xyz(1:3,i) = xyz(1:3,i) + vec + lbox/2
    do j=1,3
      if (xyz(j,i) >= lbox) xyz(j,i) = mod(xyz(j,i),lbox)
      if (xyz(j,i) < 0.00) xyz(j,i) = mod(xyz(j,i),lbox) + lbox
    enddo
  enddo
endif
!! ---------------------------------------------------------------
!! Write raster3D file, if requested.
!! ---------------------------------------------------------------
if (rendering) then
  call masker_render(.true.,unit=runit,filename=r3dfile,smooth=.true.)
  !! set raster3D "scale" (meaning windowsize) to be the long dimension
  vec = 0.0
  do i=1,nat
    xmin = min(xmin,xyz(1,i))
    ymin = min(ymin,xyz(2,i))
    zmin = min(zmin,xyz(3,i))
    xmax = max(xmax,xyz(1,i))
    ymax = max(ymax,xyz(2,i))
    zmax = max(zmax,xyz(3,i))
    vec = vec + xyz(1:3,i)
  enddo
  vec = vec/real(-nat)
  xmax = xmax - xmin
  ymax = ymax - ymin
  zmax = zmax - zmin
  scl = max(xmax,ymax,zmax)
  scl = scl + 5.0
  call masker_renderheader(vec,scl)
  !! color by shape: contact, reentrant and toroidal. Other options: "byatom"
  call masker_setcoloring("byshape")
endif
!!
icycle = 0
naccept = 0
nrg = 9999999.
xyznew = xyz
if (lbox > 0.) then
  do i=1,nat
    do j=1,3
      if (xyznew(j,i) > lbox) xyznew(j,i) = xyznew(j,i) - lbox
      if (xyznew(j,i) <= 0.00) xyznew(j,i) = xyznew(j,i) + lbox
    enddo
  enddo
endif
!! ---------------------------------------------------------------
!! Get surface and render it, if requested.
!! ---------------------------------------------------------------
bfac = 0.
if (masker_plotting()) then
  open(13,file=dotfile,status='replace',form='formatted')
endif
if (outfile==" ") then
  call masker_getms(xyznew,nat,x)
else
  call masker_getms(xyznew,nat,x,buried=bfac,SAS=sasonly)
endif
if (masker_plotting())  close(13)
if (rendering) close(runit)
if (outfile/=" ") then
  open(14,file=outfile,status='replace',form='formatted')
  do i=1,nat
    vec = xyznew(1:3,i)
    write(14,'("ATOM  ",i5,2x,a4,a3,1x,a1,i4,4x,3f8.3,2f6.2)') &
              i,atomlib(atype(i))%name(1:4), &
              resname(i)," ",resseq(i),vec,0.0,bfac(i)
  enddo
  close(14)
endif
if (allocated(bfac)) deallocate(bfac)
if (allocated(xyznew)) deallocate(xyznew)
if (allocated(xyz)) deallocate(xyz)
stop 
CONTAINS
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
end program pdbmask
