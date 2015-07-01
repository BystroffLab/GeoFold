program twomask
!! Modified version of PDBMASK2. This version calls
!! getms twice, resetting the probe radius the second time.
use masker
!! Monte Carlo simulation of spherical atoms
!! in implicit solvent. Molecular surafce area
!! defines the free energy of hydration
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
character(len=80) :: xyzfile, aline,outfile
real,dimension(:,:),allocatable :: xyz,xyznew
!!integer,dimension(:),allocatable :: atype
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,jj
real :: tt, sig=0.3,x,y,z,nrg,scl,xmin,xmax,ymin,ymax,zmin,zmax
real,dimension(3) :: vec
integer :: iargc, time
logical :: drawingatall = .true.
real :: rw2  !! second probe radius
xmin = 999.; ymin = 999.; zmin = 999.
xmax = -999.; ymax = -999.; zmax = -999.

rendering = .false.
drawing = drawingatall
!! Readcommand line
tt = 100.
ncycle = 1000
nseed = time()
jarg = iargc()
lbox = 0.00
if (jarg < 2) then
  write(*,'("Usage: xpdbmask2 coordfile output [boxsize r3dfile r1 r2]")')
  write(*,'(a,f7.2)') 'boxsize=',lbox
  write(*,'(a,f7.2)') 'proberadius r1=',rw
  write(*,'(a,f7.2)') 'proberadius r2=',rw*2
  call usage
  stop 'twomask.f90 v.14-DEC-01'
endif
call getarg(1,xyzfile)
call getarg(2,outfile)
if (outfile == "no") then
  drawingatall = .false.
  drawing = drawingatall
endif
lbox = 0.
if (jarg >=3 ) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) lbox
  if (ios/=0) stop 'Bad value for lbox'
  write(*,'("Periodic boundary: ",f8.2)') lbox
endif
if (jarg >=4 ) then
  call getarg(4,aline)
  if (rendering) then
    open(runit,file=aline,status='replace',form='formatted')
  endif
else
  if (rendering) then
    open(runit,file='junk.r3d',status='replace',form='formatted')
  endif
endif
if (jarg >=5 ) then
  call getarg(5,aline)
  read(aline,*,iostat=ios) rw
  if (ios/=0) stop 'Bad value for rw'
endif
write(*,'("Probe radius: ",f8.2)') rw
rw2 = rw*2
if (jarg >=5 ) then
  call getarg(6,aline)
  read(aline,*,iostat=ios) rw2
  if (ios/=0) stop 'Bad value for rw2'
endif
write(*,'("Probe radius 2: ",f8.2)') rw2
tt = 1.
ncycle = 1
!! Initialize masks
call initmasks
!! Read coordinates
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
if (ios/=0) stop 'Error allocating xyz'
allocate(xyznew(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
!!allocate(atype(nat),stat=ios)
!!if (ios/=0) stop 'Error allocating atype'
call allocatoms(nat)
write(*,*) nat,' Atoms allocated.'
i = 0
vec = 0.
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'Format error xyzfile'
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
    write(*,*) 'Atom ',i,jj,xyz(1:3,i)
    vec = vec + xyz(1:3,i)
  endif
enddo
vec = -vec/real(nat)
if (lbox > 0.) then
  do i=1,nat
    xyz(1:3,i) = xyz(1:3,i) + vec + lbox/2
    do j=1,3
      if (xyz(j,i) >= lbox) xyz(j,i) = mod(xyz(j,i),lbox)
      if (xyz(j,i) < 0.00) xyz(j,i) = mod(xyz(j,i),lbox) + lbox
    enddo
  enddo
endif
if (rendering) then
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
  scl = scl + 2.0
  call renderheader(vec,scl)
endif
!!
icycle = 0
naccept = 0
nrg = 9999999.
drawing = .false.

do while (icycle < ncycle)
  icycle = icycle + 1
  !! add random vectors, scaled to 3 sigma
  !do i=1,nat
  !  call ranvec(vec,sig,nseed)
  !  !! diagnostic
  !  ! write(*,*) 'Ranvec: ',vec(1:3)
  !  xyznew(1:3,i) = xyz(1:3,i) + vec(1:3)
  !enddo
  xyznew = xyz
  !! enforce box
  if (lbox > 0.) then
    do i=1,nat
      do j=1,3
        if (xyznew(j,i) > lbox) xyznew(j,i) = xyznew(j,i) - lbox
        if (xyznew(j,i) <= 0.00) xyznew(j,i) = xyznew(j,i) + lbox
      enddo
    enddo
  endif
  !! Get VDW
  ! v = 0.
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
    ! write(outfile,'(a6,i5.5,".pdb")')  "mcmask",icycle
    open(13,file=outfile,status='replace',form='formatted')
    !! must use unit=13 here.
    if (lbox > 0.) then
      x = lbox*BOXSCALE
      write(13,'("CRYST1",3f9.3,3f7.2," P 1")') x,x,x,90.,90.,90.
    endif
    do i=1,nat
      vec = xyznew(1:3,i)
      if (lbox > 0.) then
        vec = vec*BOXSCALE
      endif
      write(13,'("ATOM  ",i5,2x,a1,3x,a3,1x,a1,i4,4x,3f8.3,2f6.2)') &
                i,atomlib(atype(i))%name(1:1), &
                atomlib(atype(i))%name," ",i,vec,0., 0.
    enddo
  endif
  y = 0.
  do i=0,5
    rw = 1.5 + i*0.5
    call getms(xyznew,atype,nat,x)
    write(*,'("Total SES for rw=",f6.2,":  ",4f12.2)') rw,x,sasa,ssasa,bsasa
    y = y + x
  enddo
  ! call resetradii(rw2)  !! add rw2 to all radii
  ! call getms(xyznew,atype,nat,y)
  write(*,'("Total SES all radii:  ",f12.2)') y
  !! diagnostic
  !! write(*,*) '  MS=',x,'  VDW=',v,'  TOTAL=',x+v
  !y = x + v - nrg  ! total energy change
  !if (y < 0 ) then
  !  xyz = xyznew
  !  nrg = x + v
  !  naccept = naccept + 1
  !else
  !  if (ran(nseed) < exp(-y/tt)) then
  !    xyz = xyznew
  !    nrg = x + v
  !    naccept = naccept + 1
  !  endif
  !endif 
  !if (mod(icycle,10)==0) then
  !  write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",f9.4," total=",f12.3)') icycle,v,x,v+x
  !endif
  if (drawing) then
    close(13)
    if (rendering) close(runit)
    drawing = .false.
  endif
enddo
!! write(*,'("Fraction of moves accepted: ",f8.5)') real(naccept)/real(ncycle)
! stop 'pdbmask2 v.30-JUL-01'
stop 'pdbmask2.f90 v.12-OCT-01'
CONTAINS
!----------------------------------------------------!
subroutine resetradii(rr)
implicit none
real,intent(in) :: rr
integer :: i
do i=1,nattype
  atomlib(i)%r = atomlib(i)%r + rr
enddo
!!
end subroutine resetradii
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
!----------------------------------------------------!
end program twomask
