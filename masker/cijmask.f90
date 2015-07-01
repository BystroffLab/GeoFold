program cijmask
!! Derived from pdbmask.f90  30-mar-05
use masker
!!  Read PDB coordinates.
!! Find all contacts based on the criterion that any two atom
!! have VDW spheres within 1.0A .
!! First, call MASKER with each residue alone, save the SES values.
!! Loop over all contacting residues, call MASKER with each pair
!! of residues, save the SES values.
!! Output i,j and the difference between the SES for the paired
!! amino acids versus the sum of the SES for the aa's alone.
!! This number is the amount uf buried surface that is exposed upon 
!! unfolding, and this should be proportional to the energy of the
!! residue-residue contact.
character(len=80) :: xyzfile, aline,outfile, dockfile
real,dimension(:,:),allocatable :: xyz,xyznew
!!integer,dimension(:),allocatable :: atype
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,jj
real :: tt, sig=0.3,x,y,z,nrg,scl,xmin,xmax,ymin,ymax,zmin,zmax
real,dimension(3) :: vec
real,dimension(:,:) :: dockatom
integer,dimension(:),allocatable :: nposit,reslbl
real :: dockdist
integer :: iargc, time
logical :: drawingatall = .true., docking=.false.
xmin = 999.; ymin = 999.; zmin = 999.
xmax = -999.; ymax = -999.; zmax = -999.

rendering = .false.
smoothing = .true.
drawing = drawingatall
saying = .true.
!! Readcommand line
tt = 100.
jarg = iargc()
lbox = 0.00
rw = 1.4
if (jarg < 1) then
  write(*,'("Usage: xcijmask coordfile [proberadius]")')
  write(*,'(a,f7.2)') 'proberadius rw=',rw
  call usage
  stop 'cijmask.f90 v.13-JAN-03'
endif
call getarg(1,xyzfile)
if (jarg >=2 ) then
  call getarg(2,aline)
  read(aline,*,iostat=ios) rw
  if (ios/=0) stop 'Bad value for rw'
  write(*,'("Probe radius: ",f8.2)') rw
endif
tt = 1.
!! Initialize masks
call initmasks
!! Read coordinates
open(12,file=xyzfile,status='old',form='formatted')
nres = 0
lres = -99
i = 0
allocate(reslbl(1000),nposit(1000),stat=ios); if (ios/=0) stop 'cijmask: error allocating nposit.'
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    i = i + 1
    read(aline(23:26),*,iostat=ios) ires
    if (ires /= lres) then
      nres = nres + 1
      nposit(nres) = i
      reslbl(nres)=ires
      lres = ires
    endif
  endif
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
write(*,*) nat,' Atoms read.'
i = 0
vec = 0.
lres = -99
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
      atype(i) = -1   !! dont draw the surface
    else
      atype(i) = jj
    endif
    !! diagnostic
    !! write(*,*) 'Atom ',i,jj,xyz(1:3,i)
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
!if (rendering) then
!  !! set raster3D "scale" (meaning windowsize) to be the long dimension
!  vec = 0.0
!  do i=1,nat
!    xmin = min(xmin,xyz(1,i))
!    ymin = min(ymin,xyz(2,i))
!    zmin = min(zmin,xyz(3,i))
!    xmax = max(xmax,xyz(1,i))
!    ymax = max(ymax,xyz(2,i))
!    zmax = max(zmax,xyz(3,i))
!    vec = vec + xyz(1:3,i)
!  enddo
!  vec = vec/real(-nat)
!  xmax = xmax - xmin
!  ymax = ymax - ymin
!  zmax = zmax - zmin
!  scl = max(xmax,ymax,zmax)
!  scl = scl + 2.0
!  call renderheader(vec,scl)
!endif
!!
drawing = .false.

do ires=1,nres
  call getms(xyznew,nat,x)
enddo
do ires=1,nres-3
  do jres=ires+3,nres
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
  call getms(xyznew,nat,x)
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
stop 'cijmask.f90 v.  Mar 30 2005'
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
!----------------------------------------------------!
end program cijmask
