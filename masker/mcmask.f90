program mcmask
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
character(len=80) :: xyzfile, aline,outfile,pdbout,dfile
real,dimension(:,:),allocatable :: xyz,xyznew
integer,dimension(:),allocatable :: atype
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,nfreq,iat,maccept,mcycle,lcycle
real :: tt, sig=0.1,x,y,nrg,fa
real,dimension(3) :: vec
integer :: iargc, time

!! Readcommand line
tt = 100.
ncycle = 1000
nfreq = 100
pdbout = "mcmask"
dfile = " "
nseed = time()
jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage: xmcmask coordfile boxsize [T N freq outputfile]'
  write(*,*) 'T=',tt
  write(*,*) 'N=',ncycle
  write(*,*) 'freq= sampling frequency=',nfreq
  write(*,*) 'nseed=',nseed
  write(*,*) 'base output file name=',trim(pdbout)
  stop 'mcmask.f90 v.21-JUN-01'
endif
call getarg(1,xyzfile)
call getarg(2,aline)
read(aline,*,iostat=ios) lbox
if (jarg > 2) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) tt
  if (ios/=0) stop 'Bad temperature'
  if (tt < 0) stop 'Temperature must be > 0'
endif
if (jarg > 3) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) ncycle
  if (ios/=0) stop 'Bad ncycle'
  if (ncycle < 0) stop 'Ncycle must be > 0'
endif
if (jarg > 4) then
  call getarg(5,aline)
  read(aline,*,iostat=ios) nfreq
  if (ios/=0) stop 'Bad nfreq'
  if (nfreq < 0) stop 'nfreq must be > 0'
endif
if (jarg > 5) then
  call getarg(6,pdbout)
  dfile = trim(pdbout) // ".dist"
endif
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
rewind(12)
allocate(xyz(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
allocate(xyznew(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
allocate(atype(nat),stat=ios)
if (ios/=0) stop 'Error allocating atype'
write(*,*) nat,' Atoms allocated.'
i = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'Format error xyzfile'
    do j=1,nattype
      if (aline(14:17)==atomlib(j)%name) exit
    enddo
    atype(i) = j
    !! diagnostic
    write(*,*) 'Atom ',i,j,xyz(1:3,i)
  endif
enddo
!!
if (dfile /= " ") then
  open(66,file=dfile,form='formatted',status='replace')
endif
icycle = 0
lcycle = 0
naccept = 0
maccept = 0
nrg = 9999.
drawing = .false.
iat = 0
xyznew = xyz
do while (icycle < ncycle)
  icycle = icycle + 1
  !! add random vectors, scaled to 3 sigma
  iat = mod(iat,nat)+1    !! cycle through atoms in order
  !! do i=1,nat
    call ranvec(vec)
    !! diagnostic
    ! write(*,*) 'Ranvec: ',vec(1:3)
    xyznew(1:3,iat) = xyz(1:3,iat) + vec(1:3)
  !! enddo
  !! enforce box
  if (lbox > 0.) then
  !! do i=1,nat
    do j=1,3
      if (xyznew(j,iat) > lbox) xyznew(j,iat) = xyznew(j,iat) - lbox
      if (xyznew(j,iat) <= 0.00) xyznew(j,iat) = xyznew(j,iat) + lbox
    enddo
  !! enddo
  endif
  !! Get VDW
  v = 0.
  do i=1,nat-1
    do j=i+1,nat
      d = distmod(xyznew(1:3,i),xyznew(1:3,j),lbox)
      r = atomlib(atype(i))%r + atomlib(atype(j))%r
      v = v + vdwnrg(d,r)
      if ((mod(icycle,100)==0).and.(dfile /= " ")) then
        write(66,*) d
      endif
    enddo
  enddo
  !! get molecular surface area
  if (mod(icycle,nfreq)==0) then
    drawing = .true.
    write(outfile,'(a,i7.7,".pdb")')  trim(pdbout),icycle
    open(13,file=outfile,status='replace',form='formatted')
    !! must use unit=13 here.
    if (lbox > 0.) then
      write(13,'("CRYST1",3f9.3,3f7.2," P 1")') &
        BOXSCALE*lbox,BOXSCALE*lbox,BOXSCALE*lbox,90.,90.,90.
    endif
    do i=1,nat
      vec = xyznew(1:3,i) 
      if (lbox > 0.) call inbox(vec,lbox)
      write(13,'("ATOM  ",i5,2x,a1,3x,a4,1x,a1,i4,3x,3f8.3,2f6.2)') &
                i,atomlib(atype(i))%name(1:1), &
                atomlib(atype(i))%name," ",i,vec(1:3),0., 0.
    enddo
  endif
  call getms(xyznew,atype,nat,x)
  !! diagnostic
  ! write(*,*) 'Cycle ',icycle,'  MS=',x,'  VDW=',v
  y = x + v - nrg  ! total energy change
  if (y < 0 ) then
    xyz = xyznew
    nrg = x + v
    maccept = maccept + 1
  else
    if (ran(nseed) < exp(-y/tt)) then
      xyz = xyznew
      nrg = x + v
      maccept = maccept + 1
    endif
  endif 
  if (mod(icycle,nfreq)==0) then
    mcycle = icycle - lcycle
    fa = 100*real(maccept)/real(mcycle)
    write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",f9.4," total=",f12.3, " %accept=",f6.1)') icycle,v,x,v+x,fa
    naccept = naccept + maccept
    maccept = 0
    lcycle = icycle
  endif
  if (drawing) then
    close(13)
    drawing = .false.
  endif
enddo
write(*,'("Fraction of moves accepted: ",f8.5)') real(naccept)/real(icycle)
if (dfile /= " ") then
  close(66)
endif
stop 'mcmask v.20-JUL-01'
CONTAINS
!----------------------------------------------------!
real function vdwnrg(d,r)
real,intent(in) :: d,r
real,parameter :: eps=1.2264  !  kJ/mol
real :: x
x = (r/d)**6
y = eps*(x*x - x)
vdwnrg = y
end function vdwnrg
!----------------------------------------------------!
subroutine ranvec(vec)
!! contained in program that uses masker module
real,dimension(3),intent(out) :: vec
real :: phi,psi,x,y,z
integer :: imask
!! chose random point on the surface of the template mask
x = ran(nseed)
imask = nint(x*MAXATOM)
x = 3*sig*(ran(nseed)**3)
vec = x*mxyz(1:3,imask)
end subroutine ranvec
!----------------------------------------------------!
end program mcmask
