!! Brownian Molecular Dynamics using Masks
!! This is a canabolized version of mcmask.f90, where the surface area
!! derivatives from masker2.mod are used.
!! Current version (16-AUG-01) does not do Lengevin heat bath
!! This version should conserve momentum.
program bmdmask
use masker2
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
real,dimension(:,:),allocatable :: xyz,xyznew,msderiv,frc,vel
real,dimension(:,:),allocatable :: atomms
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,nfreq,iat,maccept,mcycle,lcycle
real :: tt, sig=0.1,x,y,nrg,fa,timestep,dd,ddmax
real,dimension(3) :: vec,vab,uab
integer :: iargc, time

!! Readcommand line
tt = 100.
timestep = 10.0 !! ps
ncycle = 1000
nfreq = 100
pdbout = "bmdmask"
dfile = " "
nseed = time()
jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage: xbmdmask coordfile boxsize [T N freq outputfile]'
  write(*,*) 'T=',tt
  write(*,*) 'N=',ncycle
  write(*,*) 'freq= sampling frequency=',nfreq
  write(*,*) 'nseed=',nseed
  write(*,*) 'base output file name=',trim(pdbout)
  stop 'bmdmask.f90 v.16-AUG-01'
endif
call getarg(1,xyzfile)
call getarg(2,aline)
read(aline,*,iostat=ios) lbox  ! global
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
call allocatoms(nat)
allocate(xyz(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
allocate(xyznew(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
! allocate(atype(nat),stat=ios)
! if (ios/=0) stop 'Error allocating atype'
allocate(msderiv(nat,nat),stat=ios)
if (ios/=0) stop 'Error allocating msderiv'
allocate(atomms(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating atomms'
allocate(vel(3,nat),frc(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating frc, vel'
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
drawing = .false.
iat = 0
xyznew = xyz
frc = 0.
vel = 0.
do while (icycle < ncycle)
  icycle = icycle + 1
  !! add random vectors, scaled by T
  !  do i=1,nat
  !    call ranvec(vec)
    !! diagnostic
    ! write(*,*) 'Ranvec: ',vec(1:3)
  !    xyznew(1:3,iat) = xyz(1:3,iat) + vec(1:3)
  !  enddo
  !! enforce box
  if (lbox > 0.) then
    do i=1,nat
      do j=1,3
        if (xyz(j,iat) > lbox) xyz(j,iat) = xyz(j,iat) - lbox
        if (xyz(j,iat) <= 0.00) xyz(j,iat) = xyz(j,iat) + lbox
      enddo
    enddo
  endif
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
      vec = xyz(1:3,i) 
      if (lbox > 0.) call inbox(vec,lbox)
      write(13,'("ATOM  ",i5,2x,a1,3x,a3,1x,a1,i4,4x,3f8.3,3f6.2)') &
                i,atomlib(atype(i))%name(1:1), &
                atomlib(atype(i))%name," ",i,vec(1:3),atomms(1:3,i)
    enddo
  endif
  !! Get surface area energy and derivatives 
  call getms(xyz,nat,x,atomms,msderiv)
  !! diagnostic
  ! write(*,*) "MSDERIV==================================="
  ! do i=1,nat
  !   write(*,*) msderiv(:,i)
  ! enddo
  !! Get VDW
  v = 0.
  do i=1,nat-1
    do j=i+1,nat
      call boxaround(xyz(1:3,i),xyz(1:3,j),lbox,vab)
      d = sqrt(dotprod(vab,vab))
      !! Write distance data for radial distrance distribution
      if ((mod(icycle,100)==0).and.(dfile /= " ")) then
        write(66,*) d
      endif
      r = atomlib(atype(i))%r + atomlib(atype(j))%r
      v = v + vdwnrg(d,r,vf)  !! vf is deriv in kJ/mol/A
      uab = vab/d
      frc(1:3,i) = frc(1:3,i) + vf*uab
      frc(1:3,j) = frc(1:3,j) - vf*uab
      !! add MS derivs to forces
      y = msderiv(i,j)+msderiv(j,i)
      frc(1:3,i) = frc(1:3,i) + y*uab
      frc(1:3,j) = frc(1:3,j) - y*uab
    enddo
  enddo
  !! diagnostic
  ! write(*,*) 'Cycle ',icycle,'  MS=',x,'  VDW=',v
  y = x + v - nrg  ! total energy change
  !! add shifts to coords
  ddmax = 0.2
  do i=1,nat
    vab = vel(1:3,i)*timestep + frc(1:3,i)*fscale*timestep*timestep/(2*atomlib(atype(i))%m)
    dd = sqrt(dotprod(vab,vab))
    if (dd > ddmax) then
      ddmax = dd
      write(*,*) "Big shift for ",i,dd
      write(*,*) "    m=",atomlib(atype(i))%m
      write(*,*) "    frc=",frc(1:3,i)
    endif
    xyznew(1:3,i) = xyz(1:3,i) + vab
    vel(1:3,i) = vel(1:3,i) + frc(1:3,i)*fscale*timestep/(atomlib(atype(i))%m)
    !! diagnostic
    ! write(*,*) "i,xyz,xyznew,vel ",i,xyz(1:3,i),xyznew(1:3,i),frc(1:3,i)
  enddo
  !! !! Check for big shifts. Scale shifts back to 0.2A max.
  !! if (ddmax > 0.2) then
    !! write(*,'("Cycle ",i9," LARGEST SHIFT = ",f9.2)') icycle, ddmax
    !! dd = 0.2/ddmax
    !! do i=1,nat
      !! vab = vel(1:3,i)*timestep + frc(1:3,i)*fscale*timestep*timestep/(2*atomlib(atype(i))%m)
      !! xyznew(1:3,i) = xyz(1:3,i) + dd*vab
      !! vel(1:3,i) = vel(1:3,i) + frc(1:3,i)*fscale*timestep/(atomlib(atype(i))%m)
    !! enddo
  !! endif
  xyz(1:3,:) = xyznew(1:3,:)
  frc = 0.
  !! if (y < 0 ) then
  !!  xyz = xyznew
    nrg = x + v
  !!   maccept = maccept + 1
  !! else
  !!   if (ran(nseed) < exp(-y/tt)) then
  !!     xyz = xyznew
  !!     nrg = x + v
  !!     maccept = maccept + 1
  !!   endif
  !! endif 
  if (mod(icycle,nfreq)==0) then
    !! mcycle = icycle - lcycle
    !! fa = 100*real(maccept)/real(mcycle)
    !! write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",f9.4," total=",f12.3, " %accept=",f6.1)') icycle,v,x,v+x,fa
    write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",f9.4," total=",f12.3 )') icycle,v,x,nrg
    !! naccept = naccept + maccept
    !! maccept = 0
    !! lcycle = icycle
  endif
  if (drawing) then
    close(13)
    drawing = .false.
  endif
enddo
!! write(*,'("Fraction of moves accepted: ",f8.5)') real(naccept)/real(icycle)
if (dfile /= " ") then
  close(66)
endif
stop 'bmdmask v.16-AUG-01'
CONTAINS
!----------------------------------------------------!
real function vdwnrg(d,r,vf)
real,intent(in) :: d,r
real,intent(out) :: vf
real,parameter :: eps=1.2264  !  kJ/mol
real,parameter :: vmax=20.
real :: x
x = (r/d)**6
y = eps*(x*x - x)
vdwnrg = y
vf = eps*((6*x/d) - (12*x*x/d))
if (abs(vf) > vmax) vf = vmax*vf/abs(vf)
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
end program bmdmask
