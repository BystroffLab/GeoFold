!! Re-canabalize to mcmask2.f90  29-NOV-01
!!-----------------------------------------
!! Brownian Molecular Dynamics using Masks
!! This is a canabolized version of mcmask.f90, where the surface area
!! derivatives from masker2.mod are used.
!! Current version (16-AUG-01) does not do Lengevin heat bath
!! This version should conserve momentum.
!!-----------------------------------------
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
program mcmask
use masker
implicit none
character(len=80) :: xyzfile, aline,outfile,pdbout,dfile,lastcycle
real,dimension(:,:),allocatable :: xyz,xyznew,msderiv,frc,vel,xyzsave
real,dimension(:,:),allocatable :: atomms
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,nfreq,iat,maccept,mcycle,lcycle
integer :: ios=0, nat
real :: tt, sig=0.8,x,y,z,nrg,fa,timestep,dd,ddmax,boltzk,kt,nava,ms,d,r,v,sumd,vf
real,dimension(3) :: vec,vab,uab,dvec
integer :: iargc, time
real(kind=kind_8),dimension(:,:),allocatable ::  shift

!! Readcommand line
tt = 370.  !! degrees kelvin
boltzk = 0.008314517   !! kJ/mol/degK
kt = boltzk*tt
timestep = 10.0 !! ps
ncycle = 1000
nfreq = 100
pdbout = "mcmask"
dfile = " "
lastcycle = "lastcycle.pdb"
! nrg = HUGE(nrg)
!nseed = time()
jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage: xmcmask coordfile boxsize [T N freq VMDoutputfile timestep lastpdb]'
  write(*,*) 'T=',tt,' kT=',kt
  write(*,*) 'N=',ncycle
  write(*,*) 'freq (sampling frequency)=',nfreq
  !write(*,*) 'nseed=',nseed
  write(*,*) 'timestep=',timestep,'*10^(-6)s'
  write(*,*) 'VMD output file name=',trim(pdbout)
  write(*,*) 'lastpdb=',trim(lastcycle)
  stop 'mcmask2.f90 v.16-MAY-06'
endif
call getarg(1,xyzfile)
cAll getarg(2,aline)
read(aline,*,iostat=ios) lbox  ! global
if (jarg > 2) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) tt
  if (ios/=0) stop 'Bad temperature'
  if (tt < 0) stop 'Temperature must be > 0'
  kt = boltzk*tt
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
if (jarg > 6) then
  call getarg(7,aline)
  read(aline,*,iostat=ios) timestep
  if (ios/=0) stop 'Bad timestep'
  if (timestep < 0) stop 'sig must be > 0'
endif
timestep = 0.000001*timestep
if (jarg > 7) then
  call getarg(8,lastcycle)
endif
!!
!! Initialize masks
drawing = .true.
saying = .false.
call initmasks
!! Read coordinates
open(12,file=xyzfile,status='old',form='formatted')
open(13,file=pdbout,status='replace',form='formatted') !! VMD output
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
allocate(xyz(3,nat),stat=ios); if (ios/=0) stop 'Error allocating xyz'
allocate(shift(3,nat),stat=ios); if (ios/=0) stop 'Error allocating shift'
allocate(xyznew(3,nat),stat=ios); if (ios/=0) stop 'Error allocating xyz'
allocate(xyzsave(3,nat),stat=ios); if (ios/=0) stop 'Error allocating xyzsave'
! allocate(atype(nat),stat=ios); if (ios/=0) stop 'Error allocating atype'
! allocate(msderiv(nat,nat),stat=ios)
! if (ios/=0) stop 'Error allocating msderiv'
! allocate(atomms(3,nat),stat=ios)
! if (ios/=0) stop 'Error allocating atomms'
! allocate(vel(3,nat),frc(3,nat),stat=ios)
! if (ios/=0) stop 'Error allocating frc, vel'
write(*,*) nat,' Atoms allocated.'
i = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'Format error xyzfile'
    do j=1,nattype  !! global
      if (aline(14:17)==atomlib(j)%name) exit
    enddo
    atype(i) = j
  endif
enddo
close(12)
!! Temporary, force all atoms to be methanes
atype(1:nat) = 15
!! Move all atoms to the box
if (lbox > 0.) then
  do i=1,nat
     do j=1,3
       if (xyz(j,i) > lbox) xyz(j,i) = xyz(j,i) - lbox*(int(xyz(j,i)/lbox))
       if (xyz(j,i) <= 0.00) xyz(j,i) = xyz(j,i) - lbox*(int(xyz(j,i)/lbox)-1)
     enddo
     !! diagnostic
     write(*,'(a,i4,a,3f8.3,a,i4)') 'Atom ',i,' xyz=',xyz(1:3,i),' type=',atype(i)
  enddo
endif
!!
if (dfile /= " ") then
  open(66,file=dfile,form='formatted',status='replace')
endif
icycle = 0
lcycle = 0
drawing = .false.
iat = 0
xyznew = xyz
xyzsave = xyz
! frc = 0.
! vel = 0.
iat = 0
do while (icycle < ncycle)
  icycle = icycle + 1
  !! add random vectors, scaled by T
  !do iat=1,nat
  !  iat = mod(iat,nat)+1
  !  call ranvec(dvec,sig,nseed)
  !  !! diagnostic
  !  ! write(*,*) 'Ranvec: ',dvec(1:3)
  !  xyznew(1:3,iat) = xyznew(1:3,iat) + dvec(1:3)
  !  ! dd = sqrt(sum(dvec*dvec))
  !enddo
  !! enforce box
  if (lbox > 0.) then
    !! do i=1,nat
      do j=1,3
        if (xyznew(j,iat) > lbox) xyznew(j,iat) = xyznew(j,iat) - lbox
        if (xyznew(j,iat) <= 0.00) xyznew(j,iat) = xyznew(j,iat) + lbox
      enddo
    !! enddo
  endif
  !! get molecular surface area
  if (mod(icycle,nfreq)==0) then
    ! drawing = .false.    !! temporary: dont draw surface, just atoms
    ! write(outfile,'(a,i8.8,".pdb")')  trim(pdbout),icycle   !! icycle up to 99999999
    ! open(13,file=outfile,status='replace',form='formatted')
    !! must use unit=13 here.
    if (lbox > 0.) then
      if (drawing) then
        write(13,'("CRYST1",3f9.3,3f7.2," P 1")') &
        BOXSCALE*lbox,BOXSCALE*lbox,BOXSCALE*lbox,90.,90.,90.
      else
        write(13,'("CRYST1",3f9.3,3f7.2," P 1")') &
        lbox,lbox,lbox,90.,90.,90.
      endif
    endif
    do i=1,nat
      vec = xyznew(1:3,i) 
      if (drawing.and.lbox > 0.) call inbox(vec,lbox)
      write(13,'("ATOM  ",i5,2x,a1,3x,a3,1x,a1,i4,4x,3f8.3,i9,f7.1)') &
                i,atomlib(atype(i))%name(1:1), &
                atomlib(atype(i))%name," ",i,vec(1:3),icycle,tt
    enddo
    write(13,'("TER",/,"ENDMDL")')
    ! if (.not.drawing) close(13)
    !close(13)
  endif
  !! Get surface area energy 
  !! real(kind=kind_8),intent(inout),dimension(3,nat),optional :: gradient

  shift = 0.  !! gradient initialized to zero
  call getms(xyznew,nat,ms,gradient=shift)  !! returns ms in kJ/mol
  !! Get VDW
  v = 0.
  do i=1,nat-1
    do j=i+1,nat
      call boxaround(xyznew(1:3,i),xyznew(1:3,j),lbox,vab)
      d = sqrt(sum(vab*vab))
      vab = vab/d
      !! Write distance data for radial distrance distribution
      if ((mod(icycle,nfreq)==0).and.(dfile /= " ")) then
        write(66,*) d
      endif
      r = atomlib(atype(i))%r + atomlib(atype(j))%r
      v = v + vdwnrg(d,r,vf)  !! vf is deriv in kJ/mol/A
      ! v = v + vdwnrg(d,r)  !! Lennard-Jones in kJ/mol
      shift(1:3,i) = shift(1:3,i) + vf*vab/2
      shift(1:3,j) = shift(1:3,j) - vf*vab/2
    enddo
  enddo
  nrg = v + ms
  !! apply shifts
  sig = tt*timestep
  do i=1,nat
    call ranvec(dvec,sig,nseed)
    xyznew(1:3,i) = xyznew(1:3,i) + dvec(1:3) + timestep*shift(1:3,i)
  enddo
  sumd = sqrt(sum((xyznew(2:3,:)-xyz(1:3,:))*(xyznew(1:3,:)-xyz(1:3,:))/real(nat)))
  !! diagnostic
  write(*,'(i9," rmsshift=",e12.4e2)') icycle, sumd
  xyz(1:3,:) = xyznew(1:3,:)
  !y = ms + v - nrg  ! total energy change, kJ
  !if (y < 0. ) then  !! accept
  !  nrg = ms + v
  !  maccept = maccept + 1
  !  xyz(1:3,:) = xyznew(1:3,:)
  !  !! diagnostic
  !  !!write(*,'(i4,3f7.3,"A dE=",e12.4e2,"  -dE/kT=",e12.4e2)') iat,dvec,y,(-y/kt)
  !else
  !  if (ran(nseed) < exp(-y/kt)) then  !! accept
  !  !! diagnostic
  !  !!write(*,'(i4,3f7.3,"A dE=",e12.4e2,"  -dE/kT=",e12.4e2,$)') iat,dvec,y,(-y/kt)
  !    nrg = ms + v
  !    maccept = maccept + 1
  !    xyz(1:3,:) = xyznew(1:3,:)
  !    !! write(*,*) ' accept'
  !  else  !! reject
  !    xyznew(1:3,:) = xyz(1:3,:)  !! reset, keep nrg
  !  !! diagnostic
  !  !!  write(*,'(i4,3f7.3,"A dE=",e12.4e2,"  -dE/kT=",e12.4e2,$)') iat,dvec,y,(-y/kt)
  !  !!   write(*,*) ' reject'
  !  endif
  !endif 
  if (mod(icycle,nfreq)==0) then
    !z = real(maccept)/(icycle-lcycle)
    !sumd = 0.
    !do j=1,nat
    !  call boxaround(xyz(1:3,j),xyzsave(1:3,j),lbox,vab)
    !  d = sqrt(sum(vab*vab))
    !  sumd = sumd + d
    !enddo
    !sumd = sumd/real(nat)
 !write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",E12.4e2," total=",E12.4e2," accept=",f6.3," shift=",f7.3 )') icycle,v,ms,nrg,z,sumd
    write(*,'("Cycle ",i8," VDW=",E12.4e2," MS=",E12.4e2," total=",E12.4e2," rms=",f9.3 )') icycle,v,ms,nrg,sumd
    !naccept = naccept + maccept
    !maccept = 0
    lcycle = icycle
    !xyzsave = xyz
  endif
enddo
!write(*,'("Fraction of moves accepted: ",f8.5)') real(naccept)/real(icycle)
!if (dfile /= " ") then
!  close(66)
!endif
close(13)
open(14,file=lastcycle,status="replace",form="formatted")
if (lbox > 0.) then
  write(14,'("CRYST1",3f9.3,3f7.2," P 1")') &
    lbox,lbox,lbox,90.,90.,90.
endif
do i=1,nat
  write(14,'("ATOM  ",i5,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3,i9,f7.1)') &
            i,atomlib(atype(i))%name(1:3), &
            atomlib(atype(i))%name(1:3)," ",i,xyz(1:3,i),0,tt
enddo
close(14)
write(0,*) 'mcmask v.16-MAY-06'
stop 
CONTAINS
!----------------------------------------------------!
real function vdwnrg(d,r,vf)
real,intent(in) :: d,r
real,optional :: vf
real,parameter :: eps=1.2264  !  kJ/mol
real,parameter :: vmax=20.
real :: x,y
x = (r/d)**6
y = 4*eps*(x*x - x)
vdwnrg = y
if (present(vf)) then
  vf = 4*eps*((6*x/d) - (12*x*x/d))
  if (abs(vf) > vmax) vf = vmax*vf/abs(vf)
endif
end function vdwnrg
!----------------------------------------------------!
! subroutine ranvec(vec)
!! MOVED TO masker.f90
!! contained in program that uses masker module
! real,dimension(3),intent(out) :: vec
! real :: phi,psi,x,y,z
! integer :: imask
! !! chose random point on the surface of the template mask
! x = ran(nseed)
! imask = nint(x*MAXATOM)
! x = 3*sig*(ran(nseed)**3)
! vec = x*mxyz(1:3,imask)
! end subroutine ranvec
!----------------------------------------------------!
end program mcmask
