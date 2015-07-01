program diffmask
!! Take the difference between two masks
!! mask1 - mask2
!! Fri Jun 15 16:07:32 EDT 2001
!! C.Bystroff
use masker
implicit none
!integer,parameter :: MAXATOM=4096
integer,parameter :: MAXPSI=360
character(len=80) :: aline,pdbfile,binfile,datfile,outfile
integer :: i,j,k,L,jarg,ios,iargc,natm
integer,parameter :: MASKSIZE=MAXATOM/KND
!real,parameter :: DTHETA=3.0   !! degree increments
integer :: MMASK=MAXATOM,NTHETA=30  !! 0-90 deg
integer(kind=KND),dimension(:,:,:),allocatable :: masklib
integer(kind=KND),dimension(MASKSIZE) :: amask
integer :: nn,ipsi,npsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit
integer :: imask1,imask2,itheta1,itheta2
integer,dimension(:),allocatable :: psiposit,nphi
real :: dpsi,dphi
real :: phi,psi,theta
real,parameter :: pi=3.1415927,radius=10.
real,parameter :: rad=pi/180.
i = 0
jarg = iargc()
if (jarg < 4) then
  write(*,*) 'Usage: xdiffmask pbdfile binfile datfile output.pdb'
  write(*,*) 'pdbfile is a set of points at radius = ',radius
  write(*,*) 'binfile is a library of masks from xbinarymask'
  write(*,*) 'datfile is the logfile from xbinarymask'
  write(*,*) 'output.pdb is a pdb file of masked points'
  stop 'diffmask.f90 v. 3-APR-01'
endif
call getarg(1,pdbfile)
call getarg(2,binfile)
call getarg(3,datfile)
call getarg(4,outfile)
!! at this point there should be exactly MAXATOM in xyz()
!! read log file from binarymask
open(11,file=datfile,status='old',form='formatted')
ipsi = 0
npsi = 0
do
  !! the format/keywords for these lines are specifically those
  !! output by xbinarymask (binarymask.f90)
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:4)=="NPSI".and.npsi==0) then
    read(aline(7:),*,iostat=ios) npsi
    if (ios/=0) stop 'Format error 1'
    if (allocated(psiposit)) deallocate(psiposit)
    if (allocated(nphi)) deallocate(nphi)
    allocate(psiposit(0:npsi),nphi(0:npsi),stat=ios)
    if (ios/=0) stop 'Problem allocating psiposit, nphi'
  endif
  if (aline(1:4)=="PSI=".and.npsi/=0) then
    read(aline(19:24),*,iostat=ios) nphi(ipsi)
    if (ios/=0) stop 'Format error 2'
    read(aline(32:38),*,iostat=ios) psiposit(ipsi)
    if (ios/=0) stop 'Format error 3'
    write(*,*) nphi(ipsi), psiposit(ipsi)
    ipsi = ipsi + 1
  endif
  if (aline(1:5)=="NMASK".and.npsi/=0) then
    read(aline(7:),*,iostat=ios) nmask
    if (ios/=0) stop 'Format error 4'
  endif
enddo
close(11)
write(*,'("NPSI =",i7)') npsi
write(*,'("NMASK=",i7)') nmask
allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
if (ios/=0) stop 'Unable to allocate memory for masklib'
!! read binary format mask library, created by xbinarymask
open(12,file=binfile,status='old',form='unformatted')
  read(12) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
close(12)
!!diagnostic
 write(*,*)  masklib(1:MASKSIZE,10,50)
!! prompt for which mask to output
amask = -1
!! do
  write(*,'("MASK1: PSI, PHI, THETA? ",$)')
  read(*,*,iostat=ios) psi,phi,theta
  if (ios/=0) exit
  psi = psi*rad
  phi = phi*rad
  dpsi = DTHETA*rad
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi > npsi) stop 'Bug in psi calculation'
  dphi = 2*pi/real(nphi(ipsi))
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  imask = psiposit(ipsi) + iphi
  itheta = nint(theta/DTHETA)
  write(*,'("ipsi, iphi, itheta, imask:",4i6)') ipsi, iphi, itheta, imask
  if (itheta <= 0 .or. itheta > NTHETA) stop 'Theta out of range'
  itheta1 = itheta; imask1 = imask
  !!
  write(*,'("MASK2: PSI, PHI, THETA? ",$)')
  read(*,*,iostat=ios) psi,phi,theta
  if (ios/=0) exit
  psi = psi*rad
  phi = phi*rad
  dpsi = DTHETA*rad
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi > npsi) stop 'Bug in psi calculation'
  dphi = 2*pi/real(nphi(ipsi))
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  imask = psiposit(ipsi) + iphi
  itheta = nint(theta/DTHETA)
  write(*,'("ipsi, iphi, itheta, imask:",4i6)') ipsi, iphi, itheta, imask
  if (itheta <= 0 .or. itheta > NTHETA) stop 'Theta out of range'
  itheta2 = itheta; imask2 = imask
  !!
  amask = iand(masklib(:,itheta1,imask1),not(masklib(:,itheta2,imask2)))
  !!
!! enddo
!! read the mask and output the dots
open(11,file=pdbfile,status='old',form='formatted')
open(13,file=outfile,status='replace',form='formatted')
write(13,'("REMARK  Mask ",i8," theta=",f8.2)') imask,theta
i = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (i > MAXATOM) stop 'Too many atoms'
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(amask(ibyte),ibit)) then
    write(13,'(a)') trim(aline)
  endif
enddo
close(11)
close(13)
end program diffmask
