program colormasks
!! color masks (vary b-value) by theta
!! mask1 - mask2
!! Fri Jun 15 16:07:32 EDT 2001
!! C.Bystroff
implicit none
integer,parameter :: MAXATOM=4096,MAXPSI=360
character(len=80) :: aline,pdbfile,binfile,datfile,outfile
integer :: i,j,k,L,jarg,ios,iargc,natm
integer,parameter :: MASKSIZE=MAXATOM/32
real :: DTHETA   !! degree increments
integer,parameter :: MMASK=MAXATOM
integer :: NTHETA
integer,dimension(:,:,:),allocatable :: masklib
integer,dimension(MASKSIZE) :: amask
integer :: nn,ipsi,npsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit
integer :: imask1,imask2,itheta1,itheta2
integer,dimension(:),allocatable :: psiposit,nphi
real :: dpsi,dphi,bb,cen(3)
real :: phi,psi,theta
real,dimension(3,MAXATOM) :: mxyz
real,parameter :: pi=3.1415927,radius=10.
real,parameter :: rad=pi/180.
DTHETA=3.75   !! degree increments
NTHETA=int(90/DTHETA)  !! 0-90 deg
i = 0
jarg = iargc()
if (jarg < 4) then
  write(*,*) 'Usage: xcolormasks pbdfile binfile datfile output.pdb'
  write(*,*) 'pdbfile is a set of points at radius = ',radius
  write(*,*) 'binfile is a library of masks from xbinarymask'
  write(*,*) 'datfile is the logfile from xbinarymask'
  write(*,*) 'output.pdb is a pdb file of masked points'
  stop 'colormasks.f90 v. 3-APR-01'
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
    DTHETA = 180./npsi
    NTHETA = int(90/DTHETA)  ! note: DTHETA should be an integer fraction of 90 deg.
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
write(*,'("NPSI =",i7)') npsi
write(*,'("NMASK=",i7)') nmask
close(11)
allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
if (ios/=0) stop 'Unable to allocate memory for masklib'
!! read the mask coordinates
open(11,file=pdbfile,status='old',form='formatted')
i = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (i > MAXATOM) stop 'Too many atoms'
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
  ! if (btest(amask(ibyte),ibit)) then
  !   write(13,'(a)') trim(aline)
  ! endif
enddo
mxyz = mxyz/10.
close(11)
!! read binary format mask library, created by xbinarymask
open(12,file=binfile,status='old',form='unformatted')
read(12,iostat=ios) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
if (ios/=0) stop 'Error reading mask library'
close(12)
!!diagnostic
!! write(*,*)  masklib(1:MASKSIZE,10,50)
!! prompt for which mask to output
amask = -1
!! do
  write(*,'("MASK1: PSI, PHI? ",$)')
  read(*,*,iostat=ios) psi,phi
  if (ios/=0) exit
  psi = psi*rad
  phi = phi*rad
  dpsi = DTHETA*rad
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi > npsi) stop 'Bug in psi calculation'
  dphi = 2*pi/real(nphi(ipsi))
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  imask = psiposit(ipsi) + iphi
bb = 0.
open(13,file=outfile,status='replace',form='formatted')
write(13,'("REMARK  Mask ",i8," psi/phi=",2f8.2)') imask,psi,phi
amask = -1
do itheta=1,NTHETA
  bb = bb + 1.
  ! itheta = nint(theta/DTHETA)
  ! write(*,'("ipsi, iphi, itheta, imask:",4i6)') ipsi, iphi, itheta, imask
  ! if (itheta <= 0 .or. itheta > NTHETA) stop 'Theta out of range'
  itheta1 = itheta; imask1 = imask
  !!
  amask = iand(amask,not(masklib(:,itheta,imask)))
  cen = 0.
  call drawsurface(13,cen,10.,itheta,amask,bb,"C",1)  !! output the sasa surface
  amask = masklib(:,itheta,imask)
  !!
enddo
close(13)
CONTAINS
!!------------------------------------------------------------
subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,iskip)
implicit none
integer,intent(in) :: iunit,iat,iskip
real,dimension(3),intent(in) :: cen
real,intent(in) :: r,bb
real :: d
character(len=1),intent(in) :: ch
integer,dimension(MASKSIZE),intent(in) :: amask
integer :: i,j,ibyte,ibit
real,dimension(3) :: avec,bvec
real,parameter :: DCUT=6.0
!
j = 0
do i=1,MAXATOM
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(amask(ibyte),ibit)) then
    j = j + 1
    if (mod(j,iskip)==0) then
      avec = cen + r*mxyz(1:3,i)
      ! if (ishow /= 0) then
      !   bvec = avec - showvec
      !   if (sqrt(dotprod(bvec,bvec)) > DCUT) cycle
      ! endif
      write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
      i,ch,iat,avec(1:3),0.,bb
    endif
  endif
enddo
end subroutine drawsurface
end program colormasks
