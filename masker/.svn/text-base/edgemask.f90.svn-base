program edgemask
!! Tue Apr  3 16:22:13 EDT 2001
!! C.Bystroff
implicit none
use masker, only : MAXATOM
integer,parameter :: MAXPSI=360
character(len=80) :: aline,pdbfile,binfile,logfile,outfile,template
character(len=80) :: edgefile
integer :: i,j,k,L,m,jarg,ios,iargc,natm
integer,parameter :: MASKSIZE=MAXATOM/32
real,parameter :: DTHETA=5.0   !! five degree increments
integer,parameter :: MMASK=MAXATOM,NTHETA=18  !! 0-90 deg
integer,dimension(MASKSIZE,NTHETA,MMASK) :: masklib
integer,dimension(MASKSIZE) :: amask,bmask,cmask,dmask,emask
integer :: nn,ipsi,npsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv
integer :: imask2,itheta2
integer,dimension(:),allocatable :: psiposit,nphi
integer,dimension(100) :: saveimask, saveitheta
real,dimension(MMASK) :: maskpsi,maskphi
real :: dpsi,dphi
real :: phi,psi,theta,x,d
real,dimension(3) :: avec,bvec,cvec
real,parameter :: pi=3.1415927,radius=10.
real,parameter :: rad=pi/180.
real :: dotprod
!! INITIALIZE
i = 0
jarg = iargc()
dpsi = DTHETA*rad
if (jarg < 4) then
  write(*,*) 'Usage: xmaskpdb masktemplate masklib mask.log output '
  write(*,*) 'masktemplate is a set of points at radius = ',radius
  write(*,*) 'binfile is a library of masks from xbinarymask'
  write(*,*) 'mask.log is the logfile from xbinarymask'
  write(*,*) 'output is a pdb file of masked SASA points'
  write(*,*) 'edgefile is a pdb file of masked edge (collar) points'
  stop 'edgemask.f90 v.17-APR-01'
endif
call getarg(1,template)
call getarg(2,binfile)
call getarg(3,logfile)
call getarg(4,outfile)
!! at this point there should be exactly MAXATOM in xyz()
!! read log file from binarymask
open(11,file=logfile,status='old',form='formatted')
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
!! read binary format mask library, created by xbinarymask
open(12,file=binfile,status='old',form='unformatted')
  read(12) masklib(1:MASKSIZE,1:NTHETA,1:NMASK)
close(12)
!!diagnostic
!!write(*,*)  masklib(1:MASKSIZE,10,50)
!! get psi and phi for each mask
do ipsi=0,npsi
  psi = ipsi*dpsi
  dphi = 2*pi/real(nphi(ipsi))
  do iphi=1,nphi(ipsi)
    phi = (iphi-1)*dphi
    imask = psiposit(ipsi) + iphi - 1
    maskpsi(imask) = psi
    maskphi(imask) = phi
    !! diagnostic
    write(*,*) imask,psi,phi
  enddo
enddo
!!
amask = -1
bmask = -1
open(11,file=template,status='old',form='formatted')
open(13,file=outfile,status='replace',form='formatted')
nm = 0
do
  write(*,'("PSI, PHI, THETA? ",$)')
  read(*,*,iostat=ios) psi,phi,theta
  if (ios/=0) exit
  psi = psi*rad
  phi = phi*rad
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi > npsi) stop 'Bug in psi calculation'
  dphi = 2*pi/real(nphi(ipsi))
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  imask = psiposit(ipsi) + iphi
  itheta = nint(theta/DTHETA)
  write(*,'("ipsi, iphi, itheta, imask:",4i6)') ipsi, iphi, itheta, imask
  if (itheta <= 0 .or. itheta > NTHETA) stop 'Theta out of range'
  !! read the mask and output the dots
  bmask = iand(masklib(:,itheta,imask),amask)
  if (any(bmask /= amask) ) then
    amask = bmask
    nm = nm + 1
    saveimask(nm) = imask
    saveitheta(nm) = itheta
  endif
  write(13,'("REMARK  Mask ",i8," theta=",f8.2)') imask,theta
enddo
j = 0
do i=1,MAXATOM
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(amask(ibyte),ibit)) then
    j = j + 1
  endif
enddo
x = real(j)/real(MAXATOM)
write(13,'("REMARK SASA = ",f8.3)') x
i = 0
cmask = 0
nn = 0
do k=1,nm
  imask = saveimask(k)
  itheta = saveitheta(k)
  bmask = iand(not(masklib(:,itheta+1,imask)),amask)
  cmask = ior(cmask,bmask)
  j = 0
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
    if (btest(cmask(ibyte),ibit)) then
      j = j + 1
    endif
  enddo
  write(13,'("REMARK edgemask ",3i8," bits")') imask,itheta,j
  if (j > 0) then
    nn = nn + 1
    saveimask(nn) = imask
    saveitheta(nn) = itheta
  endif
enddo
!! draw edgemask
write(*,'("There are ",i3," edges.")') nn
!! get verteces = intersections of edges
nv = 0
dmask = 0
emask = 0
do k=2,nn
  imask = saveimask(k)
  itheta = saveitheta(k)
  do m=1,k-1
    imask2 = saveimask(m)
    itheta2 = saveitheta(m)
    bmask = iand(not(masklib(:,itheta+1,imask)),masklib(:,itheta,imask) )
    bmask = iand(bmask,not(masklib(:,itheta2+1,imask2)) )
    bmask = iand(bmask,masklib(:,itheta2,imask2) )
    bmask = iand(cmask,bmask)  !! bmask is current verteces
    dmask = ior(dmask,bmask)  !! dmask is all verteces
    if (any(bmask /= 0)) then
      !! j = 0
      !! do i=1,MAXATOM
      !!   ibyte = (i-1)/32 + 1
      !!   ibit = mod(i-1,32)
      !!   if (btest(bmask(ibyte),ibit)) then
      !!     j = j + 1
      !!   endif
      !! enddo
      psi = maskpsi(saveimask(j))
      phi = maskphi(saveimask(j))
      avec(3) = cos(psi)
      avec(2) = sin(psi)*sin(phi)
      avec(1) = sin(psi)*cos(phi)
      psi = maskpsi(saveimask(m))
      phi = maskphi(saveimask(m))
      bvec(3) = cos(psi)
      bvec(2) = sin(psi)*sin(phi)
      bvec(1) = sin(psi)*cos(phi)
      call cros(avec,bvec,cvec)
      d = sqrt(dotprod(cvec,cvec))
      psi = acos(cvec(3)/d)
      phi = atan2(cvec(2),cvec(1))
      ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
      if (ipsi > npsi) stop 'Bug in psi calculation 2'
      dphi = 2*pi/real(nphi(ipsi))
      iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
      imask = psiposit(ipsi) + iphi
      itheta = NTHETA
      nv = 0
      emask = iand(masklib(:,itheta,imask),bmask)
      if (any(emask /= 0)) nv = nv + 1
      emask = iand(not(masklib(:,itheta,imask)),bmask)
      if (any(emask /= 0)) nv = nv + 1
      write(*,'("Vertex between edges ",i3," and ",i3," number ",i3)') m,j,nv
      write(13,'("REMARK vertex between ",3i8," bits")') m,k,nv
    endif
  enddo
enddo
rewind(11)
!! Draw all dots, labeled by chain
i = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (i > MAXATOM) stop 'Too many atoms'
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(dmask(ibyte),ibit)) then
    aline(22:22) = "D"
    write(13,'(a)') trim(aline)
  elseif (btest(cmask(ibyte),ibit)) then
    aline(22:22) = "C"
    write(13,'(a)') trim(aline)
  elseif (btest(amask(ibyte),ibit)) then
    aline(22:22) = "A"
    write(13,'(a)') trim(aline)
  endif
enddo
close(11)
close(13)
end program edgemask
