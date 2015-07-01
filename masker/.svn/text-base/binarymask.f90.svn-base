program binarymask
!!-------- MASKER part 3 --------------
!! Tue Mar  6 10:51:32 EST 2001
!! C.Bystroff
!!------------------------------------
use masker, ONLY: MAXATOM, NBIT,KND
implicit none
! integer,parameter :: MAXATOM=1024
integer,parameter :: MAXPSI=360
character(len=80) :: aline,pdbfile,binfile,datfile
integer :: i,j,k,L,jarg,ios,iargc,natm
integer,parameter :: MASKSIZE=MAXATOM/NBIT
real :: DTHETA=4.50   !! degree increments
integer,parameter :: MMASK=MAXATOM
integer :: NTHETA  !! 0-90 deg
!! integer,dimension(MASKSIZE,NTHETA,MMASK) :: masklib
integer(kind=KND),dimension(:,:,:),allocatable :: masklib
integer :: nn,ipsi,npsi,itmp,nmask,iphi,nphi
integer :: imask,ipt,npt,ibit,ibyte
real,dimension(3) :: apoint
real :: dpsi,theta,dphi,dang
real, dimension(3,MAXATOM) :: xyz
real :: phi,psi
real,parameter :: pi=3.1415927,radius=10.
real,parameter :: rad=pi/180.
!real :: dotprod ! now contained
!!=======================
!! defaults
DTHETA=4.50
NTHETA=int(90/DTHETA)  !! 0-90 deg
i = 0
jarg = iargc()
if (jarg < 3) then
  write(*,*) 'Usage: xbinarymask points.pdb maskdata datfile [dtheta]'
  write(*,*) 'MASK LIBRARY STEP 3'
  write(*,*) 'This program generates binary masks for a set'
  write(*,*) 'of evenly-spaced points on a sphere, generated'
  write(*,*) 'by xmakemask, and sorted by xsortmask.'
  write(*,*) 'Current settings:'
  write(*,*) 'MAXATOM=',MAXATOM
  write(*,*) 'MAXPSI=',MAXPSI
  write(*,*) 'DTHETA=',DTHETA
  write(*,*) 'NTHETA=',NTHETA
  write(*,*) 'KND=',KND
  write(*,*) 'MASKSIZE=',MASKSIZE,' bytes'
  write(0,*) 'binarymask.f90 v.27-NOV-01'
  stop 
endif
call getarg(1,pdbfile)
call getarg(2,binfile)
call getarg(3,datfile)
if (jarg >= 4) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) DTHETA
  if (ios/=0) stop 'Bad value for dtheta'
  NTHETA=int(90/DTHETA)  !! 0-90 deg
  write(*,*) 'DTHETA=',DTHETA
  write(*,*) 'NTHETA=',NTHETA
endif
  
open(11,file=pdbfile,status='old',form='formatted')
i = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (i > MAXATOM) stop 'Too many atoms'
  read(aline,'(30x,3f8.3)',iostat=ios) xyz(1:3,i)
  if (ios/=0) then
    i = i - 1
    cycle
  endif
enddo
close(11)
xyz = xyz/10.   !! coords from makemask are on a 10A sphere.
! npt = MAXATOM
npt = i
write(*,*) "There are ", npt, " points."
!! at this point there should be exactly MAXATOM in xyz()
dpsi = DTHETA*rad
npsi = int(pi/dpsi)  !! DPSI is the same angle as DTHETA
open(14,file=datfile,status='replace',form='formatted')
write(14,'("NPSI=",i6)') npsi
write(*,'("NPSI=",i6)') npsi
imask = 0
do ipsi=0,npsi
  psi = ipsi*dpsi
  if (ipsi == 0 .or. ipsi == npsi) then
    nphi = 1
  else
    dphi = dpsi/sin(psi)
    nphi = nint(2*pi/dphi)
  endif
  dphi = 2*pi/real(nphi)
  write(14,'("PSI=",f8.2," NPHI=",i6," Mask#=",i7)') psi/rad,nphi, imask+1
  imask = imask + nphi
enddo
nmask = imask
write(*,'("NMASK=",i7)') nmask
write(14,'("NMASK=",i7)') nmask
close(14)
allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
if (ios/=0) stop 'Error allocating masklib'
imask = 0
do ipsi=0,npsi
  psi = ipsi*dpsi
  apoint(3) = cos(psi)
  if (ipsi == 0 .or. ipsi == npsi) then
    nphi = 1
  else
    dphi = dpsi/sin(psi)
    nphi = nint(2*pi/dphi)
  endif
  dphi = 2*pi/real(nphi)
  write(*,'("PSI=",f8.2," NPHI=",i6," Mask#=",i7)') psi/rad,nphi, imask+1
  do iphi = 1,nphi
    phi = (iphi-1)*dphi
    apoint(1) = sin(psi)*cos(phi)
    apoint(2) = sin(psi)*sin(phi)
    imask = imask + 1
    if (imask > nmask) then
      write(*,*) 'BUG!! nmask exceeded. nmask=',nmask,' imask=',imask
      stop 'BUG nmask exceeded.'
    endif
    do i=1,NTHETA
      theta = i*DTHETA*rad
      !! initialize mask
      masklib(1:MASKSIZE,i,imask) = 0
      do ipt=1,npt
        dang  = acos(dotprod(apoint(1:3),xyz(1:3,ipt)))
        if (dang > theta) then
          ibyte = (ipt-1)/NBIT + 1
          ibit = mod(ipt-1,NBIT)
          masklib(ibyte,i,imask) = IBSET(masklib(ibyte,i,imask),ibit)
        endif
      enddo
      !! done making a mask
    enddo
  enddo
enddo
write(*,'("NMASK=",i7)') imask
open(12,file=binfile,status='replace',form='unformatted')
write(*,'("Writing ",a)') trim(binfile)
write(12) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
close(12)
CONTAINS
!****** PROTEAN_MATH routines for PROTEAN
!* Modified for F90. *x converted to real kind
!******************************************************************
        real function dotprod(v1,v2)
        implicit none
        real v1(3),v2(3)
        dotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
        return
        end function dotprod
!***********************************************************************
end program binarymask
