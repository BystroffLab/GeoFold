program collartest
!! Plot the estimated arc length versus the true arc length
!! for various values of theta.
use masker
!! 
character(len=80) :: aline,outfile="junk.pdb"
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,jj
real :: tt, sig=0.3,x,y,z,chi,deg,arcfrac,col
real,dimension(3) :: vec,wivec,wjvec,dvec,wvec,orig
real,dimension(3,3) :: mat
integer(kind=KND),dimension(MASKSIZE) :: amask,bmask
integer :: iargc, time, iwat, ideg, itheta,ibyte,ibit
integer :: imask,jmask,kmask,iat
logical :: drawingatall = .true.
!! real,parameter :: pi=3.14159

drawing = drawingatall
!! Readcommand line
tt = 100.
jarg = iargc()
orig = (/0.,0.,0./)
if (jarg < 1) then
  write(*,'("Usage: xcollartest point [outfile]")')
  stop 'collartest.f90 v.  Mon Jan 21 10:09:59 EST 2002'
endif
call getarg(1,aline)
read(aline,*,iostat=ios) iwat
if (jarg >= 2) then
  call getarg(2,outfile)
endif
if (ios/=0) stop 'Error reading arg 1'
call masker_plot(drawingatall)
!! Initialize masks
call masker_initmasks

wvec = mxyz(1:3,iwat)
y = sqrt(dotprod(wvec,wvec))
wvec = wvec/y
imask = getimask(wvec)
psi = acos(wvec(3))
phi = atan2(wvec(2),wvec(1))

i=0
col = 0.

write(*,'("  deg",$)')
do itheta=2,18,2
  write(*,'(f6.1,2x,$)') real(itheta)*DTHETA
enddo
write(*,*)
write(*,'(i5,$)') i
do itheta=2,18,2
  write(*,'(f8.3,$)') col
enddo
write(*,*)

open(13,file=outfile,status='replace',form='formatted')
do ideg=1,359   !! degrees
  write(*,'(i5,$)') ideg
  deg = real(ideg)
  chi = deg*pi/180.
  !! skip by 9 deg (=4.5 * 2)
  !! get collar for mask(itheta) - mask(itheta+1)
  do itheta=2,18,2
    amask = iand(not(masklib(:,itheta+1,imask)),masklib(:,itheta,imask))
    !! get first bit in mask
    !! get the reference point, first 1-bit
    do iat=1,MAXATOM
      ibyte = (iat-1)/NBIT + 1
      ibit = mod(iat-1,NBIT)
      if (btest(amask(ibyte),ibit)) exit
    enddo
    !! reference point is iat
    wivec = mxyz(1:3,iat)
    wjvec = wvec
    !! mask ij0 plane, everything to the left of 0->iat->iwat
    call cros(wivec,wjvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    !! rotate reference point deg around wvec
    call getrotS(orig,wvec,chi,mat,vec)
    wivec = mxyz(1:3,iat)
    call move_S(wivec,mat,vec)
    call cros(wivec,wjvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    kmask = getimask(dvec)
    if (ideg < 180) then
      bmask = iand(masklib(:,NTHETA,jmask),not(masklib(:,NTHETA,kmask)))
    else
      bmask = ior(masklib(:,NTHETA,jmask),not(masklib(:,NTHETA,kmask)))
    endif
    amask = iand(amask,bmask)
    j = countbits(amask)
    arcfrac = real(j)/real(collarsize(itheta,imask))
    col = arcfrac*360.
    write(*,'(f8.3,$)') col
    if (itheta==2) call drawsurface(13,orig,10.,ideg,amask,deg,"E",1)  !! output the sasa surface
  enddo
  write(*,*)
enddo
stop 'collartest.f90 v.  Mon Jan 21 10:11:28 EST 2002'
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
end program collartest
