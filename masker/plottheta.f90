program plottheta
!! For the range of theta, get the analytical SASA
!! and the estimated SASA using the mask.

use masker
!! 
implicit none
character(len=80) :: aline
integer :: i,j,itheta,imask
real :: x,x2,phi,psi,y,z,tt,theta,rms,xx,slope,interc,lastx,lastthet,dpsi
real,dimension(3) :: uvec
integer :: iargc, time, jarg, nmask
logical :: drawingatall = .false.
real,parameter :: pi=3.1415927
real,parameter :: rad=pi/180.

jarg = iargc()
if (jarg < 2) then
  write(*,*) 'plottheta phi psi'
  write(*,*) 'Uses masker module'
  write(*,*) 'Calculates the SLOPES data, slope and intercept for interpolating between'
  write(*,*) 'masks. NOTE: Choice of phi and psi should not matter.'
  stop 'plottheta.f90 v.15-DEC-01'
endif
call getarg(1,aline)
read(aline,*) phi
call getarg(2,aline)
read(aline,*) psi
z = cos(psi*rad)
y = sin(psi*rad)*sin(phi*rad)
x = sin(psi*rad)*cos(phi*rad)
uvec(1) = x
uvec(2) = y
uvec(3) = z

call masker_initmasks
imask = masker_getimask(uvec)
nmask = masker_nmask()
dpsi = masker_dpsi()
!! imask = 25
write(*,*) 'Imask = ',imask
write(*,*) 'nmask = ',nmask


lastx = 1
lastthet = 0.
do itheta=1,NTHETA
  theta = dpsi*itheta
  !! get analytical value for a unit sphere
  x2 = 0.5*(1 + cos(theta))  ! exact surface at theta (fraction of sphere)
  x = 0.
  xx = 0.
  do imask = 1,nmask
    j = masker_maskj(imask,itheta)
    y = real(j)/real(MAXATOM)
    ! write(0,'("theta, imask, fraction=",f9.4,i9,f9.5)') theta,imask,y
    x = x + y  !! estimated surface at theta
    xx = xx + y*y
  enddo
  x = x/nmask
  xx = xx/nmask
  y = xx - x*x
  rms = 0.
  if (y > 0.) rms = sqrt(y)
 !!  get Slope and intercept
  slope = (x2 - lastx)/dpsi     !! square A on a unit sphere per radian
  interc = ((x - lastx)/(x2 - lastx))*dpsi + lastthet  !! point in theta where the estimate is exactly correct
  write(*,'(i4,f10.5,f8.2,f10.6,f12.6,f7.4,2f10.5)') itheta, theta, theta/rad, x2, x, rms, slope, interc
  lastx = x2
  lastthet = theta
enddo
call masker_deallo()
end program plottheta
