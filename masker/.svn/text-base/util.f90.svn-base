program util
use masker
implicit none
integer,dimension(MASKSIZE) :: amask
real :: myr1,myrw,thet,dij,x,y,z
integer :: i,j,k

!! Initialize masks
call initmasks
!! Fri Jan 25 06:47:44 EST 2002
!! Special temporary version of masker.f90 creates ascii masklib file
!!
!! amask = -1
!! open(13,file="fullmask.mas",form='unformatted',status='new')
!! write(13) amask
!! close(13)
myr1 = 2.0
myrw = 1.5
do i=1,80
  dij = i*0.1
  thet = acos(dij/(2*(myr1+myrw)))
  x = 2*pi*myr1*myr1*(1+cos(thet))
  x = x + x
  y = getsaddle(myr1,myrw,thet,1.0)
  y = y + y
  z = x + y
  write(*,*) dij,y,x,z
enddo

end program util
