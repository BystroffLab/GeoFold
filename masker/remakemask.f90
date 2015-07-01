program makemask
!! THIS VERSION READS IN A STARTING SET OF COORDS
!! MONITORS MINIMUM DISTANCE
!! Fri Jun 15 08:50:16 EDT 2001
!!
!! This program places 'size' points on a sphere of radius 
!! 'radius'. The points are evenly spaced after
!! steepest-descent energy minimization of a repulsive
!! function.
!! CB   17-MAY-01
implicit none
real,dimension(:,:),allocatable :: dmat,xyz,dxyz
character(len=80) :: aline,outfile
integer :: i,j,nc,na,nm,nseed,TBLOCK,ncycle
real :: x,y,z,T,Tstart,Tend,tinc,eps,w(3),v(3),dd(3)
real :: epstart,epsend,dx,radius=10.,nndmax
real,parameter :: pi=3.1415927
integer :: jarg,iargc,isize,ios=0

jarg = iargc()
nseed = 1287876
TBLOCK = 1000
tinc = 0.01
epstart = 1.0
epsend = 0.1
eps = 1.0
Tstart = 20.
Tend = 1.0
ncycle = 1000
if (jarg < 2) then
  write(*,*) 'Usage: xremakemask size name [nseed ncycle epstart epsend] < inputpdbpoints'
  write(*,*) 'This program creates evenly spaced points on a sphere of radius'
  write(*,*) radius
  write(*,*) 'Starting from the coords in stdin'
  stop 'remakemask.f90 v.06-JUN-01'
endif
call getarg(1,aline)
read(aline,*,iostat=ios) isize
if (ios/=0) stop 'bad argument 1'
call getarg(2,outfile)
if (jarg >= 3) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) nseed
  if (ios/=0) stop 'bad argument 3'
endif
if (jarg >= 4) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) ncycle
  if (ios/=0) stop 'bad argument 4'
endif
if (jarg >= 5) then
  call getarg(5,aline)
  read(aline,*,iostat=ios) epstart
  if (ios/=0) stop 'bad argument 5'
endif
if (jarg >= 6) then
  call getarg(6,aline)
  read(aline,*,iostat=ios) epsend
  if (ios/=0) stop 'bad argument 6'
endif
allocate(xyz(3,isize),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
allocate(dxyz(3,isize),stat=ios)
if (ios/=0) stop 'Error allocating dxyz'
allocate(dmat(isize,isize),stat=ios)
if (ios/=0) stop 'Error allocating dmat'
dmat = 0.0
do i=1,isize; dmat(i,i) = 0.; enddo
t = Tstart
xyz = 0.
do i=1,isize
  read(*,'(a)',iostat=ios) aline
  if (ios/=0) stop 'Premature EOF (stdin)'
  read(aline(31:54),'(3f8.3)',iostat=ios) w(1:3)
  if (ios/=0) stop 'Bad PDB format'
  !  do j=1,3
  !    x = ran(nseed)-0.5
  !    w(j) = x
  !  enddo
  x = sqrt(dot_product(w,w))
  w = w/x
  xyz(1:3,i) = w
enddo
y = x
nm = 0
na = 0
nc = 0
nndmax = sqrt((4*pi)/real(isize))
write(*,'("Theoretical maximum nearest neighbor distance =",f8.5)') nndmax
dx = 99999.
do while (dx > 0.0 .and. nc < ncycle)
  z = 0.
  do i=1,isize
    !! i = int(isize*ran(nseed) + 1)
    call gradient(i,w)
    dxyz(1:3,i) = w
    x = sqrt(dot_product(w,w))
    z = max(x,z)
  enddo
  z = nndmax*eps/z
  dx = 0.
  do i=1,isize
    w = xyz(1:3,i) 
    v = w + dxyz(1:3,i)*z
    x = sqrt(dot_product(v,v))
    v = v/x
    xyz(1:3,i) = v
    w = w - v
    dx = dx + dot_product(w,w)
  enddo
  !! diagnostic
  !! write(*,'(i7,3f8.4)') 1,xyz(1:3,1)
  x = NND()
  nc = nc + 1
  write(*,'(i9," cycles NND=",f12.5," eps=",f6.3," dx= ",f12.4)') nc,x,eps,dx
  eps = (real(nc)/real(ncycle))*(epsend-epstart) + epstart
enddo
write(*,'(i9," cycles NND=",f12.3," T= ",f8.3," eps=",f6.3," accepted=",f8.4))') nc,x,T,eps,z
open(11,file=outfile,status='replace',form='formatted')
xyz = radius*xyz
do i=1,isize
  write(11,'("ATOM  ",i5,"  O   HOH ",i5,"    ",3f8.3,2f6.2)') i,i,xyz(1:3,i),0.,0.
enddo
close(11)
CONTAINS    
!----------------------------------------------------------------------------------------------!
subroutine gradient(ix,v)
integer,intent(in) :: ix
real,intent(out) :: v(3)
real :: d,x,w(3),u(3)
integer :: j,xp
real :: xscale
xp = 3
! xscale = 0.1*(nndmax/2)**(xp+1)
xscale = 0.01
v = 0.
do j=1,isize
  if (j == ix) cycle
  d = dist(ix,j)
  ! dmat is distance matrix
  dmat(j,ix) = d
  dmat(ix,j) = d
  if (d == 0.) then
    d = 0.001
  else
    ! inverse distance^3
    x = xscale/(d**xp)
    call cross(xyz(1:3,j),xyz(1:3,ix),w)
    call cross(w,xyz(1:3,ix),u)
    ! u is the tangent vector through ix, pointing away from j
    w = u*(x/sqrt(dot_product(u,u)))
    ! w is u scaled to length x (inverse d^3)
    ! x is 1/2 of the theoretical ideal distance (nndmax) when d=1/2*nndmax
    v = v + w
  endif
  ! x = xscale/(d*d)
enddo
end subroutine gradient
!----------------------------------------------------------------------------------------------!
real function NND
! minimum nearest neighbor distance 
integer :: j,k
real :: d
d = 999.
do j=1,isize-1
  do k=1,j-1
    d = min(d,dmat(j,k))
  enddo
enddo
NND = d
end function NND
!----------------------------------------------------------------------------------------------!
real function dist(ix,j)
integer,intent(in) :: ix,j
real :: d,eeps=0.001
d = (xyz(1,ix) - xyz(1,j))**2 
d = d + (xyz(2,ix) - xyz(2,j))**2 
d = d + (xyz(3,ix) - xyz(3,j))**2 
d = sqrt(d)
!! d = log(d+eeps)
!! angular difference
!! d = acos(dot_product(xyz(1:3,ix),xyz(1:3,j)))
dist = d
end function dist
!----------------------------------------------------------------------------------------------!
subroutine cross(a,b,c)
!** a X b = c, length of c  = d
real :: a(3),b(3),c(3),d
!** cross product A X B
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
!d=dsqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
end subroutine cross
!----------------------------------------------------------------------------------------------!

    


end program makemask
