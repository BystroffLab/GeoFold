program r3drotate
use vectormath
implicit none
character(len=132) :: aline
real,dimension(3,3) :: mat1,mat2,mat3
real,dimension(3) :: vec1,vec2
integer :: I,J,iargc,jarg,ios
real :: phi,psi,chi,rad,scale,newscale
real :: pi=3.1415927

rad = pi/180.

jarg = iargc()
newscale = 0.
!! write(*,*) jarg
if (jarg < 3) then
  write(*,*) "Usage xr3drot phi, psi, chi (polar angles) [x y z scale]"
  stop 'r3drotate.f90 v. 6-AUG-01'
endif

vec2 = 0.

call getarg(1,aline)
read(aline,*,iostat=ios) phi
if (ios/=0) stop 'Bad arg 1'
call getarg(2,aline)
read(aline,*,iostat=ios) psi
if (ios/=0) stop 'Bad arg 2'
call getarg(3,aline)
read(aline,*,iostat=ios) chi
if (ios/=0) stop 'Bad arg 3'
if (jarg >= 6) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) vec2(1)
  if (ios/=0) stop 'Bad arg 4'
  call getarg(5,aline)
  read(aline,*,iostat=ios) vec2(2)
  if (ios/=0) stop 'Bad arg 5'
  call getarg(6,aline)
  read(aline,*,iostat=ios) vec2(3)
  if (ios/=0) stop 'Bad arg 6'
endif
if (jarg >= 7) then
  call getarg(7,aline)
  read(aline,*,iostat=ios) newscale
  if (ios/=0) stop 'Bad arg 7'
endif


!! Copy the first 12 lines as is.
do i=1,12
  read(*,'(a)',iostat=ios) aline
  if (ios/=0) stop 'Error reading stdin'
   write(*,*) trim(aline)
enddo

!! Read the matrix and vector on the next three lines
!! NOTE: Raster3D uses the transposed rotation matrix !!
do i=1,3
  read(*,'(a)',iostat=ios) aline
  if (ios/=0) stop 'r3drotate.f90:: error reading the rotation matrix from r3d header. Check format of lines 13-15'
  read(aline,*) mat1(I,1),mat1(I,2),mat1(I,3), vec1(I)
enddo
chi = chi * rad
psi = psi * rad
phi = phi * rad

!! Rotate matrix
call getmat_S(phi,psi,chi,mat2)
call MM_S(mat1,mat2,mat3)
call move_s(vec1, mat2, vec2)

!! Output new rotation matrix
do i=1,3
  write(*,'(4f10.4)') mat3(I,1),mat3(I,2),mat3(I,3), vec1(I)
enddo

!! read line 16, translation vector and scale
read(*,*) vec1(1),vec1(2),vec1(3), scale

!! rotate translation vector, increment scale (?)
call move_s(vec1, mat2, vec2)
write(*,*) vec1(1),vec1(2),vec1(3), scale+newscale

do
  read(*,'(a)',iostat=ios) aline
  if (ios/=0) exit
  write(*,*) trim(aline)
enddo
end program r3drotate
