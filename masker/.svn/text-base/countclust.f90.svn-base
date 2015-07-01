program countclust
!! modified to report largest clusetr size too
implicit none
integer :: i,j,nat,ios,nclust,large,L
integer,dimension(:),allocatable :: flag
integer,dimension(:,:),allocatable :: edge
real,dimension(:,:),allocatable :: x
real,parameter :: r1=1.865, rw=1.4
real :: DCUT, lbox
integer :: jarg, iargc
character(len=80) :: pdbfile, aline

DCUT=2*sqrt((r1+rw)**2 - rw*rw)
lbox = 1000.

jarg = iargc()
if (jarg < 1) then
  write(*,*) 'Usage: xcountclust pdbfile [lbox]'
  write(*,*) 'This program counts the number of separate surfaces'
  write(*,*) 'that would be generated from the coordinates given'
  write(*,*) 'a constant radius ',r1,' and a probe radius ',rw
  write(*,*) 'Distance cutoff for connected atoms:',DCUT
  stop 'countclust.f90 v.08-MAR-02'
endif
call getarg(1,pdbfile)
if (jarg > 1) then
  call getarg(2,aline)
  read(aline,*,iostat=ios) lbox
  if (ios/=0) stop 'Bad argument 2 lbox'
endif
!! READ PDB FILE
open(13,file=pdbfile,form='formatted',status='old',iostat=ios)
if (ios/=0) stop 'File is missing'
i = 0
do
  read(13,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/='ATOM ') cycle
  i = i + 1
enddo
nat = i
rewind(13)
allocate(x(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating x'
i = 0
do
  read(13,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/='ATOM ') cycle
  i = i + 1
  read(aline(31:54),'(3f8.3)',iostat=ios) x(1:3,i)
enddo
close(13)
write(*,'(i8," atoms read.")') nat
!! GET edge matrix
allocate(edge(nat,nat),flag(nat),stat=ios)
if (ios/=0) stop 'Error allocating edge'
edge = 0
flag = 0
do i=1,nat
  do j=i+1,nat
    if (dist(x(1:3,i),x(1:3,j)) < DCUT) then
      edge(i,j) = 1
      edge(j,i) = 1
    endif
  enddo
enddo
!! recursive node removal
nclust = 0
do i=1,nat
  if (flag(i)==1) cycle
  flag(i)=1
  L = 1
  call removeclust(i,edge,flag,nat,L)
  nclust = nclust + 1
  large = max(large,L)
enddo
write(*,'("nclust= ",2i8)') nclust,large
if (any(flag==0)) stop 'Bug in the mix'

CONTAINS

real function dist(x,y)
real,dimension(3),intent(in) :: x,y
real :: z, q
integer :: k
z = 0.
do k=1,3
  q = abs(x(k)-y(k))
  if (q > (lbox/2.)) q = lbox - q
  z = z + q*q
enddo
if (z==0.) then
  dist = 0.
else
  dist = sqrt(z)
endif
end function dist

recursive subroutine removeclust(j,edge,flag,nat,L)
implicit none
integer,intent(in) :: j,nat
integer,intent(in),dimension(nat,nat) :: edge
integer,intent(inout),dimension(nat) :: flag
integer,intent(inout) :: L
integer :: i
do i=1,nat
  if (flag(i)==1) cycle
  if (edge(j,i)==1) then
    flag(i) = 1
    L = L + 1
    call removeclust(i,edge,flag,nat,L)
  endif
enddo
end subroutine removeclust


end program countclust
