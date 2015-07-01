program histo
implicit none
integer,dimension(100) :: hist
character(len=80) :: aline
real :: x,binsize,y,z,vol
integer :: i,j,ios=0,iskip,k
integer :: iargc,jarg
jarg = iargc()
binsize = 0.1
iskip = 0
if (jarg > 0) then
  call getarg(1,aline)
  read(aline,*,iostat=ios) binsize
  if (ios/=0) then
    write(*,*) "Usage: xhisto [binsize nskip] < data > histogram"
    stop 'Bad first argument: binsize'
  endif
endif
if (jarg > 1) then
  call getarg(2,aline)
  read(aline,*,iostat=ios) iskip
  if (ios/=0) then
    write(*,*) "Usage: xhisto [binsize nskip] < data > histogram"
    stop 'Bad 2nd argument: nskip'
  endif
endif
hist = 0
j = 0
k = 0
do
  read(*,*,iostat=ios) x
  if (ios/=0) exit
  k = k + 1
  if (k <= iskip) cycle
  if (x < 0) cycle
  i = nint(x/binsize)
  if (i == 0) cycle
  if (i > 100) cycle
  hist(i) = hist(i) + 1
  j = j + 1
enddo
write(*,*) "Done reading ",j," data points."
write(*,*) "Histogram "
y = 0
do i=1,100
  x = i*binsize
  vol = (4/3.)*3.14159*(x**3) - y
  z = hist(i)/vol
  write(*,*) x, hist(i), z
  y = (4/3.)*3.14159*(x**3)
  if (all(hist(i:100)==0)) exit
enddo


end program histo

