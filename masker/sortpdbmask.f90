program sortpdbmask
!! This program takes a PDB-formatted files of evenly-spaced points on
!! a sphere (radius='radius'), and converts them to psi and phi
!! values (latitude and longitude). Then it sorts the points
!! first in psi (latitude), by bins of size dpsi degrees, then
!! In each bin, then points are sorted in phi.
!! Tue Mar  6 10:51:32 EST 2001
!! C.Bystroff
implicit none
!! integer,parameter :: MAXATOM=2048,MAXPSI=360
integer :: maxatom, maxpsi
character(len=80) :: aline,pdbfile
integer :: i,j,k,L,jarg,ios,iargc,natm
integer :: nn,ipsi,npsi,itmp
integer,dimension(:),allocatable :: psiposit
integer,dimension(:),allocatable :: rank
real :: dpsi
real, dimension(:,:),allocatable :: xyz
real,dimension(:),allocatable :: phi,psi,sphi,spsi
real,parameter :: pi=3.1415927,radius=10.
jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage: xsortpdbmask dpsi pdbfile > outputfile'
  write(*,*) 'MASK LIBRARY STEP 2'
  write(*,*) 'This program takes the output PDB file of xmakemask'
  write(*,*) 'and sorts the points in psi and phi for input to'
  write(*,*) 'xbinarymask, which in turn creates a mask library.'
  stop 'sortpdbmask.f90 v.06-MAR-01'
endif
call getarg(1,aline)
read(aline,*,iostat=ios) dpsi
maxpsi = int(180./dpsi) + 1
if (ios/=0) stop 'Error reading argument 1'
call getarg(2,pdbfile)
open(11,file=pdbfile,form='formatted',status='old')
i = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (ios/=0) then
    i = i - 1
    cycle
  endif
enddo
rewind(11)
maxatom = i
i = 0
allocate(xyz(3,maxatom),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/="ATOM ") cycle
  i = i + 1
  if (i > maxatom) stop 'Too many atoms'
  read(aline,'(30x,3f8.3)',iostat=ios) xyz(1:3,i)
  if (ios/=0) then  !!  bad format? bug?
    i = i - 1
    cycle
  endif
enddo
close(11)
!! sort 
allocate(phi(maxatom),psi(maxatom),stat=ios)
if (ios/=0) stop 'Error allocating phi,psi'
natm = i
do i=1,natm
  psi(i) = acos(xyz(3,i)/radius)*180./pi
  phi(i) = atan2(xyz(2,i)/radius,xyz(1,i)/radius)*180./pi
enddo
!!
!! sort by psi first
!! rank(I) is the output order
allocate(rank(maxatom),stat=ios)
if (ios/=0) stop 'Error allocating rank'
rank(1) = 1
do I=2,natm
  k = I
  do while (psi(rank(k-1)) > psi(I))
    rank(k) = rank(k-1)
    k = k - 1
    if (k==1) exit
  enddo
  rank(k) = I
enddo
!! find limits of psi bins
allocate(psiposit(maxpsi),stat=ios)
if (ios/=0) stop 'Error allocating psiposit'
k = 1
npsi = 180./dpsi
do ipsi=1,npsi
  do while (psi(rank(k)) < dpsi*(ipsi-1))
    k = k + 1
  enddo
  psiposit(ipsi) = k
enddo
psiposit(npsi+1) = natm + 1
write(*,'("REMARK natm=",i8, " dpsi=",f7.2," npsi=",i5)') natm,dpsi,npsi

!! output Psi bins
do ipsi=1,npsi
  nn = psiposit(ipsi+1) - psiposit(ipsi)
  write(*,'("REMARK PSIshell ",i4,2f8.3,i4,i8)') &
    ipsi,(ipsi-1)*dpsi,ipsi*dpsi,nn,psiposit(ipsi)
enddo

!! diagnostic
!! write(*,'("Done with psi bins")')

!! for each psi bin, sort by phi
do ipsi=1,npsi
  !! write(*,*) ipsi
  do k=psiposit(ipsi)+1,psiposit(ipsi+1)-1
    j = k
    itmp = rank(k)
    if (itmp <= 0) then
      write(*,*) k,itmp
      stop
    endif
    do while (phi(rank(j-1)) > phi(itmp))
      rank(j) = rank(j-1)
      j = j - 1
      if (j == psiposit(ipsi)) exit
    enddo
    rank(j) = itmp
  enddo
enddo
!! diagnostic
!! write(*,'("Done with phi   ")')

!! Output points


do I=1,natm
  j = rank(I)
  ipsi = int(psi(j)/dpsi) + 1
  write(*,'("ATOM  ",i5,"  O   HOH ",i5,"    ",3f8.3,2f6.1)') ipsi,i,xyz(1:3,j),psi(j),phi(j)
enddo

  


end program sortpdbmask
