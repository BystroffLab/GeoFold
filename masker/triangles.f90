program triangles
!!  Get all triangles for a full viewmask.
!! Triangles are used for rendering by Raster3d
use vectormath, only: dist, cros
use masker
!!
integer,parameter :: MASKSIZE=MAXATOM/NBIT
logical(kind=1),dimension(:),allocatable :: done
integer :: ibyte, ibit, natom,i,j,k,ii,jj,kk,h,m
real :: d,dd
integer,dimension(:),allocatable :: mp
integer,dimension(:),allocatable :: tmp
integer,dimension(:,:),allocatable :: tri
integer :: iargc, jarg, ntri
character(len=200) :: outfile,outputtri,envstr
real, dimension(:,:), pointer :: mxyz
integer(kind=KND),dimension(MASKSIZE) :: vmask
logical :: maker3d=.true.

jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage: xtriangles outputR3Dfile triangles'
  write(*,*) 'Uses environment variables MASKLIB  VMASK MASKTEMPLATE  MASKDAT  ATOMLIB   '
  write(*,*) 'NBIT=',NBIT
  call getenv("MASKLIB",envstr); write(*,'(a,a)') "MASKLIB=",trim(envstr)
  call getenv("VMASK",envstr); write(*,'(a,a)') "VMASK=",trim(envstr)
  call getenv("MASKTEMPLATE",envstr); write(*,'(a,a)') "MASKTEMPLATE=",trim(envstr)
  call getenv("MASKDAT",envstr); write(*,'(a,a)') "MASKDAT=",trim(envstr)
  call getenv("ATOMLIB",envstr); write(*,'(a,a)') "ATOMLIB=",trim(envstr)
  call masker_usage()
  ! stop 'triangles.f90 v.29-AUG-01'
  stop 'triangles.f90 v.  Fri Jun 29 14:09:07 EDT 2007'
endif

call getarg(1,outfile)
call getarg(2,outputtri)
maker3d = .not.((trim(outfile)=="none").or.(trim(outfile)=="None"))

if (maker3d) open(13,file=outfile,status='replace',form='formatted')
open(14,file=outputtri,status='replace',form='formatted')

!! drawing
call masker_plot(.true.)
!! get mask coords 
call masker_initmasks
!!  get private point coordinates
call masker_getmxyz(mxyz)
!! diagnostic
! write(*,*) mxyz
!!  get private VMASK
call masker_getvmask(vmask)
!! get index of viewmask points (mp)
allocate(tmp(5000),stat=ios)
if (ios/=0) stop 'Error allocating tmp()'
j = 0
ntri = 0
do i=1,MAXATOM
  ibyte = (i-1)/NBIT + 1
  ibit = mod(i-1,NBIT)
  if (btest(vmask(ibyte),ibit)) then
    j = j + 1
    tmp(j) = i
  endif
enddo
natom = j
write(*,*) natom,' atoms in viewmask.'
allocate(mp(natom),stat=ios)
if (ios/=0) stop 'Error allocating mp()'
mp(1:natom) = tmp(1:natom)
deallocate(tmp)
allocate(tri(3,10*natom),stat=ios)
if (ios/=0) stop 'Error allocating tri()'
  i = 1
  dd = 999.
  do j=1,natom
    if (i==j) cycle
    d = dist(mxyz(1:3,mp(i)),mxyz(1:3,mp(j)))
    if (d < dd) then
      jj = j
      dd = d
    endif
  enddo
  j = jj
    !! shortest distance is i->j
    !! i->j is a side
    !! get the atom nearest to both
  !! diagnostic
  !! write(*,*) mp(i), ' is nearest to ', mp(j), dd
  k = nearest(i,j)
  ! write(*,*) mp(k), ' is nearest to ', mp(i), mp(j)
  !! i->j->k is a triangle
  !! Now recursively circle i and every neighbor of i until
  !! all points are done.
  ! call circle(i,j,k,j)
  !! all work is done in diamond()
  call diamond(i,j,k)

if (maker3d) write(13,*) 0
if (maker3d) close(13)
close(14)
call masker_deallo()
CONTAINS
!!------------------------------------------------------------
!  real function dist(x,y)
!    real,dimension(3),intent(in) :: x,y
!    real :: d
!    integer :: i
!    d = 0.
!    do i=1,3
!      d = d + (x(i)-y(i))**2
!    enddo
!    dist = sqrt(d)
!  end function dist
!!------------------------------------------------------------

  integer function nearest(i,j)
    implicit none
    !! get nearest point to i and j 
    integer,intent(in) :: i,j
    real :: dd,d
    integer :: k,kk
    dd = 999.
    kk = 0
    do k=1,natom
      if (k==i) cycle
      if (k==j) cycle
      d = dist(mxyz(1:3,mp(i)),mxyz(1:3,mp(k))) &
        + dist(mxyz(1:3,mp(j)),mxyz(1:3,mp(k)))
      if (d < dd) then
        kk = k
        dd = d
      endif
    enddo
    nearest = kk
    !! diagnostic
    !   write(*,*) mp(kk), " is nearest to ", mp(i), mp(j), dd
  end function nearest
!!------------------------------------------------------------
  recursive subroutine diamond(i,j,k)
  integer,intent(in) :: i,j,k
  integer :: ii,jj,kk
  integer :: h
  if (i==j) stop 'triangles:: diamond: i==j'
  if (i==k) stop 'triangles:: diamond: i==k'
  if (j==k) stop 'triangles:: diamond: j==k'
  ii = i; jj = j; kk = k;
  if (onthelist(ii,jj,kk)) return
  call writetri(ii,jj,kk)
  if (.not.twosides(ii,jj)) then
    h = nextnearest(ii,jj,kk)
    call diamond(ii,jj,h)
  endif
  if (.not.twosides(jj,kk)) then
    h = nextnearest(jj,kk,ii)
    call diamond(jj,kk,h)
  endif
  if (.not.twosides(kk,ii)) then
    h = nextnearest(kk,ii,jj)
    call diamond(kk,ii,h)
  endif
  end subroutine diamond
!!------------------------------------------------------------
  logical function twosides(i,j)
  integer, intent(in) :: i,j
  integer :: n,m
  twosides = .true.
  m = 0
  do n=1,ntri
    if ((any(tri(1:3,n)==i)).and.(any(tri(1:3,n)==j))) m = m + 1
  enddo
  if (m > 2) then
    write(0,*) i,j,' have ',m,' sides'
    stop 'Bug! more than two triangles to a side.'
  endif
  if (m == 2) return  !! already two sides
  twosides = .false.
  end function twosides
!!------------------------------------------------------------
  logical function onthelist(i,j,k)
  integer, intent(inout) :: i,j,k
  integer :: n
  onthelist = .true.
  call putinorder(i,j,k)
  do n=1,ntri
    if ((tri(1,n)==i).and.(tri(2,n)==j).and.(tri(3,n)==k)) return
  enddo
  ntri = ntri + 1
  tri(1:3,ntri) = (/i,j,k/)
  onthelist = .false.
  !! diagnostic
  !      write(*,*) "three distances  ",i,j,k, dist(mxyz(1:3,mp(i)),mxyz(1:3,mp(j))), &
  !                        dist(mxyz(1:3,mp(k)),mxyz(1:3,mp(j))), &
  !                        dist(mxyz(1:3,mp(i)),mxyz(1:3,mp(k)))
  end function onthelist
!!------------------------------------------------------------
  integer function nextnearest(i,j,h)
    implicit none
    !! get nearest point to i and j that is not h, and
    !! is on the opposite side of ij from h
    integer,intent(in) :: i,j,h
    real :: dd,d,x
    real,dimension(3) :: vij,vih,vik,vxh,vxk
    integer :: k,kk=0
    dd = 999.
    vij = mxyz(1:3,mp(j)) - mxyz(1:3,mp(i))
    vih = mxyz(1:3,mp(h)) - mxyz(1:3,mp(i))
    call cros(vij,vih,vxh)
    kk = 0
    do k=1,natom
      if (k==i) cycle
      if (k==j) cycle
      if (k==h) cycle
      d = dist(mxyz(1:3,mp(i)),mxyz(1:3,mp(k))) &
        + dist(mxyz(1:3,mp(j)),mxyz(1:3,mp(k)))
      if (d > dd) cycle
      vik = mxyz(1:3,mp(k)) - mxyz(1:3,mp(i))
      call cros(vij,vik,vxk)  !!  cross-product is < 0 if 
      x = sum(vxh*vxk)        !! tetrahedron i,j,k,h is obtuse.
      if (x > 0.) cycle
      if (twosides(i,k)) cycle
      if (twosides(j,k)) cycle
      kk = k
      dd = d
    enddo
    nextnearest = kk
    !! diagnostic
    if (kk==0) then
      write(0,*) 'No nextnearest point found!!!', i,j,h,kk,dd
    elseif (kk==i) then
      write(0,*) 'nextnearest point already in the set!!!', i,j,h,kk,dd
    elseif (kk==j) then
      write(0,*) 'nextnearest point already in the set!!!', i,j,h,kk,dd
    elseif (kk==h) then
      write(0,*) 'nextnearest point already in the set!!!', i,j,h,kk,dd
    else
      ! write(*,*) mp(kk), ' is nearest to ', mp(i), mp(j), dd
      !write(*,'(2i9,3f8.3)') i, mp(i), mxyz(1:3,mp(i))
      !write(*,'(2i9,3f8.3)') j, mp(j), mxyz(1:3,mp(j))
      !write(*,'(2i9,3f8.3)') kk, mp(kk), mxyz(1:3,mp(kk))
    endif
  end function nextnearest
!!------------------------------------------------------------
recursive subroutine circle(i,j,h,m)
  implicit none
  !! circle a point with triangles, then
  !! call circle() for each point on the perimeter
  !! return when a circled (done) point is reached.
  !! Circle point i, starting with side (i,j)
  !! find the nearest point /= h and return when it is == m
  !! NOTE: i /= j /= h,m, but h == m is OK.
  integer,intent(in) :: i,j,h,m
  integer,parameter :: MAXN=12
  integer,dimension(MAXN) :: nbor
  integer :: nn,ii,jj,kk,k,mm

  if (done(i)) return
  nn = 1
  nbor(nn) = j
  k = nextnearest(i,j,h)
  jj = j
  do while (k /= m) ! full circle
    if (.not.done(jj).and..not.done(k)) call writeinorder(i,jj,k)
    nn = nn + 1
    if (nn > MAXN) stop 'Too many neighbors'
    nbor(nn) = k
    kk = jj
    jj = k
    k = nextnearest(i,jj,kk)
  enddo
  if (.not.done(jj).and..not.done(k)) call writeinorder(i,jj,k)
  done(i) = .true.
  nn = nn + 1
  write(*,*) nn, " neighbors of ",i
  if (nn > MAXN) stop 'Too many neighbors'
  nbor(nn) = k
  do ii=1,nn-1
    jj = ii + 1
    mm = nbor(jj)
    write(*,*) 'Circle: ',nbor(ii),nbor(jj),i,mm
    call circle(nbor(ii),nbor(jj),i,mm)
  enddo
  end subroutine circle
!!------------------------------------------------------------
  subroutine writetri(i,j,k)
  integer, intent(in) :: i,j,k
  real :: red=0.2,green=0.3,blue=0.4
  integer :: ii,jj,kk
  ii = i; jj = j; kk = k
  !! diagnostic
  ! write(*,*) ii,jj,kk
  write(14,*) mp(ii),mp(jj),mp(kk)
  if (maker3d)   write(13,*) 1
  if (maker3d)   write(13,'(12f8.3)') mxyz(1:3,mp(ii)),mxyz(1:3,mp(jj)),mxyz(1:3,mp(kk)),red,green,blue
  end subroutine writetri
!!------------------------------------------------------------
  subroutine writeinorder(i,j,k)
  integer, intent(in) :: i,j,k
  real :: red=0.2,green=0.3,blue=0.4
  integer :: ii,jj,kk
  ii = i; jj = j; kk = k
  call putinorder(ii,jj,kk)
  ! write(*,*) ii,jj,kk
  if (maker3d)   write(13,*) 1
  if (maker3d)   write(13,'(12f8.3)') mxyz(1:3,mp(ii)),mxyz(1:3,mp(jj)),mxyz(1:3,mp(kk)),red,green,blue
  end subroutine writeinorder
!!------------------------------------------------------------
  subroutine putinorder(i,j,k)
  integer, intent(inout) :: i,j,k
  integer :: ii,jj,kk

  !! Write out the triplet sorted low to high
  if (i < j ) then
    if (i < k ) then
      ii = i
      if (j < k) then
        jj = j
        kk = k
      else
        jj = k
        kk = j
      endif
    else
      ii = k
      jj = i
      kk = j
    endif
  else
    if (j < k) then
      ii = j
      if (i < k) then
        jj = i
        kk = k
      else
        jj = k
        kk = i
      endif
    else
      ii = k
      jj = j
      kk = i
    endif
  endif
  i = ii; j = jj; k = kk
  end subroutine putinorder
    
end program triangles
