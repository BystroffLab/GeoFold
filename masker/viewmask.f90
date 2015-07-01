program viewmask
  !! Thu May  3 17:34:42 EDT 2007
  !! C.Bystroff
  use vectormath, ONLY: dist
  use masker, ONLY: MAXATOM,NBIT,KND,masker_drawsurface,masker_initmasks, masker_plot, masker_deallo,masker_getmxyz
  !!
  implicit none
  integer,parameter :: MASKSIZE=MAXATOM/NBIT
  character(len=80) :: aline,template,pdbfile,binfile,datfile,outfile
  integer :: i,j,k,L,jarg,ios,iargc,natm,npsi,nmask
  integer(kind=KND),dimension(MASKSIZE) :: amask
  integer :: nn,ipsi,itmp,itheta,iphi,imask,ibyte,ibit
  integer :: nat,vat,iat
  real :: d,dmin
  integer,dimension(:),allocatable :: psiposit,nphi
  real, dimension(:,:),pointer :: mxyz,vxyz
  real, dimension(3) :: cen
  real :: dpsi,dphi
  real :: phi,psi,theta,r
  real,parameter :: pi=3.1415927,radius=10.
  real,parameter :: rad=pi/180.
  logical :: drawingatall=.true.
  !!
  i = 0
  jarg = iargc()
  if (jarg < 5) then
    write(*,*) 'Usage: xviewmask template.pdb viewable.pdb datfile output.pdb output.vmas'
    write(*,*) 'viewmask finds a set of N points that are maximally dispersed '
    write(*,*) 'within the points defined in pdbfile.'
    write(*,*) 'The resulting mask can be used for drawing a surface.'
    write(*,*) 'MASK LIBRARY STEP 4b (after xbinarymask)'
    write(*,*) 'template.pdb is a large set of points at radius = ',radius
    write(*,*) 'viewable.pdb is a smaller set of points at radius = ',radius
    write(*,*) 'datfile is the logfile for the larger set from xbinarymask'
    write(*,*) 'output.pdb is a pdb file of masked points'
    write(*,*) 'output.vmas is a mask file of the selected points'
    ! stop 'viewmask.f90 v.28-AUG-01'
    write(0,*) 'viewmask.f90 v.  Thu May  3 17:34:27 EDT 2007'
    stop
  endif
  !! get command line arguments
  call getarg(1,template)  !! coordinatea of mask bits
  call getarg(2,pdbfile)   !! coordinates of smaller mask bits
  call getarg(3,datfile)   !! data from binarymask.f90 output.
  call getarg(4,outfile)   !! viewable points output in PDB format
  call getarg(5,binfile)   !! mask with viewable bits set.
  !! initialize masker arrays
  call masker_plot(drawingatall)
  call masker_initmasks()

  !! at this point there should be exactly MAXATOM in mxyz()
  !! read log file from binarymask
  open(11,file=datfile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'viewmask.f90:: error opening template.pdb file'
  ipsi = 0
  npsi = 0
  do
    !! the format/keywords for these lines are specifically those
    !! output by xbinarymask (binarymask.f90)
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:4)=="NPSI".and.npsi==0) then
      read(aline(7:),*,iostat=ios) npsi
      if (ios/=0) stop 'Format error 1'
      if (allocated(psiposit)) deallocate(psiposit)
      if (allocated(nphi)) deallocate(nphi)
      allocate(psiposit(0:npsi),nphi(0:npsi),stat=ios)
      if (ios/=0) stop 'Problem allocating psiposit, nphi'
    endif
    if (aline(1:4)=="PSI=".and.npsi/=0) then
      read(aline(19:24),*,iostat=ios) nphi(ipsi)
      if (ios/=0) stop 'Format error 2'
      read(aline(32:38),*,iostat=ios) psiposit(ipsi)
      if (ios/=0) stop 'Format error 3'
      write(*,*) nphi(ipsi), psiposit(ipsi)
      ipsi = ipsi + 1
    endif
    if (aline(1:5)=="NMASK".and.npsi/=0) then
      read(aline(7:),*,iostat=ios) nmask
      if (ios/=0) stop 'Format error 4'
    endif
  enddo
  close(11)
  write(*,'("NPSI =",i7)') npsi
  write(*,'("NMASK=",i7)') nmask
  !! read template of mask -------
  call masker_getmxyz(mxyz)
  mxyz = radius*mxyz
  !! diagnostic
  write(*,'("mxyz:: ",3f8.3))') mxyz(1:3,1:MAXATOM)
  nat = MAXATOM
  !write(*,*) 'Opening template PDB file'
  !open(11,file=template,status='old',form='formatted')
  !i = 0
  !do
  !  read(11,'(a)',iostat=ios) aline
  !  if (ios/=0) exit
  !  if (aline(1:5) /= "ATOM ") cycle
  !  i = i + 1
  !enddo
  !nat = i
  !! allocate and read template again -------
  !allocate(mxyz(3,nat),stat=ios)
  !if (ios/=0) stop 'Error allocating mxyz'
  !rewind(11)
  !i = 0
  !do
  !  read(11,'(a)',iostat=ios) aline
  !  if (ios/=0) exit
  !  if (aline(1:5) /= "ATOM ") cycle
  !  i = i + 1
  !  read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
  !enddo
  !! read second template for viewable dots --------
  close(11)
  open(11,file=pdbfile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'viewmask.f90:: error opening pdbfile'
  i = 0
  do
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:5) /= "ATOM ") cycle
    i = i + 1
  enddo
  vat = i
  write(*,*) 'Allocating ',vat,' viewable points.'
  allocate(vxyz(3,vat),stat=ios)
  if (ios/=0) stop 'Error allocating vxyz'
  rewind(11)
  i = 0
  do
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:5) /= "ATOM ") cycle
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) vxyz(1:3,i)
    !! write(*,'("vxyz:: ",i7,3f8.3))') i,vxyz(1:3,i)
  enddo
  close(11)
  !! Find the equivalent dot on the template mask
  amask = 0
  do i=1,vat
    dmin = 999.
    k = 0
    do j=1,nat
      d = dist(mxyz(1:3,j),vxyz(1:3,i))
      if (d < dmin) then
        k = j
        dmin = d
      endif
    enddo
        !! diagnostic
        write(*,'(2i7,2(3f8.3,2x),f9.5)') i,j,vxyz(1:3,i),mxyz(1:3,k),dmin
    ibyte = (k-1)/NBIT + 1
    ibit = mod(k-1,NBIT)
    amask(ibyte) = IBSET(amask(ibyte),ibit)
  enddo
  
  !! Write out the viewable points mask.
  open(12,file=binfile,status='replace',form='unformatted',iostat=ios)
  if (ios/=0) stop 'viewmask.f90:: error opening file for writing. Permissions?'
  write(12) amask
  close(12)
  open(13,file=outfile,status='replace',form='formatted',iostat=ios)
  if (ios/=0) stop 'viewmask.f90:: error opening file for writing. Permissions?'
  cen = 0.
  r = 1.
  iat = 1
  call masker_drawsurface(13,cen,r,iat,amask,10.,"V",scale=10.,ires=1)
  close(13)
  call masker_deallo()
  deallocate(vxyz,mxyz)
  !CONTAINS
  !!------------------------------------------------------------
!subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,scale,ires)
!
!subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,iskip)
!implicit none
!integer,intent(in) :: iunit,iat,iskip
!real,dimension(3),intent(in) :: cen
!real,intent(in) :: r,bb
!real :: d
!character(len=1),intent(in) :: ch
!integer,dimension(MASKSIZE),intent(in) :: amask
!integer :: i,j,ibyte,ibit
!real,dimension(3) :: avec,bvec
!real,parameter :: DCUT=6.0
!!
!j = 0
!do i=1,MAXATOM
!  ibyte = (i-1)/NBIT + 1
!  ibit = mod(i-1,NBIT)
!  if (btest(amask(ibyte),ibit)) then
!    j = j + 1
!    if (mod(j,iskip)==0) then
!      avec = cen + r*mxyz(1:3,i)
!      ! if (ishow /= 0) then
!      !   bvec = avec - showvec
!      !   if (sqrt(dotprod(bvec,bvec)) > DCUT) cycle
!      ! endif
!      write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
!      i,ch,j,avec(1:3),0.,bb
!    endif
!  endif
!enddo
!end subroutine drawsurface
!!==============================================================
!----------------------------------------------------------------------------------------------!
!        real function dist(x,y)
!        ! integer,intent(in) :: ix,j
!        real,intent(in) ::  x(3),y(3)
!        real :: d,eeps=0.001
!        d = (x(1) - y(1))**2 
!        d = d + (x(2) - y(2))**2 
!        d = d + (x(3) - y(3))**2 
!        d = sqrt(d)
!        dist = d
!        end function dist
!----------------------------------------------------------------------------------------------!
end program viewmask
