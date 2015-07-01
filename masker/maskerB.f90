module masker
  implicit none
  !! Module version of pdbmask.f90
  !! initmasks -- allocates memory, reads masks, gets collarsizes etc.
  !! getms -- accepts coordinates and atomtypes, returns MS
  !! -------------------------------------------------------------
  !!
  !! Fri Jun 15 19:55:28 EDT 2001
  !! C.Bystroff
  !! This program calculates the solvent accessible molecular surface
  !! of an arbitrary set of atoms in PDB format
  !!
  !! This program uses functions from the vectormath.f library.
  !!
  !! Compile using Fortran90
  !! i.e.:  pgf90 -o xpdbmask -Mextend pdbmask.f90 vectormath.f
  !!
  !! Required files: (environment variables)
  !! MASKLIB=4096.mas
  !! VMASK=4096.vmask
  !! VTRIANGLE=4096.tri
  !! MASKTEMPLATE=4096.pdb
  !! MASKDAT=4096.dat
  !! ATOMLIB=atoms.lib
  !! COUNTBITS=cbits32.bin
  !! COLLARS=collars4096.bin
  !!
  !! Command line arguments for xpdbmask:
  !! (1) Input file : inputpdb (pdb format)
  !! (2) Output file : output (pdb format, surface dots in HETATM lines)
  !! (3) Unit cell length [0.0] (default: no periodic cell)
  !! (4) Optional Output file : Raster3D triangulated surface (default: no output)
  !! Stdout : surface area with breakdown (sasa, saddle, bowl)
  !! -------------------------------------------------------------
  !*=========================================================================
  !** Copyright 2001-2002 Chris Bystroff
  !**
  !** Permission is hereby granted to use, execute, copy, distribute, modify,
  !** and distribute in modified form this software for non-profit purposes,
  !** provided this notice is retained in any and all copies.
  !**
  !** For for-profit usage, please contact bystrc@rpi.edu
  !*=========================================================================
  !! parameters
  integer,parameter :: MAXATOM=4096   !! number of points in mask. 
  integer,parameter :: KND=2          !! bytes per mask word
  integer,parameter :: NBIT=(KND*8)
  integer,parameter :: ILO=-(2**(15)),IHI=(2**(15))-1
  integer,parameter :: MASKSIZE=MAXATOM/NBIT    !! 3rd index of mask
  integer,parameter :: KC=1   !! thickness of collars (times DTHETA)
  integer,parameter :: runit=14  !! output unit for rendering
  integer,parameter :: watfac=500  !! max numbers of waters per atom
  real,parameter,private :: pi=3.1415927
  real,parameter :: radius=10.
  real,parameter,private :: rad=pi/180.
  real,parameter :: BOXSCALE=10.  !! scaling factor if periodic
  !!
  integer :: nmask,npsi,ishow=0,ntri    !! nmask= number of bisecting plane normals in MASKLIB
  integer(kind=KND),dimension(MASKSIZE),private :: emask,vmask
  integer,dimension(:),allocatable :: psiposit,nphi
  integer(kind=KND),dimension(:,:,:),allocatable,private :: masklib
  integer :: NTHETA
  real :: DTHETA   !! degree increments
  real,dimension(:),allocatable,private :: tau,slope,interc
  integer(kind=2),dimension(:,:),allocatable,private :: collarsize
  real,dimension(:),allocatable,private :: maskpsi,maskphi
  integer,dimension(:,:),allocatable,private :: tri
  ! real,dimension(3,MAXATOM) :: mxyz
  real,dimension(:,:),allocatable,private :: mxyz
  real,dimension(:,:),allocatable,private :: watvec
  real,dimension(3) :: cvec
  real :: rw=1.4  !! not private
  real :: wsphere, sasa, ssasa, bsasa
  integer,dimension(:),allocatable :: atype
  real :: dpsi 
  !! uses pi,npsi,dpsi,psiposit
  real :: maxd,taumin,r1,lbox
  real,dimension(3) :: showvec
  logical,private :: drawing
  logical :: periodic,rendering,smoothing,saying=.false.
  integer(kind=1),private,dimension(ILO:IHI) :: cbits
  type watertype
    real,dimension(3) :: xyz      !! coordinates of a probe position for a bowl
    integer,dimension(3) :: trng  !! indeces of base triangle for a water
    logical :: good               !! if this is false, don't use it for SES
  end type watertype
  type(watertype),dimension(:),allocatable,private :: water
  integer,private :: nwat  !! number of verteces
  integer,private :: iwat,jwat  
  type atomtype
    real :: r        ! united atom radius
    real :: m        ! united atom mass
    real :: w1,w2,w3 ! convert surface area to energy (sasa, saddle, bowl) kJ/mol/A^2
    character(len=4) :: name   !  element or other atom name
    !! format for atomlib is (a4,5f8.3)
  end type atomtype
  type(atomtype),dimension(:),allocatable :: atomlib
  integer :: nattype
  INTERFACE
    real function dotprod(a,b) 
      real,dimension(3) :: a, b
    end function dotprod
  end INTERFACE
CONTAINS
!=========================================================================!
  subroutine usage    !! show current environment variables
    character(len=80) :: aline
    write(*,*) "The following environment variables should be set (current setting):"
    call getenv("MASKLIB",aline); write(*,*) "MASKLIB=",trim(aline)
    call getenv("VMASK",aline); write(*,*) "VMASK=",trim(aline)
    call getenv("VTRIANGLE",aline); write(*,*) "VTRIANGLE=",trim(aline)
    call getenv("MASKTEMPLATE",aline); write(*,*) "MASKTEMPLATE=",trim(aline)
    call getenv("MASKDAT",aline); write(*,*) "MASKDAT=",trim(aline)
    call getenv("ATOMLIB",aline); write(*,*) "ATOMLIB=",trim(aline)
    call getenv("COUNTBITS",aline); write(*,*) "COUNTBITS=",trim(aline)
    call getenv("COLLARS",aline); write(*,*) "COLLARS=",trim(aline)
    call getenv("SLOPES",aline); write(*,*) "SLOPES=",trim(aline)
  end subroutine usage
!-------------------------------------------------------------------------!
  subroutine initmasks
  
  character(len=80) :: aline,pdbfile,binfile,logfile,outfile,template,cbitfile
  character(len=80) :: vmaskfile,libfile,vtriangle,collarfile,slopefile,fmt
  integer :: i,j,k,L,m,jarg,ios,iargc,natm,kcc=KC
  integer(kind=KND),dimension(MASKSIZE) :: amask,bmask,cmask,dmask
  integer :: nn,ipsi,itmp,itheta,iphi,imask,ibyte,ibit,nm,nv
  integer :: imask1,itheta1,imask2,itheta2,iat,nat,jat,kat,ires
  real :: dphi,sphere,col
  real :: phi,psi,theta,x,y,d,dd,r2
  real,dimension(3) :: avec,bvec,wxyz,showvec,uvec
  integer,parameter :: iunit=8,junit=15
  !! defaults
  DTHETA = 4.5
  NTHETA = int(90./DTHETA)
  dpsi = DTHETA*rad
  wsphere=rw*rw*4*pi
  !!
  !!  setenv  MASKLIB  4096.mas
  !!  setenv  VMASK  4096.vmask
  !!  setenv  VTRIANGLE  4096.tri
  !!  setenv  MASKTEMPLATE  4096.pdb
  !!  setenv  MASKDAT   4096.dat
  !!  setenv  ATOMLIB   atoms.lib
  !!  setenv COUNTBITS  cbits32.bin
  !!  setenv COLLARS  collars4096.bin
  !! 
  call getenv("MASKLIB",binfile)
  if (binfile=="") binfile = "4096.mas"
  call getenv("VMASK",vmaskfile)
  if (vmaskfile=="") vmaskfile = "4096.vmask"
  call getenv("VTRIANGLE",vtriangle)   !! triangle indeces
  if (vtriangle=="") vtriangle = "4096.tri"
  call getenv("MASKTEMPLATE",template)
  if (template=="") template = "4096.pdb"
  call getenv("MASKDAT",logfile)
  if (logfile=="") logfile = "4096.dat"
  call getenv("ATOMLIB",libfile)
  if (libfile=="") libfile = "atoms.lib"
  call getenv("COUNTBITS",cbitfile)
  if (cbitfile=="") cbitfile = "cbits32.bin"
  call getenv("COLLARS",collarfile)
  if (collarfile=="") collarfile = "collars4096.bin"
  call getenv("SLOPES",slopefile)
  if (slopefile=="") slopefile = "slopes4.5.dat"

  !! COUNTBITS
  !! initialize bit counting array, precalculated for KND-byte words
  !call initcountbits
  open(iunit,file=cbitfile,form="unformatted",status="old",iostat=ios)
  if (ios/=0) stop 'COUNTBITS file is missing'
  read(iunit,iostat=ios) cbits(ILO:IHI)
  if (ios/=0) stop 'masker.mod: ERROR reading COUNTBITS file.'
  close(iunit)

  !! ATOMLIB
  !! Read atom types, radii, weights, etc. from ATOMLIB file
  open(junit,file=libfile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'ATOMLIB file is missing'
  read(junit,*,iostat=ios) nattype
  if (ios/=0) stop 'masker.mod:ERROR reading ATOMLIB file. line 1.'
  allocate(atomlib(nattype),stat=ios)
  if (ios/=0) stop 'masker.mod: Error allocating atomlib'
  !! Thu Aug 30 07:39:02 EST 2001 new atomlib (from masker2.f90)
  do i=1,nattype
    read(junit,'(a4,5f8.3)',iostat=ios) atomlib(i)%name, atomlib(i)%m, atomlib(i)%r,  &
           atomlib(i)%w1,  atomlib(i)%w2,  atomlib(i)%w3
    if (ios/=0) stop 'masker.mod:ERROR reading ATOMLIB file. lines 2-'
  enddo
  close(junit)

  !! MASKDAT
  !! Read log file from binarymask . This contains indexing data for masks
  open(11,file=logfile,status='old',form='formatted')
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
      DTHETA = 180./npsi
      NTHETA = int(90./DTHETA)
      write(*,*) 'DTHETA=',DTHETA,'  NTHETA=',NTHETA
      allocate(psiposit(0:npsi),nphi(0:npsi),stat=ios)
      if (ios/=0) stop 'Problem allocating psiposit, nphi'
      !! allocate(tau(NTHETA),stat=ios)
      !! if (ios/=0) stop 'Problem allocating tau'
      !! write(*,*) 'psiposit and nphi allocated',npsi
    endif
    if (aline(1:4)=="PSI=".and.npsi/=0) then
      !!
      !! The format for the lines read from datfile is defined in
      !! binarymask.f90
      !!
      read(aline(19:24),*,iostat=ios) nphi(ipsi)
      if (ios/=0) stop 'Format error 2'
      read(aline(32:38),*,iostat=ios) psiposit(ipsi)
      if (ios/=0) stop 'Format error 3'
      !! write(*,*) nphi(ipsi), psiposit(ipsi)
      ipsi = ipsi + 1
    endif
    if (aline(1:5)=="NMASK".and.npsi/=0) then
      read(aline(7:),*,iostat=ios) nmask
      if (ios/=0) stop 'Format error 4'
    endif
  enddo
  close(11)
  !! Report the size of the mask and allocate memory accordingly
  ! write(*,'("NPSI =",i7)') npsi
  ! write(*,'("NMASK=",i7)') nmask
  allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for masklib'
  allocate(maskpsi(nmask),maskphi(nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for maskpsi,maskphi'
  allocate(collarsize(0:NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for collarsize'
  allocate(slope(NTHETA),interc(NTHETA),stat=ios)
  if (ios/=0) stop 'Problem allocating slope, interc'
  ! write(*,*) 'Memory allocated for masklib',MASKSIZE, NTHETA, nmask

  !! MASKLIB
  !! read binary format mask library, created by xbinarymask
  open(12,file=binfile,status='old',form='unformatted')
  read(12,iostat=ios) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
  if (ios/=0) stop 'initmasks: Error read mask library MASKLIB'
  close(12)
  !! temporary utility
  !! create fortran module for masklib
  !!  open(12,file="initmask.f",status="replace",form="formatted")
  !!  write(12,'("module initmask4096")')
!!
  !!  write(12,'("integer,parameter :: MAXATOM=4096   !! number of points in mask.")')
  !!  write(12,'("integer,parameter :: KND=2          !! bytes per mask word")')
  !!  write(12,'("integer,parameter :: NBIT=(KND*8)")')
  !!  write(12,'("integer,parameter :: NTHETA=",i4," !! 2nd index of mask")') NTHETA
  !!  write(12,'("integer,parameter :: NMASK=",i4," !! 3rd index of mask")') nmask
  !!  write(12,'("integer,parameter :: MASKSIZE=MAXATOM/NBIT    !! 1st index of mask")')
  !!  write(12,'("integer(kind=KND),dimension(MASKSIZE,NTHETA,NMASK) :: masklib= reshape((/ &")')
  !!  write(12,'("  ",$)')
  !!  j = 2
  !!  do imask=1,NMASK
  !!    do itheta=1,NTHETA
  !!      do i=1,MASKSIZE
  !!        if (masklib(i,itheta,imask)==0) then
  !!          k = 1
  !!        else
  !!          k = int(log10(float(abs(masklib(i,itheta,imask))))+1)
  !!          if (masklib(i,itheta,imask) < 0) k = k + 1
  !!        endif
  !!        write(fmt,'("(i",i1,",$)")') k
  !!        write(12,fmt) masklib(i,itheta,imask)
  !!        j = j + k + 1
  !!        if (i < MASKSIZE) then
  !!          write(12,'(",",$)')
  !!          if (j > 70) then
  !!            write(12,'(" &",/,"  ",$)')
  !!            j = 2
  !!          endif
  !!        endif
  !!      enddo
  !!      if (itheta < NTHETA) write(12,'(",",$)')
  !!      j = j + 1
  !!    enddo
  !!    if (imask < NMASK) write(12,'(",",$)')
  !!    j = j + 1
  !!  enddo
  !!  write(12,'("/),(/",i4,",",i4,",",i4,"/)) ")') MASKSIZE,NTHETA,NMASK
  !!  write(12,'("end module initmask4096")')
  !!  close(12)
  !!  stop
  !!!! diagnostic
  ! write(*,*) 'Done reading masklib'

  !! VMASK
  !! Get viewable dots mask . This is the mask for outputing
  !! points in PDB format
  open(12,file=vmaskfile,status='old',form='unformatted')
  read(12,iostat=ios) vmask(1:MASKSIZE)
  if (ios/=0) stop 'initmasks: Error read viewable mask'
  close(12)

  !! VTRIANGLE
  !! Get tesslation triangles for rendering. 
  !! vtriangle is the output point-indeces for a complete
  !! set of triangles on the mask template sphere. This is 
  !! produced by triangles.f90 and must be regenerated 
  !! if a new mask is used.
  open(12,file=vtriangle,status='old',form='formatted')
  ntri = 0
  do
    read(12,*,iostat=ios) i,j,k
    if (ios/=0) exit
    ntri = ntri + 1
  enddo
  allocate(tri(3,ntri),stat=ios)
  if (ios/=0) stop 'Error allocating tri()'
  ! write(*,*) ntri, ' Triangles for rendering.'
  rewind(12)
  do i=1,ntri
    read(12,*,iostat=ios) tri(1:3,i)
    if (ios/=0) stop 'Error re-reading triangles'
  enddo
  close(12)
  !! diagnostic
  ! write(*,*) 'Done reading vmask file'
  !!diagnostic
  ! write(*,*)  masklib(1:MASKSIZE,10,50)

  !! Get psi and phi for each mask
  !! ..AND..
  !! Calculate collarsize and tau for each theta and imask
  !! collarsize is the number of bits in the edge mask for angle theta
  !! tau is the angular width of the half-saddle for angle theta and r1
  !!
  !! PSI is the angle from the positive Z-axis of the normal to the bisecting plane
   do ipsi=0,npsi
     psi = ipsi*dpsi
     dphi = 2*pi/real(nphi(ipsi))
     !! PHI is the angle of the projection of the normal 
     !! on the XY plane to the positive X-axis, R-handed
     do iphi=1,nphi(ipsi)
       phi = (iphi-1)*dphi
       imask = psiposit(ipsi) + iphi - 1
    !   y = cos(psi)
    !   x = sin(psi)
    !   uvec = (/x*cos(phi),x*sin(phi),y/)
    !   uvec = -uvec
    !   imask2 = getimask(uvec)  !! inverse mask index
       maskpsi(imask) = psi
       maskphi(imask) = phi
       !! diagnostic
       !! write(*,*) imask,psi,phi
       !! get collarsize for each difference mask
     !!!!! NOTE: These lines are used to re-generate the data in the COLLARS file
     !!!!!       if that file is lost.
    !   do itheta=0,NTHETA-1
    !     theta = itheta*dpsi
    !     if (itheta == 0) then
    !       amask = -1
    !       amask = iand(amask,not(masklib(:,itheta+1,imask)))
    !     else 
    !       amask = iand(masklib(:,itheta,imask),not(masklib(:,itheta+1,imask)))  
    !     endif
    !     j = countbits(amask)
    !     collarsize(itheta,imask) = j    !! approximately the same for all imask
    !   enddo
    !   !! this number must be non-zero, but doesn't matter since tau=0.
    !   !! Thu Feb  7 14:35:30 EST 2002 
    !   !! Note comment 2 lines above. IT DOES MATTER. For the interpolation step when
    !   !! atom j is embedded, the collarsize at 90 deg is used, if theta = 90 deg.
    !   !! note: kcc is obselete (=1). removed.
    !   amask = iand(masklib(:,NTHETA,imask),masklib(:,NTHETA-1,imask2))  
    !   j = countbits(amask)
    !   collarsize(NTHETA,imask) = j   ! approx the same as for NTHETA-1
     enddo
   enddo
  
  !! COLLARS
  !! read collarsizes from file
  !! open(33,file=collarfile,status="replace",form="unformatted",iostat=ios)
  open(33,file=collarfile,status="old",form="unformatted",iostat=ios)
  if (ios/=0) stop 'masker.mod: Error opening COLLARS'
  !! write(33) collarsize(0:NTHETA,1:nmask)  !! use this to re-generate the file
  read(33,iostat=ios) collarsize(0:NTHETA,1:nmask)
  if (ios/=0) stop 'masker.f90: ERROR reading from COLLARS'
  close(33)
  !! stop  'Done saving collarsizes'
  !! NTHETA is 90deg, only possible if dd ~= 0. No saddle? 
  !! =====> NOTE: If hydrogens are used, this will have to be changed! <=====

  !! SLOPES
  !! Read the calibration data from file: slope of A^2 vs theta (radians) for a unit sphere
  !! indexed by itheta, which is incremented by dpsi (dtheta). 
  open(34,file=slopefile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'masker.mod: Error opening SLOPES'
  do itheta=1,NTHETA
    read(34,*) slope(itheta),interc(itheta)
  enddo
  close(34)
  !! NOTE: slopes and intercepts are generated by plottheta.f90
  !! which takes as input the masks and mask data used in this module
  !! to re-run plottheta, using this module, you must make a dummy
  !! SLOPES file with NTHETA pairs of reals, 2 on each line.
  !! slope() is in units of suare angstroms on a unit sphere, per radian.
  !! Mon Dec 17 12:47:27 EST 2001
  !! To convert to contact (sasa) surface gradient, multiply by (arcfrac*sphere/4pi)*dRdD
  !! where sphere is the surface area of a sphere of r1, arcfrac is the fraction of the
  !! collar that is exposed, and dRdD is the slope of radians per angstrom, which
  !! is d(theta)/d(Dij) = 
  !! ((r1+rw)**2 + (r2+rw)**2 + Dij)/(-2 Dij * (r1+rw)) - 1/(r1+rw) 
  !! Sat Feb  2 21:29:13 EST 2002

  !! MASKTEMPLATE
  !! Coordinates of mask points on a 10A sphere.
  !! These are only needed if drawing=.true.
  if (drawing) then
    allocate(mxyz(3,MAXATOM),stat=ios)
    if (ios/=0)  stop 'masker.mod: error allocating mxyz'
    open(11,file=template,status='old',form='formatted',iostat=ios)
    if (ios/=0) stop 'masker.mod: error opening MASKTEMPLATE file.'
    i = 0
    do 
      read(11,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:5) /= "ATOM ") cycle
      i = i + 1
      if (i > MAXATOM) stop 'Too many points in template file MASKTEMPLATE'
      read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
      if (ios/=0) then
        write(*,*) 'Error reading template coordinate file, at i=',i
        stop 'initmasks: Error reading template coordinate file MASKTEMPLATE'
      endif
      mxyz(1:3,i) = mxyz(1:3,i)/10.   !! template has 10A radius
    enddo
    close(11)
  endif
  !! 
  end subroutine initmasks
  !!----------------------------------------------------------------!!
  subroutine getms(xyz,nat,ms,peratom)
  integer, intent(in) :: nat
  real,dimension(3,nat),intent(in) :: xyz
  !! integer,dimension(nat),intent(in) :: atype
  real,intent(out) :: ms
  real,intent(out),dimension(nat) :: peratom
  !
  integer :: i,j,k,L,m,jarg,ios,iargc,natm
  integer(kind=KND),dimension(MASKSIZE) :: amask,bmask,cmask,dmask
  integer :: nn,ipsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv
  integer :: imask1,itheta1,imask2,itheta2,iat,jat,kat,ires
  integer,dimension(100) :: saveimask, saveitheta,saveiat
  real,dimension(100) :: savetheta, savedij, saver2
  real :: phi,psi,theta,x,y,d,dd,r2,r3
  integer :: kskip=100,kcc=KC, flag
  real :: dphi,sphere,col,sasanrg,ssasanrg,bsasanrg,xx,arcfrac
  real,dimension(3) :: avec,bvec,wxyz,showvec,veci,vecj,veck
  real :: dotprod
  !
  i = 0
  periodic = (lbox > 0.)
  if (ishow /= 0 .and.drawing) then
    write(*,'("Show surface around residue ",i4," only.")') ishow
  endif
  
  sasa = 0.0    !! total solvent accessible molecular surface area
  ssasa = 0.0   !! total "saddle" convex/concave surface area
  bsasa = 0.0   !! total "bowl" concave/concave surface area
  sasanrg = 0.0    !! total solvent accessible molecular surface area
  ssasanrg= 0.0   !! total "saddle" convex/concave surface area
  bsasanrg= 0.0   !! total "bowl" concave/concave surface area
  nwat = 0        !! Number of verteces for re-entrant surface calc.
  peratom = 0.
    
  !! CONTACT SURFACE  (VDW,  SAS)  
  wsphere=rw*rw*4*pi
  !! Loop over all atoms, get sasa, saddle and bowl surfaces
  do iat=1,nat
    nm = 0
    amask = -1                 ! init all bits = 1
    r1 = atomlib(atype(iat))%r
    do jat=1,nat
      if (iat == jat) cycle
      !! For every atom pair, calculate phi, psi, theta
      !! with iat as the central atom
      if (periodic) then ! get nearest copy of jat, return vector
        call boxaround(xyz(1:3,iat),xyz(1:3,jat),lbox,avec)
      else
        avec = xyz(1:3,jat) - xyz(1:3,iat)
      endif
      dd = sqrt(dotprod(avec,avec))
      r2 = atomlib(atype(jat))%r
      maxd = r1 + r2 + 2*rw
      if (dd > maxd) cycle    !! atom out of range
      avec = avec/dd
      !! returns imask and theta. If theta > pi/2, use 2*NTHETA-itheta and .not. it
      call getimask2(r1+rw,dd,r2+rw,theta,avec,imask,itheta)
      if (imask==0) cycle   !! bad triangle (bug?)
      if (imask==-1) then   !! i is completely embedded in j
        amask = 0
        exit
      endif
      if (theta < pi/2) then
        !! changed 17-DEC-01  nint
        !! itheta = nint((theta/rad)/DTHETA) 
        if (itheta==0) itheta = 1
        if (itheta==NTHETA) then
          amask = iand(not(masklib(:,itheta,imask)),amask)     
        else
          amask = iand(masklib(:,itheta,imask),amask)     
        endif
      else  !! jat is embedded in iat
        !! changed 17-DEC-01  nint
        !! itheta =  nint(((pi-theta)/rad)/DTHETA)    !! itheta index for embedded atom
        !if (imask==0) then   !! flag for completely embedded (also theta=pi)
        !  cycle    !! completely embedded, no surface, go to next atom
        !else
          amask = iand(not(masklib(:,itheta,imask)),amask)     !!invert mask to get theta > pi/2, embedded mask.
        !endif
      endif
      !! if (itheta <= 0) cycle    !! atom out of range (borderline)
                                !! borderline atoms should have gradient ~ #of points in itheta=1 mask.
      !! diagnostic
      !!   if (iat==1.and.jat==2) write(*,*) 'Theta: ',theta/rad
      !! diagnostic
      ! write(*,'(2i4,4f6.1,4i5)') iat,jat,psi/rad,phi/rad,theta/rad,dd,ipsi,iphi,itheta,imask
      !! if (itheta > NTHETA-1) then
      !!   write(*,*) "ERROR: itheta=",itheta," NTHETA=",NTHETA
      !!   write(*,*) " theta=",theta/rad," iat,jat,d=",iat,jat,dd
      !!   stop 'Theta out of range'  !! this would spot a bug
      !! endif
      !! remove the bits covered by atom jat
      !! diagnostic / experimental: keep all close atoms, don't use bmask 21-NOV-01
      !!
      !! Using masks to determine whether a water exists has problems.
      !! It is found that numerical errors occur when waters (tetrahedron verteces)
      !! are close to each other. Instead of numerical methods, exact calculations
      !! are now used to find which of these waters to keep, which to exclude
      !! Mon Nov 26 11:52:10 EST 2001
      !!
      ! if (any(bmask /= amask) ) then
      !   amask = bmask
        !! keep a list of close atoms, for the next part
        nm = nm + 1
        saveiat(nm) = jat
        saveimask(nm) = imask
        saveitheta(nm) = itheta   !! may be theta or (pi - theta) index.
        savetheta(nm) = theta     !! depending on whether theta > pi/2
        savedij(nm) = dd
      ! endif
    enddo
    j = countbits(amask)
    if (j<=1) cycle   ! ignore saddle, etc if nothing is exposed.
    ! enddo
    !
    !! CONTACT SURFACE
    sphere = 4*pi*atomlib(atype(iat))%r**2
    x = sphere*(real(j)/real(MAXATOM))  !! MAXATOM=4096
    sasa = sasa + x    !! x is the sasa before edge correction
    sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
    peratom(iat) = peratom(iat) + x
    if (drawing) then
      cmask = iand(amask,vmask)
      call drawsurface(13,xyz(1:3,iat),r1,iat,cmask,50.,"V",1)  !! output the sasa surface
    endif
    !
    !! TOROIDAL SURFACE  (SADDLE)
    !
    i = 0
    cmask = 0
    nn = 0
    do k=1,nm
      theta = savetheta(k)
      kat = saveiat(k)
      imask = saveimask(k)
      itheta = saveitheta(k)
      dd = savedij(k)
      r2 = atomlib(atype(kat))%r
      !! KC was added to provide the option of double-width (KC=2) 
      !! saddles. Mon Jun 25 16:57:13 EDT 2001
      !! Collarsizes of width KC are calculated. 
      !! Obselete if edges are thick enough to guarantee a point at every vertex.
      !!   Sat Jun 30 21:30:14 EDT 2001
      !! if (itheta > NTHETA-KC) then; kcc = 1; else; kcc=KC; endif
      !! diagnstic
      !! write(*,'(a,i4,i4,f8.5,i4,i7,i6)') &
      !! "i k theta itheta", iat, kat, theta, itheta, imask, collarsize(itheta,imask)

      if (theta >= pi/2) then  !! embedded, no saddle
        !!write(*,*) "embedded itheta=", itheta
        if (itheta <= 1) then
          bmask = amask
        else
          bmask = iand(masklib(:,itheta-1,imask),amask) !! this arc is for derivs and drawing
        endif
      else
        !!write(*,*) "not embedded itheta=", itheta
        if (itheta == NTHETA) then  !! imask is for inverse mask
          bmask = iand(masklib(:,itheta-1,imask),amask)
        else
          bmask = iand(not(masklib(:,itheta+1,imask)),amask)
        endif
      endif
      cmask = ior(cmask,bmask)  !! collect all edges
      if (drawing) then
        call drawsurface(13,xyz(1:3,iat),r1,kat,bmask,50.,"E",1)  !! output the edges 1 at a time
      endif
      !! Count the bits in the exposed arc
      j = countbits(bmask)
      arcfrac = real(j)/real(collarsize(itheta,imask))
      !! diagnostic
      !! if (iat==1.and.kat==2) then
      !!   write(*,*) iat,kat, j, countbits(amask), theta
      !! endif
      ! j = 0
      ! do i=1,MAXATOM
      !   ibyte = (i-1)/NBIT + 1
      !   ibit = mod(i-1,NBIT)
      !   if (btest(bmask(ibyte),ibit)) then
      !     j = j + 1
      !   endif
      ! enddo
      ! write(*,*) "btest:", j
      !! 
      !! Saddle surface:
      !! The following is an "exact approximation", meaning
      !! that the estimate is as good as the numerical estimates
      !! of the arc distance, theta, phi, psi.
      !! (equations derived using Mathematica)
      !! Mon May  7 10:34:43 EDT 2001
      !! 
      if (theta < pi/2) then  !! saddle exists
        !! theta = itheta*dpsi  !! why not use the real theta?  Mon Jun 18 10:49:46 EDT 2001
        !! SASA correction
        x = (theta - interc(itheta))*slope(itheta)*arcfrac*sphere
        sasa = sasa + x
        sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
        !! saddle
        y = (r1+rw)*sin(theta)
        col = 2*pi*(arcfrac)   !! this is the length of the exposed arc in radians
        !! diagnostic
        !! write(*,*) 'Saddle arc length: ',iat,kat,col," theta=",theta/rad
        x = col*((0.5*pi-theta)*y - rw*cos(theta))*rw
        taumin = 0.
        if (y < rw) then
          taumin = acos(y/rw)
          x = x + col*(rw*sin(taumin) - taumin*y)*rw
        elseif ( (r1+rw)**2 > (dd*dd) + (r2+rw)**2) then  !! kat is embedded
          taumin = acos(y/(r2+rw))
          x = x + col*(rw*sin(taumin) - taumin*y)*rw
        endif
        ssasa = ssasa + x
        peratom(iat) = peratom(iat) + x   !! add toroidal surface per atom
        ssasanrg = ssasanrg + x*atomlib(atype(iat))%w2
        if (drawing) then
          x = (pi/2) - theta  !! taumax
          !! x = itheta*dpsi  !! estimated theta
          !! output the saddle surface
          call drawsaddle(13,xyz(1:3,iat),xyz(1:3,kat),bmask,r1,iat,taumin,x,50.,"S")  
          if (rendering) then
            !! NOTE! this will fail if theta is close to 90 degrees (it never is)
            dmask = iand(not(masklib(:,itheta+3,imask)),amask) 
            call rendersaddle(xyz(1:3,iat),xyz(1:3,kat),r1,amask,dmask,x,taumin)
          endif
        endif
      else    !! theta > pi/2
        !! write(*,*) 'No Saddle area: ',iat,kat,col," theta=",theta/rad
        x = (interc(itheta) - pi + theta)*slope(itheta)*arcfrac*sphere
        sasa = sasa + x
        !! diagnostic 
        ! write(*,'(8x,i4," embedded correction=",f8.2)') kat,x
        ! if (abs(x) > 10.) then
        !   write(*,*) "theta ",theta," itheta ",itheta
        !   write(*,*) "interc ",interc(itheta), " slope ",slope(itheta)
        !   write(*,*) "arcfraC ",arcfrac, " sphere ",sphere
        ! endif
      endif
      !! Diagnostic: don't reduce saved atom set
      ! if (any(bmask /= 0)) then
      !   !! keep a reduced list of close atoms.
      !   nn = nn + 1    !!  note: nn <= k, so this is legal.
      !   saveiat(nn) = kat
      !   saveimask(nn) = imask
      !   saveitheta(nn) = itheta
      !   savetheta(nn) = theta
      ! endif
    enddo
    ! diagnostic: use un-reduced atom set
    nn = nm

    !!
    !! Using the reduced list of close atoms, get the intersections
    !! between the 'collars' .
    !! Each intersection corresponds to a water. Get the surface
    !! of the water that is in the spherical triangle between the three 
    !! atoms.  Subtract any part of that surface that falls below the
    !! plane of those three atoms. (iat, kat, jat)
    !! 
    nv = 0
    dmask = 0
    emask = 0
    do k=1,nn
      kat = saveiat(k)
      if (kat <= iat) cycle  !! calculate bowl iff iat is the first atom
      imask1 = saveimask(k)
      itheta1 = saveitheta(k)
      ! theta1 = savetheta(k)
      r2 = atomlib(atype(kat))%r
      do m=1,nn
        jat = saveiat(m)
        if (jat <= kat) cycle  !! calculate bowl iff kat is the second atom
        imask2 = saveimask(m)
        itheta2 = saveitheta(m)
        !! --------- Ignore collars. Just get the third side of the triangle
        !! Wed Nov 21 08:15:00 EST 2001
        if (periodic) then ! get nearest copy of jat, return vector
          dd = distmod(xyz(1:3,kat),xyz(1:3,jat),lbox)
        else
          avec = xyz(1:3,jat) - xyz(1:3,kat)
          dd = sqrt(dotprod(avec,avec))
        endif
        r3 = atomlib(atype(jat))%r
        maxd = r3 + r2 + 2*rw
        if (dd > maxd) cycle    !! third side too long 
        !! ----------
        ! theta2 = savetheta(m)
        !! diagnostic (bmask commented out)
        !   bmask = iand(not(masklib(:,itheta1+kcc,imask1)),masklib(:,itheta1,imask1) )
        !   bmask = iand(bmask,not(masklib(:,itheta2+kcc,imask2)) )
        !   bmask = iand(bmask,masklib(:,itheta2,imask2) )
        !   bmask = iand(cmask,bmask)  !! bmask is current verteces
        !!
        !!  Find/Save waters for RE_ENTRANT SURFACE (BOWL)
        !!
        !! diagnostic: forget the vertex. Just place a water
        !if (.true.) then
        !! if (any(bmask /= 0)) then !! if there is at least one vertex
        !  !! diagnostic
        !  !! write(*,*) 'Triangle:',iat,jat,kat
        !  !! Retrieve unit vector to atom k: avec
        !  !! save atom positions
          veci = xyz(1:3,iat)
          if (periodic) then
            call boxaround(veci,xyz(1:3,jat),lbox,avec)
            vecj = veci + avec
            call boxaround(veci,xyz(1:3,kat),lbox,avec)
            veck = veci + avec
          else
            vecj = xyz(1:3,jat)
            veck = xyz(1:3,kat)
          endif
        !  psi = maskpsi(imask1)
        !  phi = maskphi(imask1)
        !  avec(3) = cos(psi)
        !  avec(2) = sin(psi)*sin(phi)
        !  avec(1) = sin(psi)*cos(phi)
        !  !! Retrieve unit vector to atom m: bvec
        !  psi = maskpsi(imask2)
        !  phi = maskphi(imask2)
        !  bvec(3) = cos(psi)
        !  bvec(2) = sin(psi)*sin(phi)
        !  bvec(1) = sin(psi)*cos(phi)
        !  !! get the cross product. cvec is a normal to the 3-atom plane
        !  !! pointing in the direction of the L-handed water
        !  call cros(avec,bvec,cvec)
        !  d = sqrt(dotprod(cvec,cvec))
        !  cvec = cvec/d
        !  !! Get the mask that splits the vertex mask down the middle
        !  !! Zeros cover possible L-handed waters, Ones cover R-handed waters
        !  imask = getimask(cvec)
        !  itheta = NTHETA    !! hemisphere mask
        !  x = 0.
        !  !! First check on the R-hand side for verteces...
        !  emask = iand(masklib(:,itheta,imask),bmask)
        !  !! diagnostic: forget vecteces. just place all waters
        !  !! if (any(emask /= 0)) then
        !  if (.true.) then
        !    !! call waterbowl(iat,kat,jat,x,wxyz,xyz,nat,imask) !! also uses/returns emask, uses cvec
            call tetrahedron(veci,vecj,veck,wxyz,d,iat,jat,kat,1)
            if (d==0.0) cycle   !! fails triangle inequality
            if (exposed(wxyz,xyz,nat)) then
              nwat = nwat + 1
              if (nwat > watfac*nat) then
                write(*,*) 'Ran out of water space: iat=',iat,' limit=',watfac*nat
                stop 'not enough water'
              endif
              water(nwat)%xyz = wxyz; water(nwat)%trng(1:3) = (/iat,jat,kat/)
            endif
            !! diagno
        !    write(*,'("R ",4i4)') nwat,water(nwat)%trng(1:3)
        !    do i=1,3
        !      avec = water(nwat)%xyz - xyz(1:3,water(nwat)%trng(i))
        !      d = sqrt(dotprod(avec,avec))
        !      write(*,'(2f8.3,$)') d, atomlib(atype(water(nwat)%trng(i)))%r
        !    enddo
        !    write(*,*)
        !    ! call waterbowl(veci,vecj,veck,x,wxyz,imask,iat,jat,kat)  !! R-handed water (use cvec for hole)
        !    !! diagnostic, write water coords
        !    ! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        !    !   iat,"W",iat,wxyz(1:3),0.,10.0
        !    ! if (drawing) then  !! do this to force drawing surface for bowls
        !    !   dmask = iand(emask,vmask)  !! option to draw fewer dots
        !    !   !! dmask = emask
        !    !   call drawsurface(13,wxyz(1:3),rw,iat,dmask, 10.,"B",1) !! output the surface
        !    ! endif
        !  endif
        !  !! Now check for L-handed waters
        !  emask = iand(not(masklib(:,itheta,imask)),bmask)
        !  !! diagnostic: forget vecteces. just place all waters
        !  if (.true.) then
        !  ! if (any(emask /= 0)) then
            call tetrahedron(veci,veck,vecj,wxyz,d,iat,kat,jat,1)
            if (d==0.0) cycle   !! fails triangle inequality
            if (exposed(wxyz,xyz,nat)) then
              nwat = nwat + 1
              if (nwat > watfac*nat) then
                write(*,*) 'Ran out of water space: iat=',iat,' limit=',watfac*nat
                stop 'not enough water'
              endif
              water(nwat)%xyz = wxyz; water(nwat)%trng = (/iat,kat,jat/)
            endif
            !! diagno
        !    write(*,'("L ",4i4,$)') nwat,water(nwat)%trng(1:3)
        !    do i=1,3
        !      avec = water(nwat)%xyz - xyz(1:3,water(nwat)%trng(i))
        !      d = sqrt(dotprod(avec,avec))
        !      write(*,'(2f8.3,$)') d, atomlib(atype(water(nwat)%trng(i)))%r
        !    enddo
        !    write(*,*)
        !    ! if (x==0..or.drawing) then  
        !    !   !! call waterbowl(iat,kat,jat,y,wxyz,xyz,nat,imask) !! also uses/returns emask
        !    !   call waterbowl(veci,veck,vecj,y,wxyz,imask,iat,kat,jat)   !! L-handed water (use -cvec for hole)
        !    !   x = x + y
        !    !   !! diagnostic, write water coords
        !    !   ! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        !    !   !   iat,"W",iat,wxyz(1:3),0.,20.0
        !    !   if (drawing) then  !! do this to force drawing surface of bowls
        !    !     dmask = iand(emask,vmask)
        !    !     !! dmask = emask
        !    !     call drawsurface(13,wxyz(1:3),rw,iat,dmask, 11.,"B",1) !! output the surface
        !    !   endif
        !    ! else
        !    !   x = 2*x
        !    ! endif
        !  endif
        !  ! bsasa = bsasa + x
        !  ! bsasanrg = bsasanrg + x*(atomlib(atype(iat))%w3 + &
        !  !                    atomlib(atype(jat))%w3 + atomlib(atype(kat))%w3)
        !endif
      enddo
    enddo
  enddo   !! iat

  !! RE-ENTRANT SURFACE accounting for intersecting re-entrant surfaces.
  !! 
  !! This segment of the code (obselete given correction of bugs in tetrahedron?)
  !! was written to ensure that no waters sat too close to any atom.
  !! 
  xx = 4*rw*rw
  !! write(*,*) "nwat = ",nwat
  water(1:nwat)%good = .true.
  !!! obselete?
  ! do iwat = 1,nwat
  !   do iat=1,nat
  !     d = distmod(water(iwat)%xyz,xyz(1:3,iat),lbox)
  !     x = (rw + atomlib(atype(iat))%r)
  !     if (d + 0.001 < x) then
  !       water(iwat)%good = .false.
  !       ! write(*,'(a,4i4,a,i4,a,f6.2,a,f6.2)') "Water ",iwat,water(iwat)%trng(1:3)," too close to atom ",iat, "  d= ",d, " < ",x
  !       exit
  !     endif
  !   enddo
  ! enddo
  !!
  do iwat = 1,nwat
    if (.not.water(iwat)%good) cycle
    !! diagnostic, show water triangles. For debugging.
    !  write(*,*) iwat,water(iwat)%trng(1:3)
    !! write(*,'(4(3f8.3,1x))') water(iwat)%xyz(1:3),xyz(1:3,water(iwat)%trng(1)), &
    !!      xyz(1:3,water(iwat)%trng(2)),xyz(1:3,water(iwat)%trng(3))
    !! 
    !! trainglemask: get the mask  (amask) for the reentrant (bowl) surface
    !! at the position iwat. 
    !!
    call trianglemask(water(iwat)%xyz(1:3),xyz(1:3,water(iwat)%trng(1)), &
         xyz(1:3,water(iwat)%trng(2)),xyz(1:3,water(iwat)%trng(3)),amask) 
         !! this routine initializes amask
    !! Remove intersecting reentrant surfaces
    do jwat = 1,nwat
      if (.not.water(jwat)%good) cycle
      if (iwat == jwat) cycle
      if (periodic) then
        call boxaround(water(iwat)%xyz,water(jwat)%xyz,lbox,avec) 
      else
        avec = water(jwat)%xyz - water(iwat)%xyz
      endif
      d = dotprod(avec,avec)
      if (d < xx) then    !! xx = 4*rw*rw  = (twice the probe radius)^2
        if (d<=0.0) cycle  !! this would only happen (d=0) VERY RARELY
        d = sqrt(d)
        avec = avec/d
        imask = getimask(avec)
        call cosinerule(rw,d,rw,theta,flag); if (flag /= 0) cycle !! bug??
        itheta = nint((theta/rad)/DTHETA)
        if (itheta > NTHETA) stop 'BUG: Two waters '
        !! remove the bits in the iwat mask covered by water atom jwat
        !! diagnostic, monitor intersecting reentrants
        !! bmask = amask  !! save a copy for comparison
        amask = iand(masklib(:,itheta,imask),amask)     
        !! diagnostic, monitor intersecting reentrants
        !! if (any(amask /= bmask)) then
        !!   write(*,'("Intersection: iwat=",i7," jwat=",i7)') iwat,jwat
        !! endif
      endif
    enddo
    !! NOTE ON DERIVS:
    !! At this point, 'amask' contains the reentrant surface and 'water(iwat)' contains
    !! the indeces of the three atoms. To get the bowl derivs, get the three edgemasks
    !! for the three sides of the trianglemask. Each edge is the amount of surface
    !! buried for moving one atom toward the other two. How do we get pairwise from this??
    if (drawing) then  !! do this to force drawing surface of bowls
      !! diagnostic, draw water
      avec = water(iwat)%xyz(1:3)
      if (periodic) call inbox(avec,lbox)
      write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1,i7)') &
      iwat,"W",iwat,avec(1:3),0.,50.,nwat
      write(13,'("TER")') 
      !! dmask = iand(amask,vmask)
      ! call drawsurface(13,water(iwat)%xyz,rw,iwat,dmask, 50.,"B",1) !! output the surface
      call drawsurface(13,water(iwat)%xyz,rw,iwat,amask, 50.,"B",1) !! output the surface
    endif
    j = countbits(amask)
    x = wsphere*real(j)/real(MAXATOM)  !! area in A^2
    bsasa = bsasa + x
    !! diagnostic, writeout surface for each atoms
    !! write(*,'(f8.2,$)') x
    bsasanrg = bsasanrg + (x/3.)*(atomlib(atype(water(iwat)%trng(1)))%w3 + &
        atomlib(atype(water(iwat)%trng(2)))%w3 + atomlib(atype(water(iwat)%trng(3)))%w3)
    iat = water(iwat)%trng(1)
    peratom(iat) = peratom(iat) + x/3   !! add reentrant surface per atom
    iat = water(iwat)%trng(2)
    peratom(iat) = peratom(iat) + x/3   !! add reentrant surface per atom
    iat = water(iwat)%trng(3)
    peratom(iat) = peratom(iat) + x/3   !! add reentrant surface per atom
  enddo 
  y = sasa + ssasa + bsasa
  x = sasanrg + ssasanrg + bsasanrg
  if (drawing.and.saying) then
     write(*,'("Total VDW area     = ",f9.2)') sasa
     write(*,'("Total saddle area  = ",f9.2)') ssasa
     write(*,'("Total bowl area    = ",f9.2)') bsasa
     write(*,'("Total surface area = ",f9.2)') y
     write(*,'("Total VDW nrg      = ",f9.2," kJ")') sasanrg
     write(*,'("Total saddle nrg   = ",f9.2," kJ")') ssasanrg
     write(*,'("Total bowl nrg     = ",f9.2," kJ")') bsasanrg
     write(*,'("Total SES energy   = ",f9.2," kJ")') x
  endif
  ms = x
  !! close(13)
  end subroutine getms ! (xyz,atype,nat,ms)
  !!----------------------------------------------------------------!!
  logical function exposed(axyz,xyz,nat)
  implicit none
  integer,intent(in) :: nat
  real,dimension(3),intent(in) :: axyz
  real,dimension(3,nat),intent(in) :: xyz
  integer :: iat
  real :: x,d
  exposed = .false.
  do iat=1,nat
    d = distmod(axyz,xyz(1:3,iat),lbox)
    x = (rw + atomlib(atype(iat))%r)
    if (d < x) return   !! note: removed 0.001A buffer 14-MAR-02
  enddo
  exposed = .true.
  end function exposed
  !!----------------------------------------------------------------!!
  real function gettheta(r1,r2,dd)
  implicit none
  !! returns theta in radians for angle between sides r1 and dd
  !! for triangle w/sides length = r1,r2,dd
  !! using Heron's Rule
  !! More at: http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap06/area-2.html
  !!  (and many other web sites)
  real,intent(in) :: r1,r2,dd
  real :: s,x,ar
  s = (r1+r2+dd)/2.
  ar = s*(s-r1)*(s-r2)*(s-dd)
  if (ar < 0.) then
    gettheta = 0.
    return
  endif
  ar = sqrt(ar)
  x = 2*ar/dd
  x = asin(x/r1)
  gettheta = x
  !! diagnostic
  !! write(*,*) 'gettheta: ',r1,r2,dd,x
  end function gettheta
  !!----------------------------------------------------------------!!
  subroutine waterbowl(veci,vecj,veck,x,wxyz,imask,iat,jat,kat)
  !! Modified to use three vectors instead of three atom indeces.
  !! This makes it general. It is still assumed, though, that the 
  !! three vectors are the atoms responsible for the mask 'emask'
  !! Sun Aug 12 16:45:08 EST 2001
  !! 
  !! Returns the re-entrant concave/concave "bowl" surface
  !! x for three atoms iat,jat,kat for the water at the 
  !! position of the vertex in emask. (returned as wxyz)
  !! also uses/returns current emask
  !! emask is returned as the surface mask around the water wxyz.
  !! Current imask on calling this subroutine is imask for +cvec
  !! ..so it need not be recalcuated, unless bowl is R-handed.
  !!  Mon May  7 11:56:47 EDT 2001
  !! 
  implicit none
  integer,intent(in) :: imask,iat,jat,kat
  real,dimension(3),intent(in) :: veci,vecj,veck
  real,intent(out) :: x, wxyz(3)
  integer :: i,j,k,ibyte,ibit,iphi,ipsi,itheta,jmask
  real :: y,theta, d
  real,dimension(3) :: avec,bvec,dvec,jvec,vjj,vkk
  !! uses globals cvec, emask
  !! ---------
  !! uses current value of cvec, which is the cross product
  !! iat->jat X iat->kat.  If cvec.(iat->wat) > 0, it's R-handed
  !! So the order is iat,jat,kat. Otherwise, it's L-handed
  !! and the order is iat,kat,jat
  j = 0
  jvec = 0.

  !!  Not necessary to check chirality. Done on calling
  ! do i=1,MAXATOM
  !  ibyte = (i-1)/NBIT + 1
  !  ibit = mod(i-1,NBIT)
  !  if (btest(emask(ibyte),ibit)) then
  !    j = i
  !    !! For more accurate placement of the water, we average all
  !    !! vertex points. This may slow things down. How much?
  !    !! Sat Jun 30 21:32:57 EDT 2001
  !    jvec = mxyz(1:3,i) 
  !    exit
  !    !! alternatively, we may compute the water position from the tetrahedron
  !    !! dimensions...
  !    !! If we are calculating derivatives, then we should do this because
  !    !! the mask need not be re-calculated 4 times. Whether it is R-handed 
  !    !! or L-handed only matters for plotting.
  !  endif
  !enddo
  !if (j==0) then
  !  wxyz = 0.
  !  x = 0.
  !  return
  !endif
  !d = sqrt(dotprod(jvec,jvec))
  !jvec = jvec/d     !! normalized, averaged water vector
  !! diagnostic
  ! write(*,*) 'emask not=0'
  !! wxyz is water position in absolute coords (A)
  !wxyz = (r1+rw)*jvec   !! current value of r1 set in calling routine
  !d = dotprod(cvec,wxyz)
  !! diagnostic
  !! avec = veci + wxyz       ! w
  !! write(*,'(3i5,4f9.3)') iat,jat,kat,avec(1:3),d
  if (periodic) then ! get nearest copy of jat, kat
    call boxaround(veci,vecj,lbox,bvec) ! bvec = i->j
    vjj = veci + bvec ! jat
    call boxaround(veci,veck,lbox,bvec) ! bvec = i->k
    vkk = veci + bvec ! kat
  else
    vjj = vecj
    vkk = veck
  endif
  ! if (d > 0.) then !! R-handed. Use iat,jat,kat and -cvec  !! NEW: It's always R-handed
    !! get water position using RH tetrahedron , replace d
    call tetrahedron(veci,vjj,vkk,wxyz,d,iat,jat,kat,1)
    if (drawing) then  !! write the water atoms
      avec = wxyz
      if (periodic) call inbox(avec,lbox)
      write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1,i7)') iat,"W",jat,avec(1:3),0.,d,kat
      write(13,'("TER")') 
    endif
    d = dotprod(cvec,wxyz)    !! overwrites d from tetrahedron
    !!
    if (abs(d) < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
      theta = acos(d/rw)
      itheta = nint((theta/rad)/DTHETA)
      if (d > 0) then  !! waterbowl called using i,j,k. use cvec. better yet, use global 'imask'
        jmask = imask  !!  NOTE: now using local variable jmask. 
      else             !! waterbowl called using i,j,k. use -cvec. 
        dvec = -1*cvec
        jmask = getimask(dvec)        ! dont need dvec after this?
      endif
      emask = masklib(:,itheta,jmask)
    else
      emask = -1    !! all ones
    endif
    !! mask ijw plane =================
    avec = veci - wxyz       ! w->i
    bvec = vjj - wxyz        ! w->j
    !! diagnostic
    ! write(*,'("a b ",6f9.4)') avec, bvec
    call cros(avec,bvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
    !! mask jkw plane =================
    avec = bvec          !  w->j
    bvec = vkk - wxyz    ! w->k
    call cros(avec,bvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
    !! mask kiw plane =================
    avec = bvec         ! w->k
    bvec = veci - wxyz  ! w->i
    call cros(avec,bvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
  ! else !! L-handed. Use iat,kat,jat and +cvec (current imask)
  !   !! get water position using RH tetrahedron
  !   call tetrahedron(veci,vjj,vkk,wxyz,d,iat,jat,kat,-1)
  !   !!
  !   if (drawing) then  !! write the water atoms
  !      write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') iat,"W",iat,wxyz(1:3),-1.,d
  !   endif
  !   if (d < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
  !     theta = acos(d/rw)
  !     itheta = nint((theta/rad)/DTHETA)
  !     jmask = imask  !!  NOTE: now using local variable jmask. 
  !     emask = masklib(:,itheta,jmask)
  !   else
  !     emask = -1
  !   endif
  !   !! mask ikw plane =================
  !   avec = veci - wxyz  !  w->i
  !   bvec = vkk - wxyz   ! w->k
  !   call cros(avec,bvec,dvec)
  !   dvec = dvec/sqrt(dotprod(dvec,dvec))
  !   jmask = getimask(dvec)
  !   emask = iand(masklib(:,NTHETA,jmask),emask)
  !   !! mask kjw plane =================
  !   avec = bvec          ! w->k
  !   bvec = vjj - wxyz    ! w->j
  !   call cros(avec,bvec,dvec)
  !   dvec = dvec/sqrt(dotprod(dvec,dvec))
  !   jmask = getimask(dvec)
  !   emask = iand(masklib(:,NTHETA,jmask),emask)
  !   !! mask jiw plane =================
  !   avec = bvec         ! w->j
  !   bvec = veci - wxyz  ! w->i
  !   call cros(avec,bvec,dvec)
  !   dvec = dvec/sqrt(dotprod(dvec,dvec))
  !   jmask = getimask(dvec)
  !   emask = iand(masklib(:,NTHETA,jmask),emask)
  ! endif
  j = countbits(emask)
  !! do i=1,MAXATOM
  !!   ibyte = (i-1)/NBIT + 1
  !!   ibit = mod(i-1,NBIT)
  !!   if (btest(emask(ibyte),ibit)) j = j + 1
  !! enddo
  x = wsphere*real(j)/real(MAXATOM)  !! area in A^2
  end subroutine waterbowl
  !!------------------------------------------------------------
  subroutine trianglemask(wxyz,ixyz,jxyz,kxyz,amask)
  !! The tetrahedron stored in awat is assumed to be right-handed.
  implicit none
  real,dimension(3),intent(in) :: wxyz,ixyz,jxyz,kxyz
  integer(kind=KND),dimension(MASKSIZE),intent(out) :: amask
  real, dimension(3) :: wivec,wjvec,wkvec,dvec
  real :: y
  integer :: jmask
  
  amask = -1     !! all ones
  if (periodic) then
    call boxaround(wxyz,ixyz,lbox,wivec) ! w->i
    call boxaround(wxyz,jxyz,lbox,wjvec) ! w->j
    call boxaround(wxyz,kxyz,lbox,wkvec) ! w->k
  else
    wivec = wxyz - ixyz
    wjvec = wxyz - jxyz
    wkvec = wxyz - kxyz
  endif

  !! mask ijw plane =================
  call cros(wivec,wjvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  amask = iand(masklib(:,NTHETA,jmask),amask)

  !! mask jkw plane =================
  call cros(wjvec,wkvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  amask = iand(masklib(:,NTHETA,jmask),amask)

  !! mask kiw plane =================
  call cros(wkvec,wivec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  amask = iand(masklib(:,NTHETA,jmask),amask)

  !! return spherical triangle mask
  end subroutine trianglemask
  !!------------------------------------------------------------
  integer function getimask(uvec)
  implicit none
  !! return the mask index for a unit vector
  !! uses pi,npsi,dpsi,psiposit
  !!   from the main program
  real,intent(in) :: uvec(3)
  integer :: iphi,ipsi
  real :: phi,psi,dphi
  psi = acos(uvec(3))
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi<=0) then
    if (uvec(3) < 0.) ipsi = npsi
    phi = 0.
  elseif (ipsi > npsi) then
    stop 'Bug in psi calculation: getimask'
  else
    phi = atan2(uvec(2),uvec(1))
  endif
  dphi = 2*pi/real(nphi(ipsi))
  !! diagnostic
  ! write(*,'(a,3f6.2,2f6.1,2i5,f6.3)') 'getimask: uvec, psi, phi, ipsi, iphi, dphi',uvec, psi, phi, ipsi, iphi, dphi
  ! if (dphi == 0.0) stop 'BUGBUGBUGBUGBUG'
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  getimask = psiposit(ipsi) + iphi
  end function getimask
  !!------------------------------------------------------------
  subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,iskip)
  !! rendering added using triangles
  !! Wed Aug 29 10:39:53 EST 2001
  implicit none
  integer,intent(in) :: iunit,iat,iskip
  real,dimension(3),intent(in) :: cen
  real,intent(in) :: r,bb
  real :: d
  character(len=1),intent(in) :: ch
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer :: i,j,ibyte,ibit
  real,dimension(3) :: avec,bvec,rgb
  real,parameter :: DCUT=6.0
  !
  j = 0
  do i=1,MAXATOM
    ibyte = (i-1)/NBIT + 1
    ibit = mod(i-1,NBIT)
    if (btest(amask(ibyte),ibit)) then
      j = j + 1
      if (mod(j,iskip)==0) then
        avec = cen + r*mxyz(1:3,i)
        if (ishow /= 0) then
          bvec = avec - showvec
          if (sqrt(dotprod(bvec,bvec)) > DCUT) cycle
        endif
        if (periodic) call inbox(avec,lbox)
        write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        i,ch,iat,avec(1:3),0.,bb
        write(iunit,'("TER")') 
      endif
    endif
  enddo
  if (rendering) then
    if (ch == "B") then      !! bowl
      rgb = (/0.5,0.5,0.0/)  !! yellow?
      call rendersphere(cen,r,amask,rgb)
    elseif (ch == "V") then  !! SASA
      rgb = (/0.0,1.0,0.0/)  !! green
      call rendersphere(cen,r,amask,rgb)
    elseif (ch == "E") then  !! edges
      rgb = (/0.0,0.0,1.0/)  !! blue
      return
    else
      rgb = (/1.0,1.0,1.0/)  !! white
      return
    endif
  endif
  end subroutine drawsurface
  !!==========================================================
  subroutine rendersphere(cen,r,amask,rgb)
  implicit none
  real,dimension(3),intent(in) :: cen,rgb
  real,intent(in) :: r
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer :: i,j,ibyte,ibit
  real,dimension(3) :: avec,bvec,cvec
  real,dimension(3) :: anor,bnor,cnor
  !! diagnostic
  ! write(*,*) "RENDERSPHERE",cen

  !! draw all triangles that have all bits set
  do j=1,ntri
    i = tri(1,j)
    if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
    avec = cen + r*mxyz(1:3,i)
    anor = mxyz(1:3,i)
    i = tri(2,j)
    if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
    bvec = cen + r*mxyz(1:3,i)
    bnor = mxyz(1:3,i)
    i = tri(3,j)
    if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
    cvec = cen + r*mxyz(1:3,i)
    cnor = mxyz(1:3,i)
    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
    write(runit,*) 1
    write(runit,'(12f8.3)') avec,bvec,cvec,rgb
    if (smoothing) then
      write(runit,*) 7
      write(runit,'(9f8.3)') anor,bnor,cnor
    endif
  enddo
  end subroutine rendersphere
  !!==========================================================
  !!    call rendersaddle(xyz(1:3,iat),xyz(1:3,kat),r1,amask,dmask,x,taumin)
  subroutine rendersaddle(va1,va2,r,amask,bmask,thet,taumin)
  implicit none
  real,dimension(3),intent(in) :: va1,va2
  real,intent(in) :: r,thet,taumin
  real :: dtau1,dtau2,thet1,thet2,mtheta
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask,bmask  
  !! amask is atommask, bmask is atomedgemask
  integer :: i,j,ibyte,ibit,k,itau,nn,ii,ntau
  integer,dimension(3) :: sides
  real,dimension(3,3) :: trv,trn
  real,dimension(3) :: avec,bvec,cvec,vab,uab,vp
  real,dimension(3) :: anor,bnor,cnor,rgb
  real,dimension(3) :: v1,v2,v3,vt,vw,vec
  real,dimension(3) :: w1,w2,w3,wt,ww,wec,wp
  real,parameter :: spcng=0.25   !! angstrom spacing of saddle dots, reset to 0.35
  real :: chi,x,mat(3,3),dtau
  ntau = nint(rw*(thet - taumin)/spcng)   ! number of steps in tau
  !! diagnostic
  !! write(*,*) "RENDERSADDLE",va1,ntau
  if (ntau <= 0 ) return
  dtau = (thet - taumin)/real(ntau)    ! step size in tau
  if (periodic) then                   ! get nearest copy of j
    call boxaround(va1,va2,lbox,vab) ! vab = a->b
  else
    vab = va2 - va1
  endif
  uab = vab/sqrt(dotprod(vab,vab))
  !!
  mtheta = (pi/2) - thet
  rgb = (/0.5,0.0,0.5/)  !! magenta?
  !! draw strips for all triangles that have
  !! a side on the edge
  do j=1,ntri    !! loop over all triangles
    nn = 0
    do k=1,3     !! check to see if any of the vectices are in the edgemask
      i = tri(k,j)
      if (btest(bmask((i-1)/NBIT + 1),mod(i-1,NBIT))) then
        nn = nn + 1
      endif
    enddo
    if (nn == 0) cycle  !! no points in edgemask
    nn = 0
    do k=1,3     !! find triangles that have 2 points in the atommask
      i = tri(k,j)
      if (btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT))) then
        nn = nn + 1
        sides(nn) = i
      endif
    enddo
    if (nn < 2) cycle   !! 0 or 1 vertex in the mask, forget it.
    if (nn > 2) cycle   !! if nn=3, whole triangle is already drawn.
    !! write(*,*) "found an edge-side",tri(1:3,j)
    !! a side of a triangle (adjacent points) were found in the edgemask
    !! draw a strip of surface starting from this side
    ii = sides(1)    ! first point
    vw = mxyz(1:3,ii) ! direction of mask point
    call cros(uab,vw,v2)    ! v2 perpendicular to abw plane
    vp = va1 + v2           ! vp point above a
    call getrotS(va1,vp,mtheta,mat,vec)  ! matrix to rotate uab to vw
    call rotate_S(mat,uab,vw)          ! this vw is exactly mtheta from uab
    vw = (r+rw)*vw          ! vw = water position relative to a
    v1 = va1 + vw           ! location of water
    !! call cros(vab,vw,v2)     ! v2 = perpendicular to abw
    call cros(vab,v2,vt)     ! perpendicular to vab in abw plane
    x = sqrt(dotprod(vt,vt))  ! length of vt
    vt = vt/x           ! normalized vt
    x = -dotprod(mxyz(1:3,ii),vt)  ! both are unit vectors
    vt = rw*vt           ! vt, length rw
    thet1 = acos(x)  ! cos(theta). vt always points away from vw
    dtau1 = (thet1 - taumin)/real(ntau)    ! step size in tau
    v2 = v1 - v2             ! point below water
    !! 
    ii = sides(2)    ! second point
    ww = mxyz(1:3,ii) ! direction of mask point
    call cros(uab,ww,w2)    ! v2 perpendicular to abw plane
    vp = va1 + w2           ! vp point above a
    call getrotS(va1,vp,mtheta,mat,vec)  ! matrix to rotate uab to vw
    call rotate_S(mat,uab,ww)          ! this vw is exactly mtheta from uab
    ww = (r+rw)*ww          ! vw = water position relative to a
    w1 = va1 + ww           ! location of water
    !! call cros(vab,ww,w2)     ! v2 = perpendicular to abw
    call cros(vab,w2,wt)     ! perpendicular to vab in abw plane
    x = sqrt(dotprod(wt,wt))  ! length of vt
    wt = wt/x           ! normalized wt
    x = -dotprod(mxyz(1:3,ii),wt)  ! both are unit vectors
    wt = rw*wt           ! vt, length rw
    thet2 = acos(x)  ! cos(theta). vt always points away from vw
    dtau2 = (thet2 - taumin)/real(ntau)    ! step size in tau
    w2 = w1 - w2             ! point below water
    nn = 0
    !!
    !! first point in first triangle
    chi = taumin 
    call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
    v3 = v1 + vt                     ! starting position
    call move_S(v3,mat,vec)          ! rotation
    nn = mod(nn,3)+1
    trv(1:3,nn) = v3
    trn(1:3,nn) = v1 - v3
    !!
    !! second point in first triangle
    call getrotS(w1,w2,chi,mat,vec)  ! matrix and vector for rotation
    w3 = w1 + wt                     ! starting position
    call move_S(w3,mat,vec)          ! rotation
    nn = mod(nn,3)+1
    trv(1:3,nn) = w3
    trn(1:3,nn) = w1 - w3
    !! diagnostic
    !! write(*,*) 'dtau1, dtau2, thet1, thet2 ', dtau1, dtau2, thet1, thet2
    !! vp = va1 + mxyz(:,sides(1)); wp = va1 + mxyz(:,sides(2))
    !! write(*,'(a,6f8.3)') "mxyz a,b ",vp,wp
    do itau=1,ntau           ! steps separated by spcng
      chi = taumin + itau*dtau
      !! rotate first point in side
      call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
      v3 = v1 + vt                     ! starting position
      call move_S(v3,mat,vec)          ! rotation
      nn = mod(nn,3)+1
      trv(1:3,nn) = v3
      trn(1:3,nn) = v1 - v3
      avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
      !! write a triangle
      write(runit,*) 1
      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,0.5,0.0,0.5
      if (smoothing) then
        write(runit,*) 7
        write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
      endif
      !! rotate second point in side
      chi = taumin + itau*dtau
      call getrotS(w1,w2,chi,mat,vec)  ! matrix and vector for rotation
      w3 = w1 + wt                     ! starting position
      call move_S(w3,mat,vec)          ! rotation
      nn = mod(nn,3)+1
      trv(1:3,nn) = w3
      trn(1:3,nn) = w1 - w3
      avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
      !! write a triangle
      write(runit,*) 1
      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,0.5,0.5,0.0
      if (smoothing) then
        write(runit,*) 7
        write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
      endif
    enddo
    vp = va1 + r*mxyz(:,sides(1))
    nn = mod(nn,3)+1
    trv(1:3,nn) = vp
    trn(1:3,nn) = mxyz(:,sides(1))
    avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
    !! write a triangle
    write(runit,*) 1
    write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,0.0,1.0,0.0
    if (smoothing) then
      write(runit,*) 7
      write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
    endif
    vp = va1 + r*mxyz(:,sides(2))
    nn = mod(nn,3)+1
    trv(1:3,nn) = vp
    trn(1:3,nn) = mxyz(:,sides(1))
    avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
    !! write a triangle
    write(runit,*) 1
    write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,0.0,1.0,0.0
    if (smoothing) then
      write(runit,*) 7
      write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
    endif
    !! write(*,'(a,6f8.3)') "last a,b ",trv(:,mod(nn+1,3)+1),trv(:,nn)
  enddo
  end subroutine rendersaddle
  !!==========================================================
  subroutine inbox(avec,L)
  !! Put coords avec in a cubic box of size L
  !! then scale by BOXSCALE (parameter)
  real,dimension(3),intent(inout) :: avec
  real,intent(in) :: L
  integer :: i
  do i=1,3
    if (avec(i) < 0.) avec(i) = avec(i) + L
    if (avec(i) > L) avec(i) = avec(i) - L
  enddo
  avec = BOXSCALE*avec
  end subroutine inbox
  !!==========================================================
  subroutine threeinbox(avec,bvec,cvec,L)
  !! Put coords of triangle in a cubic box of size L
  !! then scale by BOXSCALE (parameter)
  real,dimension(3),intent(inout) :: avec,bvec,cvec
  real,intent(in) :: L
  real,dimension(3) :: dvec
  integer :: i
  dvec = (avec + bvec + cvec)/3.0
  do i=1,3
    if (dvec(i) < 0.) then
      avec(i) = avec(i) + L; bvec(i) = bvec(i) + L;cvec(i) = cvec(i) + L
    endif
    if (dvec(i) > L) then
      avec(i) = avec(i) - L; bvec(i) = bvec(i) - L;cvec(i) = cvec(i) - L
    endif
  enddo
  end subroutine threeinbox
  !!==========================================================
  subroutine scalethreeinbox(avec,bvec,cvec,L)
  !! Put coords of triangle in a cubic box of size L
  !! then scale by BOXSCALE (parameter)
  real,dimension(3),intent(inout) :: avec,bvec,cvec
  real,intent(in) :: L
  real,dimension(3) :: dvec
  integer :: i
  call threeinbox(avec,bvec,cvec,L)
  avec = BOXSCALE*avec
  bvec = BOXSCALE*bvec
  cvec = BOXSCALE*cvec
  end subroutine scalethreeinbox
  !!==========================================================
  real function distmod(avec,bvec,L)
  !! project bvec into the cube os size L around avec
  !! and report that distance
  real,dimension(3),intent(in) :: avec,bvec
  real,intent(in) :: L
  real,dimension(3) :: cvec
  integer :: i,j
  real :: dd

  cvec = bvec - avec
  if (periodic) then
    do i=1,3
      if (cvec(i) > (L/2.)) cvec(i) = cvec(i) - L
      if (cvec(i) <= -(L/2.)) cvec(i) = cvec(i) + L
    enddo
  endif
  dd = sqrt(dotprod(cvec,cvec))
  distmod = dd
  end function distmod
  !!==========================================================
  subroutine boxaround(avec,bvec,L,cvec)
  real,dimension(3),intent(in) :: avec,bvec
  real,dimension(3),intent(out) :: cvec
  real,intent(in) :: L    ! global variable lbox
  integer :: i
  cvec = bvec - avec
  do i=1,3
    if (cvec(i) > (L/2.)) cvec(i) = cvec(i) - L
    if (cvec(i) <= -(L/2.)) cvec(i) = cvec(i) + L
  enddo
  end subroutine boxaround
  !!==========================================================
  subroutine drawsaddle(iunit,avec,bvec,smask,r1,iat,taumin,taumax,bb,ch)
  implicit none
  !! output the saddle surface as dots
  !! Wed Jun 20 14:35:42 EDT 2001
  real,parameter :: spcng=0.25   !! angstrom spacing of saddle dots, reset to 0.35
  character(len=1),intent(in) :: ch
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: smask
  real,intent(in) :: taumin,r1,bb,taumax,avec(3),bvec(3)
  integer,intent(in) :: iat,iunit
  integer :: i,j,itau,ntau,ibit,ibyte
  real,dimension(3) :: v1,v2,v3,vt,vw,vab,vec
  real :: chi,dtau,x,mat(3,3)
  !!
  dtau = 0.
  ntau = nint(rw*(taumax - taumin)/spcng)   ! number of steps in tau
  if (ntau > 0) dtau = (taumax - taumin)/real(ntau)    ! step size in tau
  if (periodic) then                   ! get nearest copy of j
    call boxaround(avec,bvec,lbox,vab) ! vab = a->b
  else
    vab = bvec - avec
  endif
  !! diagnostic
  ! write(*,*) "DRAWSADDLE:",vab,ntau,dtau,taumin,taumax
  if (ntau == 0) return
  do i=1,MAXATOM
    ibyte = (i-1)/NBIT + 1
    ibit = mod(i-1,NBIT)
    if (btest(smask(ibyte),ibit)) then
      j = j + 1
      vw = (r1+rw)*mxyz(1:3,i) ! direction of water (rw, mxyz are global)
      v1 = avec + vw           ! location of water
      call cros(vab,vw,v2)     ! v2 = perpendicular to abw
      call cros(vab,v2,vt)     ! perpendicular to vab in abw plane
      x = sqrt(dotprod(vt,vt))  ! length of vt
      vt = (rw/x)*vt           ! normalized vt, length rw
      v2 = v1 - v2             ! point below water
      do itau=0,ntau-1         ! steps separated by spcng
        chi = taumin + itau*dtau
        call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
        v3 = v1 + vt                     ! starting position
        call move_S(v3,mat,vec)          ! rotation
        if (periodic) call inbox(v3,lbox)
        write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        i,ch,iat,v3(1:3),0.,bb
        write(iunit,'("TER")') 
      enddo
    endif
  enddo
  end subroutine drawsaddle
  !!==========================================================
  subroutine renderheader(cen,scale)
  real,dimension(3),intent(in) :: cen
  real,intent(in) :: scale
  write(runit,*) 'Masker v.30-AUG-01'
  write(runit,*) '20 20     tiles in x,y'
  write(runit,*) '16 16     pixels (x,y) per tile'
  write(runit,*) '4         anti-aliasing level 4; 3x3->2x2'
  write(runit,*) '0 0 0     black background'
  write(runit,*) 'F         no shadows cast'
  write(runit,*) '25        Phong power'
  write(runit,*) '0.25      secondary light contribution'
  write(runit,*) '0.05      ambient light contribution'
  write(runit,*) '0.25      specular reflection component'
  write(runit,*) '4.0       eye position'
  write(runit,*) '1 1 1     main light source position (from over right shoulder)'
  write(runit,*) '1 0 0 0   view matrix describing input coordinate transformation'
  write(runit,*) '0 1 0 0'
  write(runit,*) '0 0 1 0'
  write(runit,*) cen, scale
  write(runit,*) '3         mixed objects'
  write(runit,*) '*        (free format triangle and plane descriptors)'
  write(runit,*) '*        (free format sphere descriptors)'
  write(runit,*) '*        (free format cylinder descriptors)'
  end subroutine renderheader
  !!==========================================================
  subroutine tetrahedron(veci,vecj,veck,wxyz,d,iat,jat,kat,hand)
  implicit none
  !! NOTE!!: THIS ROUTINE FAILS IF THE TETRAHEDRON INEQUALITY IS VIOLATED.
  !! (...it's like the triangle inequality...)
  !! ===> NOW checks tetrahedron inequality, returns d=0 if it fails.Tue Nov 20 09:16:25 EST 2001
  !!
  !! Now returns d=relative distance from ijk plane (<0 if hand=-1). Mon Oct  1 15:16:55 EDT 2001
  !!---------------------
  !! Using three atoms and a set of three distances, 
  !! get the apex of a tetrahedron based on the three atoms
  !! and having the three distances to the apex. The tetrahedron
  !! is assumed to be R-handed. To change handedness, switch two atoms.
  !! Mon Aug 13 10:58:06 EST 2001
  integer,intent(in) :: iat,jat,kat,hand
  real,dimension(3),intent(in) :: veci,vecj,veck
  real,dimension(3),intent(out) :: wxyz
  real,intent(out) :: d
  real,dimension(3) :: vec,wvec
  real,dimension(3,3) :: mat
  real :: area,dij,djk,dik,diw,djw,dkw,ajik,ajiw,akiw,x,y,yp,z,dd
  integer :: flag
  !!
  !! write(*,*) "In tetrahedron: veci=",veci
  !! write(*,*) "In tetrahedron: vecj=",vecj
  !! write(*,*) "In tetrahedron: veck=",veck
  vec = vecj - veci
  dij = sqrt(dotprod(vec,vec))
  vec = veck - vecj
  djk = sqrt(dotprod(vec,vec))
  vec = veck - veci
  dik = sqrt(dotprod(vec,vec))
  diw = atomlib(atype(iat))%r + rw
  djw = atomlib(atype(jat))%r + rw
  dkw = atomlib(atype(kat))%r + rw
  !! write(*,*) "In tetrahedron: dij,djk,dik,diw,djw,dkw=",dij,djk,dik,diw,djw,dkw
  d = 0.0  !! this is a flag indicating that tetrahedron rule failed
  call cosinerule(dij,dik,djk,ajik,flag); if (flag /= 0) return
  call cosinerule(dij,diw,djw,ajiw,flag); if (flag /= 0) return
  call cosinerule(dik,diw,dkw,akiw,flag); if (flag /= 0) return
  !! 
  ! write(*,*) "dij,dik,djk ",dij,dik,djk
  ! write(*,*) "ajiw=",ajiw*180./3.14159," akiw=",akiw*180./3.14159," ajik=",ajik*180./3.14159
  x = diw*cos(ajiw)
  yp = diw*cos(akiw)
  y = yp/sin(ajik) - x/tan(ajik)
  d = sqrt(diw*diw - x*x - y*y)
  if (hand < 0) d = -d
  wxyz = (/x,y,d/)
  !  write(*,*) "In tetrahedron: wxyz=",wxyz
  !! diagnostic : distance to water
  !  vec = wxyz ;  dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - (/dij,0.,0./) ;  dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - (/dik*cos(ajik),dik*sin(ajik),0./); dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5)') dd
  !  write(*,*) dik*cos(ajik),dik*sin(ajik),0.
  
  !! Get frame
  call getframe(veci,vecj,veck,mat,vec)
  call move_S(wxyz,mat,veci)
  !!
  !  vec = wxyz - veci; dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - vecj; dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - veck; dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5)') dd
  !  write(*,*) "water pos = ",wxyz
  !  mat = transpose(mat)
  !  vec = 0.
  !  wvec = vecj - veci
  !  call move_S(wvec,mat,vec)
  !  write(*,*) "rotated vecj = ",wvec
  !  wvec = veck - veci
  !  call move_S(wvec,mat,vec)
  !  write(*,*) "rotated veck = ",wvec
  
  end subroutine tetrahedron
  !!----------------------------------------------------------------!!
  subroutine allocatoms(nat)
  !! This routine may be necessary if private module arrays
  !! are to be allocated.
  !! Here the private derived type 'atype' is allocated and
  !! an array of water locations (tetrahedron verteces) is allocated.
  !! Since we don't know how many waters will be needed, we allocate 
  !! 'nat', one per atom. That should be more than enough.
  integer,intent(in) :: nat
  integer :: ios=0
  if (allocated(atype)) deallocate(atype)
  allocate(atype(nat),stat=ios)
  if (ios/=0) stop 'Error allocating atype'
  if (allocated(water)) deallocate(water)
  allocate(water(watfac*nat),stat=ios)   !! this is more than enough. what should it really be?
  if (ios/=0) stop 'Error allocating water'
  end subroutine allocatoms
  !!----------------------------------------------------------------!!
  ! subroutine initcountbits
  ! implicit none
  !! Initialize the array 'cbits' which is the number of
  !! bits set to one in the byte argument
  !! integer,dimension(-128:127) :: cbits
  ! byte :: ib
  ! integer :: j,k
  ! do ib=-128,127
  !   k = 0
  !   do j=0,7
  !    if (btest(ib,j)) then
  !       k = k + 1
  !    endif
  !   enddo
  !   cbits(ib) = k
  ! enddo
  ! end subroutine initcountbits
  !!----------------------------------------------------------------!!
  integer function countbits(amask)
  implicit none
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer :: i,j,ib
  if (KND==2) then
    j = 0
    do i=1,MASKSIZE,8
      j = j + cbits(amask(i)) + cbits(amask(i+1)) +  &
              cbits(amask(i+2)) + cbits(amask(i+3)) +  &
              cbits(amask(i+4)) + cbits(amask(i+5)) +  &
              cbits(amask(i+6)) + cbits(amask(i+7))
    enddo
  elseif (KND==4) then
    j = 0
    do i=1,MASKSIZE
      do ib=0,31
        if (btest(amask(i),ib)) j = j + 1
      enddo
    enddo
  else
    stop 'countbits doesn t understand KND setting!!! should be 2 or 4'
  endif
  countbits = j
  end function countbits
  !!----------------------------------------------------------------!!
  integer function maskj(imask,itheta)
  implicit none
  integer,intent(in) :: imask,itheta
  integer(kind=KND),dimension(MASKSIZE) :: amask
  amask = masklib(:,itheta,imask)
  maskj = countbits(amask)
  end function maskj
  !!----------------------------------------------------------------!!
  !! NOTE: This routine (wordbits) gives a compilation WARNING because it is called
  !! using a 4-byte (normal) integer, but the dummy argument
  !! is a 4*byte array. It seems to work. Is there a legal way
  !! to do this? The old way would be 'equivalence'
  !!---------------------------------------------------------------!!
  ! integer function wordbits(w)
  ! implicit none
  ! integer(kind=1),dimension(4),intent(in) :: w
  ! integer :: i,j
  ! wordbits = cbits(w(1)) + cbits(w(2)) + cbits(w(3)) + cbits(w(4))
  ! end function wordbits
  !!----------------------------------------------------------------!!

  !----------------------------------------------------!
  subroutine ranvec(vec,sig,nseed)   !! used by mcmask
  implicit none
  !! contained in program that uses masker module
  real,dimension(3),intent(out) :: vec
  integer,intent(inout) :: nseed
  real,intent(in) :: sig
  real :: phi,psi,x,y,z
  integer :: imask
  !! chose random point on the surface of the template mask
  x = ran(nseed)
  imask = nint(x*MAXATOM)
  x = sig*ran(nseed)
  vec = x*mxyz(1:3,imask)
  end subroutine ranvec
  !----------------------------------------------------!
  subroutine getimask2(side1,side2,sideopp,thet,uvec,imask,itheta)
  implicit none
  real,intent(in) :: side1,side2,sideopp
  real,dimension(3),intent(in) :: uvec
  real,intent(out) :: thet
  integer,intent(out) :: imask,itheta
  integer :: flag
  !! logical :: embedded

  call cosinerule(side1,side2,sideopp,thet,flag)
  if (flag == 1) then
    if (thet == 0.0) then  !! bad triangle. sideopp too small
      imask = 0
    else
      imask = -1   !! sideopp too big, embedded
    endif
    return
  elseif  (thet > pi/2) then  !! partly embedded, get mask index in the opposite direction.
    itheta =  nint(((pi-thet)/rad)/DTHETA)    !! itheta index for embedded atom
    imask = getimask(-uvec)  !! uvec must be a unit vector
  else      !! normal triangle, not embedded
    itheta = nint((thet/rad)/DTHETA) 
    if (itheta == NTHETA) then  
      imask = getimask(-uvec)  !! close to being embedded. Use inverse mask for edges.
    else
      imask = getimask(uvec)  !! uvec must be a unit vector
    endif
  endif
  end subroutine getimask2
  !----------------------------------------------------!
  ! this routine is just a utility. Dont use it for saddle calculation  !!
  ! real function getsaddle(r1,rw,theta,arcfrac)
  ! real,intent(in) :: r1,rw,theta,arcfrac
  ! real :: x,y,col,taumin
  ! y = (r1+rw)*sin(theta)
  ! col = 2*pi*(arcfrac)   !! this is the length of the exposed arc in radians
  ! x = col*((0.5*pi-theta)*y - rw*cos(theta))
  ! taumin = 0.
  ! if (y < rw) then
  !    taumin = acos(y/rw)
  !    x = x + col*(rw*sin(taumin) - taumin*y)
  ! endif
  ! getsaddle = x
  ! end function getsaddle

  end module masker
