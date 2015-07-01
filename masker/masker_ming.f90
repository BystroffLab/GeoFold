  !! Module version of pdbmask.f90
  !! masker_initmasks -- allocates memory, reads masks, gets collarsizes etc.
  !! masker_getms -- accepts coordinates and atomtypes, returns MS
  !!
  !! This module sets the public and private variables used by the subroutines in masker.f90
  !! Some of the parameter setting cannot be changed. runit and watfac are the exceptions.
  !!
  !! Tue Jun 11 13:24:45 EDT 2002 C.Bystroff
  !! This version has been validated against MSMS. Any changes to this 
  !! algorithm should be validated by running the three-atom trajectory 
  !! (three_mask.csh) and the multi-size atom sets (doall.csh).
  !!
  !! Fri Jun 15 19:55:28 EDT 2001 C.Bystroff
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
  !!
  !! =====Modification=====
  !! 2006/12/10 Yao-ming Huang
  !!   1. Increase NSAVE to 200 

  !** ====================================================================
  !** Copyright 2001-2002 Chris Bystroff
  !**
  !** Permission is hereby granted to use, execute, copy, distribute, modify,
  !** and distribute in modified form this software for non-profit purposes,
  !** provided this notice is retained in any and all copies.
  !**
  !** For for-profit usage, please contact bystrc@rpi.edu
  !** ====================================================================

  module masker
  use vectormath
  implicit none
  !!
  !! integer,parameter :: kind_8=selected_real_kind(P=8)
  private
  integer,parameter :: MAXATOM=4096   !! number of points in mask. 
  integer,parameter :: KND=2          !! bytes per mask word
  integer,parameter :: NBIT=(KND*8)
  integer,parameter :: ILO=-(2**(15)),IHI=(2**(15))-1
  integer,parameter :: MASKSIZE=MAXATOM/NBIT    !! 3rd index of mask
  integer,parameter :: KC=1   !! thickness of collars (times DTHETA)
  integer,parameter :: runit=14  !! output unit for rendering
  integer,parameter :: watfac=40  !! max numbers of waters per atom 
  integer,parameter :: NSAVE=200  !! this is the number of nearest neighbor atoms to save
  real,parameter,private :: pi=3.1415927410125732421875
  real,parameter :: radius=10.
  real,parameter,private :: rad=pi/180.0
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
  integer,dimension(:),allocatable,target :: atype
  integer,dimension(:),allocatable,private :: nposit_mask
  integer,private :: nres_mask
  real :: dpsi 
  !! uses pi,npsi,dpsi,psiposit
  real :: maxd,taumin,r1,lbox=0.0
  real,dimension(3) :: showvec
  logical :: drawing=.false.
  logical :: periodic,rendering=.false.,smoothing=.false.,saying=.true.
  integer(kind=1),private,dimension(ILO:IHI) :: cbits
  type watertype
    real,dimension(3) :: xyz      !! coordinates of a probe position for a bowl
    integer,dimension(3) :: trng  !! indeces of base triangle for a water
    logical :: good               !! if this is false, don't use it for SES
  end type watertype
  type(watertype),dimension(:),allocatable,private :: water
  integer,private :: nwat  !! number of verteces
  integer,private,parameter :: surfacetension=1,byatom=2,byshape=3,bycharge=4
  integer,private :: iwat,jwat ,coloring=bycharge
  type atomtype
    real :: r        ! united atom radius
    real :: m        ! united atom mass
    real :: w1,w2,w3 ! convert surface area to energy (sasa, saddle, bowl) kJ/mol/A^2
    character(len=5) :: name   !  element or other atom name
  end type atomtype
  type(atomtype),dimension(:),allocatable :: atomlib
  integer :: nattype
  public masker_getms,masker_initmasks,masker_allocatoms,masker_getatypeptr,masker_deallo
  !!INTERFACE
  !!  real function dotprod(a,b) 
  !!    real,dimension(3) :: a, b
  !!  end function dotprod
  !!end INTERFACE
CONTAINS
!=========================================================================!
  subroutine usage    !! show current environment variables
    character(len=80) :: aline
    write(*,*) "The following environment variables should be set:"
    write(*,*) "set  ______=   <current setting>"
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
  subroutine usehmmstrprm(n,npos)
    !!-------------------------
    !! Call this only if using HMMSTRPRM module combined with MASKER module
    !! this allocates a private copy of nposit for MASKER
    !!-------------------------
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: npos
    integer :: ios
    if (allocated(nposit_mask)) stop 'masker.mod - ERROR: please do not call usehmmstrprm twice'
    allocate(nposit_mask(n+1),stat=ios); if (ios/=0) stop 'masker.mod - ERROR: allocating nposit_mask()'
    nposit_mask(1:n+1) = npos(1:n+1)
    nres_mask = n
  end subroutine usehmmstrprm
!-------------------------------------------------------------------------!
  integer function getresidue(iat)
    implicit none
    integer,intent(in) :: iat  !! atom number
    integer :: i,j
    getresidue = 0
    if (allocated(nposit_mask)) then !! this is a peptide
      i = 1
      do  while (nposit_mask(i) <= iat)
        i = i + 1
        if (i > nres_mask+1) return
      enddo
      getresidue = i - 1
    endif
  end function getresidue
!-------------------------------------------------------------------------!
  subroutine initmc
  call random_seed()
  end subroutine initmc
!-------------------------------------------------------------------------!
  subroutine masker_initmasks
  !! ======================== INITMASKS =================================
  !! This routine reads the masker data files and sets numerous 
  !! constants. It must be called (only once!) before calling 'masker_getms'.
  !! ====================================================================
  
  character(len=80) :: aline,pdbfile,binfile,logfile,outfile,template,cbitfile
  character(len=80) :: vmaskfile,forcefile,vtriangle,collarfile,slopefile,fmt
  character(len=1000) :: design_dir
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
  !! ====================================================================
  !! GET THE NAMES of the masker data files from the environment variable
  !! settings. This allows the user to define a complete path for
  !! each file. If no setting is found, then a default filename is used.
  !! NOTE: this filename has no directory path, so if you accept the
  !! default settings then you must be running
  !! the program in the directory where the files are kept!! 
  !! ANOTHER NOTE: an alternative to using the environment variables
  !! is to make a symbolic link (ln -s) in the current directory
  !! using the default settings as the link names  i.e.
  !! ln -s /home/database/masker/4096.mas  4096.mas
  !!
  !! Default settings:
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
  !! Warning:
  !! 'getenv' is probably NOT standard fortran90, but it works with pgf90.
  !!
  call getenv("DESIGN_HOME", design_dir)
  if (design_dir==" ") stop 'masker_initmasks - ERROR: environment variable DESIGN_HOME needed'
  call getenv("MASKLIB",binfile)
  if (binfile=="") binfile = trim(design_dir) // "/" // "4096.mas"
  call getenv("VMASK",vmaskfile)
  if (vmaskfile=="") vmaskfile = trim(design_dir) // "/" // "4096.vmask"
  call getenv("VTRIANGLE",vtriangle)   !! triangle indeces
  if (vtriangle=="") vtriangle = trim(design_dir) // "/" // "4096.tri"
  call getenv("MASKTEMPLATE",template)
  if (template=="") template = trim(design_dir) // "/" // "4096.pdb"
  call getenv("MASKDAT",logfile)
  if (logfile=="") logfile = trim(design_dir) // "/" // "4096.dat"
  call getenv("FORCEFIELD",forcefile)
  if (forcefile=="") forcefile = trim(design_dir) // "/" // "forcefield.prm"
  call getenv("COUNTBITS",cbitfile)
  if (cbitfile=="") cbitfile = trim(design_dir) // "/" // "cbits32.bin"
  call getenv("COLLARS",collarfile)
  if (collarfile=="") collarfile = trim(design_dir) // "/" // "collars4096.bin"
  call getenv("SLOPES",slopefile)
  if (slopefile=="") slopefile = trim(design_dir) // "/" // "slopes4.5.dat"

  !! ====================================================================
  !! COUNTBITS
  !! initialize bit counting array, precalculated for KND-byte words (KND=2)
  !call initcountbits
  open(iunit,file=cbitfile,form="unformatted",status="old",iostat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: missing COUNTBITS file'
  read(iunit,iostat=ios) cbits(ILO:IHI)
  if (ios/=0) stop 'masker_initmasks - ERROR: reading COUNTBITS file'
  close(iunit)

  !! ====================================================================
  !! ATOMLIB
  !! Read atom types, radii, weights, etc. from ATOMLIB file
  !! Surface tensions in kJ/A^2 are stored in this file.
  !! For example, source of surface tension for CH4 was
  !! A.P.Lyubartsev , O.Forrisdahl and A.Laaksonen
  !! J. Chem. Phys., v.108(1), pp.227-233, (1998) 
  !! This paper reports various solvation studies using water models and methane models.
  !! Numbers range from 8.5 to 15.35 kJ/Mol
  !! methane diameter is 3.73A, so the surfae area is 43.7 A^2
  !! We can pick 11 kJ/Mol as a median solvation free energy for methane. Then
  !! the surface tension is 11 kJ/Mol / 43.7A^2 = 0.252 kJ/mol/A^2
  !!---------------------------------------------------------------
  open(junit,file=forcefile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: missing FORCEFIELD file'
  aline=" "
  do while (aline(1:22)/="!-----Atom types-----!")
    read(junit,'(a)',iostat=ios) aline
  enddo
  read(junit,*,iostat=ios) aline, nattype
  if (ios/=0) stop 'masker_initmasks - ERROR: reading FORCEFIELD file(1)'
  if (allocated(atomlib)) stop 'masker_initmasks - ERROR: please do not call masker_initmasks twice'
  allocate(atomlib(nattype),stat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating atomlib()'
  !! Thu Aug 30 07:39:02 EST 2001 new atomlib (from masker2.f90)
  do while (aline(1:2)/="!#")
    read(junit,'(a)',iostat=ios) aline
  enddo
  do i=1,nattype
    read(junit,*,iostat=ios) aline, atomlib(i)%name, atomlib(i)%m, atomlib(i)%r,&
           atomlib(i)%w1,  atomlib(i)%w2,  atomlib(i)%w3
    if (ios/=0) stop 'masker_initmasks - ERROR: reading FORCEFIELD file(2)'
    !change units
    atomlib(i)%r=atomlib(i)%r*10.0
    atomlib(i)%w1=atomlib(i)%w1/100.0;atomlib(i)%w2=atomlib(i)%w2/100.0;atomlib(i)%w3=atomlib(i)%w3/100.0
  enddo
  close(junit)
  !! ====================================================================
  !! MASKDAT
  !! Read log file from binarymask.f90 . This contains indexing data for masks
  !! allocate memory based on these numbers.
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
      if (ios/=0) stop 'masker_initmasks - Format error 1'
      DTHETA = 180./npsi
      NTHETA = int(90./DTHETA)
      !! write(*,*) 'DTHETA=',DTHETA,'  NTHETA=',NTHETA
      allocate(psiposit(0:npsi),nphi(0:npsi),stat=ios)
      if (ios/=0) stop 'masker_initmasks - ERROR: allocating psiposit(), nphi()'
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
      if (ios/=0) stop 'masker_initmasks - Format error 2'
      read(aline(32:38),*,iostat=ios) psiposit(ipsi)
      if (ios/=0) stop 'masker_initmasks - Format error 3'
      !! write(*,*) nphi(ipsi), psiposit(ipsi)
      ipsi = ipsi + 1
    endif
    if (aline(1:5)=="NMASK".and.npsi/=0) then
      read(aline(7:),*,iostat=ios) nmask
      if (ios/=0) stop 'masker_initmasks - Format error 4'
    endif
  enddo
  close(11)
  !! Report the size of the mask and allocate memory accordingly
  ! write(*,'("NPSI =",i7)') npsi
  !! write(*,'("NMASK=",i7)') nmask
  allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating masklib()'
  allocate(maskpsi(nmask),maskphi(nmask),stat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating maskpsi(), maskphi()'
  allocate(collarsize(0:NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating collarsize()'
  allocate(slope(NTHETA),interc(NTHETA),stat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating slope(), interc()'
  ! write(*,*) 'Memory allocated for masklib',MASKSIZE, NTHETA, nmask

  !! ====================================================================
  !! MASKLIB
  !! read binary format mask library, created by binarymask.f90
  open(12,file=binfile,status='old',form='unformatted')
  read(12,iostat=ios) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
  if (ios/=0) stop 'masker_initmasks - ERROR: reading MASKLIB library'
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
  !!!! Diagnostic
  ! write(*,*) 'Done reading masklib'

  !! VMASK
  !! Get viewable dots mask . This is the mask for outputing
  !! a regular subset of dots in PDB format
  open(12,file=vmaskfile,status='old',form='unformatted')
  read(12,iostat=ios) vmask(1:MASKSIZE)
  if (ios/=0) stop 'masker_initmasks - ERROR: reading VMASK file'
  close(12)

  !! ====================================================================
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
  if (ios/=0) stop 'masker_initmasks - ERROR: allocating tri()'
  ! write(*,*) ntri, ' Triangles for rendering.'
  rewind(12)
  do i=1,ntri
    read(12,*,iostat=ios) tri(1:3,i)
    if (ios/=0) stop 'masker_initmasks - ERROR: re-reading triangles'
  enddo
  close(12)
  !! Diagnostic
  ! write(*,*) 'Done reading vmask file'
  !!Diagnostic
  ! write(*,*)  masklib(1:MASKSIZE,10,50)

  !! ====================================================================
  !! CALCULATE POLAR COORDINATES PSI AND PHI FOR EACH MASK
  !! ..AND..
  !! Calculate collarsize and tau for each theta and imask
  !! collarsize == the number of bits in the edge mask for angle theta
  !! tau        == the angular width of the half-saddle for angle theta and r1
  !!
  !! PSI   ==  the angle from the positive Z-axis of the normal to the bisecting plane (latitude)
  !! PHI   ==  the angle from the positive X-axis of the projection on the XY-plane of the normal to the bisecting plane (longitude)
  !!
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
       !! Diagnostic
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
  
  !! ====================================================================
  !! COLLARS
  !! read collarsizes from file. These are the precalculated sizes of edgemasks.
  !! open(33,file=collarfile,status="replace",form="unformatted",iostat=ios)
  open(33,file=collarfile,status="old",form="unformatted",iostat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: missing COLLARS file'
  !! write(33) collarsize(0:NTHETA,1:nmask)  !! use this to re-generate the file
  read(33,iostat=ios) collarsize(0:NTHETA,1:nmask)
  if (ios/=0) stop 'masker_initmasks - ERROR: reading COLLARS file'
  close(33)
  !! stop  'Done saving collarsizes'
  !! NTHETA is 90deg, only possible if dd ~= 0. No saddle? 
  !! =====> NOTE: If hydrogens are used, this will have to be changed! <=====

  !! ====================================================================
  !! SLOPES
  !! Read the calibration data from file: slope of A^2 vs theta (radians) for a unit sphere
  !! indexed by itheta, which is incremented by dpsi (dtheta). 
  !! These slopes are used for interpolating linearly between masks.
  open(34,file=slopefile,status='old',form='formatted',iostat=ios)
  if (ios/=0) stop 'masker_initmasks - ERROR: missing SLOPES file'
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

  !! ====================================================================
  !! MASKTEMPLATE
  !! Coordinates of mask points on a 10A sphere.
  !! These are only needed if drawing=.true.
  if (drawing) then
    allocate(mxyz(3,MAXATOM),stat=ios)
    if (ios/=0)  stop 'masker_initmasks - ERROR: allocating mxyz()'
    open(11,file=template,status='old',form='formatted',iostat=ios)
    if (ios/=0) stop 'masker_initmasks - ERROR: missing MASKTEMPLATE file'
    i = 0
    do 
      read(11,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:5) /= "ATOM ") cycle
      i = i + 1
      if (i > MAXATOM) stop 'masker_initmasks - ERROR: Too many points in MASKTEMPLATE file'
      read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
      if (ios/=0) then
        write(*,*) 'ERROR: reading template coordinate file, at i=',i
        stop 'masker_initmasks - ERROR: reading template coordinate from MASKTEMPLATE file'
      endif
      mxyz(1:3,i) = mxyz(1:3,i)/10.   !! template has 10A radius
    enddo
    close(11)
  endif
  !! 
  end subroutine masker_initmasks
  !!----------------------------------------------------------------!!
  subroutine masker_getms(xyz,nat,ms,buried,gradient,start,SAS)
  !! ======================== GETMS ====================================
  !! This is the subroutine that you call to get the SES, which is
  !! returned in 'ms'. nat is the number of atoms. atype are their atom 
  !! types, and xyz are the coordinates in orthogonal angstroms.
  !! ===================================================================
  !! Added gradient to return the shift vector for each atom
  !! Sun Oct 26 09:11:31 EST 2003  CB
  implicit none
  integer,intent(in) :: nat
  real,intent(out) :: ms
  logical,intent(in) :: SAS
  real,dimension(3,nat),intent(in) :: xyz
  real,dimension(nat),intent(out),optional :: buried
  real(kind=kind_8),dimension(3,nat),intent(inout),optional :: gradient
  integer,intent(in),optional :: start
  !
  integer :: i,j,jj,k,L,m,jarg,ios,iargc,natm,ii,kwat,istart
  integer(kind=KND),dimension(MASKSIZE) :: amask,bmask,cmask,dmask
  integer :: nn,ipsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv
  integer :: imask1,itheta1,imask2,itheta2,iat,jat,kat,ires,jtri,jxtri,jcap,jarc
  integer,dimension(NSAVE) :: saveimask, saveitheta,saveiat
  real,dimension(NSAVE) :: savetheta, savedij, saver2
  real :: phi,psi,theta,x,y,d,dd,r2,r3
  real,parameter :: delt=0.01
  integer :: kskip=100,kcc=KC, flag
  real :: dphi,sphere,col,sasanrg,ssasanrg,bsasanrg,xx,arcfrac,abur,acap
  real :: dsasa,c2,dij,th,da,asc,sdlarea=0.,scale=10.
  real,dimension(3) :: avec,bvec,wxyz,showvec,veci,vecj,veck,iv,jv,kv,wv,dvec
  real,dimension(3) :: vij,cij,vik,viw
  real,dimension(3) :: saddlecolor,bowlcolor
  !!real :: dotprod
  logical :: sendburied,sendgrad
  real,parameter :: eps=0.05  !! sync this
  real :: dik,djk,diw,djw,dkw,ddij,ddjk,ddik
  !
  saddlecolor = (/0.5,1.0,0.0/)
  bowlcolor = (/1.0,0.8,0.0/)
  sendburied = present(buried)  !! if an array is passed, use it to send back
                                !! contact surfaces
  !! if (sendburied) write(*,*) "Sendburied = TRUE"
  sendgrad = present(gradient)  !! if an gradient is passed, get shifts, send back
  istart=1
  if(present(start)) istart=start
  i = 0
  periodic = (lbox > 0.)
  !!====== WARNING!!!  periodic has been turned back on when using sendgrad. Will it work??? =====
  !!  Tue May 16 15:50:41 EDT 2006
  if (sendgrad) periodic = .false.
  if (ishow /= 0 .and.drawing) then
    write(*,'("Show surface around residue ",i4," only.")') ishow
  endif
  
  sasa = 0.0      !! total solvent accessible molecular surface area
  ssasa = 0.0     !! total "saddle" convex/concave surface area
  bsasa = 0.0     !! total "bowl" concave/concave surface area
  sasanrg = 0.0   !! total solvent accessible molecular surface area
  ssasanrg= 0.0   !! total "saddle" convex/concave surface area
  bsasanrg= 0.0   !! total "bowl" concave/concave surface area
  nwat = 0        !! Number of verteces for re-entrant surface calc.
  if (sendgrad) then
    ! if (allocated(gradient)) then
    !! should be initialized in calling program
    !  gradient = 0.  !! initialize shift vectors to zero
    ! else
    !   write(0,'("MASKER::masker_getms: bug, gradient not allocated.")')
    !   stop 
    ! endif
  endif

 
  !! ====================================================================
  !! CONTACT SURFACE  (VDW,  SAS)  
  wsphere=rw*rw*4*pi
  !! Loop over all atoms, get sasa, saddle and bowl surfaces
  if (sendburied) buried = 0.
  do iat=istart,nat
    if (atype(iat)<=0) cycle   !! atoms that are flagged with negative numbers or unknown type
                               !! are simply ignored.
    if (atype(iat)>nattype) cycle
    if (atomlib(atype(iat))%name(1:1)=="H") cycle   !! ignore hydrogens
    if (atomlib(atype(iat))%name(1:1)=="D") cycle   !! ignore dummy atoms
    if (xyz(1,iat)>998.0) cycle  !! coordinates not assigned
    nm = 0
    amask = -1                 ! init all bits = 1
    bmask = -1
    r1 = atomlib(atype(iat))%r
    do jat=1,nat
      if (atype(jat)<=0) cycle
      if (atype(jat)>nattype) cycle
      if (atomlib(atype(jat))%name(1:1)=="H") cycle   !! ignore hydrogens
      if (atomlib(atype(jat))%name(1:1)=="D") cycle   !! ignore dummy atoms
      if (xyz(1,jat)>998.0) cycle
      if (iat == jat) cycle
      !! For every atom pair, calculate phi, psi, theta
      !! with iat as the central atom
      if (periodic) then ! get nearest copy of jat, return vector
        call boxaround(xyz(1:3,iat),xyz(1:3,jat),lbox,avec)
      else
        avec = xyz(1:3,jat) - xyz(1:3,iat)
      endif
      r2 = atomlib(atype(jat))%r
      maxd = r1 + r2 + 2*rw
      if (farapart(avec,maxd,dd)) cycle
      !! dd = sqrt(dotprod(avec,avec))
      !! if (dd > maxd) cycle    !! atom out of range
      if (dd > 0.) then
        !! write(*,*) '>>>>>>>>>>>>>>> dd=',dd,avec(1:3)
        avec = avec/dd
        !! returns imask and theta. If theta > pi/2, use 2*NTHETA-itheta and .not. it
        call getimask2(r1+rw,dd,r2+rw,theta,avec,imask,itheta)
      else
        write(*,*) '############## dd=',dd,avec(1:3)
        write(*,'("MASKER ERROR: iat,jat=",2i5," xyz=",6f8.2)') iat,jat,xyz(1:3,iat),xyz(1:3,jat)
        stop 'masker_getms - ERROR: bad coordinates'  !! bad coordinates, two atoms on top of each other. Ignore the second one.
      endif
      if (imask==0) cycle   !! bad triangle (bug?)
      if (imask==-1) then   !! i is completely embedded in j
        amask = 0
        exit          !! nevermind the other atoms j
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
          !!invert mask to get theta > pi/2, embedded mask.
          amask = iand(not(masklib(:,itheta,imask)),amask)     
        !endif
      endif
      !! if (itheta <= 0) cycle    !! atom out of range (borderline)
                                   !! borderline atoms should have gradient ~ #of points in itheta=1 mask.
      !! Diagnostic
      !!   if (iat==1.and.jat==2) write(*,*) 'Theta: ',theta/rad
      !! Diagnostic
      ! write(*,'(2i4,4f6.1,4i5)') iat,jat,psi/rad,phi/rad,theta/rad,dd,ipsi,iphi,itheta,imask
      !! if (itheta > NTHETA-1) then
      !!   write(*,*) "ERROR: itheta=",itheta," NTHETA=",NTHETA
      !!   write(*,*) " theta=",theta/rad," iat,jat,d=",iat,jat,dd
      !!   stop 'Theta out of range'  !! this would spot a bug
      !! endif
      !! remove the bits covered by atom jat
      !! Diagnostic / experimental: keep all close atoms, don't use bmask 21-NOV-01
      !!
      !! Using masks to determine whether a water exists has problems.
      !! It is found that numerical errors occur when waters (tetrahedron verteces)
      !! are close to each other. Instead of numerical methods, exact calculations
      !! are now used to find which of these waters to keep, which to exclude
      !! Mon Nov 26 11:52:10 EST 2001
      !!
      if (any(bmask /= amask) ) then
        !! if any changes occurred in the mask...
        bmask = amask
      endif      !! moved forward in BUG test. and kept.
      !! BUG test: does checking for changes exclude neighbors that
      !! are needed for reentrant surface ?? 
      !! Apparently yes. Keep all close atoms.
      !!-----------
        !! keep a list of close atoms, for the next part
        !! sort the list by theta, descending. Lowest theta first
        do k=nm,1,-1
          if (savetheta(k)<theta) exit
          saveiat(k+1) = saveiat(k)
          saveimask(k+1) = saveimask(k)
          saveitheta(k+1) = saveitheta(k)
          savetheta(k+1) = savetheta(k)
          savedij(k+1) = savedij(k)
        enddo
        k = k + 1
        saveiat(k) = jat
        saveimask(k) = imask
        saveitheta(k) = itheta   !! may be theta or (pi - theta) index.
        savetheta(k) = theta     !! depending on whether theta > pi/2
        savedij(k) = dd
        nm = nm + 1
        if (nm >= NSAVE-1) stop 'masker_getms - ERROR: NSAVE must be increased. Is the probe radius too big?'
      !! endif  !! endif moved backward in bug test. left there.
    enddo
    j = countbits(amask)
    if (j < 1) cycle   ! ignore saddle, etc if nothing is exposed.
    ! enddo
    !
    !! CONTACT SURFACE
    !! This should be modified to add an edge correction for each exposed edge. (done)
    sphere = 4*pi*atomlib(atype(iat))%r**2
    x = (real(j)/real(MAXATOM))
    if (sendburied) buried(iat) = x  !! 'buried' is really "exposed"
    x = sphere*x
    sasa = sasa + x    !! x is the sasa before edge correction
    sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
    if (drawing) then
      cmask = iand(amask,vmask)
      !! Diagnostic
      !! write(*,*) 'Drawing contact surface for atom ',iat
      ires = getresidue(iat)
      call drawsurface(13,xyz(1:3,iat),r1,iat,cmask,50.,"V",scale,ires)  !! output the sasa surface
    endif
    if (sendburied.and..not.sendgrad) cycle
    if (SAS) cycle 
    !
    !! ====================================================================
    !! TOROIDAL SURFACE  (SADDLE)
    !
    i = 0
    cmask = 0
    emask = 0
    dmask = amask     !! dmask is the original amask with all of the edges.
                      !! edges are removed from amask as we go along.
    nn = 0
    do k=1,nm
      theta = savetheta(k)   !! sorted low to high theta
      kat = saveiat(k)
      imask = saveimask(k)
      itheta = saveitheta(k)
      dd = savedij(k)
      r2 = atomlib(atype(kat))%r
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
          emask = iand(masklib(:,itheta-1,imask),dmask)
          !! remove the atomedgemask from the atommask to avoid overcounting edges.
          !! Thu Jun  6 06:42:29 EDT 2002
          amask = iand(amask,not(bmask))
        else
          bmask = iand(not(masklib(:,itheta+1,imask)),amask)
          emask = iand(not(masklib(:,itheta+1,imask)),dmask)
          !! remove the atomedgemask from the atommask to avoid overcounting edges.
          !! Thu Jun  6 06:42:29 EDT 2002
          amask = iand(amask,not(bmask))
        endif
      endif
      if (drawing) then
        cmask = ior(cmask,bmask)  !! collect all edges
        call drawsurface(13,xyz(1:3,iat),r1,kat,bmask,50.,"E",scale)  !! output the edges 1 at a time
      endif
      !! Count the bits in the exposed arc
      j = countbits(bmask)
      jj = countbits(emask)
      jarc = nint(real(j+jj)/2.)
      if (jarc < 1) cycle  !! no saddle, no bowl. Ignore this atom from here on
      !! write(0,'("Edgemask contains ",i9)') jarc
      arcfrac = real(jarc)/real(collarsize(itheta,imask))
      !! 
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
        !! saddle correction term jun-6 = (1.+0.04*exp(1.-arcfrac))
        ! col = 2*pi*arcfrac*(1.+0.02*exp(1.-arcfrac))
        !! 
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
        sdlarea = x
        !! write(0,'("Saddlearea ",2(i4,1x,a4),2f9.3)') iat,atomlib(atype(iat))%name,kat,atomlib(atype(kat))%name,sdlarea,arcfrac
        ssasa = ssasa + x
        ssasanrg = ssasanrg + x*atomlib(atype(iat))%w2
        if (drawing) then
          th = (pi/2) - theta  !! taumax
          !! x = itheta*dpsi  !! estimated theta
          !! output the saddle surface
          call drawsaddle(13,xyz(1:3,iat),xyz(1:3,kat),bmask,r1,iat,taumin,th,50.,"S",scale)  
          if (rendering) then
            !! NOTE! this will fail if theta is close to 90 degrees (it never is)
            !! dmask = iand(not(masklib(:,itheta+3,imask)),amask)   !! dont reassign dmask
            call rendersaddle(xyz(1:3,iat),xyz(1:3,kat),r1,amask,dmask,th,taumin,saddlecolor,iat,kat)
          endif
        endif
        !! if calculating derivs
        if (sendgrad) then
          !!======= toroidal surface derivative =============
          !! this is numerical and should eventually be replaced by
          !! a semi- analytical solution.
          r2 = atomlib(atype(kat))%r
          th = (itheta+1)*dpsi  !! theta for next shell
          !! th = theta + delt  !! small change in theta
          !! if (th > pi/2) th=pi/2
          y = (r1+rw)*sin(th)
          if (y > r2+rw) then
            c2 = 0.
          else
            c2 = sqrt((r2+rw)**2 - (y)**2)/y   !! This will fail if r1 >> r2
          endif
          !! fixed sign error  in next line 9-jan-03
          dij = -((r1+rw)*cos(th) + (r2+rw)*c2)   !! distance assoc w dpsi shift
          dsasa = -(real(jarc)/real(MAXATOM))*sphere*atomlib(atype(iat))%w1  !! kJ/mol for dpsi shift
                     !!  jarc = points in edgemask
          !! msderiv(iat,kat) = msderiv(iat,kat) + dsasa/dij
          xx = dsasa/(dij*2)       !!   divided by two because 1/2 is applied to each atom.
          !write(0,'("Contact shift: ",2i4,f9.4)') iat,kat,xx
          call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,kat),xx)
          call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,iat),xx)
                          !! kJ/mol.A, negative dij is decrease in distance
          !! get saddle for dij-shift, same arc length
          xx = col*((0.5*pi-th)*y - rw*cos(th))
          taumin = 0.
          if (y < rw) then
            taumin = acos(y/rw)
            xx = xx + col*(rw*sin(taumin) - taumin*y)
          endif
          !! xx is the half-saddle shifted by dij, so (xx - x) is the change in A^2,
          !! so (xx - x)*atomlib(atype(iat))%w2/dij is the energy deriv.
          !! msderiv(iat,kat) = msderiv(iat,kat) + (xx - x)*atomlib(atype(iat))%w2/dij
          xx = (xx-sdlarea)*atomlib(atype(iat))%w2/(dij*2)    !! kJ/mol/A
          ! write(0,'("Toroidal shift: ",2i4,f9.4)') iat,kat,xx
          !! write(0,'("Toroidal derivs: ",7f9.2)') xx,sdlarea,dij
          call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,kat),xx)
          call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,iat),xx)
        endif
      else    !! theta > pi/2
        !! write(*,*) 'No Saddle area: ',iat,kat,col," theta=",theta/rad
        x = (interc(itheta) - pi + theta)*slope(itheta)*arcfrac*sphere
        sasa = sasa + x
        sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
        ! write(*,'(8x,i4," embedded correction=",f8.2)') kat,x
        ! if (abs(x) > 10.) then
        !   write(*,*) "theta ",theta," itheta ",itheta
        !   write(*,*) "interc ",interc(itheta), " slope ",slope(itheta)
        !   write(*,*) "arcfraC ",arcfrac, " sphere ",sphere
        ! endif
        !!----------------------------------------------------------------------
        !! NOTE: currently there is no gradient calculated for embedded surfaces.
        !!----------------------------------------------------------------------
      endif   !! end if saddle exists
      !! end if calculating derivs
      !! create re-reduced atom set from atoms that have saddle
      !! if (any(bmask /= 0)) then
      !  nn = nn + 1    !!  note: nn <= k, so this is legal.
      !  saveiat(nn) = kat
      !  saveimask(nn) = imask
      !  saveitheta(nn) = itheta
      !  savetheta(nn) = theta
      !  savedij(nn) = dd
      !! endif
    enddo
    nn = nm
    !! ====================================================================
    !! REENTRANT SURFACE (bowl)
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
      if (kat <= iat) cycle  !! calculate bowl if iat is the first atom
      imask1 = saveimask(k)
      itheta1 = saveitheta(k)
      ! theta1 = savetheta(k)
      r2 = atomlib(atype(kat))%r
      do m=1,nn
        jat = saveiat(m)
        if (jat <= kat) cycle  !! calculate bowl if kat is the second atom
        imask2 = saveimask(m)
        itheta2 = saveitheta(m)
        !! --------- Ignore collars. Just get the third side of the triangle
        !! Wed Nov 21 08:15:00 EST 2001
        ! revised for speed, jun-7
        ! if (periodic) then ! get nearest copy of jat, return vector
        !   dd = distmod(xyz(1:3,kat),xyz(1:3,jat),lbox)
        ! else
        !   avec = xyz(1:3,jat) - xyz(1:3,kat)
        !   dd = sqrt(dotprod(avec,avec))
        ! endif
        if (periodic) then ! get nearest copy of jat, return vector
          call boxaround(xyz(1:3,kat),xyz(1:3,jat),lbox,avec)
        else
          avec = xyz(1:3,jat) - xyz(1:3,kat)
        endif
        r3 = atomlib(atype(jat))%r
        maxd = r3 + r2 + 2*rw
        if (farapart(avec,maxd,dd)) cycle
        !! if (dd > maxd) cycle    !! third side too long 
        !! ----------
        ! theta2 = savetheta(m)
        !   bmask = iand(not(masklib(:,itheta1+kcc,imask1)),masklib(:,itheta1,imask1) )
        !   bmask = iand(bmask,not(masklib(:,itheta2+kcc,imask2)) )
        !   bmask = iand(bmask,masklib(:,itheta2,imask2) )
        !   bmask = iand(cmask,bmask)  !! bmask is current verteces
        !!
        !!  Find/Save waters for RE_ENTRANT SURFACE (BOWL)
        !!
        !if (.true.) then
        !! if (any(bmask /= 0)) then !! if there is at least one vertex
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
        !  !! if (any(emask /= 0)) then
        !  if (.true.) then
        !    !! call waterbowl(iat,kat,jat,x,wxyz,xyz,nat,imask) !! also uses/returns emask, uses cvec
        call tetrahedron(veci,vecj,veck,wxyz,d,iat,jat,kat,1)
        if (d==0.0) cycle   !! fails triangle inequality
        if (exposed(wxyz,xyz,nat)) then
          nwat = nwat + 1
          if (nwat > watfac*nat) then
            write(*,*) 'Ran out of water space: iat=',iat,' limit=',watfac*nat
            stop 'masker_getms - ERROR: not enough water(1). reset watfac'
          endif
          water(nwat)%xyz = wxyz; water(nwat)%trng(1:3) = (/iat,jat,kat/)
        endif
        !    write(*,'("R ",4i4)') nwat,water(nwat)%trng(1:3)
        !    do i=1,3
        !      avec = water(nwat)%xyz - xyz(1:3,water(nwat)%trng(i))
        !      d = sqrt(dotprod(avec,avec))
        !      write(*,'(2f8.3,$)') d, atomlib(atype(water(nwat)%trng(i)))%r
        !    enddo
        !    write(*,*)
        !    ! call waterbowl(veci,vecj,veck,x,wxyz,imask,iat,jat,kat)  !! R-handed water (use cvec for hole)
        !    !! , write water coords
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
        !  !! : forget vecteces. just place all waters
        !  if (.true.) then
        !  ! if (any(emask /= 0)) then
        call tetrahedron(veci,veck,vecj,wxyz,d,iat,kat,jat,1)
        if (d==0.0) cycle   !! fails triangle inequality
        if (exposed(wxyz,xyz,nat)) then
          nwat = nwat + 1
          if (nwat > watfac*nat) then
            write(*,*) 'Ran out of water space: iat=',iat,' limit=',watfac*nat
            stop 'masker_getms - ERROR: not enough water(2). reset watfac'
          endif
          water(nwat)%xyz = wxyz; water(nwat)%trng = (/iat,kat,jat/)
        endif
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
        !    !   !! , write water coords
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
  enddo  !! iat=
  ms = sasanrg
  if (sendburied.and..not.sendgrad) return
  if (SAS) return 

  !! RE-ENTRANT SURFACE accounting for intersecting re-entrant surfaces.
  !! 
  !! This segment of the code (obselete given correction of bugs in tetrahedron?)
  !! was written to ensure that no waters sat too close to any atom.
  !! === This code has since been commented out. ===
  !! 
  sphere = 4*pi*rw**2
  ! xx = 4*rw*rw
  xx = 2*rw
  !! write(*,*) "nwat = ",nwat
  water(1:nwat)%good = .true.
  !!! obselete?
  !! Including this loop does not fix the bug of missing waters.
  !do iwat = 1,nwat
  !  do iat=1,nat
  !    d = distmod(water(iwat)%xyz,xyz(1:3,iat),lbox)
  !    x = (rw + atomlib(atype(iat))%r)
  !    if (d + 0.001 < x) then
  !      water(iwat)%good = .false.
  !      write(*,'(a,4i4,a,i4,a,f6.2,a,f6.2)') "Water ",iwat,water(iwat)%trng(1:3)," too close to atom ",iat, "  d= ",d, " < ",x
  !      exit
  !    endif
  !  enddo
  !enddo
  !!
  do iwat = 1,nwat
    ! if (.not.water(iwat)%good) cycle
    !! show water triangles. For debugging.
    !  write(*,*) iwat,water(iwat)%trng(1:3)
    !! write(*,'(4(3f8.3,1x))') water(iwat)%xyz(1:3),xyz(1:3,water(iwat)%trng(1)), &
    !!      xyz(1:3,water(iwat)%trng(2)),xyz(1:3,water(iwat)%trng(3))
    !! 
    !! trainglemask: get the mask  (amask) for the reentrant (bowl) surface
    !! at the position iwat. 
    !!
    iat = water(iwat)%trng(1)
    jat = water(iwat)%trng(2)
    kat = water(iwat)%trng(3)
    if(.not.(iat.ge.istart.or.jat.ge.istart.or.kat.ge.istart))cycle 
    ! write(0,'("WATER: ",4i5)') iwat,iat,jat,kat
    call trianglemask(water(iwat)%xyz(1:3),xyz(1:3,iat),xyz(1:3,jat),xyz(1:3,kat),amask) 
         !! this routine initializes amask
    !! count bits in complete triangle (before removing intersections)
    jtri = countbits(amask)
    !!
    !! Remove intersecting reentrant surfaces. If any part of any bowl
    !! surface is closer to a different water, such as the left-hand mirror
    !! water on the same 3-atom base, then remove those mask points.
    !!
    nn = 0; jcap=0
    kwat = 0
    do jwat = 1,nwat
      ! if (.not.water(jwat)%good) cycle
      if (iwat == jwat) cycle
      !! ------ look for waters with at least two atoms in common ------
      j = 0
      do iat=1,3
        do jat=1,3
          if (water(iwat)%trng(iat)==water(jwat)%trng(jat)) j = j + 1
        enddo
      enddo
      if (j < 2 ) cycle   !! ----- comment out for testing ------
      if (periodic) then  !! NOTE: periodic should be .false. when sendgrad in .true.
        call boxaround(water(iwat)%xyz,water(jwat)%xyz,lbox,avec) 
      else
        avec = water(jwat)%xyz - water(iwat)%xyz
      endif
      xx = 2*rw
      if (farapart(avec,xx,d)) cycle
      ! d = dotprod(avec,avec)
      ! if (d < xx) then    !! xx = 4*rw*rw  = (twice the probe radius)^2
      !   if (d<=0.0) then
      !     write(*,'("Warning: two waters EXACTLY on top of each other.")')
      !     cycle  !! this would only happen (d=0) VERY RARELY
      !   endif
      !   d = sqrt(d)
      if (d > 0.) then
        avec = avec/d
        imask = getimask(avec)
        ! write(*,*) '>>>>>>>>>>>>>>> d=',d,avec(1:3),imask
      else
        write(*,*) '>>>>>>>>>>>>>>> water-water distance d=',d,' vec=',avec(1:3)
        cycle
      endif
      call cosinerule(rw,d,rw,theta,flag); if (flag /= 0) cycle !! bug??
      itheta = nint((theta/rad)/DTHETA)
      if (itheta > NTHETA) stop 'masker_getms - BUG: Two waters. check the calculations'
      !! write(*,*) 'itheta=',itheta,' imask=',imask
      if (itheta>0) then
        amask = iand(masklib(:,itheta,imask),amask)     
        nn = nn + 1
        kwat = jwat
        jcap =  countbits(not(masklib(:,itheta,imask)))
      endif
      !!----------------------------------------
      !! At this point we have two triangles with two atoms in common.
      !! A toroid exists between the two triangles. We can calculate and 
      !! draw it here, instead of above. This would account for
      !! partial, but not full, saddles.
      !! Mon Jan 19 23:00:47 EST 2004
      !!----------------------------------------
      !!
      !! monitor intersecting reentrants
      !! if (any(amask /= bmask)) then
      !!   write(*,'("Intersection: iwat=",i7," jwat=",i7)') iwat,jwat
      !!  endif
      ! endif
    enddo
    jxtri = countbits(amask)
    ! write(0,'("jtri, jxtri, jcap= ",3i8," nn=",i4," fract=",f8.2)') jtri, jxtri, jcap, nn,real(jxtri)/real(jtri)
    ddij=0.; ddjk=0.; ddik=0.; abur=0.
    !! get triangle derivative
    if (sendgrad) then
      !!
      iat = water(iwat)%trng(1)
      jat = water(iwat)%trng(2)
      kat = water(iwat)%trng(3)
      abur = real(jxtri)/real(jtri)
      ! dvec = xyz(1:3,jat) - xyz(1:3,iat)  !! ij
      ! dij = sqrt(dotprod(dvec,dvec))
      ! dvec = xyz(1:3,kat) - xyz(1:3,iat)  !! ik
      ! dik = sqrt(dotprod(dvec,dvec))
      ! dvec = xyz(1:3,kat) - xyz(1:3,jat)  !! jk
      ! djk = sqrt(dotprod(dvec,dvec))
      ! dvec = water(iwat)%xyz(1:3) - xyz(1:3,iat)  !! iw
      ! diw = sqrt(dotprod(dvec,dvec))
      ! dvec = water(iwat)%xyz(1:3) - xyz(1:3,jat)  !! jw
      ! djw = sqrt(dotprod(dvec,dvec))
      ! dvec = water(iwat)%xyz(1:3) - xyz(1:3,kat)  !! kw
      ! dkw = sqrt(dotprod(dvec,dvec))
      !! deriv of spherical triangle w/respect to each pair
      !! for example, ddij is the change in spher. tri. w/change in ij distance
      !! multiplied by w3, in kJ/A^2 we get kJ/A 
      !! call strianglederiv(rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik)  !! A^2/A
      call strianglederiv(rw,xyz(1:3,iat),xyz(1:3,jat),xyz(1:3,kat),water(iwat)%xyz(1:3),ddij,ddjk,ddik)  
      ddij = ddij
      ddjk = ddjk
      ddik = ddik
      !! write(0,'("============ striangle:",3i5,3f8.3,$)') iat,jat,kat,ddij,ddjk,ddik
      !! 
      !! w3 > 0 means increased area --> increased energy, 
      !! so, if dA/dDij > 0., then dDij should be < 0.
      !! if xx > 0, then dDij < 0., so 
      !! xx = dA/dDij times w3
      !! 
      xx = ddij*atomlib(atype(iat))%w3/2                    !! dTRI/di->j
      call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,jat),xx)
      xx = ddij*atomlib(atype(jat))%w3/2                    !! dTRI/dj->i
      call addgrad(gradient(1,jat),xyz(1,jat),xyz(1,iat),xx)
      xx = ddjk*atomlib(atype(jat))%w3/2                    !! dTRI/dj->k
      call addgrad(gradient(1,jat),xyz(1,jat),xyz(1,kat),xx)
      xx = ddjk*atomlib(atype(kat))%w3/2                    !! dTRI/dk->j
      call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,jat),xx)
      xx = ddik*atomlib(atype(iat))%w3/2                    !! dTRI/di->k
      call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,kat),xx)
      xx = ddik*atomlib(atype(kat))%w3/2                    !! dTRI/dk->i
      call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,iat),xx)
    endif
    ddij=0.; ddjk=0.; ddik=0.; abur=0.
    if (jxtri /= jtri) then  !! then there was an intersecting surface (5th body)
      if (sendgrad) then
        iat = water(iwat)%trng(1)
        jat = water(iwat)%trng(2)
        kat = water(iwat)%trng(3)
        !! area of buried portian of spherical cap
        ! abur = real(jtri - jxtri)/real(jtri)  !! fraction buried
        abur = real(jtri - jxtri)/real(jcap)  !! fraction of cap buried
        dvec = xyz(1:3,jat) - xyz(1:3,iat)  !! ij
        vij = dvec
        dij = sqrt(dotprod(dvec,dvec))
        dvec = xyz(1:3,kat) - xyz(1:3,iat)  !! ik
        vik = dvec
        dik = sqrt(dotprod(dvec,dvec))
        dvec = xyz(1:3,kat) - xyz(1:3,jat)  !! jk
        djk = sqrt(dotprod(dvec,dvec))
        dvec = water(iwat)%xyz(1:3) - xyz(1:3,iat)  !! iw
        viw = dvec
        diw = sqrt(dotprod(dvec,dvec))
        dvec = water(iwat)%xyz(1:3) - xyz(1:3,jat)  !! jw
        djw = sqrt(dotprod(dvec,dvec))
        dvec = water(iwat)%xyz(1:3) - xyz(1:3,kat)  !! kw
        dkw = sqrt(dotprod(dvec,dvec))
        call cros(vij,vik,cij)
        cij = cij/sqrt(dotprod(cij,cij))
        d = (rw-abs(dotprod(viw,cij)))
        acap = 2*pi*rw*d
        d = d/rw
        ! abur = ((jtri - jxtri)/real(MAXATOM))*sphere
        ! abur = abur/acap
        ! write(0,*) "Fract Cap height=",d," fraction of cap buried=",abur, " Cap surface area",acap
        !! get the derivative of the excluded cap
        !! The distances dij,dik etc are already calculated. above.
        call dCapdDij(abur,rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik)
        !! 
        !!write(0,'("============ dCapdDij:",3i8,f10.2,3f10.4)') iat,jat,kat,xx,ddij,ddjk,ddik
        !! if w3 > 0, then increased area means increased energy, so 
        !! if dA/dDij > 0, then an increase in distance means increased energy,
        !! so the vector should move i toward j, xx should be > 0.
        !!  But dCapdDij() is the change in the amount of surface lost due to intersection.
        !! So dCapdDij > 0 means a negative contribution to the area. Therefore
        !! xx = -dA/dDij times w3.  We divide this vector by two because half
        !! 
        xx = -ddij*atomlib(atype(iat))%w3/2                    !! dCAP/di->j
        call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,jat),xx)
        xx = -ddij*atomlib(atype(jat))%w3/2                    !! dCAP/dj->i
        call addgrad(gradient(1,jat),xyz(1,jat),xyz(1,iat),xx)
        xx = -ddjk*atomlib(atype(jat))%w3/2                    !! dCAP/dj->k
        call addgrad(gradient(1,jat),xyz(1,jat),xyz(1,kat),xx)
        xx = -ddjk*atomlib(atype(kat))%w3/2                    !! dCAP/dk->j
        call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,jat),xx)
        xx = -ddik*atomlib(atype(iat))%w3/2                    !! dCAP/di->k
        call addgrad(gradient(1,iat),xyz(1,iat),xyz(1,kat),xx)
        xx = -ddik*atomlib(atype(kat))%w3/2                    !! dCAP/dk->i
        call addgrad(gradient(1,kat),xyz(1,kat),xyz(1,iat),xx)
      endif
    endif
    !! if (sendgrad) write(0,'(2i7,f8.5,3f8.3)') jtri,jxtri,abur,ddij,ddjk,ddik  !! finishes the line.
    if (jxtri < 1) cycle
    x = wsphere*real(jxtri)/real(MAXATOM)  !! area in A^2
    asc = wsphere*real(jtri-jxtri)/real(MAXATOM)   !! area of buried spherical cap.
    bsasa = bsasa + x
    !! writeout surface for each atoms
    ! write(*,'("bsasa = ",f8.2,$)') x
    bsasanrg = bsasanrg + (x/3.)*(atomlib(atype(water(iwat)%trng(1)))%w3 + &
        atomlib(atype(water(iwat)%trng(2)))%w3 + atomlib(atype(water(iwat)%trng(3)))%w3)
    if (drawing) then  !! do this to force drawing surface of bowls
      !! , draw water
      avec = water(iwat)%xyz(1:3)
      iat = water(iwat)%trng(1)
      jat = water(iwat)%trng(2)
      kat = water(iwat)%trng(3)
      if (periodic) call inbox(avec,lbox)
      write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1,i7)') &
      iwat,"W",iwat,scale*avec(1:3),0.,50.,nwat
      write(13,'("TER")') 
      dmask = iand(amask,vmask)
      !! 
      !! write(*,*) "Water number ",iwat
      ! call drawsurface(13,water(iwat)%xyz,rw,iwat,dmask, 50.,"B",1) !! output the surface
      call drawsurface(13,water(iwat)%xyz,rw,iwat,amask, 50.,"B",scale) !! output the surface
      call renderstriangle(rw,avec,xyz,nat,iat,jat,kat,amask,bowlcolor,intrs=(jxtri/=jtri),kwat=kwat)
    endif
    !! NOTE ON DERIVS: (old idea)
    !! At this point, 'amask' contains the reentrant surface and 'water(iwat)' contains
    !! the indeces of the three atoms. To get the bowl derivs, get the three edgemasks
    !! for the three sides of the trianglemask. Each edge is the amount of surface
    !! buried for moving one atom toward the other two. How do we get pairwise from this??
    !! get probe_edge_masks----------------------
    !! NEW IDEA:
    !!  The derivatives for the bowls can be calculated analytically, using the coordinates
    !!  of the tetrahedron, and the formula for the surface of a spherical
    !!  triangle (R^2(A+B+C-pi)) . If there is a hole in the bowl, then its size is
    !!  the spherical cap (surface = 2pi.R.h) possibly truncated (spherical frustoms omitted!)
    !!  The gradient of the cap is calculated, then multiplied by the fraction truncated using
    !!  mask numbers. No need to use probe_edge_masks.
    !!  Sat Dec  6 15:51:14 UTC 2003

    !!===== NOTE: move the if (sendgrad) section to inside the water loop. Call
    !! subroutine strianglederiv(rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik)
    !! and subroutine dCapdDij(rw,dij,dik,djk,diw,djw,dkw,ddij,ddjk,ddik) and multiply the latter
    !! results by the ratio of countbits(amask) before and after intersections.
    !!===== 31-dec-03

  enddo 

  !! ====================================================================
  !! REPORT the surfcae area breakdown, if desired (meaning if draw.and.saying)
  y = sasa + ssasa + bsasa
  x = sasanrg + ssasanrg + bsasanrg
  if (saying) then
     ! write(*,'("Total VDW area     = ",f9.2)') sasa
     ! write(*,'("Total saddle area  = ",f9.2)') ssasa
     ! write(*,'("Total bowl area    = ",f9.2)') bsasa
     ! write(*,'("Total surface area = ",f9.2)') y
     ! write(*,'("Total VDW nrg      = ",f9.2," kJ")') sasanrg
     ! write(*,'("Total saddle nrg   = ",f9.2," kJ")') ssasanrg
     ! write(*,'("Total bowl nrg     = ",f9.2," kJ")') bsasanrg
     ! write(*,'("Total SES energy   = ",f9.2," kJ")') x
  endif
  ms = x
  !! close(13)
  end subroutine masker_getms ! (xyz,atype,nat,ms)
  !!----------------------------------------------------------------!!
  subroutine addgrad(grad,at1,at2,x)
    implicit none
    real(kind=kind_8),dimension(3),intent(inout) :: grad
    real,dimension(3),intent(in) :: at1,at2
    real,intent(in) :: x
    real,dimension(3) :: v,u
    !! add a vector to grad in the direction at1-at2 with magnitude x
    v = at2 - at1
    call unitvec(v,u)
    v = u*x
    !! 
    !! write(0,'("addgrad: ",4f9.3,6f8.2)') x,v(1:3),at1(1:3),at2(1:3)
    grad = grad + v
  end subroutine addgrad
  !!----------------------------------------------------------------!!
  logical function farapart(avec,d,dd)
  implicit none
  real,dimension(3),intent(in) :: avec
  real,intent(in) :: d
  real,intent(inout) :: dd
  farapart = .true.
  if (abs(avec(1))>d) return
  if (abs(avec(2))>d) return
  if (abs(avec(3))>d) return
  dd = dotprod(avec,avec)
  if (dd > d*d) return
  dd = sqrt(dd)
  farapart = .false.
  end function farapart
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
    if (atype(iat)<=0) cycle
    if (atype(iat)>nattype) cycle
    if (atomlib(atype(iat))%name(1:1)=="H") cycle   !! ignore hydrogens
    if (atomlib(atype(iat))%name(1:1)=="D") cycle   !! ignore dummy atoms
    d = distmod(axyz,xyz(1:3,iat),lbox)
    x = (rw + atomlib(atype(iat))%r) - 0.01  !! replacing 0.001A buffer 4-JUN-02
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
  !! write(*,*) 'gettheta: ',r1,r2,dd,x
  end function gettheta
  !!----------------------------------------------------------------!!
  subroutine waterbowl(veci,vecj,veck,x,wxyz,imask,iat,jat,kat)
  !! >>>>>>>>> obselete ? <<<<<<<<<<<<
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
  ! write(*,*) 'emask not=0'
  !! wxyz is water position in absolute coords (A)
  !wxyz = (r1+rw)*jvec   !! current value of r1 set in calling routine
  !d = dotprod(cvec,wxyz)
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
    !! Diagnostic
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
    real, dimension(3) :: wivec,wjvec,wkvec,dvec,evec,fvec
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
    evec = dvec
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    amask = iand(masklib(:,NTHETA,jmask),amask)
  
    !! mask kiw plane =================
    call cros(wkvec,wivec,dvec)
    fvec = dvec
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    amask = iand(masklib(:,NTHETA,jmask),amask)
  
    !!==== check to see if the planes are nearly parallel
    y = netangle(dvec,evec)
    !! 
    !! write(*,*) "ij versus jk angle =",y/rad
    if (y > 2*dpsi) return
    y = netangle(dvec,fvec)
    !! 
    !! write(*,*) "ij versus ik angle =",y/rad
    if (y > 2*dpsi) return
    y = netangle(evec,fvec)
    !! 
    !! write(*,*) "ik versus jk angle =",y/rad
    if (y > 2*dpsi) return
    call masktrianglecircle(wivec,wjvec,wkvec,amask)
    !! return spherical triangle mask
  end subroutine trianglemask
  !!------------------------------------------------------------
  real function netangle(dvec,evec)
    real,dimension(3),intent(in) :: dvec,evec
    real,dimension(3) :: avec,bvec
    real :: x
    avec = dvec/sqrt(dotprod(dvec,dvec))
    bvec = evec/sqrt(dotprod(evec,evec))
    x = acos(dotprod(avec,bvec))
    if (x > pi/2) x = pi - x
    netangle = x
  end function netangle
  !!------------------------------------------------------------
  subroutine masktrianglecircle(wivec,wjvec,wkvec,amask)
    !! Thu Jan  1 15:39:27 EST 2004
    !! mask a circle in the direction of the average wx vector
    real,dimension(3),intent(in) :: wivec,wjvec,wkvec
    integer(kind=KND),dimension(MASKSIZE),intent(inout) :: amask
    real,dimension(3) :: avec,bvec
    real :: y
    integer :: jmask,itheta
    !!
    bvec = -(wivec + wjvec + wkvec )
    avec = bvec/sqrt(dotprod(bvec,bvec))
    y = max(netangle(avec,wivec),netangle(avec,wjvec),netangle(avec,wkvec))
    itheta = nint((y/rad)/DTHETA)
    jmask = getimask(avec)
    amask = iand(not(masklib(:,itheta,jmask)),amask)
    !!write(*,*) "In masktrianglecircle. maxtheta=",y," itheta=",itheta
  end subroutine masktrianglecircle
  !!------------------------------------------------------------
  subroutine expandedtrianglemask(wxyz,ixyz,jxyz,kxyz,amask,iv,jv,kv)
  !! The tetrahedron stored in awat is assumed to be right-handed.
  !!-----------------------
  !! Reentrant surface is expanded by 0.1A in all direcitions without moving the location
  !! of the water probe. This is done by getting (e.g.) the cross-product of the w->i
  !! vector and the i->j vector and moving by 0.1 (tangent to the water sphere and away from
  !! atom k perpendiculr to the ijw plane.) The three new atoms are used to get the
  !! bowl surface, and that mask is passed back. We don't care about the mask itself,
  !! just its area. In the calling function (masker_getms), the difference between this mask
  !! and the un-perturbed reentrant mask area is the derivative to be applied to the
  !! three atoms equally.
  implicit none
  real,dimension(3),intent(in) :: wxyz,ixyz,jxyz,kxyz
  real,dimension(3),intent(out) :: iv,jv,kv
  integer(kind=KND),dimension(MASKSIZE),intent(out) :: amask
  real, dimension(3) :: wivec,wjvec,wkvec,dvec
  real :: y
  real,parameter :: eps=0.05  !! sync this
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

  iv = ixyz
  jv = jxyz
  kv = kxyz
  !! move i and j perpendicular to ijw plane =================
  call cros(wivec,wjvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  iv = iv + eps*dvec
  jv = jv + eps*dvec
  !! move j and k perpendicular to jkw plane =================
  call cros(wjvec,wkvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  kv = kv + eps*dvec
  jv = jv + eps*dvec
  !! 
  !! write(0,'("Expandedtriangle: i j k=",9f9.4)') iv,jv,kv
  !! move k and i perpendicular to kiw plane =================
  call cros(wkvec,wivec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  kv = kv + eps*dvec
  iv = iv + eps*dvec

  !! Get new sides of tetrahedron
  if (periodic) then
    call boxaround(wxyz,iv,lbox,wivec) ! w->i
    call boxaround(wxyz,jv,lbox,wjvec) ! w->j
    call boxaround(wxyz,kv,lbox,wkvec) ! w->k
  else
    wivec = wxyz - iv
    wjvec = wxyz - jv
    wkvec = wxyz - kv
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
  !! 
  !! write(0,'("Expandedtriangle:",3f9.4)') dvec
  jmask = getimask(dvec)
  amask = iand(masklib(:,NTHETA,jmask),amask)

  !! mask kiw plane =================
  call cros(wkvec,wivec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  amask = iand(masklib(:,NTHETA,jmask),amask)

  !! return spherical triangle mask for expanded atoms
  end subroutine expandedtrianglemask
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
    stop 'masker.mod - Bug: psi calculation: getimask'
  else
    phi = atan2(uvec(2),uvec(1))
  endif
  !! Diagnostic
  !  write(*,'(a,3f6.2,2f6.1,2i5,f6.3)') 'getimask: uvec, psi, phi, ipsi, iphi, dphi',uvec, psi, phi, ipsi, iphi
  dphi = 2*pi/real(nphi(ipsi))
  ! if (dphi == 0.0) stop 'BUGBUGBUGBUGBUG'
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  getimask = psiposit(ipsi) + iphi
  end function getimask
  !!------------------------------------------------------------
  subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,scale,ires)
  !! rendering added using triangles
  !! Wed Aug 29 10:39:53 EST 2001
  implicit none
  integer,intent(in) :: iunit,iat  !  ,iskip
  real,dimension(3),intent(in) :: cen
  real,intent(in) :: r,bb,scale
  integer,intent(in),optional :: ires
  real :: d
  character(len=1),intent(in) :: ch
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer :: i,j,ibyte,ibit,irs
  real,dimension(3) :: avec,bvec,rgb
  real,parameter :: DCUT=6.0
  !
  j = 0
  irs = 0
  if (present(ires)) irs = ires
  do i=1,MAXATOM
    ibyte = (i-1)/NBIT + 1
    ibit = mod(i-1,NBIT)
    if (btest(amask(ibyte),ibit)) then
      ! j = j + 1
      ! if (mod(j,iskip)==0) then
        avec = cen + r*mxyz(1:3,i)
        if (ishow /= 0) then
          bvec = avec - showvec
          if (sqrt(dotprod(bvec,bvec)) > DCUT) cycle
        endif
        if (periodic) call inbox(avec,lbox)
        !! write(iunit,'("MODEL ",i5)') iat
        write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        i,ch,iat,scale*avec(1:3),0.,bb
        !! write(iunit,'("ENDMDL")') 
      ! endif
    endif
  enddo
  if (rendering) then
    if (ch == "B") then      !! bowl
      rgb = (/0.5,0.5,0.0/)  !! yellow?
      !! call rendersphere(cen,r,amask,rgb)  !! replaced by renderstriangle()
    elseif (ch == "V") then  !! SASA
      rgb = (/0.0,1.0,0.0/)  !! green
      call rendersphere(cen,r,amask,rgb,iat,irs)  !! no longer uses amask
    elseif (ch == "E") then  !! edges
      !! dont render edges
      rgb = (/0.0,0.0,1.0/)  !! blue
      return
    else
      !! dont render 
      rgb = (/1.0,1.0,1.0/)  !! white
      return
    endif
  endif
  end subroutine drawsurface
  !!==========================================================
  subroutine rendersphere(cen,r,amask,rgb,iat,ires)
  implicit none
  real,dimension(3),intent(in) :: cen
  real,dimension(3),intent(in) :: rgb
  real,intent(in) :: r
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer,intent(in) :: iat,ires
  integer :: i,j,ibyte,ibit
  real,dimension(3) :: avec,bvec,cvec,rgbtype
  real,dimension(3) :: anor,bnor,cnor
  logical :: cheating=.true.  !! use raster3D sphere object instead
  !! Diagnostic
  ! write(*,*) "RENDERSPHERE",cen

  rgbtype = rgb
  !! easy way, just draw the whole sphere as a smooth ball
  if (cheating) then
    !! color based on atom type
    call surfacecolor(iat,ires,shape=1,rgb=rgbtype)
    write(runit,*) 2
    write(runit,'(9f8.3)') cen(1:3),r,rgbtype
    return
  endif
  !!
  !! draw all triangles that have all bits set
  !! Changed: draw a triangle if any bit is set
  do j=1,ntri
    !if ((btest(amask((tri(1,j)-1)/NBIT + 1),mod(tri(1,j)-1,NBIT))).or.  &
    !    (btest(amask((tri(2,j)-1)/NBIT + 1),mod(tri(2,j)-1,NBIT))).or.  &
    !    (btest(amask((tri(3,j)-1)/NBIT + 1),mod(tri(3,j)-1,NBIT)))) then
      i = tri(1,j)
      !!  if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
      avec = cen + r*mxyz(1:3,i)
      anor = mxyz(1:3,i)
      i = tri(2,j)
      !!  if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
      bvec = cen + r*mxyz(1:3,i)
      bnor = mxyz(1:3,i)
      i = tri(3,j)
      !!  if (.not.(btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT)))) cycle
      cvec = cen + r*mxyz(1:3,i)
      cnor = mxyz(1:3,i)
      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
      write(runit,*) 1
      write(runit,'(12f8.3)') avec,bvec,cvec,rgb
      if (smoothing) then
        write(runit,*) 7
        write(runit,'(9f8.3)') anor,bnor,cnor
      endif
    !endif
  enddo
  end subroutine rendersphere
  !!==========================================================
  subroutine renderstriangle(r,wat,xyz,nat,iat,jat,kat,amask,rgb,intrs,kwat)
  implicit none
  integer,intent(in) :: nat,iat,jat,kat
  real,dimension(3,nat),intent(in) :: xyz
  real,dimension(3),intent(in) :: rgb,wat
  real,intent(in) :: r
  logical,intent(in) :: intrs
  !! integer :: iat,jat,kat
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask
  integer,intent(in) :: kwat
  integer :: i,j,k,nvert,ivert,jvert,kvert,ip,ipx,itau,flg,nout,im,np
  real,dimension(3) :: vij,vjk,vik,vwi,vwj,vwk,vv,vp,vio,cen,rgb2,rgb3,vijk
  real,dimension(3) :: wij,wjk,wik,vstart,vend,vi,vj,vk
  real,dimension(3) :: avec,bvec,cvec,vec,kv1,kv2,kv3
  real,dimension(3) :: anor,bnor,cnor,rgbi,rgbj,rgbk
  real,dimension(3,3) :: vert,norm,matij,matik,matjk,rgbv
  real,dimension(3,3,2) :: matxcng
  real :: dtauj, dtauk, dtaui,wsphere
  logical :: nomask=.true.
  integer :: nsect=10
  wsphere=r*r*4*pi
  !!
  rgbv(1:3,1) = rgb; rgbv(1:3,2) = rgb;rgbv(1:3,3) = rgb;
  vi = xyz(1:3,iat)
  vj = xyz(1:3,jat)
  vk = xyz(1:3,kat)
  call surfacecolor(iat,0,shape=3,rgb=rgbi)
  call surfacecolor(jat,0,shape=3,rgb=rgbj)
  call surfacecolor(kat,0,shape=3,rgb=rgbk)
  if (intrs) then
    nsect=50
    if (kwat/=0) then
      kv1 = xyz(1:3,water(kwat)%trng(1))
      kv2 = xyz(1:3,water(kwat)%trng(2))
      kv3 = xyz(1:3,water(kwat)%trng(3))
    else
      kv1 = vi
      kv2 = vj
      kv3 = vk
    endif
  else
    nsect=10
  endif
  !! get water normals
  rgb2 = (/0.,.8,0.8/)
  rgb3 = (/0.4,.2,0.9/)
  cen = wat
  vwi = vi - wat
  vwj = vj - wat
  vwk = vk - wat
  vij = vj - vi
  vik = vk - vi
  !! normals to the tetrahedron planes, pointing outwards
  call cros(vwi,vwj,wij)  ;  wij = wij/sqrt(dotprod(wij,wij))
  call cros(vwj,vwk,wjk)  ;  wjk = wjk/sqrt(dotprod(wjk,wjk))
  call cros(vwk,vwi,wik)  ;  wik = wik/sqrt(dotprod(wik,wik))
  call cros(vij,vik,vijk)  ;  vijk = vijk/sqrt(dotprod(vijk,vijk))
  if (intrs.and.kwat/=0) then
    vv = kv2 - kv1
    vp = kv3 - kv1
    call cros(vv,vp,vijk)  ;  vijk = vijk/sqrt(dotprod(vijk,vijk))
  endif

  if (nomask) then  !! don't use edgemask
    nout = 0
    np = 0
    vwi = vwi/sqrt(dotprod(vwi,vwi))   !! unit vector from w to i
    vwj = vwj/sqrt(dotprod(vwj,vwj))   !! unit vector from w to j
    vwk = vwk/sqrt(dotprod(vwk,vwk))   !! unit vector from w to k
    dtauj = acos(dotprod(vwi,vwj))/real(nsect)  !! angle increment for arc from vwi to vwj
    dtauk = acos(dotprod(vwi,vwk))/real(nsect)  !! angle increment for arc from vwi to vwk
    dtaui = acos(dotprod(vwj,vwk))/real(nsect)  !! angle increment for arc from vwj to vwk
    !! write(0,'("dtauj dtauk: ",2f9.4)') dtauj,dtauk
    vv = wat + wij
    call getrotS(wat,vv,dtauj,matij,vec)    !! incremental rotation around normal to wij plane, i to j
    vv = wat - wik
    call getrotS(wat,vv,dtauk,matik,vec)    !! incremental rotation around normal to wik plane, i to k
    vv = wat - wjk
    call getrotS(wat,vv,dtaui,matjk,vec)    !! incremental rotation around normal to wjk plane, j to k
    vstart = -vwi
    vend = -vwk
    i = 1
    im = 1
    matxcng = 0.
    do j=1,3; matxcng(j,j,:) = 1.; enddo
    matxcng(1:3,1:3,im) = matik(1:3,1:3)
    im = mod(im,2)+1
    do itau=1,nsect
      !! set up a row of triangles going from i to k
      norm(1:3,i) = vstart
      vert(1:3,i) = wat - r*norm(1:3,i)
      call interpolatecolor(rgbi,rgbj,rgbk,wij,wjk,wik,vstart,vwi,vwj,vwk,rgbv(1:3,i))
      !! move start of row
      call rotate_S(matij,vstart,vv)
      vstart = vv
      i = mod(i,3) + 1
      norm(1:3,i) = vstart
      vert(1:3,i) = wat - r*norm(1:3,i)
      call interpolatecolor(rgbi,rgbj,rgbk,wij,wjk,wik,vstart,vwi,vwj,vwk,rgbv(1:3,i))
      !! move end of row
      call rotate_S(matjk,vend,vv)
      vend = vv
      !! get angle across row
      if (itau<nsect) then
        dtauj = acos(dotprod(vstart,vend))/real(nsect-itau)
        call cros(vstart,vend,vv)  
        vv = wat + vv
        call getrotS(wat,vv,dtauj,matxcng(1,1,im),vec)
      endif
      nvert = 2*(nsect-itau)+1
      do ivert=1,nvert
        np = np + 1
        i = mod(i,3) + 1
        j = mod(i,3) + 1
        !! rotate point j to make point i (j exists already)
        im = mod(im,2)+1    !! switch matrix
        call rotate_S(matxcng(1,1,im),norm(1:3,j),vv)
        !! write(0,'("Rotating ",3f9.4," to ",3f9.4,f9.2," degrees.")') norm(1:3,j),vv,dtauk/rad
        norm(1:3,i) = vv
        vert(1:3,i) = wat - r*norm(1:3,i)
        call interpolatecolor(rgbi,rgbj,rgbk,wij,wjk,wik,vv,vwi,vwj,vwk,rgbv(1:3,i))
        if (intrs) then
          !! project excluded point back to ijk plane? someday?...
          !! avec = vert(1:3,1); bvec=vert(1:3,2); cvec=vert(1:3,3)
          !! anor = norm(1:3,1); bnor=norm(1:3,2); cnor=norm(1:3,3)
          flg = 0
          do k=1,3
            vv = vert(1:3,k)  - kv1
            if (dotprod(vv,vijk) < 0.) flg = flg + 1
            !!write(0,'(a10,i3,7f8.3)') "dotprod", k,vv,vijk,dotprod(vv,vijk)
          enddo
          if (flg>0) cycle
          !! call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgb)
          nout = nout + 1
          call renderonetriangle(vert(1,1),vert(1,2),vert(1,3),norm(1,1),norm(1,2),norm(1,3),rgbv)
        else
          nout = nout + 1
          call renderonetriangle(vert(1,1),vert(1,2),vert(1,3),norm(1,1),norm(1,2),norm(1,3),rgbv)
        endif
      enddo
    enddo
    if (nout==0) stop "masker.mod - BUG: entire reentrant surface is intersecting: renderstriangle"
    !! write(0,'("renderstrianlge: fract of triangle rendered=",f8.2)') (real(nout)/real(np))
    return   
  endif
  !! using mask
  do j=1,ntri
    nvert = 0
    if ((btest(amask((tri(1,j)-1)/NBIT + 1),mod(tri(1,j)-1,NBIT)))) nvert = nvert + 1
    if ((btest(amask((tri(2,j)-1)/NBIT + 1),mod(tri(2,j)-1,NBIT)))) nvert = nvert + 1
    if ((btest(amask((tri(3,j)-1)/NBIT + 1),mod(tri(3,j)-1,NBIT)))) nvert = nvert + 1
    i = tri(1,j)
    avec = cen + r*mxyz(1:3,i)
    anor = -mxyz(1:3,i)
    i = tri(2,j)
    bvec = cen + r*mxyz(1:3,i)
    bnor = -mxyz(1:3,i)
    i = tri(3,j)
    cvec = cen + r*mxyz(1:3,i)
    cnor = -mxyz(1:3,i)
    if (nvert==3) then
      call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgbv)
    elseif (nvert==2) then
      !! if (nvert==2) cycle
      do i=1,3
        if ((btest(amask((tri(i,j)-1)/NBIT + 1),mod(tri(i,j)-1,NBIT)))) exit
      enddo
      ivert = i
      jvert = mod(ivert,3) + 1
      avec = cen + r*mxyz(1:3,tri(jvert,j))
      anor = -mxyz(1:3,tri(jvert,j))
      kvert = mod(jvert,3) + 1
      bvec = cen + r*mxyz(1:3,tri(kvert,j))
      bnor = -mxyz(1:3,tri(kvert,j))
      !! how many planes does it cross?
      ip = 0; ipx=0
      if (dotprod(vij,mxyz(1:3,tri(ivert,j))) >= 0.) then; ip = ip + 1 ;ipx=ipx+1; endif
      if (dotprod(vjk,mxyz(1:3,tri(ivert,j))) >= 0.) then; ip = ip + 1 ;ipx=ipx+2; endif
      if (dotprod(vik,mxyz(1:3,tri(ivert,j))) >= 0.) then; ip = ip + 1 ;ipx=ipx+4; endif
      if (ip == 3) cycle  !! impossible unless a bug
      if (ip == 0) then  ! inside striangle, plot it. Might be intersecting...
        cvec = cen + r*mxyz(1:3,tri(ivert,j))
        cnor = -mxyz(1:3,tri(ivert,j))
        call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgbv)
      elseif (ip == 1) then  ! crosses one plane
        !! what plane does the triangle cross?
        k = 0
        if     (dotprod(vij,mxyz(1:3,tri(ivert,j))) >= 0.) then !! ivert is outside wij plane
          vp = vij
          k = 1
        elseif (dotprod(vjk,mxyz(1:3,tri(ivert,j))) >= 0.) then !! ivert is outside wjk plane
          vp = vjk
          k = 2
        elseif (dotprod(vik,mxyz(1:3,tri(ivert,j))) >= 0.) then !! ivert is outside wik plane
          vp = vik
          k = 3
        endif
        !! get new vertex for ivert-jvert
        call cros(mxyz(1:3,tri(jvert,j)),mxyz(1:3,tri(ivert,j)),vio)
        !  write(0,'("io:",2i5,5x,6f9.4)') jvert,ivert, vio(1:3),vp(1:3)
        call cros(vp,vio,vv)
        vv = vv/sqrt(dotprod(vv,vv))
        ! write(0,'("I:",2i5,3f9.4,5x,3f9.4)') k,ivert,mxyz(1:3,tri(ivert,j)), vv(1:3)
        cvec = cen + r*vv
        cnor = -vv
        call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgbv)
        !! get new vertex for ivert-kvert
        call cros(mxyz(1:3,tri(kvert,j)),mxyz(1:3,tri(ivert,j)),vio)
        call cros(vp,vio,vv)
        vv = vv/sqrt(dotprod(vv,vv))
        ! write(0,'("K:",2i5,3f9.4,5x,3f9.4)') k,kvert,mxyz(1:3,tri(kvert,j)), vv(1:3)
        avec = cen + r*vv
        anor = -vv
       ! call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgb2)
      elseif (ip==2) then   !! cross two planes
        !! just replace ivert with corner of spherical triangle.
        if (ipx==3) then
          vv = vwj/sqrt(dotprod(vwj,vwj))
        elseif (ipx==5) then
          vv = vwi/sqrt(dotprod(vwi,vwi))
        elseif (ipx==6) then
          vv = vwk/sqrt(dotprod(vwk,vwk))
        endif
        cvec = cen +r*vv
        cnor = -vv
        call renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgbv)
      endif
      !! if the missing vertex does not cross one of the planes, omit the triangle
    elseif (nvert==1) then
      !! fill on this missing code to get truncated triangles with only
      !! one vertex in the box
    else
    endif
  enddo
  end subroutine renderstriangle
  !!==========================================================
  subroutine interpolatecolor(rgbi,rgbj,rgbk,wij,wjk,wik,vv,vi,vj,vk,rgbv)
    !! chose an interpolated color according to the position of vv
    !! in the spherical triangle made by vwi,vwj,vwk
    real,dimension(3),intent(in) :: rgbi,rgbj,rgbk,wij,wjk,wik,vv,vi,vj,vk
    real,dimension(3),intent(out) :: rgbv
    real :: di,dj,dk,xi,xj,xk,dsum
    di = dotprod(vv,wjk)/dotprod(vi,wjk)
    dj = dotprod(vv,wik)/dotprod(vj,wik)
    dk = dotprod(vv,wij)/dotprod(vk,wij)
    dsum = (di+dj+dk)
    xi = di/dsum
    xj = dj/dsum
    xk = dk/dsum
    rgbv = xi*rgbi + xj*rgbj + xk*rgbk
  end subroutine interpolatecolor
  !!==========================================================
  !! NOTE: calls to this routine have been changed to pass a (3,3) real.
  !!  Wed Nov  9 11:20:33 EST 2005
  subroutine renderonetriangle(avec,bvec,cvec,anor,bnor,cnor,rgb)
    implicit none
    real,dimension(3),intent(inout) :: avec,bvec,cvec
    real,dimension(3),intent(in) :: anor,bnor,cnor
    real,dimension(3,3),intent(in) :: rgb
    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
    write(runit,*) 1
    write(runit,'(12f8.3)') avec,bvec,cvec,rgb(1:3,1)
    if (smoothing) then
      write(runit,*) 7
      write(runit,'(9f8.3)') anor,bnor,cnor
      write(runit,*) 17
      write(runit,'(9f7.3)') rgb(1:3,1),rgb(1:3,2),rgb(1:3,3)
    endif
  end subroutine renderonetriangle
  !!==========================================================
  !!    call rendersaddle(xyz(1:3,iat),xyz(1:3,kat),r1,amask,dmask,x,taumin)
  subroutine rendersaddle(va1,va2,r1,amask,bmask,thet,taumin,rgb,iatm,jatm)
  implicit none
  real,dimension(3),intent(in) :: va1,va2
  real,dimension(3),intent(inout) :: rgb
  real,intent(in) :: r1,thet,taumin
  real :: mtheta
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: amask,bmask  
  !! amask is atommask, bmask is atomedgemask
  integer,intent(in),optional :: iatm,jatm
  integer :: i,j,ibyte,ibit,k,itau,nn,ii,ntau,ires,jres
  integer,dimension(3) :: sides
  real,dimension(3,4) :: trv,trn,rgbv
  real,dimension(3) :: avec,bvec,cvec,vab,uab,vp
  !! real,dimension(3) :: anor,bnor,cnor
  !! real,dimension(3) :: v1,v2,v3,vt,vw,vec
  real,dimension(3) :: vv,vw,vv1,va,vw1,vwa,v3,vt,vec,rgbi,rgbj
  !! real,dimension(3) :: w1,w2,w3,wt,ww,wec,wp
  real,parameter :: spcng=0.25   !! angstrom spacing of saddle dots, reset to 0.35
  real :: chi,x,mat(3,3),dtau,y
  ntau = nint(rw*(thet - taumin)/spcng)   ! number of steps in tau
  !--- coloring
  ires = getresidue(iatm)
  call surfacecolor(iatm,ires,shape=2,rgb=rgbi)
  jres = getresidue(jatm)
  call surfacecolor(jatm,jres,shape=2,rgb=rgbj)
  !---
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
  !! rgb = (/0.5,0.0,0.5/)  !! magenta?
  !! draw strips for all triangles that have
  !! a side on the edge
  !!-------------------------------------------------------
  !! Instead of using edgemask and triangles, just draw the complete saddle.
  !! It is opaque anyway and the reentrant surface will cover it where
  !! it is not.
  !! --------- CB Tue Jan  6 09:36:02 EST 2004
  ! do j=1,ntri    !! loop over all triangles
  !   nn = 0
  !   do k=1,3     !! check to see if any of the vectices are in the edgemask
  !     i = tri(k,j)
  !     if (btest(bmask((i-1)/NBIT + 1),mod(i-1,NBIT))) then
  !       nn = nn + 1
  !     endif
  !   enddo
  !   if (nn == 0) cycle  !! no points in edgemask
  !   nn = 0
  !   do k=1,3     !! find triangles that have 2 points in the atommask
  !     i = tri(k,j)
  !     if (btest(amask((i-1)/NBIT + 1),mod(i-1,NBIT))) then
  !       nn = nn + 1
  !       sides(nn) = i
  !     endif
  !   enddo
  !   if (nn < 2) cycle   !! 0 or 1 vertex in the mask, forget it.
  !   if (nn > 2) cycle   !! if nn=3, whole triangle is already drawn.
  !   !! write(*,*) "found an edge-side",tri(1:3,j)
  !   !! a side of a triangle (adjacent points) were found in the edgemask
  !   !! draw a strip of surface starting from this side
  !   ii = sides(1)    ! first point
  !  vw = mxyz(1:3,ii) ! direction of mask point
  !!------------------------------------------------------------------------------
  !! new method, draw full saddle
  call cros(uab,mxyz(1:3,1),vw)    ! vw is perpendicular to vab
  vw = vw/sqrt(dotprod(vw,vw))   ! vw is unit vector
  vwa = (/0.,0.,0./)      ! origin
  call getrotS(vwa,vw,mtheta,mat,vec)  ! mat is mtheta rotation, vec is 0
  v3 = uab
  call move_S(v3,mat,vec)          ! v3 is unit vec in direction of start water
  vw1 = va1 + (r1+rw)*v3   ! position of water, to start
  vv1 = va1 + v3*r1               ! starting saddle point
  call cros(uab,v3,vwa)    ! 
  va = vw1 + vwa           ! vw1 and va1 are starting axis of saddle arc (tau)
  do itau=1,ntau           ! steps separated by spcng
    vv = vv1
    vw = vw1
    chi = (itau-1)*dtau
    !--- color triangles 
    if (coloring /= 0 .and. coloring /= byshape) then
      y = (itau-1)*dtau/thet
      rgb = (1-y)*rgbi + y*(rgbj + rgbi)/2.
      rgbv(1:3,1) = (1-y)*rgbi + y*(rgbj + rgbi)/2.
      rgbv(1:3,3) = rgbv(1:3,1)
      y = itau*dtau/thet
      rgbv(1:3,2) = (1-y)*rgbi + y*(rgbj + rgbi)/2.
      rgbv(1:3,4) = rgbv(1:3,2)
    endif
    !---
    call getrotS(vw,va,chi,mat,vec)  ! matrix to rotate vv around saddle axis
    call move_S(vv,mat,vec)          ! v3 is unit vec in direction of start water
    trv(1:3,1) = vv                 ! vertex
    trn(1:3,1) = vw - vv            ! normal vector
    vv = vv1
    chi = itau*dtau
    call getrotS(vw,va,chi,mat,vec)  ! matrix to rotate vv around saddle axis
    call move_S(vv,mat,vec)          ! v3 is unit vec in direction of start water
    trv(1:3,2) = vv                 ! vertex
    trn(1:3,2) = vw - vv            ! normal vector
    chi = 20*rad
    call getrotS(va1,va2,chi,mat,vec)  ! matrix to rotate around va1-va2 axis
    do j=20,360,20
      call move_S(vw,mat,vec)          ! move water
      vv = trv(1:3,1)
      call move_S(vv,mat,vec)          ! move vertex 1 --> 3
      trv(1:3,3) = vv                 ! vertex
      trn(1:3,3) = vw - vv            ! normal vector
      vv = trv(1:3,2)
      call move_S(vv,mat,vec)          ! move vertex 2 --> 4
      trv(1:3,4) = vv                 ! vertex
      trn(1:3,4) = vw - vv            ! normal vector
      !! draw two triangles
      avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
      !! write a triangle
      write(runit,*) 1
      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,rgb(1:3)
      if (smoothing) then
        write(runit,*) 7
        write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
        if (coloring /= 0 .and. coloring /= byshape) then
          write(runit,*) 17
          write(runit,'(9f8.3)') rgbv(1:3,1),rgbv(1:3,2),rgbv(1:3,3)
        endif
      endif
      avec = trv(1:3,2); bvec = trv(1:3,3); cvec = trv(1:3,4)
      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
      !! write a triangle
      write(runit,*) 1
      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,rgb(1:3)
      if (smoothing) then
        write(runit,*) 7
        write(runit,'(9f8.3)') trn(1:3,2),trn(1:3,3),trn(1:3,4)
        if (coloring /= 0 .and. coloring /= byshape) then
          write(runit,*) 17
          write(runit,'(9f8.3)') rgbv(1:3,2),rgbv(1:3,3),rgbv(1:3,4)
        endif
      endif
      trv(1:3,1) = trv(1:3,3)         !! copy 3 --> 1
      trn(1:3,1) = trn(1:3,3)
      trv(1:3,2) = trv(1:3,4)         !! copy 4 --> 2
      trn(1:3,2) = trn(1:3,4)
    enddo  !! j
  enddo  !! itau
  !!------------------------------------------------------------------------------
  !    
!
!    call cros(uab,vw,v2)    ! v2 perpendicular to abw plane
!    vp = va1 + v2           ! vp point above a
!    call getrotS(va1,vp,mtheta,mat,vec)  ! matrix to rotate uab to vw
!    call rotate_S(mat,uab,vw)          ! this vw is exactly mtheta from uab
!    vw = (r+rw)*vw          ! vw = water position relative to a
!    v1 = va1 + vw           ! location of water
!    !! call cros(vab,vw,v2)     ! v2 = perpendicular to abw
!    call cros(vab,v2,vt)     ! perpendicular to vab in abw plane
!    x = sqrt(dotprod(vt,vt))  ! length of vt
!    vt = vt/x           ! normalized vt
!    x = -dotprod(mxyz(1:3,ii),vt)  ! both are unit vectors
!    vt = rw*vt           ! vt, length rw
!    thet1 = acos(x)  ! cos(theta). vt always points away from vw
!    dtau1 = (thet1 - taumin)/real(ntau)    ! step size in tau
!    v2 = v1 - v2             ! point below water
!    !! 
!    ii = sides(2)    ! second point
!    ww = mxyz(1:3,ii) ! direction of mask point
!    call cros(uab,ww,w2)    ! v2 perpendicular to abw plane
!    vp = va1 + w2           ! vp point above a
!    call getrotS(va1,vp,mtheta,mat,vec)  ! matrix to rotate uab to vw
!    call rotate_S(mat,uab,ww)          ! this vw is exactly mtheta from uab
!    ww = (r+rw)*ww          ! vw = water position relative to a
!    w1 = va1 + ww           ! location of water
!    !! call cros(vab,ww,w2)     ! v2 = perpendicular to abw
!    call cros(vab,w2,wt)     ! perpendicular to vab in abw plane
!    x = sqrt(dotprod(wt,wt))  ! length of vt
!    wt = wt/x           ! normalized wt
!    x = -dotprod(mxyz(1:3,ii),wt)  ! both are unit vectors
!    wt = rw*wt           ! vt, length rw
!    thet2 = acos(x)  ! cos(theta). vt always points away from vw
!    dtau2 = (thet2 - taumin)/real(ntau)    ! step size in tau
!    w2 = w1 - w2             ! point below water
!    nn = 0
!    !!
!    !! first point in first triangle
!    chi = taumin 
!    call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
!    v3 = v1 + vt                     ! starting position
!    call move_S(v3,mat,vec)          ! rotation
!    nn = mod(nn,3)+1
!    trv(1:3,nn) = v3
!    trn(1:3,nn) = v1 - v3
!    !!
!    !! second point in first triangle
!    call getrotS(w1,w2,chi,mat,vec)  ! matrix and vector for rotation
!    w3 = w1 + wt                     ! starting position
!    call move_S(w3,mat,vec)          ! rotation
!    nn = mod(nn,3)+1
!    trv(1:3,nn) = w3
!    trn(1:3,nn) = w1 - w3
!    !! Diagnostic
!    !! write(*,*) 'dtau1, dtau2, thet1, thet2 ', dtau1, dtau2, thet1, thet2
!    !! vp = va1 + mxyz(:,sides(1)); wp = va1 + mxyz(:,sides(2))
!    !! write(*,'(a,6f8.3)') "mxyz a,b ",vp,wp
!    do itau=1,ntau           ! steps separated by spcng
!      chi = taumin + itau*dtau
!      !! rotate first point in side
!      call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
!      v3 = v1 + vt                     ! starting position
!      call move_S(v3,mat,vec)          ! rotation
!      nn = mod(nn,3)+1
!      trv(1:3,nn) = v3
!      trn(1:3,nn) = v1 - v3
!      avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
!      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
!      !! write a triangle
!      write(runit,*) 1
!      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,rgb(1:3)
!      if (smoothing) then
!        write(runit,*) 7
!        write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
!      endif
!      !! rotate second point in side
!      chi = taumin + itau*dtau
!      call getrotS(w1,w2,chi,mat,vec)  ! matrix and vector for rotation
!      w3 = w1 + wt                     ! starting position
!      call move_S(w3,mat,vec)          ! rotation
!      nn = mod(nn,3)+1
!      trv(1:3,nn) = w3
!      trn(1:3,nn) = w1 - w3
!      avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
!      if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
!      !! write a triangle
!      write(runit,*) 1
!      write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,rgb(1:3)
!      if (smoothing) then
!        write(runit,*) 7
!        write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
!      endif
!    enddo
!    vp = va1 + r*mxyz(:,sides(1))
!    nn = mod(nn,3)+1
!    trv(1:3,nn) = vp
!    trn(1:3,nn) = mxyz(:,sides(1))
!    avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
!    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
!    !! write a triangle
!    write(runit,*) 1
!    write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,rgb(1:3)
!    if (smoothing) then
!      write(runit,*) 7
!      write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
!    endif
!    vp = va1 + r*mxyz(:,sides(2))
!    nn = mod(nn,3)+1
!    trv(1:3,nn) = vp
!    trn(1:3,nn) = mxyz(:,sides(1))
!    avec = trv(1:3,1); bvec = trv(1:3,2); cvec = trv(1:3,3)
!    if (periodic) call scalethreeinbox(avec,bvec,cvec,lbox)
!!    !! write a triangle
!    write(runit,*) 1
!    write(runit,'(9f8.3,3f6.2)') avec,bvec,cvec,0.0,1.0,0.0
!    if (smoothing) then
!      write(runit,*) 7
!      write(runit,'(9f8.3)') trn(1:3,1),trn(1:3,2),trn(1:3,3)
!    endif
!    !! write(*,'(a,6f8.3)') "last a,b ",trv(:,mod(nn+1,3)+1),trv(:,nn)
!!  enddo
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
  subroutine drawsaddle(iunit,avec,bvec,smask,r1,iat,taumin,taumax,bb,ch,scale)
  implicit none
  !! output the saddle surface as dots
  !! Wed Jun 20 14:35:42 EDT 2001
  real,parameter :: spcng=0.25   !! angstrom spacing of saddle dots, reset to 0.35
  character(len=1),intent(in) :: ch
  integer(kind=KND),dimension(MASKSIZE),intent(in) :: smask
  real,intent(in) :: taumin,r1,bb,taumax,avec(3),bvec(3),scale
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
  !! Diagnostic
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
        !! write(iunit,'("MODEL ",i5)')  iat
        write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        i,ch,iat,scale*v3(1:3),0.,bb
        !! write(iunit,'("ENDMDL")') 
      enddo
    endif
  enddo
  end subroutine drawsaddle
  !!==========================================================
  subroutine renderheader(cen,scale)
  !! This is the header data for a Raster3D file (*.r3d)
  real,dimension(3),intent(in) :: cen
  real,intent(in) :: scale
  write(runit,*) 'Masker v. 5-JAN-04'
  write(runit,*) '10 10     tiles in x,y'
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
  ! write(*,*) "In tetrahedron: veci=",veci
  ! write(*,*) "In tetrahedron: vecj=",vecj
  ! write(*,*) "In tetrahedron: veck=",veck
  vec = vecj - veci
  dij = sqrt(dotprod(vec,vec))
  vec = veck - vecj
  djk = sqrt(dotprod(vec,vec))
  vec = veck - veci
  dik = sqrt(dotprod(vec,vec))
  !! 
  ! write(*,*) "In tetrahedron: ",iat,jat,kat,atype(iat),atype(jat),atype(kat)
  diw = atomlib(atype(iat))%r + rw
  djw = atomlib(atype(jat))%r + rw
  dkw = atomlib(atype(kat))%r + rw
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
  yp = diw*diw - x*x - y*y
  if (yp < 0.) return  !! angles violate triangle inequality (prism inequality??) Jul  2 15:04:29 EDT 2006
  !!--------------------
  !! NOTE: It may be possible that the cosinerule gives valid angles for all four triangles,
  !! but this is not sufficient condition that the four triangles constitute a tetrahedron.
  !! The sum of any two angles having a common vertex must be greater than the third angle
  !! at the same vertex. I'm not sure this condition is always met. --cb
  !!--------------------
  d = sqrt(yp)
  if (hand < 0) d = -d
  wxyz = (/x,y,d/)
  !  write(*,*) "In tetrahedron: wxyz=",wxyz
  !! Diagnostic : distance to water
  !  vec = wxyz ;  dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - (/dij,0.,0./) ;  dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5,$)') dd
  !  vec = wxyz - (/dik*cos(ajik),dik*sin(ajik),0./); dd=sqrt(dotprod(vec,vec)); write(*,'(f8.5)') dd
  !  write(*,*) dik*cos(ajik),dik*sin(ajik),0.
  
  !! Get frame of reference
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
  subroutine masker_allocatoms(nat)
  !! This routine may be necessary if private module arrays
  !! are to be allocated.
  !! Here the private derived type 'atype' is allocated and
  !! an array of water locations (tetrahedron verteces) is allocated.
  !! Since we don't know how many waters will be needed, we allocate 
  !! 'nat', one per atom. That should be more than enough.
  integer,intent(in) :: nat
  integer :: ios=0
  if (allocated(atype)) deallocate(atype)
  allocate(atype(nat),stat=ios); if (ios/=0) stop 'masker_allocatoms - ERROR: allocating atype()'
  atype = 1  !! default to carbon
  if (allocated(water)) deallocate(water)
  allocate(water(watfac*nat),stat=ios)   ; if (ios/=0) stop 'masker_allocatoms - ERROR: allocating water()'
  !! this is more than enough. what should it really be?
  end subroutine masker_allocatoms
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
    stop 'countbits does not understand KND setting!!! should be 2 or 4'
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
  !! real ,intrinsic :: ran
  integer :: imask
  !! chose random point on the surface of the template mask
  !! x = ran(nseed)
  !write(0,*) "in ranvec"
  call random_number(x)
  imask = nint(x*MAXATOM)
  call random_number(x)
  !! x = sig*ran(nseed)
  x = sig*x
  !write(0,*) "in ranvec",mxyz(1:3,imask)
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
    if (itheta == 0) imask = 0  !! almost completely embedded, ignore tiny error. Dec 8 2006 cb
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
  subroutine setcoloring(str)
    !! integer,private,parameter :: surfacetension=1,byatom=2,byshape=3,bycharge=4
    implicit none
    character(len=*) :: str
    if (str(1:6)=="surface") then
      coloring = surfacetension
    elseif (str(1:4)=="atom") then
      coloring = byatom
    elseif (str(1:6)=="byatom") then
      coloring = byatom
    elseif (str(1:8)=="bycharge") then
      coloring = bycharge
    elseif (str(1:6)=="charge") then
      coloring = bycharge
    elseif (str(1:5)=="shape") then
      coloring = byshape
    elseif (str(1:7)=="byshape") then
      coloring = byshape
    else
      coloring = 0
    endif
  end subroutine setcoloring
  !----------------------------------------------------!
  subroutine surfacecolor(iat,ires,shape,rgb)
    implicit none
    integer,intent(in) :: iat,shape,ires
    real :: tension !! surface tension from ATOMLIB
    real,dimension(3),intent(out) :: rgb
    real,dimension(3) :: polarrgb=(/0.,1.,0./),nprgb=(/1.,1.,1./)
    real,parameter :: mintension=-1., maxtension=1.
    real :: x,y
    !--- function for convertign surface tension to  a RGB color
    !-- interpolate from -1 to 1, -1 is polar, therefore greenish, +1 is nonpolar, therefore whitish
    select case (coloring)
      case (surfacetension)
        tension = atomlib(atype(iat))%w1
        x = max(tension,mintension)
        x = min(x,maxtension)
        y = (x - mintension)/(maxtension-mintension)
        rgb = (1-y)*polarrgb + y*nprgb
      case (byshape)
        if (shape==1) then !! VDW
          rgb = (/0.,1.,0./)
        elseif (shape==2) then !! saddle
          rgb = (/0.5,1.,0./)
        elseif (shape==3) then  !! bowl
          rgb = (/0.8,0.8,0./)
        else
          rgb = (/1.0,1.0,1./)
        endif
      case (bycharge)
        !if (ires==0) then  !! if there is no residue number, use surface tension
          x = atomlib(atype(iat))%w1
        !else
        !  x = atomcharge(iat,ires)
        !endif
        if (x < 0.) then  !! color reddish
            x = -x
          rgb = (/1.,1-x,1-x/)
        else      !! color bluish
          rgb = (/1-x,1-x,1./)
        endif
      case (byatom)
        select case(atomlib(atype(iat))%name(1:1))
          case ("C")
            rgb = (/0.8,0.8,0.8/)
          case ("O") !! oxygen
            rgb = (/1.,0.1,0.1/)
          case ("N")  !! nitrogen
            rgb = (/0.2,0.3,0.8/)
          case ("S")  !! sulfur
            rgb = (/0.7,0.7,0.0/)
          case default
            rgb = (/1.,1.,1./)
        end select
    end select 
  end subroutine surfacecolor
  !----------------------------------------------------!
  subroutine masker_getatypeptr(atype_ptr)
    implicit none
    integer, dimension(:), pointer::atype_ptr

    atype_ptr=>atype
    return
  end subroutine masker_getatypeptr
  !----------------------------------------------------!
  subroutine masker_deallo
    implicit none

    if (allocated(nposit_mask)) deallocate(nposit_mask)
    if (allocated(atomlib)) deallocate(atomlib)
    if (allocated(psiposit)) deallocate(psiposit)
    if (allocated(nphi)) deallocate(nphi)
    if (allocated(masklib)) deallocate(masklib)
    if (allocated(maskpsi)) deallocate(maskpsi)
    if (allocated(maskphi)) deallocate(maskphi)
    if (allocated(collarsize)) deallocate(collarsize)
    if (allocated(slope)) deallocate(slope)
    if (allocated(interc)) deallocate(interc)
    if (allocated(tri)) deallocate(tri)
    if (allocated(mxyz)) deallocate(mxyz)
    if (allocated(atype)) deallocate(atype)
    if (allocated(water)) deallocate(water)
    if (allocated(tau)) deallocate(tau)
    return
  end subroutine masker_deallo
  !----------------------------------------------------!

end module masker
