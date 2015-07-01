module masker2
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
  !! Required files: masktemplate (t2048.pdb)
  !!                 masklib (t2048.mas)
  !!                 mask.log (t2048.log)
  !! Input file : inputpdb (pdb format)
  !! Output file : output (pdb format, surface dots in HETATM lines)
  !! Stdout : surface area with breakdown (sasa, saddle, bowl)
  !! -------------------------------------------------------------
  !*=========================================================================
  !** Copyright 2001 Chris Bystroff
  !**
  !** Permission is hereby granted to use, execute, copy, distribute, modify,
  !** and distribute in modified form this software for non-profit purposes,
  !** provided this notice is retained in any and all copies.
  !**
  !** For for-profit usage, please contact bystrc@rpi.edu
  !*=========================================================================
  !!  CHANGES IN MASKER2.F90 versus MASKER.F90
  !!========= atomic surface areas =============
  !! Current version returns the atomic surface areas in three parts
  !!====== Pairwise numerical derivatives of MS =========
  ! To get the pairwise derivative of the surface w/respect to distance ij,
  ! we shift the distance by d = cos(theta) - cos(theta+dtheta), (ij omitted)
  ! which of course is smaller for smaller theta.
  ! Then the derivative is dA/d, where dA = d(sasa)*wv + d(saddle)*ws + d(bowl)*wb
  ! Note that this is actually in energy terms, not surface area, since the
  ! units of w are kJ/A^2. Thus we get kJ/A, energy vs distance.
  ! d, d(sasa) and d(saddle) are pre-calculated for each theta.
  ! d(sasa) is summed by one additional masking step, masking out all 
  ! but the points occluded by a shift of d. This value multiplied by the 
  ! edgemask ratio is = d(sasa)
  ! d(saddle) is the difference between saddles of theta and theta+dtheta 
  ! using the same arc length. Care is taken at the point where taumin 
  ! becomes a non-zero angle.
  ! d(bowl) is not pre-calculated (how??). Instead, the distance ij is changed 
  ! by d and the bowl is recalculated with the same ik and kj distances and
  ! the same orientation.
  ! CB  Sat Aug 11 07:37:18 EST 2001
  !!
  integer,parameter :: MAXATOM=4096,MAXPSI=360
  integer,parameter :: MASKSIZE=MAXATOM/32
  real :: DTHETA   !! degree increments
  integer,parameter :: MMASK=6600
  integer :: NTHETA
  integer,parameter :: KC=1   !! thickness of collars (times DTHETA)
  integer :: nmask,npsi,ishow=0
  integer,dimension(MASKSIZE) :: emask,vmask
  integer,dimension(:),allocatable :: psiposit,nphi
  integer,dimension(:,:,:),allocatable :: masklib
  real,dimension(:),allocatable :: tau
  integer(kind=2),dimension(:,:),allocatable :: collarsize
  real,dimension(:),allocatable :: maskpsi,maskphi
  real,dimension(3,MAXATOM) :: mxyz
  real,dimension(3) :: cvec
  real,parameter :: pi=3.1415927,radius=10.
  real,parameter :: rad=pi/180.,rw=1.4
  real,parameter :: wsphere=rw*rw*4*pi
  real,parameter :: fscale=0.01   !! converts kJ/mol/A to g*A/mol/ps
  real,parameter :: BOXSCALE=10.  !! scaling factor if periodic
  integer,dimension(:),allocatable :: atype
  real :: dpsi 
  !! uses pi,npsi,dpsi,psiposit
  real :: maxd,taumin,r1,lbox
  real,dimension(3) :: showvec
  logical :: drawing, periodic
  type atomtype
    real :: r        ! united atom radius
    real :: m        ! united atom mass
    real :: w1,w2,w3 ! convert surface area to energy (sasa, saddle, bowl) kJ/mol/A^2
    character(len=4) :: name   !  element or other atom name 
    !! format for atomlib is (a4,4f8.3) 
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
  subroutine initmasks
  
  character(len=80) :: aline,pdbfile,binfile,logfile,outfile,template
  character(len=80) :: vmaskfile,libfile
  integer :: i,j,k,L,m,jarg,ios,iargc,natm,kcc=KC
  integer,dimension(MASKSIZE) :: amask,bmask,cmask,dmask
  integer :: nn,ipsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv
  integer :: imask1,itheta1,imask2,itheta2,iat,nat,jat,kat,ires
  real :: dphi,sasa,sphere,ssasa,bsasa,col
  real :: phi,psi,theta,x,y,d,dd,r2
  real,dimension(3) :: avec,bvec,wxyz,showvec
  !! defaults
  DTHETA = 4.5
  NTHETA = int(90./DTHETA)
  dpsi = DTHETA*rad   !! stepsize in radians
  !! 
  call getenv("MASKLIB",binfile)
  call getenv("VMASK",vmaskfile)
  call getenv("MASKTEMPLATE",template)
  call getenv("MASKDAT",logfile)
  call getenv("ATOMLIB",libfile)
  !!
  open(15,file=libfile,status='old',form='formatted')
  read(15,*) nattype
  allocate(atomlib(nattype),stat=ios)
  if (ios/=0) stop 'masker.mod: Error allocating atomlib'
  do i=1,nattype
    read(15,'(a4,5f8.3)') atomlib(i)%name, atomlib(i)%m, atomlib(i)%r,  &
           atomlib(i)%w1,  atomlib(i)%w2,  atomlib(i)%w3
  enddo
  close(15)
  !! at this point there should be exactly MAXATOM in xyz()
  !! read log file from binarymask
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
      allocate(tau(NTHETA),stat=ios)
      if (ios/=0) stop 'Problem allocating tau'
      write(*,*) 'psiposit and nphi allocated',npsi
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
  write(*,'("NPSI =",i7)') npsi
  write(*,'("NMASK=",i7)') nmask
  allocate(masklib(MASKSIZE,NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for masklib'
  allocate(maskpsi(nmask),maskphi(nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for maskpsi,maskphi'
  allocate(collarsize(0:NTHETA,nmask),stat=ios)
  if (ios/=0) stop 'Unable to allocate memory for collarsize'
   write(*,*) 'Memory allocated for masklib',MASKSIZE, NTHETA, nmask
  !! read binary format mask library, created by xbinarymask
  open(12,file=binfile,status='old',form='unformatted')
  read(12,iostat=ios) masklib(1:MASKSIZE,1:NTHETA,1:nmask)
  if (ios/=0) stop 'initmasks: Error read mask library'
  close(12)
  !! diagnostic
  write(*,*) 'Done reading masklib'
  !! get viewable dots mask
  open(12,file=vmaskfile,status='old',form='unformatted')
  read(12,iostat=ios) vmask(1:MASKSIZE)
  if (ios/=0) stop 'initmasks: Error read viewable mask'
  close(12)
  !! diagnostic
  write(*,*) 'Done reading vmask file'
  ! write(*,*)  masklib(1:MASKSIZE,10,50)
  !! get psi and phi for each mask
  !! AND
  !! Calculate collarsize and tau for each theta and imask
  !! collarsize is the number of bits in the edge mask for angle theta
  !! tau is the angular width of the half-saddle for angle theta and r1
  !!
  do ipsi=0,npsi
    psi = ipsi*dpsi
    dphi = 2*pi/real(nphi(ipsi))
    do iphi=1,nphi(ipsi)
      phi = (iphi-1)*dphi
      imask = psiposit(ipsi) + iphi - 1
      maskpsi(imask) = psi
      maskphi(imask) = phi
      !! write(*,*) imask,psi,phi
      !! get collarsize for each difference mask
      do itheta=0,NTHETA-1
        !! NOTE: KC and kcc refer to the edge thickness.
        !! Thickness can be fixed at 1 is there is always a point at the vertex
        !! of intersecting edges.  --cb
        !! if (itheta > NTHETA-KC) then; kcc = 1; else; kcc=KC; endif
        theta = itheta*dpsi
        !! If the depth of the saddle is greater than one water radius
        !! then two saddles connect, and the integration limits are zero to 90-theta.
        !! Othewise, calculate the lower integration limit (taumin).
        if (imask==1) then   !! only do this for one mask. Same for all masks.
          if (itheta > 0) then
            if ((r1+rw)*sin(theta) >= rw) then
              tau(itheta) = (pi/2.) - theta
            else
              x = acos((r1+rw)*sin(theta)/rw)
              tau(itheta) = (pi/2.) - theta - x
            endif
          endif
        endif
        if (itheta > 0) then
          amask = iand(masklib(:,itheta,imask),not(masklib(:,itheta+kcc,imask)))
        else
          amask = -1
          amask = iand(amask,not(masklib(:,itheta+kcc,imask)))
        endif
        j = 0
        do i=1,MAXATOM
          ibyte = (i-1)/32 + 1
          ibit = mod(i-1,32)
          if (btest(amask(ibyte),ibit)) then
            j = j + 1
          endif
        enddo
        !! KC (global) is the thickness of the collar (kcc is local)
        !! NO LONGER ASSUMING ALL IMASK HAVE SAME SET OF COLLAR SIZES.
        collarsize(itheta,imask) = j    !! approximately the same for all imask
        !! write(*,'("theta tau size ",2f8.2,f8.0)') &
        !!   theta/rad,tau(itheta)/rad,collarsize(itheta)
      enddo
      !! this number must be non-zero, but doesn't matter since tau=0.
      collarsize(NTHETA,imask) = 1  
    enddo
  enddo
  !! NTHETA is 90deg, only possible if dd ~= 0. No saddle? 
  !! =====> NOTE: If hydrogens are used, this will have to be changed! <=====
  tau(NTHETA) = 0.0
  !! diagnostic
   write(*,*) 'Done getting phi/psi. Getting collars.'
  !!
  ! write(*,*) 'Opening template PDB file'
  open(11,file=template,status='old',form='formatted')
  i = 0
  do 
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    if (aline(1:5) /= "ATOM ") cycle
    i = i + 1
    if (i > MAXATOM) stop 'Too many points in template file'
    read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
    if (ios/=0) then
      write(*,*) 'Error reading template coordinate file, at i=',i
      stop 'initmasks: Error reading template coordinate file'
    endif
    mxyz(1:3,i) = mxyz(1:3,i)/10.   !! template has 10A radius
  enddo
  close(11)
  !! 
  end subroutine initmasks
  !!----------------------------------------------------------------!!
  subroutine allocatoms(nat)
  integer,intent(in) :: nat
  integer :: ios=0
  allocate(atype(nat),stat=ios)
  if (ios/=0) stop 'Error allocating atype'
  end subroutine allocatoms
  !!----------------------------------------------------------------!!
  subroutine getms(xyz,nat,ms,atomms,msderiv)
  integer, intent(in) :: nat
  real,dimension(3,nat),intent(in) :: xyz
  real,dimension(3,nat),intent(out) :: atomms !! 3 types of surface
  real,dimension(nat,nat),intent(out) :: msderiv !! pairwise deriv E vs d.
  real,intent(out) :: ms
  !
  integer :: i,j,k,L,m,jarg,ios,iargc,natm
  integer,dimension(MASKSIZE) :: amask,bmask,cmask,dmask
  integer :: nn,ipsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv,n
  integer :: imask1,itheta1,imask2,itheta2,iat,jat,kat,ires
  integer,dimension(100) :: saveimask, saveitheta,saveiat
  real,dimension(100) :: savetheta
  real :: phi,psi,theta,x,y,d,dd,r2,dsasa,th,xx,rat,c2
  integer :: kskip=100,kcc=KC
  real :: dphi,sasa,sphere,ssasa,bsasa,col,dij
  real,dimension(3) :: avec,bvec,wxyz,showvec
  real,dimension(3) :: veci,vecj,veck,vech,uvec,vvec
  real :: dotprod
  !
  i = 0
  sasa = 0.0    !! total solvent accessible molecular surface area
  ssasa = 0.0   !! total "saddle" convex/concave surface area
  bsasa = 0.0   !! total "bowl" concave/concave surface area
  !! diagnostic?? or permanent, or switchable?
  atomms = 0.   !! initialize atomic ms to 0.
  msderiv = 0.   !! initialize derivatives to 0.
  periodic = (lbox > 0.)
  if (ishow /= 0 .and.drawing) then
    write(*,'("Show surface around residue ",i4," only.")') ishow
  endif
  
  sasa = 0.0    !! total solvent accessible molecular surface area
  ssasa = 0.0   !! total "saddle" convex/concave surface area
  bsasa = 0.0   !! total "bowl" concave/concave surface area
    
  !! Loop over all atoms, get sasa, saddle and bowl surfaces
  do iat=1,nat
    nm = 0
    amask = -1                 ! init all bits = 1
    bmask = -1
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
      !! write(*,*) iat,jat,dd,maxd
      avec = avec/dd
      imask = getimask(avec)  !! avec must be a unit vector
      !   psi = acos(avec(3))
      !   phi = atan2(avec(2),avec(1))
      theta = gettheta(r1+rw,r2+rw,dd)   !! Heron's rule
      if (theta <= 0.) cycle
      itheta = nint((theta/rad)/DTHETA)
      ! if (iat==1.and.jat==2) write(*,*) 'Theta: ',theta/rad
      ! write(*,'(2i4,4f6.1,4i5)') iat,jat,psi/rad,phi/rad,theta/rad,dd,ipsi,iphi,itheta,imask
      if (itheta <= 0) cycle    !! atom out of range (borderline)
              !! borderline atoms should have gradient ~ #of points in itheta=1 mask.
      if (itheta > NTHETA) stop 'Theta out of range'  !! this would spot a bug
      !! remove the bits covered by atom jat
      bmask = iand(masklib(:,itheta,imask),amask)     !!vector argument works!!
      if (any(bmask /= amask) ) then
        amask = bmask
        !! keep a list of close atoms, for the next part
        nm = nm + 1
        saveiat(nm) = jat
        saveimask(nm) = imask
        saveitheta(nm) = itheta
        savetheta(nm) = theta
      endif
    enddo
    j = 0
    do i=1,MAXATOM
      ibyte = (i-1)/32 + 1
      ibit = mod(i-1,32)
      if (btest(amask(ibyte),ibit)) then
        j = j + 1
      endif
    enddo
    !! Calculate convex/convex molecular surface for this atom
    sphere = 4*pi*atomlib(atype(iat))%r**2   !! sphere is used later in this loop
    x = sphere*(real(j)/real(MAXATOM))
    !! Tue Aug  7 10:19:36 EST 2001
    !! diagnostic?
    atomms(1,iat) = x
    sasa = sasa + x*atomlib(atype(iat))%w1
    if (drawing) then
      cmask = iand(amask,vmask)
      call drawsurface(13,xyz(1:3,iat),r1,iat,cmask,50.,"V",1)  !! output the sasa surface
    endif
    i = 0
    cmask = 0
    nn = 0
    do k=1,nm
      kat = saveiat(k)
      imask = saveimask(k)
      itheta = saveitheta(k)
      theta = savetheta(k)
      !! KC was added to provide the option of double-width (KC=2) 
      !! saddles. Mon Jun 25 16:57:13 EDT 2001
      !! Collarsizes of width KC are calculated. 
      !! Obselete if edges are thick enough to guarantee a point at every vertex.
      !!   Sat Jun 30 21:30:14 EDT 2001
      !! if (itheta > NTHETA-KC) then; kcc = 1; else; kcc=KC; endif
      bmask = iand(not(masklib(:,itheta+kcc,imask)),amask)
      cmask = ior(cmask,bmask)
      j = 0
      do i=1,MAXATOM
        ibyte = (i-1)/32 + 1
        ibit = mod(i-1,32)
        if (btest(bmask(ibyte),ibit)) then
          j = j + 1
        endif
      enddo
      !! 
      !! Saddle surface:
      !! The following is an "exact approximation", meaning
      !! that the estimate is as good as the numerical estimates
      !! of the arc distance, theta, phi, psi.
      !! (equations derived using Mathematica)
      !! Mon May  7 10:34:43 EDT 2001
      !! 
      theta = itheta*dpsi  !! why not use the real theta?  Mon Jun 18 10:49:46 EDT 2001
          !! Sat Aug 11 08:08:19 EST 2001 The real theta gives a less realistic
          !! PMF than the discretized theta. (saw-toothed, rather than terraced.)
      y = (r1+rw)*sin(theta)
      rat = (real(j)/real(collarsize(itheta,imask)))
      col = 2*pi*rat
      x = col*((0.5*pi-theta)*y - rw*cos(theta))
      taumin = 0.
      if (y < rw) then
        taumin = acos(y/rw)
        x = x + col*(rw*sin(taumin) - taumin*y)
      endif
      !! diagnostic?
      atomms(2,iat) = atomms(2,iat) + x  !! area
      ssasa = ssasa + x*atomlib(atype(iat))%w2   !! energy
      if (drawing) then
        y = (pi/2) - theta
        call drawsaddle(13,xyz(1:3,iat),xyz(1:3,kat),bmask,r1,iat,taumin,y,50.,"S")  
        !! output the saddle surface
      endif
      !! if calculating derivs
        r2 = atomlib(atype(kat))%r
        th = (itheta+1)*dpsi  !! theta for next shell
        y = (r1+rw)*sin(th)
        c2 = sqrt((r2+rw)**2 - (y)**2)/y   !! This will fail if r1 >> r2
        dij = (r1+rw)*cos(th) + (r2+rw)*c2   !! distance assoc w dpsi shift
        dsasa = -(real(j)/real(MAXATOM))*sphere*atomlib(atype(iat))%w1  !! kJ/mol SASA lost
                   !!  j = points in edgemask
        msderiv(iat,kat) = msderiv(iat,kat) + dsasa/dij  
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
        msderiv(iat,kat) = msderiv(iat,kat) + (xx - x)*atomlib(atype(iat))%w2/dij
      !! end if calculating derivs
      if (any(bmask /= 0)) then
        !! keep a reduced list of close atoms.
        nn = nn + 1    !!  note: nn <= k, so this is legal.
        saveiat(nn) = kat
        saveimask(nn) = imask
        saveitheta(nn) = itheta
        savetheta(nn) = theta
      endif
    enddo
    !!
    !! Using the reduced list of close atoms, get the intersections
    !! between the 'collars' .
    !! Each intersection corresponds to a water. Get the surface
    !! of the water that is in the spherical triangle between the three 
    !! atoms.  Subtract any part of that surface that falls below the
    !! plane of those three atoms. (iat, kat, jat)
    !! 
    !! Here we use "double-width" collars to be sure of finding a vertex
    !! point.  It was found that single width callors sometimes missed
    !! a vertex. However, using double-width collars in the drawsaddle routine
    !! makes an ugly saddle.  So KC=1    Mon Jun 25 16:53:38 EDT 2001
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
      do m=1,nn
        jat = saveiat(m)
        if (jat <= kat) cycle  !! calculate bowl iff kat is the second atom
        imask2 = saveimask(m)
        itheta2 = saveitheta(m)
        ! theta2 = savetheta(m)
        !! write(*,*) imask1,itheta1,imask2,itheta2
        !! check the distance jat -> kat
          if (periodic) then ! get nearest copy of jat, return vector
            call boxaround(xyz(1:3,jat),xyz(1:3,kat),lbox,avec)
          else
            avec = xyz(1:3,jat) - xyz(1:3,kat)
          endif
          dd = dotprod(avec,avec)
          maxd =  atomlib(atype(jat))%r +  atomlib(atype(kat))%r + 2*rw
          if (dd > maxd*maxd) cycle    !! atom out of range, triangle inequality
        bmask = iand(not(masklib(:,itheta1+kcc,imask1)),masklib(:,itheta1,imask1) )
        bmask = iand(bmask,not(masklib(:,itheta2+kcc,imask2)) )
        bmask = iand(bmask,masklib(:,itheta2,imask2) )
        bmask = iand(cmask,bmask)  !! bmask is current verteces
        !! diagnostic
        !! dmask = ior(dmask,bmask)  !! dmask is all verteces (used only for display)
        !! call drawsurface(13,xyz(1:3,iat),r1,iat,dmask, -10.,"D",1) !! output the surface
        !!
        if (any(bmask /= 0)) then !! if there is at least one vertex
          !! save atom positions
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
          !! diagnostic
          ! write(*,*) 'Vertex:',iat,jat,kat
          !! Retrieve unit vector to atom k: avec
          psi = maskpsi(imask1)
          phi = maskphi(imask1)
          avec(3) = cos(psi)
          avec(2) = sin(psi)*sin(phi)
          avec(1) = sin(psi)*cos(phi)
          !! Retrieve unit vector to atom m: bvec
          psi = maskpsi(imask2)
          phi = maskphi(imask2)
          bvec(3) = cos(psi)
          bvec(2) = sin(psi)*sin(phi)
          bvec(1) = sin(psi)*cos(phi)
          !! get the cross product. This is a normal to the 3-atom plane
          call cros(avec,bvec,cvec)
          d = sqrt(dotprod(cvec,cvec))
          cvec = cvec/d
          !! Get the mask that splits the vertex mask down the middle
          imask = getimask(cvec)
          itheta = NTHETA    !! hemisphere mask
          x = 0. !! bowl area
          n = 0  !! number of verteces
          !! First check one side for verteces...
          emask = iand(masklib(:,itheta,imask),bmask)
          if (any(emask /= 0)) then
            n = 1
            !! diagnostic
            !! write(*,*) 'Vertex 1:',iat,kat,jat
            !! calculate surface, returning surface mask and water coords
            call waterbowl(veci,veck,vecj,x,wxyz,imask,iat,kat,jat) 
                 !! also uses/returns emask, uses cvec
            !! diagnostic, write water coords
            !! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
            !!   iat,"W",iat,wxyz(1:3),0.,10.0
            if (drawing) then  !! do this to force drawing surface for bowls
              ! dmask = iand(emask,vmask)
              dmask = emask
              call drawsurface(13,wxyz(1:3),rw,iat,dmask, 10.,"B",1) !! output the surface
            endif
          endif
          !! ...then check the other
          emask = iand(not(masklib(:,itheta,imask)),bmask)
          if (any(emask /= 0)) then
            n = n + 1
            !! diagnostic
            ! write(*,*) 'Vertex 2:',iat,jat,kat
            if (x==0..or.drawing) then  
              !! diagnostic
              !!write(*,*) 'Vertex 2:',iat,kat,jat
              call waterbowl(veci,veck,vecj,y,wxyz,imask,iat,kat,jat) 
                 !! also uses/returns emask, uses cvec
              x = x + y
              !! diagnostic, write water coords
              !! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
              !! iat,"W",iat,wxyz(1:3),0.,20.0
              if (drawing) then  !! do this to force drawing surface of bowls
                ! dmask = iand(emask,vmask)
                dmask = emask
                call drawsurface(13,wxyz(1:3),rw,iat,dmask, 10.,"B",1) !! output the surface
              endif
            else
              x = 2*x
            endif
          endif
          !! If calculating bowl derivs
          dij = -0.1   !! fixed shift 0.1A
          vvec = vecj - veci
          call unitvec(vvec,uvec)  !! unit vector from i to j
          vech = veci - uvec*dij  !! position of i shifted dij toward j
          !!    write(*,*) 'Deriv 1:',iat,kat,jat,vech
          call waterbowl(vech,veck,vecj,xx,wxyz,imask,iat,kat,jat) 
          msderiv(iat,jat) = msderiv(iat,jat) + &
            (xx*n - x)*(atomlib(atype(iat))%w3+atomlib(atype(jat))%w3)/(2*dij)
          vvec = veck - vecj
          call unitvec(vvec,uvec)  !! unit vector from j to k
          vech = vecj - uvec*dij  !! position of j shifted dij toward k
          !!    write(*,*) 'Deriv 3:',iat,kat,jat,vech
          call waterbowl(veci,veck,vech,xx,wxyz,imask,iat,kat,jat) 
          msderiv(jat,kat) = msderiv(jat,kat) + &
            (xx*n - x)*(atomlib(atype(jat))%w3+atomlib(atype(kat))%w3)/(2*dij)
          vvec = veci - veck
          call unitvec(vvec,uvec)  !! unit vector from k to i
          vech = veck - uvec*dij  !! position of k shifted dij toward i
          !!    write(*,*) 'Deriv 2:',iat,kat,jat,vech
          call waterbowl(veci,vech,vecj,xx,wxyz,imask,iat,kat,jat) 
          msderiv(kat,iat) = msderiv(kat,iat) + &
            (xx*n - x)*(atomlib(atype(kat))%w3+atomlib(atype(iat))%w3)/(2*dij)
          !! end if calculating bowl derivs
          !! diagnostic??
          !! For lack of a better way, the bowl surface is divided evenly between the 3 atoms
          atomms(3,iat) = atomms(3,iat) + x/3.
          atomms(3,jat) = atomms(3,jat) + x/3.
          atomms(3,kat) = atomms(3,kat) + x/3.
          bsasa = bsasa + x*(atomlib(atype(iat))%w3 + &
                             atomlib(atype(jat))%w3 + &
                             atomlib(atype(kat))%w3)/3.
        endif
      enddo
    enddo
  enddo
  x = sasa + ssasa + bsasa
  if (drawing) then
     write(*,'("Weighted VDW area = ",f9.2)') sasa
     write(*,'("Weighted saddle area  = ",f9.2)') ssasa
     write(*,'("Weighted bowl area    = ",f9.2)') bsasa
     write(*,'("Weighted surface area = ",f9.2)') x
  endif
  ms = x
  close(13)
  end subroutine getms ! (xyz,nat,ms,atomms,msderiv)
  !!----------------------------------------------------------------!!
  real function gettheta(r1,r2,dd)
  !! returns theta in radians for angle between r1 and dd
  !! for triangle w/sides r1,r2,dd
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
  !! three vectors are the atoms responsible for the make 'emask'
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
  !!  !!   diagnostic
  !! write(*,*) 'Waterbowl start', iat, jat, kat
  !! write(*,*) iat, veci
  !! write(*,*) jat, vecj
  !! write(*,*) kat, veck
  j = 0
  jvec = 0.
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
    if (btest(emask(ibyte),ibit)) then
      j = i
      !! For more accurate placement of the water, we average all
      !! vertex points. This may slow things down. How much?
      !! Sat Jun 30 21:32:57 EDT 2001
      jvec = jvec + mxyz(:,i) 
      exit
      !! alternatively, we may compute the water position from the tetrahedron
      !! dimensions...
      !! If we are calculating derivatives, then we should do this because
      !! the mask need not be re-calculated 4 times. Whether it is R-handed 
      !! or L-handed only matters for plotting.
    endif
  enddo
  if (j==0) then
    wxyz = 0.
    x = 0.
    return
  endif
  d = sqrt(dotprod(jvec,jvec))
  jvec = jvec/d     !! normalized, averaged water vector
  !! diagnostic
  ! write(*,*) 'emask not=0'
  !! wxyz is water position in absolute coords (A)
  wxyz = (r1+rw)*jvec   !! current value of r1 set in calling routine
  d = dotprod(cvec,wxyz)
  !! diagnostic
  avec = veci + wxyz       ! w
  write(*,'(3i5,4f9.3)') iat,jat,kat,avec(1:3),d
  if (periodic) then ! get nearest copy of jat, kat
    call boxaround(veci,vecj,lbox,bvec) ! bvec = i->j
    vjj = veci + bvec ! jat
    call boxaround(veci,veck,lbox,bvec) ! bvec = i->k
    vkk = veci + bvec ! kat
  else
    vjj = vecj
    vkk = veck
  endif
  if (d > 0.) then !! R-handed. Use iat,jat,kat and -cvec
    !! get water position using RH tetrahedron , replace d
    call tetrahedron(veci,vjj,vkk,wxyz,d,iat,jat,kat)
    if (drawing) then  !! write the water atoms
       write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') iat,"W",iat,wxyz(1:3),0.,0.
    endif
    !!
    if (d < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
      theta = acos(d/rw)
      itheta = nint((theta/rad)/DTHETA)
      dvec = -1*cvec
      jmask = getimask(dvec)        ! dont need dvec after this?
      emask = masklib(:,itheta,jmask)
    else
      emask = -1
    endif
    !! mask ijw plane =================
    !! diagnostic
    write(*,*) veci, vjj
    avec = veci - wxyz       ! w->i
    bvec = vjj - wxyz        ! w->j
    call cros(avec,bvec,dvec)
    y = sqrt(dotprod(dvec,dvec))
    dvec = dvec/y
    jmask = getimask(dvec)
    !! write(*,'(5x,3f7.3,i7)') dvec(1:3),jmask
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
  else !! L-handed. Use iat,kat,jat and +cvec (current imask)
    !! get water position using RH tetrahedron
    call tetrahedron(veci,vkk,vjj,wxyz,d,iat,kat,jat)
    !!
    if (drawing) then  !! write the water atoms
       write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') iat,"W",iat,wxyz(1:3),0.,-1.
    endif
    if (d < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
      theta = acos(d/rw)
      itheta = nint((theta/rad)/DTHETA)
      jmask = imask  !!  NOTE: now using local variable imask. 
      emask = masklib(:,itheta,jmask)
    else
      emask = -1
    endif
    !! mask ikw plane =================
    avec = veci - wxyz  !  w->i
    bvec = vkk - wxyz   ! w->k
    call cros(avec,bvec,dvec)
    dvec = dvec/sqrt(dotprod(dvec,dvec))
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
    !! mask kjw plane =================
    avec = bvec          ! w->k
    bvec = vjj - wxyz    ! w->j
    call cros(avec,bvec,dvec)
    dvec = dvec/sqrt(dotprod(dvec,dvec))
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
    !! mask jiw plane =================
    avec = bvec         ! w->j
    bvec = veci - wxyz  ! w->i
    call cros(avec,bvec,dvec)
    dvec = dvec/sqrt(dotprod(dvec,dvec))
    jmask = getimask(dvec)
    emask = iand(masklib(:,NTHETA,jmask),emask)
  endif
  j = 0
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
    if (btest(emask(ibyte),ibit)) j = j + 1
  enddo
  x = wsphere*real(j)/real(MAXATOM)  !! area in A^2
  end subroutine waterbowl
  !!------------------------------------------------------------
  integer function getimask(uvec)
  !! return the mask index for a unit vector
  !! uses pi,npsi,dpsi,psiposit
  !!   from the main program
  real,intent(in) :: uvec(3)
  integer :: iphi,ipsi
  real :: phi,psi,dphi
  psi = acos(uvec(3))
  ipsi = mod(int((psi/dpsi)+npsi+0.5),npsi)
  if (ipsi > npsi) stop 'Bug in psi calculation: getimask'
  if (ipsi==0) then
    if (uvec(3) < 0.) ipsi = npsi
    phi = 0.
  else
    phi = atan2(uvec(2),uvec(1))
  endif
  dphi = 2*pi/real(nphi(ipsi))
  iphi = mod(int((phi/dphi)+0.5) + nphi(ipsi),nphi(ipsi))
  getimask = psiposit(ipsi) + iphi
  end function getimask
  !!------------------------------------------------------------
  subroutine drawsurface(iunit,cen,r,iat,amask,bb,ch,iskip)
  implicit none
  integer,intent(in) :: iunit,iat,iskip
  real,dimension(3),intent(in) :: cen
  real,intent(in) :: r,bb
  real :: d
  character(len=1),intent(in) :: ch
  integer,dimension(MASKSIZE),intent(in) :: amask
  integer :: i,j,ibyte,ibit
  real,dimension(3) :: avec,bvec
  real,parameter :: DCUT=6.0
  !
  j = 0
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
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
      endif
    endif
  enddo
  end subroutine drawsurface
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
  real function distmod(avec,bvec,L)
  !! project bvec into the cube os size L around avec
  !! and report that distance
  real,dimension(3),intent(in) :: avec,bvec
  real,intent(in) :: L
  real,dimension(3) :: cvec
  integer :: i,j
  real :: dd

  cvec = bvec - avec
  do i=1,3
    if (cvec(i) > (L/2.)) cvec(i) = cvec(i) - L
    if (cvec(i) <= -(L/2.)) cvec(i) = cvec(i) + L
  enddo
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
  subroutine drawsaddle(iunit,avec,bvec,smask,r1,iat,taumin,thet,bb,ch)
  !! output the saddle surface as dots
  !! Wed Jun 20 14:35:42 EDT 2001
  real,parameter :: spcng=0.35   !! angstrom spacing of saddle dots
  character(len=1),intent(in) :: ch
  integer,dimension(MASKSIZE),intent(in) :: smask
  real,intent(in) :: taumin,r1,bb,thet,avec(3),bvec(3)
  integer,intent(in) :: iat,iunit
  integer :: i,j,itau,ntau,ibit,ibyte
  real,dimension(3) :: v1,v2,v3,vt,vw,vab,vec
  real :: chi,dtau,x,mat(3,3)
  !!
  ntau = nint(rw*(thet - taumin)/spcng)   ! number of steps in tau
  dtau = (thet - taumin)/real(ntau)    ! step size in tau
  if (periodic) then                   ! get nearest copy of j
    call boxaround(avec,bvec,lbox,vab) ! vab = a->b
  else
    vab = bvec - avec
  endif
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
    if (btest(smask(ibyte),ibit)) then
      j = j + 1
      vw = (r1+rw)*mxyz(1:3,i) ! direction of water (rw, mxyz are global)
      v1 = avec + vw           ! location of water
      call cros(vab,vw,v2)     ! v2 = perpendicular to abw
      call cros(vab,v2,vt)     ! perpendicular to vab in abw plane
      x = sqrt(dotprod(vt,vt))  ! length of vt
      vt = (rw/x)*vt           ! normalized vt, length rw
      v2 = v1 - v2             ! point below water
      do itau=0,ntau           ! steps separated by spcng
        chi = taumin + itau*dtau
        call getrotS(v1,v2,chi,mat,vec)  ! matrix and vector for rotation
        v3 = v1 + vt                     ! starting position
        call move_S(v3,mat,vec)          ! rotation
        if (periodic) call inbox(v3,lbox)
        write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
        i,ch,iat,v3(1:3),0.,bb
      enddo
    endif
  enddo
  end subroutine drawsaddle
  !!==========================================================
  subroutine tetrahedron(veci,vecj,veck,wxyz,d,iat,jat,kat)
  !! Using three atoms and a set of three distances, 
  !! get the apex of a tetrahedron based on the three atoms
  !! and having the three distances to the apex. The tetrahedron
  !! is assumed to be R-handed. To change handedness, switch two atoms.
  !! Mon Aug 13 10:58:06 EST 2001
  integer,intent(in) :: iat,jat,kat
  real,dimension(3),intent(in) :: veci,vecj,veck
  real,dimension(3),intent(out) :: wxyz
  real,intent(out) :: d
  real,dimension(3) :: vec,wvec
  real,dimension(3,3) :: mat
  real :: area,dij,djk,dik,diw,djw,dkw,ajik,ajiw,akiw,x,y,yp,z
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
  call heronsrule(dij,dik,djk,area,ajik)
  call heronsrule(dij,diw,djw,area,ajiw)
  call heronsrule(dik,diw,dkw,area,akiw)
  !! 
  x = diw*cos(ajiw)
  yp = diw*cos(akiw)
  y = yp/sin(ajik) - x/tan(ajik)
  d = sqrt(diw*diw - x*x - y*y)
  wxyz = (/x,y,d/)
  !! write(*,*) "In tetrahedron: wxyz=",x,y,d,wxyz
  !! Get frame
  call getframe(veci,vecj,veck,mat,vec)
  call move_S(wxyz,mat,vec)
  
  !! diagnostic
  !! write(*,*) "In tetrahedron: angs=",ajik,ajiw,akiw
  !! write(*,*) "In tetrahedron: wxyz=",wxyz
  end subroutine tetrahedron
  end module masker2
