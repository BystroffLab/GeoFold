program pdbmask
!! Special use version: trajectory included
!! -------------------------------------------------------------
!!
!! Tue Apr  3 16:22:13 EDT 2001
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
implicit none
integer,parameter :: MAXATOM=4096,MAXPSI=360
character(len=80) :: aline,pdbfile,binfile,logfile,outfile,template
character(len=80) :: edgefile
integer :: i,j,k,L,m,jarg,ios,iargc,natm
integer,parameter :: MASKSIZE=MAXATOM/32
real,parameter :: DTHETA=2.5   !! five degree increments
integer,parameter :: MMASK=6600,NTHETA=36  !! 0-90 deg
integer,dimension(MASKSIZE,NTHETA,MMASK) :: masklib
integer,dimension(MASKSIZE) :: amask,bmask,cmask,dmask,emask
integer :: nn,ipsi,npsi,itmp,nmask,itheta,iphi,imask,ibyte,ibit,nm,nv
integer :: imask1,itheta1,imask2,itheta2,iat,nat,jat,kat,ishow,ires
integer,dimension(:),allocatable :: psiposit,nphi
integer,dimension(100) :: saveimask, saveitheta,saveiat
real,dimension(NTHETA) :: collarsize,tau
real,dimension(MMASK) :: maskpsi,maskphi
real,dimension(3,MAXATOM) :: mxyz
real,dimension(3,80000) :: xyz
real :: dpsi,dphi,sasa,sphere,ssasa,wsphere,bsasa,col
real :: phi,psi,theta,x,y,d,dd,maxd,r1,r2,rw,taumin,zz
real,dimension(3) :: avec,bvec,cvec,wxyz,showvec
real,parameter :: pi=3.1415927,radius=10.
real,parameter :: rad=pi/180.
real :: dotprod
!! INITIALIZE
i = 0
jarg = iargc()
dpsi = DTHETA*rad
rw = 1.3
r1 = 2.0 
r2 = 2.0
maxd = r1 + r2 + 2*rw
sphere = 4*pi*r1**2
wsphere = 4*pi*rw**2
sasa = 0.0    !! total solvent accessible molecular surface area
ssasa = 0.0   !! total "saddle" convex/concave surface area
bsasa = 0.0   !! total "bowl" concave/concave surface area
!!
outfile = "junk.pdb"
if (jarg < 3) then
  write(*,*) 'Usage: xpdbmask_traj masktemplate masklib mask.log output '
  write(*,*) 'MASK LIBRARY STEP 4'
  write(*,*) 'This program takes a set of binary masks and '
  write(*,*) 'calculates the molecular surface fo the input coordinates.'
  write(*,*) 
  write(*,*) 'masktemplate is a set of points at radius = ',radius
  write(*,*) 'binfile is a library of masks from xbinarymask'
  write(*,*) 'mask.log is the logfile from xbinarymask'
  write(*,*) 'inputpdb is a pdb file of masked SASA points'
  write(*,*) 'output is a pdb file of masked edge (collar) points'
  write(*,*) 'MAXATOM=',MAXATOM
  write(*,*) 'DTHETA=',DTHETA
  write(*,*) 'MMASK=',MMASK
  stop 'pdbmask.f90 v.18-MAY-01'
endif
call getarg(1,template)
call getarg(2,binfile)
call getarg(3,logfile)
if (jarg >= 4) then
  call getarg(4,outfile)
endif
!!call getarg(4,pdbfile)
!!call getarg(5,outfile)
ishow = 0
if (jarg > 5) then
  call getarg(6,aline)
  read(aline,*,iostat=ios) ishow
  !! ishow is a residue to draw the surface around
  write(*,'("Show surface around residue ",i4," only.")') ishow
endif
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
    allocate(psiposit(0:npsi),nphi(0:npsi),stat=ios)
    if (ios/=0) stop 'Problem allocating psiposit, nphi'
  endif
  if (aline(1:4)=="PSI=".and.npsi/=0) then
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
!! read binary format mask library, created by xbinarymask
open(12,file=binfile,status='old',form='unformatted')
  read(12) masklib(1:MASKSIZE,1:NTHETA,1:NMASK)
close(12)
!!diagnostic
!!write(*,*)  masklib(1:MASKSIZE,10,50)
!! get psi and phi for each mask
do ipsi=0,npsi
  psi = ipsi*dpsi
  dphi = 2*pi/real(nphi(ipsi))
  do iphi=1,nphi(ipsi)
    phi = (iphi-1)*dphi
    imask = psiposit(ipsi) + iphi - 1
    maskpsi(imask) = psi
    maskphi(imask) = phi
    !! diagnostic
    !! write(*,*) imask,psi,phi
  enddo
enddo
!!
!! Calculate collarsize and tau for each theta
!! collarsize is the number of bits in the edge mask for angle theta
!! tau is the angular width of the half-saddle for angle theta and r1
ipsi = 1
iphi = 1
imask = psiposit(ipsi) + iphi - 1
do itheta=1,NTHETA-1
  theta = itheta*dpsi
  !! If the depth of the saddle is greater than one water radius
  !! then two saddles connect, and the integration limits are zero to 90-theta.
  !! Othewise, calculate the lower integration limit (taumin).
  if ((r1+rw)*sin(theta) >= rw) then
    tau(itheta) = (pi/2.) - theta
  else
    x = acos((r1+rw)*sin(theta)/rw)
    tau(itheta) = (pi/2.) - theta - x
  endif
  amask = iand(masklib(:,itheta,imask),not(masklib(:,itheta+1,imask)))
  j = 0
  do i=1,MAXATOM
    ibyte = (i-1)/32 + 1
    ibit = mod(i-1,32)
    if (btest(amask(ibyte),ibit)) then
      j = j + 1
    endif
  enddo
  collarsize(itheta) = j    !! approximately the same for all imask
  !! diagnostic
  write(*,'("theta tau size ",2f8.2,f8.0)') &
    theta/rad,tau(itheta)/rad,collarsize(itheta)
enddo
!! NTHETA is 90deg, only possible if dd ~= 0.
tau(NTHETA) = 0.0
!! this number must be non-zero, but doesn't matter since tau=0.
collarsize(NTHETA) = 1  
!!
open(11,file=template,status='old',form='formatted')
i = 0
do 
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5) /= "ATOM ") cycle
  i = i + 1
  if (i > MAXATOM) stop 'Too many points in template file'
  read(aline(31:54),'(3f8.3)',iostat=ios) mxyz(1:3,i)
  mxyz(1:3,i) = mxyz(1:3,i)/10.
enddo
close(11)
!! 
open(13,file=outfile,status='replace',form='formatted')
!! open(14,file=pdbfile,status='old',form='formatted')
iat = 0
showvec = 0.
nat = 2
xyz(1:3,1:nat) = 0.0
!do
!  read(14,'(a)',iostat=ios) aline
!  if (ios/=0) exit
!  if (aline(1:5)=="ATOM ") then
!    iat = iat + 1
!    aline(22:22) = " "
!    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,iat)
!    read(aline(23:26),*,iostat=ios) ires
!    write(13,'(a)') trim(aline)
!    if ((ires==ishow).and.(aline(13:16)==" CA ")) showvec = xyz(1:3,iat)
!  elseif (aline(1:4)=="TER ") then
!    exit
!  else
!    write(13,'(a)') trim(aline)
!    cycle
!  endif
!enddo
!nat = iat
!write(*,'(i9," atoms read.")') nat
  
zz = 0.2
DO WHILE (zz < maxd)
sasa = 0; ssasa = 0; bsasa = 0.
zz = zz + 0.1
xyz(3,2) = zz
!! Loop over all atoms, get sasa, saddle and bowl surfaces
do iat=1,nat
  nm = 0
  amask = -1                 ! init all bits = 1
  bmask = -1
  do jat=1,nat
    if (iat == jat) cycle
    !! For every atom pair, calculate phi, psi, theta
    !! with iat as the central atom
    avec = xyz(1:3,jat) - xyz(1:3,iat)
    dd = sqrt(dotprod(avec,avec))
    if (dd > maxd) cycle    !! atom out of range
    avec = avec/dd
    imask = getimask(avec)  !! avec must be a unit vector
    !   psi = acos(avec(3))
    !   phi = atan2(avec(2),avec(1))
    theta = gettheta(r1+rw,r2+rw,dd)   !! Heron's rule
    if (theta <= 0.) cycle
    itheta = nint((theta/rad)/DTHETA)
    !! diagnostic
    !! write(*,'(2i4,4f6.1,4i5)') iat,jat,psi/rad,phi/rad,theta/rad,dd,ipsi,iphi,itheta,imask
    if (itheta <= 0) cycle    !! atom out of range (borderline)
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
  x = sphere*(real(j)/real(MAXATOM))
  sasa = sasa + x
  call drawsurface(13,xyz(1:3,iat),r1,iat,amask,50.,"C",8)  !! output the sasa surface
  i = 0
  cmask = 0
  nn = 0
  do k=1,nm
    kat = saveiat(k)
    imask = saveimask(k)
    itheta = saveitheta(k)
    bmask = iand(not(masklib(:,itheta+1,imask)),amask)
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
    theta = itheta*dpsi
    y = (r1+rw)*sin(theta)
    col = 2*pi*(real(j)/collarsize(itheta))
    x = col*((0.5*pi-theta)*y - rw*cos(theta))
    if (y < rw) then
      taumin = acos(y/rw)
      x = x + col*(rw*sin(taumin) - taumin*y)
    endif
    ssasa = ssasa + x
    if (any(bmask /= 0)) then
      !! keep a reduced list of close atoms.
      nn = nn + 1    !!  note: nn <= k, so this is legal.
      saveiat(nn) = kat
      saveimask(nn) = imask
      saveitheta(nn) = itheta
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
  nv = 0
  dmask = 0
  emask = 0
  do k=1,nn
    kat = saveiat(k)
    if (kat <= iat) cycle  !! calculate bowl iff iat is the first atom
    imask1 = saveimask(k)
    itheta1 = saveitheta(k)
    do m=1,nn
      jat = saveiat(m)
      if (jat <= kat) cycle  !! calculate bowl iff kat is the second atom
      imask2 = saveimask(m)
      itheta2 = saveitheta(m)
      !! diagnostic
      !! write(*,*) imask1,itheta1,imask2,itheta2
      bmask = iand(not(masklib(:,itheta1+1,imask1)),masklib(:,itheta1,imask1) )
      bmask = iand(bmask,not(masklib(:,itheta2+1,imask2)) )
      bmask = iand(bmask,masklib(:,itheta2,imask2) )
      bmask = iand(cmask,bmask)  !! bmask is current verteces
      !! dmask = ior(dmask,bmask)  !! dmask is all verteces (used only for display)
      if (any(bmask /= 0)) then !! if there is at least one vertex
        !! Retreive unit vector to atom k: avec
        psi = maskpsi(imask1)
        phi = maskphi(imask1)
        avec(3) = cos(psi)
        avec(2) = sin(psi)*sin(phi)
        avec(1) = sin(psi)*cos(phi)
        !! Retreive unit vector to atom m: bvec
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
        nv = 0             !! if no vertex, nv = 0
        x = 0.
        !! First check one side for verteces...
        emask = iand(masklib(:,itheta,imask),bmask)
        if (any(emask /= 0)) then
          nv = 1       !! if only R-hand surface, nv=1
          !! calculate surface, returning surface mask and water coords
          call waterbowl(iat,kat,jat,x,wxyz) !! also uses/returns emask, uses cvec
          !! diagnostic, write water coords
          !! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
          !!   iat,"W",iat,wxyz(1:3),0.,10.0
          call drawsurface(13,wxyz(1:3),rw,iat,emask,-20.,"E",8) !! output the surface
        endif
        !! ...then check the other
        emask = iand(not(masklib(:,itheta,imask)),bmask)
        if (any(emask /= 0)) then
          !! if (nv == 1) then
          if (.false.) then  !! do this to force drawing surface
            x = 2*x
          else
            call waterbowl(iat,kat,jat,y,wxyz) !! also uses/returns emask
            x = x + y
            !! diagnostic, write water coords
            !! write(13,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
            !! iat,"W",iat,wxyz(1:3),0.,20.0
            call drawsurface(13,wxyz(1:3),rw,iat,emask,20.,"E",8) !! output the surface
          endif
          nv = 2  !! useful flag?
        endif
        bsasa = bsasa + x
      endif
    enddo
  enddo
enddo
x = sasa + ssasa + bsasa
!write(*,'("Total VDW area = ",f9.2)') sasa
!write(*,'("Total saddle area  = ",f9.2)') ssasa
!write(*,'("Total bowl area    = ",f9.2)') bsasa
!write(*,'("Total surface area = ",f9.2)') x
write(*,*) zz,sasa,ssasa,bsasa,x
!
ENDDO

close(13)
!!close(14)
CONTAINS
!!----------------------------------------------------------------!!
real function gettheta(r1,r2,dd)
!! returns theta in radians for angle between r1 and dd
!! for triangle w/sides r1,r2,dd
!! using Heron's Rule
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
gettheta = asin(x/r1)
end function gettheta
!!----------------------------------------------------------------!!
subroutine waterbowl(iat,jat,kat,x,wxyz) 
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
integer,intent(in) :: iat,jat,kat
real,intent(out) :: x, wxyz(3)
integer :: i,j,k,ibyte,ibit,iphi,ipsi,itheta,jmask
real :: y,theta
real,dimension(3) :: avec,bvec,dvec
!! ---------
!! uses current value of cvec, which is the cross product
!! iat->jat X iat->kat.  If cvec.(iat->wat) > 0, it's R-handed
!! So the order is iat,jat,kat. Otherwise, it's L-handed
!! and the order is iat,kat,jat
j = 0
do i=1,MAXATOM
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(emask(ibyte),ibit)) then
    j = i
    exit
  endif
enddo
if (j==0) then
  wxyz = 0.
  x = 0.
  return
endif
!! wxyz is water position in absolute coords (A)
wxyz = (r1+rw)*mxyz(:,j)
d = dotprod(cvec,wxyz)
!! diagnostic
!! write(*,'(3i5,f9.3)') iat,jat,kat,d
if (d > 0.) then !! R-handed. Use iat,jat,kat and -cvec
  if (d < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
    theta = acos(d/rw)
    itheta = nint((theta/rad)/DTHETA)
    dvec = -1*cvec
    jmask = getimask(dvec)
    emask = masklib(:,itheta,jmask)
  else
    emask = -1
  endif
  !! mask ijw plane
  avec = -wxyz
  bvec = xyz(:,jat) - (wxyz + xyz(:,iat))
  call cros(avec,bvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  !! write(*,'(5x,3f7.3,i7)') dvec(1:3),jmask
  emask = iand(masklib(:,NTHETA,jmask),emask)
  !! mask jkw plane
  avec = bvec
  bvec = xyz(:,kat) - (wxyz + xyz(:,iat))
  call cros(avec,bvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  emask = iand(masklib(:,NTHETA,jmask),emask)
  !! mask kiw plane
  avec = bvec
  bvec = -wxyz
  call cros(avec,bvec,dvec)
  y = sqrt(dotprod(dvec,dvec))
  dvec = dvec/y
  jmask = getimask(dvec)
  emask = iand(masklib(:,NTHETA,jmask),emask)
else !! L-handed. Use iat,kat,jat and +cvec (current imask)
  if (-d < rw) then  !! there's a hole in the bowl, get mask (don't reset global imask)
    theta = acos(-d/rw)
    itheta = nint((theta/rad)/DTHETA)
    jmask = imask  !!  NOTE: using global variable imask. DEPENDS ON CALLING SEQUENCE!!
    emask = masklib(:,itheta,jmask)
  else
    emask = -1
  endif
  !! mask ikw plane
  avec = -wxyz
  bvec = xyz(:,kat) - (wxyz + xyz(:,iat))
  call cros(avec,bvec,dvec)
  dvec = dvec/sqrt(dotprod(dvec,dvec))
  jmask = getimask(dvec)
  emask = iand(masklib(:,NTHETA,jmask),emask)
  !! mask kjw plane
  avec = bvec
  bvec = xyz(:,jat) - (wxyz + xyz(:,iat))
  call cros(avec,bvec,dvec)
  dvec = dvec/sqrt(dotprod(dvec,dvec))
  jmask = getimask(dvec)
  emask = iand(masklib(:,NTHETA,jmask),emask)
  !! mask jiw plane
  avec = bvec
  bvec = -wxyz
  call cros(avec,bvec,dvec)
  dvec = dvec/sqrt(dotprod(dvec,dvec))
  jmask = getimask(dvec)
  !! write(*,'(5x,3f7.3,i7)') dvec(1:3),jmask
  emask = iand(masklib(:,NTHETA,jmask),emask)
endif
j = 0
do i=1,MAXATOM
  ibyte = (i-1)/32 + 1
  ibit = mod(i-1,32)
  if (btest(emask(ibyte),ibit)) then
    j = j + 1
  endif
enddo
x = wsphere*real(j)/real(MAXATOM)
wxyz = xyz(1:3,iat) + wxyz
end subroutine waterbowl
!!------------------------------------------------------------
integer function getimask(uvec)
!! return the mask index for a unit vector
!! uses pi,npsi,dpsi,psiposit
!!   from the main program
real,intent(in) :: uvec(3)
integer :: iphi,ipsi,imask
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
!!
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
      write(iunit,'("HETATM",i5,"  O   HOH ",a1,i4,"    ",3f8.3,2f6.1)') &
      i,ch,iat,avec(1:3),0.,bb
    endif
  endif
enddo
end subroutine drawsurface
!!==========================================================
end program pdbmask
