!!============= included in masker.f90 ==================!
  subroutine getms_d(xyz,atype,nat,ms,peratom,msderiv)
  integer, intent(in) :: nat
  real,dimension(3,nat),intent(in) :: xyz
  real,dimension(3,nat),intent(out) :: peratom !! 3 types of surface
  real,dimension(nat,nat),intent(out) :: msderiv
  integer,dimension(nat),intent(in) :: atype
  real,intent(out) :: ms
  !
  integer :: i,j,jj,k,L,m,jarg,ios,iargc,natm
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
    
  !! CONTACT SURFACE  (VDW,  SAS)  
  wsphere=rw*rw*4*pi
  !! Loop over all atoms, get sasa, saddle and bowl surfaces
  do iat=1,nat
    nm = 0
    amask = -1                 ! init all bits = 1
    bmask = -1                 ! init all bits = 1
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
      !! dd = sqrt(dotprod(avec,avec))
      r2 = atomlib(atype(jat))%r
      maxd = r1 + r2 + 2*rw
      if (farapart(avec,maxd,dd)) cycle
      !! if (dd > maxd) cycle    !! atom out of range
      avec = avec/dd
      !! returns imask and theta. If theta > pi/2, use 2*NTHETA-itheta and .not. it
      call getimask2(r1+rw,dd,r2+rw,theta,avec,imask,itheta)
      if (imask==0) cycle   !! bad triangle (bug?)
      if (imask==-1) then   !! i is completely embedded in j
        amask = 0
        exit
      endif
      if (theta < pi/2) then
        if (itheta==0) itheta = 1
        if (itheta==NTHETA) then
          amask = iand(not(masklib(:,itheta,imask)),amask)     
        else
          amask = iand(masklib(:,itheta,imask),amask)     
        endif
      else  !! jat is embedded in iat
        !!invert mask to get theta > pi/2, embedded mask.
        amask = iand(not(masklib(:,itheta,imask)),amask)     
      endif
      if (any(bmask /= amask) ) then
        !! if any changes occurred in the mask...
        bmask = amask
        !! keep a list of close atoms, for the next part
        !! sort the list ny theta, descending. Lowest theta first
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
      endif
    enddo
    j = countbits(amask)
    if (j < 1) cycle   ! ignore saddle, etc if nothing is exposed.
    !
    !! CONTACT SURFACE
    sphere = 4*pi*atomlib(atype(iat))%r**2
    x = sphere*(real(j)/real(MAXATOM))
    sasa = sasa + x    !! x is the sasa before edge correction
    sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
    if (drawing) then
      cmask = iand(amask,vmask)
      call drawsurface(13,xyz(1:3,iat),r1,iat,cmask,50.,"V",1)  !! output the sasa surface
    endif
    !
    !! TOROIDAL SURFACE  (SADDLE)
    !
    i = 0
    cmask = 0
    emask = 0
    dmask = amask
    nn = 0
    do k=1,nm
      theta = savetheta(k)
      kat = saveiat(k)
      imask = saveimask(k)
      itheta = saveitheta(k)
      ! dd = savedij(k)
      r2 = atomlib(atype(kat))%r
      if (theta >= pi/2) then  !! embedded, no saddle
        !!write(*,*) "embedded itheta=", itheta
        if (itheta <= 1) then
          bmask = amask
        else
          bmask = iand(masklib(:,itheta-1,imask),amask) !! this arc is for derivs and drawing
        endif
      else
        if (itheta == NTHETA) then  !! imask is for inverse mask
          bmask = iand(masklib(:,itheta-1,imask),amask)
          emask = iand(masklib(:,itheta-1,imask),dmask)
          amask = iand(amask,not(bmask))
        else
          bmask = iand(not(masklib(:,itheta+1,imask)),amask)
          emask = iand(not(masklib(:,itheta+1,imask)),dmask)
          amask = iand(amask,not(bmask))
        endif
      endif
      if (drawing) then
        cmask = ior(cmask,bmask)  !! collect all edges
        call drawsurface(13,xyz(1:3,iat),r1,kat,bmask,50.,"E",1)  !! output the edges 1 at a time
      endif
      !! Count the bits in the exposed arc
      j = countbits(bmask)
      if (j < 1) cycle
      jj = countbits(emask)
      j = nint(real(j+jj)/2.)
      arcfrac = real(j)/real(collarsize(itheta,imask))
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
        !! CONTACT derivative is the slope of the contact surface with theta
        !!         times the exposed arc, over the distance. dS/dD = (dS/dQ)/(dD/dQ)
        x = 
        !! SASA correction
        x = (theta - interc(itheta))*slope(itheta)*arcfrac*sphere
        sasa = sasa + x
        sasanrg = sasanrg + x*atomlib(atype(iat))%w1  !! kJ/mol
        !! saddle
        y = (r1+rw)*sin(theta)
        col = 2*pi*(arcfrac)   !! this is the length of the exposed arc in radians
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
        ssasanrg = ssasanrg + x*atomlib(atype(iat))%w2
        if (drawing) then
          x = (pi/2) - theta  !! taumax
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
      endif
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
          call boxaround(xyz(1:3,kat),xyz(1:3,jat),lbox,avec)
        else
          avec = xyz(1:3,jat) - xyz(1:3,kat)
        endif
        r3 = atomlib(atype(jat))%r
        maxd = r3 + r2 + 2*rw
        if (farapart(avec,maxd,dd)) cycle
        !! if (dd > maxd) cycle    !! third side too long 
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
      enddo
    enddo
  enddo

  !! RE-ENTRANT SURFACE accounting for intersecting re-entrant surfaces.
  !! 
  !! This segment of the code (obselete given correction of bugs in tetrahedron?)
  !! was written to ensure that no waters sat too close to any atom.
  !! 
  xx = 4*rw*rw
  water(1:nwat)%good = .true.
  do iwat = 1,nwat
    if (.not.water(iwat)%good) cycle
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
        amask = iand(masklib(:,itheta,imask),amask)     
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
    bsasanrg = bsasanrg + (x/3.)*(atomlib(atype(water(iwat)%trng(1)))%w3 + &
        atomlib(atype(water(iwat)%trng(2)))%w3 + atomlib(atype(water(iwat)%trng(3)))%w3)
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
  end subroutine getms_d ! (xyz,atype,nat,ms)
