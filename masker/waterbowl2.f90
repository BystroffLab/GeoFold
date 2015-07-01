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
