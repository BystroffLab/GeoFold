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
