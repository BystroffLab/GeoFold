module karpluslazaridis
  use masker

  private
  !Parameters according to (Lazaridis & Karplus, 1999)
  character (len=4),dimension(17),parameter :: atomtypes = (/"C   ","CR  ",&
  "CH1E","CH2E","CH3E","CR1E","NH1 ","NR  ","NH2 ","NH3 ","NC2 ","N   ",&
  "OH1 ","O   ","OC  ","S   ","SH1E"/)
  !angstroms cubed
  real,dimension(17),parameter :: volumes = (/14.7,8.3,23.7,22.4,30.0,4.4,4.4,&
    11.2,11.2,11.2,0.0,10.8,10.8,10.8,14.7,21.4/)
  !kcal/mol --> need J or kJ?  It looks like we want J
  real,dimension(17),parameter :: Grefs = (/0.000,-3723.760,-782.408,1556.448,&
    4556.376,238.488,-24894.800,-15982.880,-22802.800,-83680.000,-41840.000,&
    -4184.000,-24769.280,-22300.720,-41840.000,-13556.160,-8577.200/) ! J/mol
  real,dimension(17),parameter :: Gfrees = (/0.00,-5857.60,-1046.00,2175.68,&
    6276.00,334.72,-37237.60,-16736.00,-32635.20,-83680.00,-41840.00,-6485.20,&
    -28032.80,-24476.40,-41840.00,-17154.40,-11296.80/) !J/mol
  real,dimension(17),parameter :: Hrefs = (/0.000,9288.480,3665.184,-2552.240,&
    -7443.336,-4071.032,-37902.856,-19472.336,-37773.152,-104600.000,-50208.000,&
    -5230.000,-38760.576,-24212.808,-50208.000,-18723.400,-18723.400/) !J/mol
  !angstroms
  real,dimension(17),parameter :: CorrLengths = (/3.5,3.5,3.5,3.5,3.5,3.5,3.5,&
    3.5,3.5,6.0,6.0,3.5,3.5,3.5,6.0,3.5,3.5/)
  real,parameter :: pi=3.1415927410125732421875
  real,dimension(17),parameter :: Vdws = (/2.1,2.1,2.1,2.1,2.1,2.1,1.6,1.6,1.6,&
    1.6,1.6,1.6,1.6,1.6,1.6,1.89,1.89/)

  public: getKplType, kpl_gsolv

contains
  integer function getKplType(aline) result(type)
    implicit none
    character(len=*),intent(in) :: aline
    character(len=3) :: res
    character(len=4) :: atom
    integer,dimension(20,14) :: types
    character(len=3),dimension(20) :: residues
    character(len=4),dimension(20,14) :: atomnames
    integer :: i,j

    !initialize residues
    residues = (/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',&
                 'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'/)

    !initialize types, There are 17 possibilities for the KPL model
    types = 0
    types(1,1:5) = (/7,3,1,14,5/)
    types(2,1:6) = (/7,3,1,14,4,17/)
    types(3,1:8) = (/7,3,1,14,4,1,15,15/)
    types(4,1:9) = (/7,3,1,14,4,4,1,15,15/)
    types(5,1:11) = (/7,3,1,14,4,2,6,6,6,6,6/)
    types(6,1:4) = (/7,3,1,14/)
    types(7,1:10) = (/7,3,1,14,4,2,7,6,6,7/)
    types(8,1:8) = (/7,3,1,14,3,4,5,5/)
    types(9,1:9) = (/7,3,1,14,4,4,4,4,10/)
    types(10,1:8) = (/7,3,1,14,4,3,5,5/)
    types(11,1:8) = (/7,3,1,14,4,4,16,5/)
    types(12,1:8) = (/7,3,1,14,4,1,14,9/)
    types(13,1:7) = (/12,3,1,14,4,4,4/)
    types(14,1:9) = (/7,3,1,14,4,4,1,14,9/)
    types(15,1:11) = (/7,3,1,14,4,4,4,7,1,11,11/)
    types(16,1:6) = (/7,3,1,14,4,13/)
    types(17,1:7) = (/7,3,1,14,4,13,5/)
    types(18,1:7) = (/7,3,1,14,3,5,5/)
    types(19,1:14) = (/7,3,1,14,4,2,6,2,7,2,6,6,6,6/)
    types(20,1:12) = (/7,3,1,14,4,2,6,6,6,6,2,13/)

    !initialize atomnames
    atomnames = "    "
    atomnames(1,1:5) = (/'N   ','CA  ','C   ','O   ','CB  '/)
    atomnames(2,1:6) = (/'N   ','CA  ','C   ','O   ', 'CB  ','SG  '/)
    atomnames(3,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','OD1 ','OD2 '/)
    atomnames(4,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                         'OE1 ','OE2 '/)
    atomnames(5,1:11) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                          'CD2 ','CE1 ','CE2 ','CZ  '/)
    atomnames(6,1:4) = (/'N   ','CA  ','C   ','O   '/)
    atomnames(7,1:10) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','ND1 ',&
                          'CD2 ','CE1 ','NE2 '/)
    atomnames(8,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG1 ','CG2 ','CD1 '/)
    atomnames(9,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                         'CE  ','NZ  '/)
    atomnames(10,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ','CD2 '/)
    atomnames(11,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','SD  ','CE  '/)
    atomnames(12,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','OD1 ','ND2 '/)
    atomnames(13,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  '/)
    atomnames(14,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                          'OE1 ','NE2 '/)
    atomnames(15,1:11) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                           'NE  ','CZ  ','NH1 ','NH2 '/)
    atomnames(16,1:6) = (/'N   ','CA  ','C   ','O   ', 'CB  ','OG  '/)
    atomnames(17,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','OG1 ','CG2 '/)
    atomnames(18,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG1 ','CG2 '/)
    atomnames(19,1:14) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                           'CD2 ','NE1 ','CE2 ','CE3 ','CZ2 ','CZ3 ','CH2 '/)
    atomnames(20,1:12) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                           'CD2 ','CE1 ','CE2 ','CZ  ','OH  '/)

   read(aline(13:16),*) atom
   read(aline(18:20),*) res
   do i = 1, 21
     if(i == 21) then !Unrecognized amino acid
       type = 0
       exit
     endif
     if(residues(i) == res) exit
   enddo
   do j = 1, 15
     if(j == 15 .or. i == 21) then !Unrecognized atomname
       type = 0
       exit
     endif
     if(atomnames(i,j) == atom) exit
   enddo
   if(i /= 21 .and. j /= 15) type = types(i,j)
  end function getKplType

  real function get_dist(atom1,atom2) result(dist)
    implicit none
    real,dimension(3),intent(in) :: atom1,atom2
    real,dimension(3) :: temp

    temp = atom1-atom2
    temp = temp*temp
    dist = sqrt(temp(1)+temp(2)+temp(3))
  end function get_dist

  real function kpl_gsolv (xyz,nat,start) result(gsolv)
    ! Calculates the solvation energy according to the Karplus-Lazaridis solvation
    ! model.  It follows a really annoying equation that I don't want to put here...
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat),intent(in) :: xyz
    integer, intent(in) :: start
    integer :: iatom,jatom
    real :: rij
    real :: gref,gfree,corrL,vdw,Vj,alpha,xi,dist

    gsolv = 0.

    do iatom = 1,nat
      gsolv = gsolv + Grefs(atype(iatom))
      do jatom = 1, nat
        if(iatom == jatom) cycle
        if(atype(iatom) <= 0 .or. atype(jatom) <= 0) cycle
        dist = get_dist(xyz(:,iatom),xyz(:,jatom))
        alpha = (2*Gfrees(atype(iatom)))/(sqrt(pi*CorrLengths(atype(iatom))))
        xi = (dist-Vdws(atype(iatom)))/CorrLengths(atype(iatom))
        gsolv = gsolv - ((alpha*exp(-2*xi))/(4*pi*dist*dist))*volumes(atype(jatom))
  end function kpl_gsolv

end module karpluslazaridis
