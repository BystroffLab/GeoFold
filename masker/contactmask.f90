program contactmask
!! Formerly unfold_energy
!! CB  Tue Apr 15 19:48:03 EDT 2008
use masker
implicit none
!!
! This program reads a PDB file and calculates the energy
! between each pair of residues in the sequence.
! The energy function is a linear combination of an Hbond
! term (Hbond) and a solvation term (SAS) .
!! 
character(len=200) :: xyzfile, aline,outfile
real,dimension(:,:),allocatable :: xyz,xyznew
integer,dimension(:),allocatable :: saveatype
integer :: i,j,k,nseed,jj,n,nn,isize,jsize,ios=0
real :: tt, sig=0.3,x,y,z,nrg,scl,xmin,xmax,ymin,ymax,zmin,zmax
real :: pairwise,gsolv,sasa
real,dimension(:),allocatable :: residuesasa,buried,residueGsolv
real,dimension(3) :: vec
integer :: iargc, jarg, time, ounit=22, nhet, nat, nres, ires, jres, nunk
logical :: drawingatall = .false., verbose=.false.
character(len=11) :: lastresstring
character :: altloc
integer,dimension(:),allocatable :: nposit
! integer, dimension(:), pointer::atype  ! public
!! set default values and constants --------------------------------------------
xmin = 999.; ymin = 999.; zmin = 999.
xmax = -999.; ymax = -999.; zmax = -999.
nunk=0
!! rendering = .false.
!! smoothing = .false.
!! saying = .false.
lastresstring="           "
!! Read command line ----------------------------------
nseed = time()
jarg = iargc()
!! lbox = 0.00  !! masker global
if (jarg < 1) then
  write(*,'("Usage: xcontactmask PDBfile [output proberadius]")')
  stop 'contactmask.f90 v.  Sat Jan 21 10:19:42 EST 2006'
endif
call getarg(1,xyzfile)
if (jarg >=2 ) then
  call getarg(2,outfile)
  ounit = 22
  open(ounit,file=outfile,status='replace',form='formatted',iostat=ios)
  if (ios/=0) stop 'contactmask:: error attempting to open file. wrong permissions?'
else
  ounit = 6
endif
if (jarg >=3 ) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) rw
  if (ios/=0) stop 'Bad value for rw (Probe radius)'
  write(*,'("Probe radius: ",f8.2)') rw
endif
!! Initialize masks
call masker_plot(drawingatall)
call masker_initmasks()
! call masker_getatypeptr(atype)
!!-------- read atom coordinates from file, just count them
open(12,file=xyzfile,status='old',form='formatted')
i = 0
j = 0
nres = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5)/='ATOM '.and.aline(1:6)/='HETATM') cycle
  if (aline(17:17)/=" ") then
    write(*,'("WARNING: alt loc indicator detected!!")')
    if (altloc==" ") altloc = aline(17:17)
    if (aline(17:17)==altloc) then
      aline(17:17)=" "
    endif
  else
    altloc = " "
  endif
  if (aline(17:17)/=" ") cycle
  if (aline(1:5)=='ATOM ') i = i + 1
  if (aline(1:6)=='HETATM') j = j + 1
  if (lastresstring/=aline(17:27)) then
    nres = nres + 1
    lastresstring = aline(17:27)
    ! write(*,'("Residue ",i6," ",a6," ",a11)') nres,aline(1:6),lastresstring
  endif
enddo
nat = i
nhet = j
!! diagnostic
write(*,*) 'Number of atoms:', nat
write(*,*) 'Number of hetatoms:', nhet
write(*,*) 'Number of residues:', nres
nat = nat + nhet
!!-------- allocate arrays
allocate(xyz(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
allocate(xyznew(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyznew'
allocate(saveatype(nat),stat=ios)
if (ios/=0) stop 'Error allocating saveatype'
allocate(nposit(nres+1),stat=ios)
if (ios/=0) stop 'Error allocating nposit'
allocate(residuesasa(nres),residueGsolv(nres),stat=ios)
if (ios/=0) stop 'Error allocating residuesasa'
allocate(buried(nat),stat=ios)
if (ios/=0) stop 'Error allocating buried'
!!--------  allocate atype() and other arrays within MASKER
call masker_allocatoms(nat)
write(*,*) 'Total atoms allocated:', nat
!!-------- rewind, read atom coordinates from file
rewind(12)
i = 0
ires = 0
vec = 0.
lastresstring = "           "
altloc = " "
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    if (aline(17:17)/=" ") then
      if (altloc==" ") altloc = aline(17:17)
      if (aline(17:17)==altloc) then
        aline(17:17)=" "
      endif
    else
      altloc = " "
    endif
    if (aline(17:17)/=" ") cycle
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'Format error xyzfile'
    if (lastresstring/=aline(17:27)) then
      ires = ires + 1
      nposit(ires) = i
      lastresstring = aline(17:27)
    endif
    jj = 0
    !!!This part may need to be modified for WZS
    if(wzs) then
      jj = getatype(aline)
    else
      do j=1,nattype
        if (aline(14:17)==atomlib(j)%name(1:4)) then  !! allocated in MASKER::initmasks
          jj = j
          exit
        endif
      enddo
    endif
    if (jj==0) then
      if (verbose) write(*,'("Unknown atom type: ",a4,"  this atom will be ignored!")') aline(14:17)
      saveatype(i) = -1   !! dont use this atom
      nunk = nunk+1
    else
      saveatype(i) = jj
    endif
    !!!
    !! diagnostic
    !! write(*,*) 'Atom ',i,jj,xyz(1:3,i)
    vec = vec + xyz(1:3,i)
  endif
enddo
nposit(nres+1) = i+1
if (any(nposit(1:nres+1)==0)) stop 'BUG BUG BUG: nposit == 0'
if (nunk>0) then
  write(*,'(i9," Unknown atom types. These atom will be ignored!",$)')  nunk
  write(*,'("(set contactmask.f90 to verbose to see atom names.)")')
endif
!!-------- mov atoms to the center of mass
vec = -vec/real(nat)
do i=1,nat
  xyz(1:3,i) = xyz(1:3,i) + vec
enddo
!!-------- get SASA for individual residues
do ires=1,nres
  isize = nposit(ires+1) - nposit(ires)
  xyznew(1:3,1:isize) = xyz(1:3,nposit(ires):nposit(ires+1)-1)
  atype(1:isize) = saveatype(nposit(ires):nposit(ires+1)-1)
  !write(0,'("Calculating ",i9,": size=",i5," ",$)') ires,isize
  !write(0,*) xyznew(1:3,1:isize)
  !write(0,*) atype(1:isize) 
  call masker_getms(xyznew,isize,x,buried=buried,SAS=.true.,ses=sasa)
  !! ignore x, uses global varable sasa
  residuesasa(ires) = sasa
  residueGsolv(ires) = x  !!! Need this for solvation modeling!
  !write(*,'(f9.3)') sasa
enddo
!!-------- get SASA for residue pairs and subtract from the sum of individual residue SASA.
do ires=1,nres-1
  isize = nposit(ires+1) - nposit(ires)
  xyznew(1:3,1:isize) = xyz(1:3,nposit(ires):nposit(ires+1)-1)
  atype(1:isize) = saveatype(nposit(ires):nposit(ires+1)-1)
  do jres=ires+1,nres
    jsize = nposit(jres+1) - nposit(jres)
    xyznew(1:3,isize+1:isize+jsize) = xyz(1:3,nposit(jres):nposit(jres)+jsize-1)
    atype(isize+1:isize+jsize) = saveatype(nposit(jres):nposit(jres)+jsize-1)
    jsize = jsize + isize
    call masker_getms(xyznew,jsize,x,buried=buried,SAS=.true.,ses=sasa)
    pairwise = residuesasa(ires) + residuesasa(jres) - sasa
    !! This is the energy of unfolding, not folding.
    gsolv = residueGsolv(ires) + residueGsolv(jres) - x
    if (abs(pairwise) > 0.01) then
      write(ounit,'(2i6,2f9.3)') ires,jres,pairwise,gsolv
    endif
  enddo
enddo
deallocate(xyz,stat=ios)
deallocate(xyznew,stat=ios)
deallocate(saveatype,stat=ios)
deallocate(nposit,stat=ios)
deallocate(residuesasa,stat=ios)
deallocate(buried,stat=ios)
call masker_deallo()
if (ounit/=6) close(ounit)
stop 'contactmask.f90  v.  Sat Jan 21 13:25:39 EST 2006'
CONTAINS

!For WZS modeling
integer function getatype(aline)
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

  !initialize types
  types = 0
  types(1,1:5) = (/21,8,11,27,13/) 
  types(2,1:6) = (/21,8,11,27,10,29/) 
  types(3,1:8) = (/21,8,11,27,13,1,23,23/)
  types(4,1:9) = (/21,8,11,27,13,13,1,23,23/)
  types(5,1:11) = (/21,8,11,27,13,15,15,15,15,15,15/)
  types(6,1:4) = (/21,8,11,27/)
  types(7,1:10) = (/21,8,11,27,13,2,17,17,2,17/)
  types(8,1:8) = (/21,8,11,27,13,13,13,13/)
  types(9,1:9) = (/21,8,11,27,13,13,13,3,18/)
  types(10,1:8) = (/21,8,11,27,13,13,13,13/)
  types(11,1:8) = (/21,8,11,27,13,10,28,10/)
  types(12,1:8) = (/21,8,11,27,13,6,25,20/)
  types(13,1:7) = (/21,8,11,27,7,7,7/)
  types(14,1:9) = (/21,8,11,27,13,13,6,25,20/)
  types(15,1:11) = (/21,8,11,27,13,13,4,19,4,19,19/)
  types(16,1:6) = (/21,8,11,27,5,24/)
  types(17,1:7) = (/21,8,11,27,5,24,13/)
  types(18,1:7) = (/21,8,11,27,13,13,13/)
  types(19,1:14) = (/21,8,11,27,13,12,12,14,22,14,14,14,14,14/)
  types(20,1:12) = (/21,8,11,27,13,30,30,30,30,30,9,26/)

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
    if(i == 21) then
      getatype = 0
      exit
    endif
    if(residues(i) == res) exit
  enddo
  do j = 1, 15
    if(j == 15 .or. i == 21) then
      getatype = 0
      exit
    endif
    if(atomnames(i,j) == atom) exit
  enddo
  if(i /= 21 .and. j /= 15) getatype = types(i,j)
end function getatype
  
!----------------------------------------------------!
end program contactmask
