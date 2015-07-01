program pdb2cij
!! This program creates a contact map from a PDB file.
!! Output, a list on contacts in the format: i  j  1.00
!! The numbering in the PDB file is used as-is.
!! The definition of residue-residue distance is the distance between
!! beta-carbons (alpha-carbons for glycine).
!! REVISED: 24-APR-05 to use a linked list. is it faster?
!! Seems faster. Adding grid algorithm. Is it even faster?
!! CB.     Sun Apr 24 09:49:37 EDT 2005
!!----------------------------------------------------------------------
character(len=80) :: filename,outfile,aline
!! real, dimension(3,1000) :: atoms
!! integer, dimension(1000,1000) :: Cij
type atomtype
  real,dimension(3) :: xyz
  integer,dimension(3) :: box
  integer :: ires
  type(atomtype),pointer :: next
end type atomtype
type(atomtype),pointer :: atomroot, atoms
type cijtype
  integer :: i,j,cij
  type (cijtype),pointer :: next
end type
type (cijtype),pointer :: Cij,Cijroot
integer :: L,M,ios
real :: dcut=8.0
integer :: bcut

!! filename = "1aaj.pdb"
!! outfile = "1aaj.cij"
jarg = iargc()
if (jarg < 1) then
  write(*,'("Usage: xpdb2cij pdbfile [dcut=",f6.2,"]> cijfile")') dcut
  write(*,*) 'This program creates a list of contacts from coordinates.'
  write(*,*) 'The cutoff is ',dcut,' between C-beta atoms (or C-alpha for Gly)'
  stop 'pdb2cij.f90 v.  Wed May  7 16:18:53 EDT 2008'
endif
!! atoms = 0.
  
call getarg(1,filename)
if (jarg >= 2) then
  call getarg(2,aline)
  read(aline,*,iostat=ios) dcut
  if (ios/=0) stop 'pdb2cij: bad value for argument 2, dcut. Must be a float'
endif
bcut = int(dcut)+1


allocate(atomroot,stat=ios)
if (ios/=0) stop 'pdb2cij: error allocating atomroot.'
atoms => atomroot
call readpdb(filename,L,atoms)
allocate(Cijroot,stat=ios)
if (ios/=0) stop 'pdb2cij: error allocating a pointer.'
atoms => atomroot
Cij => Cijroot
call calculate(L,Cij,atoms)
Cij => Cijroot
call writeCij(L,Cij)
atoms => atomroot
Cij => Cijroot
call cleanup(Cij,atoms)

CONTAINS   ! end of program, start of subroutines

!!--------------------------------------------------------------!!
subroutine readpdb(myfile,M,atoms)
character(len=*),intent(in) :: myfile
type(atomtype),pointer :: atoms
integer,intent(out) :: M
character(len=80) :: aline
integer :: ios, L

open(11,file=myfile,status="old",form="formatted")
L = 0
do
  !  read coordinates from the file
  read(11,'(a)',iostat=ios) aline
  !! diagnostic
  !! write(*,*) trim(aline)
  if (ios /= 0) exit
  if (aline(1:4) /= "ATOM") cycle
  if (aline(18:20) == "GLY") then
    if (aline(13:16) /= " CA ") cycle
  else
    if (aline(13:16) /= " CB ") cycle
  endif
  L = L + 1
  allocate(atoms%next,stat=ios)
  if (ios/=0) stop 'pdb2cij: failed to allocate atom pointer.'
  atoms => atoms%next
  read(aline(23:26),*) atoms%ires
  read(aline(31:54),*) atoms%xyz(1:3)
  atoms%box(1:3) = int(atoms%xyz(1:3))
  nullify(atoms%next)
  !! diagnostic
  !! write(*,*) atoms%ires,atoms%xyz(1:3),atoms%box(1:3)
enddo
close(11)
M = L  !! last residue read

end subroutine readpdb
!!--------------------------------------------------------------!!
subroutine calculate(L,Cij,atoms)
integer,intent(in) :: L
type(atomtype),pointer :: atoms
type(atomtype),pointer :: jatoms
integer :: I,J,ios,k
type (cijtype),pointer :: Cij
real :: d,dcut2
dcut2 = dcut*dcut

do while (associated(atoms%next))
  atoms => atoms%next
  if (atoms%xyz(1)==0.) cycle   !! this should not happen. much.
  jatoms => atoms
  do while (associated(jatoms%next))
    jatoms => jatoms%next
    if (any(abs(jatoms%box-atoms%box)>bcut)) continue
    d = 0
    do k=1,3
      x = (atoms%xyz(k)-jatoms%xyz(k))
      d = d + x*x
    enddo
    if (d > dcut2) cycle
!      do I=1,L-1
!        if (atoms(1,I)==0.) cycle
!        do J=I+1,L
!          if (atoms(1,J)==0.) cycle
!          d = (atoms(1,I)-atoms(1,J))**2  &
!            + (atoms(2,I)-atoms(2,J))**2  &
!            + (atoms(3,I)-atoms(3,J))**2
!          d = sqrt(d)
!          if (d < dcut) then
    allocate(Cij%next,stat=ios)
    if (ios/=0) stop 'pdb2cij.f90: error allocating a pointer to Cij'
    Cij => Cij%next
    Cij%i = atoms%ires
    Cij%j = jatoms%ires
    nullify(Cij%next)
  enddo
enddo

end subroutine calculate
!!--------------------------------------------------------------!!

subroutine writeCij(L,Cij)
!! character(len=*),intent(in) :: myfile
integer,intent(in) :: L
integer :: I,J
type (cijtype),pointer :: Cij

! open(12,file=myfile,status="new")
do 
  if (.not.associated(Cij%next))  exit
  Cij => Cij%next
  write(*,'(2i6,f5.2)') Cij%i, Cij%j, 1.0
enddo
!do I=1,L
!  do J=1,L
!    if (Cij(I,J) == 1) then
!      ! write(12,'(2i4,f5.2)') I,J,1.0
!      write(*,'(2i4,f5.2)') I,J,1.0
!    endif
!    !! write(12,"(i2,$)") Cij(I,J)
!  enddo
!  !! write(12,*) 
!enddo
! close(12)
end subroutine writeCij
!!--------------------------------------------------------------!!

subroutine cleanup(Cij,atoms)
  type (cijtype),pointer :: Cij, cptr
  type(atomtype),pointer :: atoms, jatoms
  do 
    if (.not.associated(Cij%next))  exit
    cptr => Cij%next
    deallocate(Cij)
    Cij => cptr
  enddo
  deallocate(Cij)
  do while (associated(atoms%next))
    jatoms => atoms%next
    deallocate(atoms)
    atoms => jatoms
  enddo
  deallocate(atoms)
end subroutine cleanup


!!--------------------------------------------------------------!!
end program pdb2cij
