module seams_pdbtools
  use seams_sequences
  use seams_utils
  !--------------------------------------------------------------------------
  ! Structures used in the pdbtools module to read the PDB and CIJ files
  !--------------------------------------------------------------------------
  ! Structure for PDB info
  type atomtype
    real,dimension(3)     :: xyz
    integer,dimension(3)    :: box
    integer         :: ires
    type(atomtype),pointer  :: next
  end type atomtype

  ! Structure for create cij (contact lists)
  type cijtype
    integer :: i,j,cij
    type (cijtype),pointer  :: next
  end type

CONTAINS
!-seams_pdbtools-------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  ! Function for PDB files
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ! getBetaResiduesPdb: Returns the beta residues of a protein
  !         Uses an external program which run the "stride" program
  !-------------------------------------------------------------------------------
  subroutine getBetaResiduesPdb (pdbFilename, betasArray)
    implicit none
    character (*), intent (in)      :: pdbFilename
    integer, allocatable, intent (out)  :: betasArray (:)
    character (len=200), allocatable  :: outputArray (:)
    integer               :: n, i
    real :: x
    character (len=500) :: command, strideDir, homedir,tmpdir

    call random_seed()
    ! Create the stridePath where stride program reside
    CALL get_environment_variable("GDIR", homedir)
    if (homedir=="")  then
            write(0,'("WARNING: GDIR environment variable not set!")')
    else
      strideDir = trim (homedir)//"/seams"
    endif

    CALL get_environment_variable("TMPDIR",tmpdir)
    if(tmpdir == "") tmpdir = "./tmp"
    command = trim (strideDir)//"/seams_stride.py " // trim(pdbFilename) // " "//trim(tmpDir)//"/"
    write(0,*) command
    call execute (command, outputArray)
    n = size (outputArray)

    allocate (betasArray (n))
    do i=1, n
      betasArray (i)  = str2Int (outputArray (i))
    enddo
  endsubroutine

  !-------------------------------------------------------------------------------
  ! Given the input  (pdbFilename) returns the contact map matrix (contactMatrix)
  !-------------------------------------------------------------------------------
  subroutine getContactMapPdb (pdbFilename, contactMatrix)
    implicit none
    character(*), intent (in)      :: pdbFilename
    integer, allocatable, intent (out) :: contactMatrix (:,:)
    integer, allocatable         :: contactsArray (:,:)
    integer                :: contact (2)
    integer                :: sizePdb, sizeContacts, i

    sizePdb = getSizePdb (pdbFilename)

    allocate (contactMatrix (sizePdb, sizePdb))
    contactMatrix = 0

    call getContactsPdb (pdbFilename, contactsArray, 8.0, 3) !sequence_int2d array

    sizeContacts = lengthSeq (contactsArray)

    do i=1, sizeContacts
      contact = getAtSeq(contactsArray, i)
      contactMatrix(contact(2),contact(1)) = 1
    enddo
  endsubroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read cotacts (AA1, AA2) from a contacts file (.cij)
  ! Contacts are written in two arrays
  !------------------------------------------------------------------------------
  subroutine getContactMapFromCij (cijFilename, contactMatrix)
    implicit none
    character (len=:),allocatable,intent(in)        :: cijFilename
    integer, allocatable, intent (out):: contactMatrix (:,:)
    integer, allocatable        :: contactsArray (:,:)
    character (len=50)          :: strLine
    integer               :: i,aax,aay,reason,sizeMatrix,contact (2)
    real                :: enr

    call createSeq (contactsArray)

    ! Read the file
    open (10, file=cijFilename)
    do
      read (10, "(A)", iostat=reason) strLine
      if (reason /=0) exit

      if (index (strLine, "!") > 0) cycle

      read (strLine,'(2i6,f5.2)', iostat=reason) aax, aay, enr

      call appendSeq (contactsArray, (/aax, aay/))
    enddo
        close (10)

    ! Create the contact map matrix
    sizeMatrix = maxval (contactsArray)

    allocate (contactMatrix (sizeMatrix, sizeMatrix))
    contactMatrix = 0

    do i=1, lengthSeq (contactsArray)
      contact = getAtSeq (contactsArray, i)
      contactMatrix (contact(2),contact(1)) = 1
    enddo
  endsubroutine
  !------------------------------------------------------------------------------
  ! getContactsPdb: calculateContacts the contacts map from a PDB file
  !------------------------------------------------------------------------------
  subroutine getContactsPdb (pdbFilename, contactsArray, minCut, minSeparation)
    implicit none
    character(*), intent (in)            :: pdbFilename
    real, intent (in)                 :: minCut !=8.0
    integer, intent (in)                :: minSeparation != 3
    integer, allocatable, dimension (:,:), intent (out) :: contactsArray
    integer  :: L, M, ios, bcut, outputUnit=6   != file unit to write results

    type(atomtype),pointer    :: atomroot, atoms
    type (cijtype),pointer    :: Cij,Cijroot

    allocate (atomroot,stat=ios)
    if (ios/=0) stop 'pdb2cij: error allocating atomroot.'

    atoms => atomroot
    call readpdb (pdbFilename, L, atoms)

    allocate (Cijroot, stat=ios)
    if (ios/=0) stop 'pdb2cij: error allocating a pointer.'

    atoms => atomroot
    Cij => Cijroot
    call calculateContacts (L, Cij, atoms, minCut, minSeparation)

    ! Write the contact sequence to an Array
    Cij => Cijroot
    call createSeq (contactsArray)
    do
      if (.not.associated(Cij%next))  exit
      Cij => Cij%next
      call appendSeq (contactsArray, (/Cij%i, Cij%j/))
    enddo

    atoms => atomroot
    Cij => Cijroot
    call cleanup(Cij,atoms)
  endsubroutine
  !-------------------------------------------------------------------------------
  ! Calculate the contacts of the protein in the "atoms" structure
  !-------------------------------------------------------------------------------
  subroutine calculateContacts (L, Cij, atoms, dcut, separation )
    implicit none
    integer,intent(in) :: L
    real, intent (in) :: dcut
    type(atomtype),pointer :: atoms
    type(atomtype),pointer :: jatoms
    integer :: I,J,ios,k, sep, separation, c
    type (cijtype),pointer :: Cij
    real :: d,dcut2, bcut, x
    dcut2 = dcut*dcut

    bcut = int(dcut)+1
    c=0
    do while (associated(atoms%next))
      !!c = c + 1
      !!print *, c
      atoms => atoms%next
      if (atoms%xyz(1)==0.) cycle !! this should not happen. much.
      jatoms => atoms

      sep = 1
      do while (associated(jatoms%next) .and. sep < separation)
      jatoms => jatoms%next
      sep = sep + 1
      enddo

      do while (associated(jatoms%next))
      jatoms => jatoms%next
      if (any(abs(jatoms%box-atoms%box)>bcut)) continue
      d = 0
      do k=1,3
        x = (atoms%xyz(k)-jatoms%xyz(k))
        d = d + x*x
      enddo

      if (d > dcut2) cycle

      allocate(Cij%next,stat=ios)
      if (ios/=0) stop 'pdb2cij.f90: error allocating a pointer to Cij'
      Cij => Cij%next
      Cij%i = atoms%ires
      Cij%j = jatoms%ires
      nullify(Cij%next)
      enddo
    enddo
  endsubroutine
  !
  !------------------------------------------------------------------------------
  ! readpdb: Read the PDB and renumber the residue residures starting to 1
  !------------------------------------------------------------------------------
  subroutine readpdb (pdbFilename, M, atoms)
    implicit none
    character(*),intent(in) :: pdbFilename
    integer,intent(out) :: M
    type(atomtype),pointer :: atoms
    character(len=80) :: aline
    integer :: ios, L, residueBase, residueNumber
    logical :: firstAtom = .true.  ! To renumbering residues starting at 1

    open(11,file=pdbFilename,status="old",form="formatted")
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

      if (firstAtom) then
        read(aline(23:26),*) residueNumber
        residueBase = residueNumber - 1
        firstAtom = .false.
      endif
      if(aline(27:27)==" ") read(aline(23:26),*) residueNumber

            !! We cannot use residue numbr because the file might have multiple
            !! chains. Instead, just increment a counter.
      !! atoms%ires = residueNumber - residueBase
      atoms%ires = L
            !! Note that chain ID is lost.
            !! Thu Jun 19 12:10:06 EDT 2014
      read(aline(31:54),*) atoms%xyz(1:3)
      atoms%box(1:3) = int(atoms%xyz(1:3))
      nullify(atoms%next)
      !! diagnostic
      !! write(*,*) atoms%ires,atoms%xyz(1:3),atoms%box(1:3)
    enddo
    close(11)
    M = L  !! last residue read
  endsubroutine

  !------------------------------------------------------------------------------
  ! getSizePdb: Get the number of residues of a PDB file
  !------------------------------------------------------------------------------
  function getSizePdb (pdbFilename) result (N)
    implicit none
    character(len=*),intent(in) :: pdbFilename
    character(len=80)     :: aline
    integer           :: N, ios

    open(11,file=pdbFilename,status="old",form="formatted")
    N = 0
    do
      !  read coordinates from the file
      read(11,'(a)',iostat=ios) aline
      if (ios /= 0) exit

      if (aline(1:4) /= "ATOM") cycle
      if (aline(13:16) /= " CA ") cycle

      N = N + 1
    enddo
    close(11)
  endfunction

  !--------------------------------------------------------------!!
  ! writeCij: write the residues "i", "j" in a array (file)
  !!--------------------------------------------------------------!!

  subroutine writeCij (L, Cij, minCut, minSeparation, outputUnit)
    implicit none
    integer,intent(in) :: L
    integer :: I,J,ios
    type (cijtype),pointer :: Cij
    real, intent (in)       :: minCut !=8.0
    integer, intent (in)      :: minSeparation != 3
    character (len=50)  :: output
    integer, optional :: outputUnit

    write (outputUnit, '((a) (f6.2))') "# Distance:", minCut
    write (outputUnit, '((a),(i6))'), "# Separation:", minSeparation

    do
      if (.not.associated(Cij%next))  exit
      Cij => Cij%next
      !write(*,'(2i6,f5.2)') Cij%i, Cij%j, 1.0
      write (output,'(I6,(A),I6,(A),F5.2)') Cij%i, " ", Cij%j, " ", 1.0
      write (outputUnit, FMT='(a)') adjustl (output)
      if (ios/=0) stop 'BUG allcoords'
    enddo
    close (outputUnit)
  end subroutine writeCij

  !--------------------------------------------------------------!!
  ! Free memory of this two structures used in getContactsPdb
  !--------------------------------------------------------------!!

  subroutine cleanup (Cij, atoms)
    implicit none
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
endmodule
