! L.Garreta, C. Bystroff  2014
!--------------------------------------------------------------------
! pdb2seams.f90
!--------------------------------------------------------------------

program pdb2seams
  use seams_main
  use seams_sequences
  use seams_graph
  use seams_utils

  implicit none
    character (len=100)       :: pdbFilename, contactsFilename
    type(sequence_node), pointer    :: seamsSeq, graphSeq, cyclesSeq
    integer, allocatable        :: betaSheetContactsSeq (:,:)
    integer             :: numberAminos
    real, allocatable             :: energyMatrix (:,:)
    type (barrel_type), allocatable :: barrelsArray (:)

  if (iargc() < 1) then
    print *,"USAGE : pdb2seams <input pdb filename> [input energies file] "
    print *,"Identifies beta barrels in PDB files. Writes barrel data, including all seams    "
    print *,"all residues H-bonding across a seam, and all buttons.   "
    print *,"INPUT: A protein structure in PDB format"
    STOP 'pdb2seams.f90 v.  Wed Jun 18 22:19:04 EDT 2014'
  endif

  call getarg (1, pdbFilename)
  numberAminos = getSizePdb (pdbFilename)
    !! diagnostic
    write(*,*) "Size of PDB file = ",numberAminos

  
  ! Get bsheets, convert to a graph, get all cycles, convert to barrel array
  call getBetaSheetsPdb (pdbFilename, seamsSeq)
  call getBetaSheetGraph (seamsSeq, betaSheetContactsSeq, graphSeq)
  call getAllCyclesGraphSeq (graphSeq, cyclesSeq)
  !stops at printGraphSeq???
  call printGraphSeq(graphSeq)
  ! write (0,*) "we got this far?"
  call createBarrels (cyclesSeq, seamsSeq, barrelsArray)  ! Create barrels struct from seams and cycles info
  
  ! write (0,*) "barrels created?"
    !! diagnostic

  ! Calculate and assign energies if energies file given
  if (iargc() == 2) then
    call getarg (2, contactsFilename)
    ! write(0,*)"contactsFilename = ",contactsFilename
    !! assigning energy to barrels not working?
    call assignEnergyToBarrels (barrelsArray, contactsFilename, numberAminos, energyMatrix)
    ! write(0,*) "energy assigned"
    call getButtonsSeams (barrelsArray, energyMatrix) ! Get residues making contact with seam's bsheets
  endif

  ! write number of barrels, barrels, seams, buttons
  call writeBarrels (barrelsArray)

CONTAINS
!---------------------------------------------------------------------------
! Create a barrel structure from a graph structure (cycles and the seams)
!---------------------------------------------------------------------------
subroutine createBarrels (cyclesSeq, seamsSeq, barrelsArray)
  implicit none
  type(sequence_node), pointer    :: cyclesSeq, seamsSeq
  type (barrel_type), allocatable, intent (out) :: barrelsArray (:)
  type(seam_data), target     :: cycleReg, seam ! cycle data type record 
  type(barrel_type), target     :: barrel

  integer        :: cycleSize, nCycles, i, j, id

  nCycles = lengthSeq (cyclesSeq)
  allocate (barrelsArray(nCycles))

  do i=1, nCycles
    cycleReg = getAtSeq(cyclesSeq, i)
    cycleSize = cycleReg%n
    do j=1, cycleSize
      id = cycleReg%x(1, j)
      seam = getAtSeq(seamsSeq, id)
      seam%id = id
      seam%nButtons = 0
      barrel%nSeams = cycleReg%n
      barrel%seams (j) = seam
    enddo
    barrelsArray (i) = barrel
  enddo
endsubroutine

!---------------------------------------------------------------------------
! Print the entire sequence to the screen (or a file)
!---------------------------------------------------------------------------
subroutine writeBarrels(barrelsArray)
  implicit none
  type (barrel_type), allocatable, intent (in)  :: barrelsArray(:)
  integer        :: nBarrels, nSeams, i, j, k, nButtons
  type(seam_data), target     :: seam ! cycle data type record 
  type(barrel_type), target     :: barrel 
  type(button_type), target     :: button 

  nBarrels = size(barrelsArray)

  write (*,"(A)") "# ========== Results from PDB2SEAMS =========="
  write (*,"(A)") "# NBARRELS  <Number of Barrels>"
  write (*,"(A)") "# BARREL :  <id>  <Number of Seams>"
  write (*,"(A)") "# SEAM :    <id>  <Number of Buttons> <Total Energy> <x1-beta1> <x2-beta1> <x1-beta2> <x2-beta2>"
  write (*,"(A)") "# BUTTON :  <id>  <Energy-Beta1> <Energy-beta2>"

  write (*,'("NBARRELS  ",I5)') nBarrels
  do i=1, nBarrels
    barrel = barrelsArray(i)
    nSeams = barrel%nSeams
    write (*, '("BARREL ",I5," ",I5)') i, nSeams
    do j=1, nSeams
      seam = barrel%seams(j)
      write (*, '("SEAM    ",I5,I5,f10.3,4I5)')  &
             seam%id, seam%nButtons, seam%energy, seam%segments
      nButtons = barrel%seams(j)%nButtons
      do k=1, nButtons
        button = barrel%seams(j)%buttons(k)
        write (*, '("BUTTON  ", I5, 2f10.3)') &
           button%residue, button%energyBeta1, button%energyBeta2
      enddo
    enddo
  enddo
endsubroutine
endprogram
