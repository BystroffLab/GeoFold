module seams_main
  use seams_pdbtools
  use seams_sequences
  use seams_graph
  use seams_utils
  use seam_debug

  integer, parameter :: MIN_CONTACTS_REGION=15  !15 !! Minimal number of contacts for regions
  integer, parameter :: MIN_CONTACTS_BETAS=5    !5 !! Minimal number of beta contacts
  integer, parameter :: MINIMUM_OVERLAP=3 ! Minimum overlap between beta strands in contact to be a seam
  integer :: nResidues

  !! Structure for barrels
  type :: barrel_type
    integer :: nSeams  ! number of barrels_array
    type (seam_data) :: seams (100) ! Array of seams
  endtype

CONTAINS
!-seams_main-------------------------------------------------------------------
  !------------------------------------------------------------------------------
  ! seams_main:
  !      Functions of the module seams_main
  !      Given a contact map (contactMatrix), it detects the regions
  !      (betas, helix, coils), and remove the ones that do not have contiguous
  !      contacts greater than a threshold
  !------------------------------------------------------------------------------
  ! Return True/False if the protein (pdbFilename) is a beta barrel
  ! - It uses the graph property that if the number of edges connecting the
  ! is equal to the number of vertices, then this graph has a complete cycle
  !----------------------------------------------------------------------------
  subroutine getBetaSheetsPdb (pdbFilename, betaSheetSeq)
    implicit none
      character (*), intent (in)      :: pdbFilename        ! PDB
      type(sequence_node), pointer      :: betaSheetSeq
      integer, allocatable          :: contactMatrix (:,:), betaResidues (:)
      ! integer               :: nResidues, minContacts !! Minimum number of contacts per region
      integer               :: minContacts !! Minimum number of contacts per region

    minContacts = MIN_CONTACTS_REGION   ! Minimal number of contacts for beta sheet
    nResidues = getSizePdb(pdbFilename)

    call getContactMapPdb (pdbFilename, contactMatrix) ! Get the contact map from the PDB
    !call writeMatrix (contactMatrix, "in.mat")           ! for debugging
    ! call dwriteMatrix(contactMatrix,nResidues,"/home/walcob/Desktop/getContactMapPdb.cij")
    call getBetaResiduesPdb  (pdbFilename, betaResidues) ! Get the residues in betas (according to the stride program)
    ! call dwriteMatrix(contactMatrix,nResidues,"/home/walcob/Desktop/getBetaResiduesPdb.cij")
    ! Testing....
    call bridgeBulges (contactMatrix)
    call detectBetaSheets (contactMatrix, betaResidues, minContacts, betaSheetSeq)
    ! call dwriteMatrix(contactMatrix,nResidues,"/home/walcob/Desktop/detectBetaSheets.cij")
    !call writeMatrix (contactMatrix, "out.mat")          ! for debugging
    ! call bridgeBulges (contactMatrix)
    ! call dwriteMatrix(contactMatrix,nResidues,"/home/walcob/Desktop/bridgeBulges.cij")
    ! call detectBetaSheets (contactMatrix, betaResidues, minContacts, betaSheetSeq)
    ! call dwriteMatrix(contactMatrix,nResidues,"/home/walcob/Desktop/detectBetaSheets2.cij")
    !call writeMatrix (contactMatrix, "bulge.mat")          ! for debugging
  endsubroutine getBetaSheetsPdb

  !--------------------------------------------------------------------
  ! Get the Total energy (seam + buttons) for a seam "iSeam"
  ! Returns the total energy and the side with low energy (beta1 or beta2)
  !--------------------------------------------------------------------
  function getEnergySeam (barrelsArray, energyMatrix, iBarrel, iSeam, side, iflags) result(energySeam)
    implicit none
      type (barrel_type), intent (in)    :: barrelsArray (:)
      real, intent (in)                  :: energyMatrix (:,:)
      integer, intent (in)               :: iBarrel, iSeam
      integer, intent (out)              :: side
      character, intent (in), optional   :: iflags       ! Unfolded residues not taken into account
      real                               :: energySeam, energyBeta1, energyBeta2
      integer                            :: nContacts, k, contact (2), c1, c2
      type(seam_data), target            :: seam
      integer                            :: nButtons, residue, x1,x2,y1,y2

    seam = barrelsArray(iBarrel)%seams(iSeam)
    nContacts = seam%n
    energySeam = 0.0
    do k=1, nContacts
      contact = seam%x (:, k)
      c1 = contact (1)
      c2 = contact (2)
      if (present (iflags) .and. iflags(c1:c1) == '.' .or. iflags (c2:c2) == '.') cycle
      energySeam = energySeam + energyMatrix (contact(1), contact(2))
    enddo

    !! Get energy of seam's buttons
    y1 = seam%segments (1) ! Beta 1
    y2 = seam%segments (2) ! Beta 1
    x1 = seam%segments (3) ! Beta 2
    x2 = seam%segments (4) ! Beta 2

    energyBeta1 = 0
    energyBeta2 = 0
    nButtons = seam%nButtons
    do k=1, nButtons
      residue = seam%buttons (k)%residue
      if (present (iflags)) then
        if (iflags (residue:residue) == '.') cycle !! broken contact
        energyBeta1 = energyBeta1 + getEnergyButtonBeta (residue, y1, y2, energyMatrix, iflags)
        energyBeta2 = energyBeta2 + getEnergyButtonBeta (residue, x1, x2, energyMatrix, iflags)
      else
        energyBeta1 = energyBeta1 + getEnergyButtonBeta (residue, y1, y2, energyMatrix)
        energyBeta2 = energyBeta2 + getEnergyButtonBeta (residue, x1, x2, energyMatrix)
      endif
    enddo

    energySeam = energySeam + MINVAL ((/energyBeta1, energyBeta2/))
    if (energyBeta1 < energyBeta2) then
      side = -1
    else
      side = 1
    endif
  endfunction
  !------------------------------------------------------------------------------
  ! Calculate the energy for the residue with the beta (coords x1, x2)
  ! It uses the energy matrix, and an optional iflags array are with  marked broken contacts
  !------------------------------------------------------------------------------
  function getEnergyButtonBeta (residue, x1, x2, energyMatrix, iflags) result (energy)
    implicit none
    integer, intent (in)              :: residue, x1, x2 ! coordinates beta
    real, intent (in)                 :: energyMatrix (:,:)
    character, intent (in), optional  :: iflags              ! Unfolded residues not taken into account
    integer, allocatable              :: contactsResidueSeq (:)! residues forming contacts with input "residue"
    real                              :: energy
    integer                           :: contact, i, j

    if (present (iflags)) then
      print *, ">>>", iflags
      stop
    endif
    call getContactsResidue (residue, energyMatrix, contactsResidueSeq)
    energy = 0
    do i=1, lengthSeq (contactsResidueSeq)
      contact = getAtSeq (contactsResidueSeq, i)   ! get contact

      if (contact < x1 .or. contact > x2) cycle  ! not in beta
      if (present (iflags)) then
        if (iflags (contact:contact) == '.') cycle  ! Broken
      endif

      energy = energy + energyMatrix (residue, contact)
    enddo
  endfunction
  !------------------------------------------------------------------------------
  ! Get all contacts for an residue acid
  ! OPTIMIZE by reading or receiving directly the contact list
  !------------------------------------------------------------------------------
  subroutine getContactsResidue (residue, energyMatrix, contactsSeq)
    integer, intent (in)              :: residue
    real, intent (in)           :: energyMatrix (:,:)
    integer, allocatable, intent (out)  :: contactsSeq (:)

    n = size (energyMatrix, 1)
    call createSeq (contactsSeq)
    do i=1, n
      if (energyMatrix (residue, i) /= 0 ) call appendSeq (contactsSeq, i)
    enddo
  endsubroutine

  !-----------------------------------------------------------------------
  ! Get contacts's residues making contact with the seams's bsheets (buttons)
  !-----------------------------------------------------------------------
  subroutine getButtonsSeams (barrelsArray, energyMatrix)
    implicit none
      type (barrel_type), intent (inout) :: barrelsArray (:)
      real, intent (in)           :: energyMatrix (:,:)
      integer                     :: nResidues, nButtons, i,j,k,x1,x2,y1,y2
      type(barrel_type), target   :: barrel
      type(seam_data), target     :: seam
      type(button_type), target   :: button
      integer, allocatable        :: contactsResidueSeq (:), buttonsSeq (:)
      integer, allocatable        :: contactsBeta1Seq (:), contactsBeta2Seq (:)
      real                        :: energyB1, energyB2

    nResidues = size (energyMatrix, 1) ! number of residues
    do i=1, size (barrelsArray)   ! number of barrels
      barrel = barrelsArray (i)
      do j=1, barrel%nSeams    !number of Seams
        seam = barrel%seams (j)
        call createSeq (buttonsSeq)
        y1 = seam%segments(1)
        y2 = seam%segments(2)
        x1 = seam%segments(3)
        x2 = seam%segments(4)

        !! All residues....optimize to only the ones forming contacts
        nButtons = 0
        do k=1, nResidues
          if (x1<=k.and.k<=x2 .or. y1<=k.and.k<=y2) cycle ! Belongs to seam

          energyB1 =  getEnergyButtonBeta (k,y1,y2, energyMatrix)
          if (energyB1 == 0) cycle ! No button

          energyB2 =  getEnergyButtonBeta (k,x1,x2, energyMatrix)
          if (energyB2 == 0) cycle ! No button

          !print *, ">>> Buttons: ", i, j, k

          button%residue = k
          button%energyBeta1 = energyB1
          button%energyBeta2 = energyB2

          nButtons = nButtons + 1
          seam%buttons (nButtons) = button

          barrelsArray (i)%seams(j)%buttons (nButtons) = button

          call appendSeq (buttonsSeq, k)
        enddo
        barrelsArray (i)%seams(j)%nButtons = nButtons
      enddo
     enddo
    endsubroutine

  !------------------------------------------------------------------------------
  ! Assign energy to beta sheets using the energies in the input file "energyFilename"
  ! Also, calculate the delimiting points of each beta sheet to be used by geofold
  !------------------------------------------------------------------------------
  subroutine assignEnergyToBarrels (barrelsArray, energyFilename, nResidues, energyMatrix)
    implicit none
      type (barrel_type), intent (inout) :: barrelsArray (:)
      character (*), intent (in)      :: energyFilename
      integer, intent (in)            :: nResidues
      real, allocatable, intent (out) :: energyMatrix (:,:)
      real                            :: energy
      integer                         :: nContacts, nSeams, k, i, j, contact (2)
      type(seam_data), target         :: seam

    call loadEnergyToMatrix (energyFilename, energyMatrix, nResidues)

    do k=1, size (barrelsArray)
      nSeams = barrelsArray(k)%nSeams
      do i=1, nSeams
        seam = barrelsArray(k)%seams(i)
        nContacts = seam%n
        energy = 0.0
        do j=1, nContacts
          contact = seam%x (:, j)
          energy = energy + energyMatrix (contact(1), contact(2))
        enddo
        barrelsArray (k)%seams (i)%energy = energy
      enddo
    enddo
  endsubroutine

  !------------------------------------------------------------------------------
  ! Load contacts energies to a matrix
  ! Read the energies from file and create and fill a matrix with them
  !------------------------------------------------------------------------------
  subroutine loadEnergyToMatrix (energyFilename, energyMatrix, nResidues)
    implicit none
      character (*), intent (in)    :: energyFilename
      real, allocatable, intent (out)       :: energyMatrix (:,:)
      integer, intent (in)        :: nResidues
      character (len=50)        :: aline
      integer             :: aa1, aa2, ios
      real                :: energy, x

    allocate (energyMatrix (nResidues, nResidues))
    energyMatrix = 0
    open(99,file=energyFilename,status="old",form="formatted")
    do
      read(99,'(a)',iostat=ios) aline
      if (ios /= 0) exit
      if (index (aline, "#") /= 0) cycle ! a comment

      !read(aline, '(2i6,2f9.3)',iostat=ios) aa1, aa2, energy, x
      read(aline, '(2i6,2f9.3)',iostat=ios) aa1,aa2,x,energy
            if (ios/=0) then
               write(0,*) trim(aline)
               stop
            endif

      energyMatrix (aa1, aa2) = energy
      energyMatrix (aa2, aa1) = energy
    enddo
    close(99)
  endsubroutine
  !------------------------------------------------------------------------------
  ! Given an input sequence of beta sheets (betaSheetSeq) it finds the contacts
  ! that connects each beta sheet with the others (betaSheetContactsSeq)
  ! This informations can be used to detect if there is a cycle
  ! NOTE: There is a fast checking if connected, using only limit points of the region
  !------------------------------------------------------------------------------
  subroutine getBetaSheetContacts (betaSheetSeq, betaSheetContactsSeq)
    implicit none
    type(sequence_node), pointer, intent (in) :: betaSheetSeq      ! Sequence of beta sheets each one with is set of contacts
    integer, allocatable, intent (out)      :: betaSheetContactsSeq (:,:) ! Contacts (x,y) which conneact the beta sheets
    integer                   :: n, m, i, k
    type(seam_data), target         :: r1, r2

    call createSeq (betaSheetContactsSeq)

    k = lengthSeqG (betaSheetSeq)
    do n=1, k
      r1 = getAtSeq (betaSheetSeq, n)
      do m=1, k
        r2 = getAtSeq (betaSheetSeq, m)

        ! Check comparing same regions and contacts visited
        if (n==m)                   cycle
        if (existsSeq (betaSheetContactsSeq, (/n,m/)))   cycle
        if (existsSeq (betaSheetContactsSeq, (/m,n/)))   cycle

        if (isConnectedRegionsFast (r1, r2)) then
          call appendSeq (betaSheetContactsSeq, (/n, m/))
        endif
      enddo
    enddo
  endsubroutine

  ! Create an undirected graph from the betasheets
    ! edges are duplicated !! n-->m and m-->n
  subroutine getBetaSheetGraph (betaSheetSeq, betaSheetContactsSeq, graphSeq)
    implicit none
    type(sequence_node), pointer, intent (in) :: betaSheetSeq    ! Sequence of beta sheets each one with is set of contacts
    integer, allocatable, intent (out)      :: betaSheetContactsSeq (:,:) ! Contacts (x,y) which connect the beta sheets
    type(sequence_node), pointer, intent (in) :: graphSeq      ! Undirected Graph
    integer                   :: n, m, i, k
    type(seam_data), target         :: r1, r2

    call createSeq (betaSheetContactsSeq)
    call createGraphSeq (graphSeq)

    k = lengthSeqG (betaSheetSeq)
    do n=1, k
      call addVertexGraphSeq (graphSeq, n)
      r1 = getAtSeq (betaSheetSeq, n)
      do m=1, k
        r2 = getAtSeq (betaSheetSeq, m)

        ! Check comparing same regions and contacts visited
        if (n==m)                   cycle
        !if (existsSeq (betaSheetContactsSeq, (/n,m/)))   cycle
        !if (existsSeq (betaSheetContactsSeq, (/m,n/)))   cycle

        if (isConnectedRegionsFast (r1, r2)) then
          call appendSeq (betaSheetContactsSeq, (/n, m/))
          call addEdgeGraphSeq (graphSeq, n, m )
        endif
      enddo
    enddo
  endsubroutine

  !------------------------------------------------------------------------------
  ! isConnectedRegion:
  ! - Return True/False if input regions "reg1" and "reg2" are connected
  ! - Connected means, there is almost one contact shared with the regions
  !------------------------------------------------------------------------------
   function isConnectedRegionsFast (reg1, reg2) result (value)
        !! Are two seams connected by a vertex?
    implicit none
    type(seam_data), target    :: reg1, reg2
    logical         :: value
    integer         :: x1, x2, y1, y2, w1, w2, h1, h2, last
    
    value = .false.

    y1   = reg1%segments(1)
    last = reg1%segments(2)
    w1   = last - y1 + 1

    x1   = reg1%segments(3)
    last = reg1%segments(4)
    h1   = last - x1 + 1

    y2   = reg2%segments(1)
    last = reg2%segments(2)
    w2   = last - y2 + 1

    x2   = reg2%segments(3)
    last = reg2%segments(4)
    h2 = last - x2 + 1

    if (y1+w1>y2 .and. y2+w2>y1) then
        if(getOverlap(reg1%segments(1),reg1%segments(2),reg2%segments(1),reg2%segments(2)) >= MINIMUM_OVERLAP) value = .true.
    elseif(y1+w1>x2 .and. x2+h2>y1) then
        if(getOverlap(reg1%segments(1),reg1%segments(2),reg2%segments(3),reg2%segments(4)) >= MINIMUM_OVERLAP) value = .true.    
    elseif(x1+h1>y2 .and. y2+w2>x1) then
        if(getOverlap(reg1%segments(3),reg1%segments(4),reg2%segments(1),reg2%segments(2)) >= MINIMUM_OVERLAP) value = .true.
    elseif(x1+h1>x2 .and. x2+h2>x1) then
        if(getOverlap(reg1%segments(3),reg1%segments(4),reg2%segments(3),reg2%segments(4)) >= MINIMUM_OVERLAP) value = .true.
    else
      value = .false.
    endif
    endfunction
    
integer function getOverlap(x1,x2,y1,y2) result(overlap)
    implicit none
    integer, intent(in) :: x1,x2,y1,y2
    
    if(x1 < y1) then
        if(x2 < y2) then
            overlap = x2 - y1
        else
            overlap = y2 - y1
        endif
    else
        if(y2 < x2) then
            overlap = y2 - x1
        else
            overlap = x2 - x1
        endif
    endif
endfunction
      
  function isConnectedRegions (reg1, reg2) result (value)
    type(seam_data), target    :: reg1, reg2
    logical         :: value
    integer         :: nReg1, nReg2, i, x, y
    integer, allocatable  :: r2y (:), r2x (:)

    nReg1 = reg1%n

    value = .false.
    do i = 1, nReg1
      y = reg1%x (1, i)
      x = reg1%x (2, i)

      nReg2 = reg2%n
      allocate (r2y (nReg2))
      allocate (r2x (nReg2))

      r2y = reg2%x (1, 1:nReg2)
      r2x = reg2%x (2, 1:nReg2)

      if (any (y==r2y) .or. any(y==r2x) .or. any (x==r2x) .or. any (x==r2y)) then
        value = .true.
        return
      endif

      deallocate (r2y)
      deallocate (r2x)
    enddo
  endfunction

  !----------------------------------------------------------------------------
  ! detectBetaSheets
  ! Detect all the regions from the contactMatrix which has more contacts
  ! than "minContacts". The other regions are cleared (zeroes)
  !----------------------------------------------------------------------------
  subroutine detectBetaSheets (contactMatrix, betaResidues, minContacts, betaSheetSeq)
    implicit none
      integer, intent (inout)    :: contactMatrix (:,:)        ! Matrix NxN of contacts
      integer,allocatable,intent(in) :: betaResidues (:)           ! info of the residues (x,y) in beta structures
      integer, intent (in)       :: minContacts            ! Minimal number of contacts for bsheet
      type (sequence_node), pointer  :: betaSheetSeq     ! Final sequence of found bsheets
      integer, allocatable       :: contacts(:,:), visitedSeq(:,:)   ! Secuence of integer 2DHandle with seq2I
      integer, allocatable       :: betaContactsSeq (:,:)      ! Beta sheet contacts ((x,y),(..),..)
      integer            :: minX, minY, maxX, maxY
      integer            :: n, i, j, sh (2)
      integer,save :: count
      character(len=:),allocatable :: base,filename,str_count
      character(len=6000) :: buffer

    count = 0
    call createSeq (visitedSeq)
    call createSeq (betaSheetSeq)
    base = "/home/walcob/Desktop/detectBetaSheets_"
    call itoa(count,str_count)
    filename = base//str_count//".cij"
    ! call dwriteMatrix(contactMatrix,nResidues,filename)
    n = nrows (contactMatrix)
    do i=1, n
      do j=1, n
        if (existsSeq (visitedSeq, (/j, i/)) .or. contactMatrix (j, i) == 0) cycle
        call getContacts (contactMatrix, (/j, i/), contacts)
        call extendsSeq (visitedSeq, contacts)

                !! Erase contact regions smaller than minContacts MIN_CONTACTS_REGION
        if (lengthSeq (contacts) < minContacts) then
          call clearRegion (contactMatrix, contacts)
          !!!!!!BDW DEBUGGING!!!!!
          count = count + 1
          call itoa(count,str_count)
          filename = base//str_count//".cij"
          ! call dwriteMatrix(contactMatrix,nResidues,filename)c      
          cycle
        endif
        !call selectBetas (betaResidues, contacts, contactMatrix, betaContactsSeq)
        !if ( isEmptySeq (betaContactsSeq)) cycle
        !call appendSeq (betaSheetSeq, betaContactsSeq)

        call selectBetas (betaResidues, contacts, contactMatrix, betaContactsSeq)
        if (lengthSeq (betaContactsSeq) < MIN_CONTACTS_BETAS) then
          call clearRegion (contactMatrix, contacts)
          cycle
        endif

        call appendSeq (betaSheetSeq, betaContactsSeq)
       enddo
    enddo
  endsubroutine



  !-------------------------------------------------------------------------------
  ! Select only beta contacts from contact sequence using the info in the beta sequence
  ! Also clear the no-beta contacts in the contactMatrix
  !-------------------------------------------------------------------------------
  subroutine selectBetas (betaResidues, contactSeq, contactMatrix, betaContactsSeq)
    integer, allocatable, intent (in) :: betaResidues (:)
    integer, allocatable, intent (in) :: contactSeq(:,:)
    integer, intent (inout)       :: contactMatrix (:,:)
    integer, allocatable, intent (out)          :: betaContactsSeq (:,:)
    integer               :: n, i, x, y, contact (2)
        integer, parameter                  :: betaextend=1
!maybe increase betaextend to 2?
    call createSeq (betaContactsSeq)
    n = lengthSeq (contactSeq)
    do i=1, n
      contact = getAtSeq (contactSeq, i)
      x = contact (1)
      y = contact (2)
      ! if (any (betaResidues .eq. x) .and. (any (betaResidues .eq. y))) then
            !! consider +/- betaextend residue from the sheet to be beta contacts
      if (any (abs(betaResidues-x)<=betaextend) .and. &
               (any (abs(betaResidues-y)<=betaextend))) then
        call appendSeq (betaContactsSeq, contact)
      else
        contactMatrix (x,y) = 0
      endif
     enddo
  endsubroutine

  !----------------------------------------------------------------------------
  ! getNeighbors
  ! Get the neighbors of of the point (x,y), the surrounding poings
  !----------------------------------------------------------------------------
  subroutine getNeighbors (contactMatrix, point, neighbors)
    implicit none
    integer, dimension (:,:), intent (in) :: contactMatrix
    integer, allocatable, dimension (:,:),intent (inout)   :: neighbors
    integer, intent (in)                   :: point(2)
    integer                          :: n, i, j, x, y

    !print *, "Get neighbors of", point(1), point (2)
    n = nrows (contactMatrix)
    call createSeq (neighbors)
    x = point (1)
    y = point (2)

    do i=x-1, x+1
      do j=y-1, y+1
        if (i<1 .or. j <1 .or. i == n+1 .or. j == n+1) cycle
        if (i /=x  .or. j /= y) call appendSeq (neighbors, (/i,j/))
      enddo
    enddo
  endsubroutine
  !----------------------------------------------------------------------------
  ! Return the sequence of friend of a i,j point
  ! The contacts correspond to all interconected points to the i.j point
  !----------------------------------------------------------------------------
  subroutine getContacts (contactMatrix, point, contacts)
    implicit none
    integer, dimension (:,:), intent (in) :: contactMatrix
    integer, intent (in)          :: point (2)
    integer, allocatable, dimension (:,:), intent (out) :: contacts
    integer, allocatable, dimension (:,:) :: pointQueue, neighborSeq
    integer                 :: tmpPoint(2), x, y, n, i, j, neighbor(2)
    character (len=15)            :: text

    call createSeq (contacts)
    call createSeq (pointQueue)
    call appendSeq (contacts, point)
    call appendSeq (pointQueue, point)
    do while ( .not. isEmptySeq (pointQueue) )
      tmpPoint = popSeq (pointQueue,1)
      call getNeighbors (contactMatrix, tmpPoint, neighborSeq)
      do while ( .not. isEmptySeq (neighborSeq) )
        neighbor = popSeq (neighborSeq, 1)
        x = neighbor (1)
        y = neighbor (2)
        if ( contactMatrix (x, y) == 0) cycle
        if ( .not. existsSeq (contacts, neighbor) ) then
          call appendSeq (pointQueue, neighbor)
          call appendSeq (contacts, (/x,y/))
        endif
       enddo
     enddo
  endsubroutine
  !----------------------------------------------------------------------------
  ! clearRegion:
  ! Crear the contactMatrix region (assign 0) correspondig to the points of input sequence
  !----------------------------------------------------------------------------
  subroutine clearRegion (contactMatrix, pointSeq)
    implicit none
    integer, dimension (:,:), intent (inout):: contactMatrix
    integer, allocatable, dimension (:,:), intent (inout) :: pointSeq
    integer                 :: minContacts, point (2)
    integer, allocatable, dimension (:,:) :: contacts ! Handle with seq2I
    do while ( .not. isEmptySeq (pointSeq) )
      point = popSeq (pointSeq, 1)
      contactMatrix (point(1), point (2)) = 0
    enddo
  endsubroutine

  !----------------------------------------------------------------------------
  ! writeRegion
  !----------------------------------------------------------------------------
  subroutine writeRegion (contacts)
    integer, dimension (:,:), intent (in):: contacts
    integer :: n, rows
    character (len=10) :: formatStr

    n = size (contacts)
    rows = nrows (contacts)

    !formatStr = pasteS (str(n), "(I5)")

    !write(*, formatStr) (a(1:2,j), j=1,rows)
  endsubroutine

  !----------------------------------------------------------------------------
  ! bridgeBulge  -- add contacts wherever surrounding contacts exist.
  !----------------------------------------------------------------------------
    subroutine bridgeBulges(contacts)
      integer,dimension(:,:),intent(inout) :: contacts
      integer :: i,j, dim(2)
      integer,parameter :: bulgesize=2

        dim = shape (contacts)

        do i=1, dim(1)-bulgesize-1
            do j=1, dim(2)-bulgesize-1
                if (contacts(i,j)/=0.and. &
                    contacts(i+bulgesize+1,j)/=0.and. &
                    any(contacts(i+1:i+bulgesize,j)==0)) then
                        contacts(i+1:i+bulgesize,j) = 1
                elseif (contacts(i,j)/=0.and. &
                    contacts(i,j+bulgesize+1)/=0.and. &
                    any(contacts(i,j+1:j+bulgesize)==0)) then
                        contacts(i,j+1:j+bulgesize) = 1
                endif
            enddo
        enddo
    endsubroutine bridgeBulges

endmodule
