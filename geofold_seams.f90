!------------------------------------------------------------------------
! Goal:  Loads "seams" from seams files
!   Seams files contains info of protein's beta barrels_array
!   Each beta barrel contains a list of seams (2 bstrands)
!   A seam contains its number (pos in contact map) , energy, and segments
! Author: Luis Garreta
! Date: Jun/06/2013
!------------------------------------------------------------------------
module geofold_seams
  use geofold_global
  ! Move for seams
  type :: seammove_type
    integer    :: barrel
    integer    :: seam
    real       :: energy
    integer    :: side
  endtype

  !call test() !To check the module

CONTAINS
!----------------------------------------------------------------------------
! Dummy function to call pivots function
!----------------------------------------------------------------------------
!function geofold_seams_inseam (f, ires, jres ) result (inseam)
! implicit none
!   type(intermediate), POINTER :: f
!   integer, intent (in)        :: ires, jres
!   logical                     :: inseam
!    inseam = geofold_pivots_inseam(f, barrels_array, ires, jres )
!end function geofold_seams_inseam
!----------------------------------------------------------------------------
! Only to test the main function of the module
!----------------------------------------------------------------------------
subroutine test (filename)
  character (len=*)                 :: filename
  integer, allocatable              :: allButtons (:)
    integer                           :: i, iseam, nButtons

    allocate(allButtons(200),stat=ios) ; if (ios/=0) stop 'geofold_seams:: error allocating allbuttons.'

  ! Test read and write a barrel structure from .pdb and .sas files
  call geofold_seams_read (filename)
  call geofold_seams_write ()

  ! Test get all buttons
    do i=1, size (barrels_array)
      !! diagnostic
      write(*,*) "size (barrels_array)=" , size (barrels_array), &
           " i=",i," nSeams=",barrels_array(i)%nSeams
      do iseam=1, barrels_array(i)%nSeams
      call getAllButtons (barrels_array(i), iseam, allButtons, nButtons)
      print *, allButtons (1:nButtons)
    enddo
  enddo
    if(allocated(allButtons)) deallocate(allButtons)
endsubroutine test
!----------------------------------------------------------------------------
! Return the points which delimit the region (A point: minX, minY, maxX, maxY)
! The file as the format " N : X1 Y1 X2 Y2 ...."
! N: Number of values in the list of coordinates
! X1, Y1: coordintates
!----------------------------------------------------------------------------
! Modified: C.Bystroff Mon Aug 12 17:42:18 EDT 2013
!   Strictly keyworded input.
!----------------------------------------------------------------------------
subroutine geofold_seams_read(seamsFilename)
  implicit none
  character (*), intent(in)     :: seamsFilename
  type (seam_type)                             :: seam
  type (button_type)                           :: button
  character (len=500)                          :: aline, idStr*9
  character (len=MAXRES)                       :: flag
  integer :: nBarrels, nSeams, nButtons, id, k, i, j, ios, iunit, ires

    if (allocated(barrels_array)) then
      deallocate(barrels_array)
    endif
    iunit = pickunit(99)
  open(iunit,file=seamsFilename,status="old",form="formatted")
    i=0; j=0; k=0; ios=0; nSeams=0; nButtons=0;
  do
    read(iunit,'(a)',iostat=ios) aline
    if (ios /= 0) exit
    if (aline(1:1)=="#") cycle ! a comment
        if (aline(1:9)=="NBARRELS ") then   !! 1
      read(aline(9:),*) nBarrels ! read line ">>>>>>" nbb
      if (.not. allocated(barrels_array)) allocate (barrels_array(nBarrels))
          write(0,*) nBarrels," barrels allocated"
          if (nBarrels==0) exit
          k = 0  !! barrel number
          i = 0  !! seam number
        elseif (aline(1:7)=="BARREL ") then  !! 2
          !if (nSeams/=0.and.i/=nSeams) write(*,*) "WARNING: wrong number of seams read for barrel",&
          !                             k," Expected nSeams=",nSeams," got",i
          !if (i/=0.and.j/=0.and.nButtons/=0.and.j/=nButtons) &
          !      write(*,*) "WARNING: wrong number of buttons read for barrel",&
          !                 k," seam ",i," Expected",nButtons," got",j
          k = k + 1
          i = 0
          j = 0  !! button number
        read(aline(8:), *) id, nSeams
          if (k/=id) write(*,*) "WARNING: seam id number out of order.",k,id
      barrels_array(k)%nSeams = nSeams
        elseif (aline(1:5)=="SEAM ") then  !! 3
          !if (nButtons/=0.and.j/=nButtons) write(*,*) "WARNING: wrong number of buttons for barrel",&
          !                             k," seam ",i," Expected",nButtons," got",j
          i = i + 1
          j = 0
      read(aline(6:),*,iostat=ios) barrels_array(k)%seams(i)%id, barrels_array(k)%seams(i)%nButtons, &
                            barrels_array(k)%seams(i)%energy, barrels_array(k)%seams(i)%segments
      if(ios /= 0) write(0,*) aline(6:)
      nButtons = barrels_array(k)%seams(i)%nButtons
        elseif (aline(1:7)=="BUTTON ") then  !! 4
          j = j + 1
      read(aline(7:),*) barrels_array(k)%seams(i)%buttons(j)%residue, &
                            barrels_array(k)%seams(i)%buttons(j)%energyBeta1, &
                            barrels_array(k)%seams(i)%buttons(j)%energyBeta2
        elseif (aline(1:7)=="U1FLAG ") then  !! 4
           flag = trim(adjustl(aline(8:)))
           do ires=1,geofold_nres
              barrels_array(k)%seams(i)%u1flag(ires) = flag(ires:ires)
           enddo
        elseif (aline(1:7)=="U2FLAG ") then  !! 4
           flag = trim(adjustl(aline(8:)))
           do ires=1,geofold_nres
              barrels_array(k)%seams(i)%u2flag(ires) = flag(ires:ires)
           enddo
        elseif (aline(1:7)=="ENERGY ") then  !! This energy replaces the one in the SEAM line
      read(aline(8:),*) barrels_array(k)%seams(i)%energy
        else
          cycle
    endif
  enddo
    !if (nSeams/=0.and.i/=nSeams) write(*,*) "WARNING: wrong number of seams read for a barrel",&
    !                                   k," Expected nSeams=",nSeams," got",i
    !if (nButtons/=0.and.j/=nButtons) write(*,*) "WARNING: wrong number of buttons read for a barrel",&
    !                                   k," seam ",i," Expected",nButtons," got",j
  close(iunit)
endsubroutine
!-------------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! Get all buttons for one seam
!-----------------------------------------------------------------------
subroutine getAllButtons (barrel, iseam, allButtons, nButtons)
  implicit none
  type(barrel_type), intent(in)     :: barrel
    integer, intent(in)                   :: iseam
  integer, intent (out)                 :: allButtons (:), nButtons
  type(seam_type)                 :: seam ! cycle data type record
  type(button_type), target     :: button
  integer                       :: i, nSeams

    allButtons = 0
    nButtons = 0
  nSeams = barrel%nSeams
    i = abs(iseam)
    if (i==0) return
    if (i>nSeams) then
      write(*,*) 'iseam =',iseam,' nSeams=',nSeams
      stop 'geofold_seams:: Error in subroutine getAllButtons: iseam out of range'
    endif
  seam = barrel%seams(i)
  nButtons = seam%nButtons
    allButtons(1:nButtons) = seam%buttons(1:nButtons)%residue
endsubroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
function isContactInSeam (contact, seam) result (isInSeam)
  integer, intent (in)                 :: contact (2)
  type (seam_type), intent (in)        :: seam
  logical                              :: isInSeam
  integer               :: xres, yres, xMin, xMax, yMin, yMax

  xMin=seam%segments(1)
  xMax=seam%segments(2)
  yMin=seam%segments(3)
  yMax=seam%segments(4)

  isInSeam = .false.
  yres = contact (1)
  xres = contact (2)
  if (xMin<=xres.and.xres<=xMax .or. yMin<=yres.and.yres<=yMax) then
    isInSeam = .true.
    write (*, "(A4, 7I3)") ,">>> ", seam%id, xres, yres, xMin, xMax, yMin, yMax
  endif

endfunction
!--------------------------------------------------------------------
! Get the Total energy (seam + buttons) for a seam "iSeam"
! Returns the total energy and the side with low energy (beta1 or beta2)
!--------------------------------------------------------------------
!function getEnergySeam (contactsFilename, nResidues, iBarrel, iSeam, side, iflag) result(energySeam)
! implicit none
!     character (*), intent (in)       :: contactsFilename
!   integer                          :: nResidues, iBarrel, iSeam
!   integer, intent (out)            :: side
!
!   !CHARACTER, dimension(1600), intent (in) :: iflags     ! flags that define the intermediate
!   character, intent (in), optional :: iflag (:)       ! Unfolded residue not taken into account
!     !type (intermediate), POINTER         :: f
!   real                             :: energySeam, energyBeta1, energyBeta2
!   real, pointer, save          :: energyMatrix (:,:)
!   integer                          :: k, contact (2), c1, c2
!   type(seam_type), target          :: seam
!   integer                          :: nButtons, residue, x1,x2,y1,y2
!        logical,save          :: firstpass=.true.
!
! if (firstpass) then
!   if (associated(energyMatrix)) deallocate(energyMatrix)
!      nullify(energyMatrix)
!      allocate (energyMatrix (nResidues, nResidues))
!      call loadContactsToEnergyMatrix (contactsFilename, energyMatrix, nResidues)
!    else
!       !! diagnostic
!       !write(*,*) "In getEnergySeam, not firstpass. energyMatrix(10,10:20)=",&
!       !           energyMatrix(10,10:20)
!    endif
!    firstpass = .false.
!
!    !! barrels_array is global
! seam = barrels_array(iBarrel)%seams(iSeam)
! energySeam = 0.0
! do k=1, seam%n
!   contact = seam%x(1:2, k)
!   c1 = contact(1)
!   c2 = contact(2)
!   if (iflag(c1) == '.' .or. iflag(c2) == '.') cycle
!   energySeam = energySeam + energyMatrix(contact(1), contact(2))
! enddo
!
! !! Get energy of seam's buttons
! y1 = seam%segments (1) ! Beta 1
! y2 = seam%segments (2) ! Beta 1
! x1 = seam%segments (3) ! Beta 2
! x2 = seam%segments (4) ! Beta 2
!
! energyBeta1 = 0
! energyBeta2 = 0
! nButtons = seam%nButtons
! do k=1, nButtons
!   residue = seam%buttons(k)%residue
!   if (iflag (residue) == '.') cycle !! broken contact
!   energyBeta1 = energyBeta1 + getEnergyButtonBeta (residue, y1, y2, energyMatrix, iflag)
!   energyBeta2 = energyBeta2 + getEnergyButtonBeta (residue, x1, x2, energyMatrix, iflag)
! enddo
!
! energySeam = energySeam + MINVAL((/energyBeta1, energyBeta2/))
! if (energyBeta1 < energyBeta2) then
!   side = -1
! else
!   side = 1
! endif
!endfunction
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Calculate the energy for the residue with the beta (coords x1, x2)
! It uses the energy matrix, and an optional iflags array are with  marked broken contacts
!------------------------------------------------------------------------------
!function getEnergyButtonBeta (residue, x1, x2, energyMatrix, f) result (energy)
function getEnergyButtonBeta (residue, x1, x2, energyMatrix, iflag) result (energy)
  implicit none
  integer, intent (in)              :: residue, x1, x2 ! coordinates beta
  real, intent (in)           :: energyMatrix (:,:)
  character, intent (in), optional :: iflag (:)       ! Unfolded residue not taken into account
  !type (intermediate), POINTER, optional         :: f
  integer, allocatable              :: allContacts (:)! residue forming contacts with input "residue"
  real                              :: energy
  integer                           :: contact, i, j

  call getContactsResidue (residue, energyMatrix, allContacts)
  energy = 0
  do i=1, size (allContacts)
    contact = allContacts (i)   ! get contact

    if (contact < x1 .or. contact > x2) cycle  ! not in beta
    if (present (iflag) .and. iflag (contact) == '.') cycle  ! Broken

    energy = energy + energyMatrix (residue, contact)
  enddo
endfunction

!------------------------------------------------------------------------------
! Get all contacts for an residue acid
! OPTIMIZE by reading or receiving directly the contact list
!------------------------------------------------------------------------------
subroutine getContactsResidue (residue, energyMatrix, allContacts)
  integer, intent (in)               :: residue
  real, intent (in)            :: energyMatrix (:,:)
  integer, allocatable, intent (out) :: allContacts (:)
  integer                            :: allContactsTmp (500), nContacts

  n = size (energyMatrix, 1)
  !call createSeq (contactsSeq)
  nContacts=0
  do i=1, n
    if (energyMatrix (residue, i) /= 0 ) then
      nContacts = nContacts + 1
      allContactsTmp  (nContacts) = i
    endif
    !if (energyMatrix (residue, i) /= 0 ) call appendSeq (contactsSeq, i)
  enddo

  allocate (allContacts (nContacts))
  allContacts (1:nContacts) = allContactsTmp(1:nContacts)
endsubroutine

!------------------------------------------------------------------------------
! Load contacts energies to a matrix
! Read the energies from file and create and fill a matrix with them
!------------------------------------------------------------------------------
subroutine loadContactsToEnergyMatrix (contactsFilename, energyMatrix, nResidues)
  implicit none
    character (*), intent (in)    :: contactsFilename
    real, pointer           :: energyMatrix (:,:)
    integer, intent (in)        :: nResidues
    character (len=50)        :: aline
    integer             :: aa1, aa2, ios,iunit
    real                :: energy, x

  !if (associated(energyMatrix)) deallocate(energyMatrix)
    !iunit = pickunit(99)
    !allocate (energyMatrix (nResidues, nResidues))
  energyMatrix = 0
  open(iunit,file=contactsFilename,status="old",form="formatted")
  do
    read(iunit,'(a)',iostat=ios) aline
    if (ios /= 0) exit
    if (index (aline, "#") /= 0) cycle ! a comment

    read(aline, *), aa1, aa2, energy, x

    energyMatrix (aa1, aa2) = energy
    energyMatrix (aa2, aa1) = energy
  enddo
  close(iunit)
endsubroutine

!----------------------------------------------------------------------------
! Print the all the beta barrels from the barrels_array
! Only used to test visually the geofold_seams_read
!----------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Print the barrel information to a file
!---------------------------------------------------------------------------
subroutine geofold_seams_write(ounit)
    !! old name writeBarrels
  implicit none
    integer,optional :: ounit
  integer        :: nBarrels, nSeams, i, j, k, nButtons
  type(seam_type), pointer      :: seam ! cycle data type record
  type(barrel_type), pointer      :: barrel
  type(button_type), pointer      :: button
    integer :: iunit=6, ires

    if (present(ounit)) iunit=ounit
    if (.not.allocated(barrels_array)) &
    stop 'geofold_seams.f90:: ERROR in geofold_seams_write -- global barrels_array not allocated.'
  nBarrels = size(barrels_array)

  write (iunit,"(A)") "# Output from geofold_seams.f90 :: geofold_seams_write()"
  write (iunit,"(A)") "# NBARRELS  <Number of Barrels>"
  write (iunit,"(A)") "# BARREL :  <id>  <Number of Seams>"
  write (iunit,"(A)") &
    "# SEAM :    <id>  <Number of Buttons> <Total Energy> <x1-beta1> <x2-beta1> <x1-beta2> <x2-beta2>"
  write (iunit,"(A)") "# BUTTON :  <id>  <Energy-Beta1> <Energy-beta2>"
  write (iunit,"(A,I5)") "NBARRELS ", nBarrels
  do i=1, nBarrels
    barrel => barrels_array(i)
    nSeams = barrel%nSeams
    write(iunit, "(A,I5,A,I5)") "BARREL   ", i, " ", nSeams
    do j=1, nSeams
      seam => barrel%seams(j)
      write(iunit, "(A,I5,I5,f10.3,4I5)"), "SEAM     ", &
             j, seam%nButtons, seam%energy, seam%segments
            !! Note: %id is not used, but j would be useful.
      nButtons = seam%nButtons
      do k=1, nButtons
        button => seam%buttons (k)
        write(iunit, "(A, I5, 2f10.3)") &
          "BUTTON   ", button%residue, button%energyBeta1, button%energyBeta2
      enddo
        write(iunit, "(A,1600a1)") &
        "U1FLAG ", (seam%u1flag(ires),ires=1,geofold_nres)
        write(iunit, "(A,1600a1)") &
        "U2FLAG ", (seam%u2flag(ires),ires=1,geofold_nres)
    enddo
  enddo
endsubroutine geofold_seams_write
!!====================================================================================
!---------------------------------------------------------------------------------
! Read a contact list file. Only for testing
!---------------------------------------------------------------------------------
!subroutine read_contacts(contactsFilename, contacts_array)
! implicit none
! character (*), intent(in)     :: contactsFilename
! integer, allocatable, intent (out):: contacts_array (:,:)
! character (len=500)         :: aline
!    integer                           :: aax,aay,reason
!    real                              :: enr
!    character (len=50)                :: strLine
! integer                                  :: contactsTmp (2,2000), nContacts
!
! ! Read the file
! nContacts = 0
! open (10, file=contactsFilename)
! do
!   read (10, "(A)", iostat=reason) strLine
!   if (reason /=0) exit
!   if (index (strLine, "#") > 0) cycle
!
!   nContacts = nContacts + 1
!   read (strLine,*, iostat=reason) aax, aay, enr
!
!   contactsTmp (:, nContacts) = (/aax, aay/)
! enddo
! close(10)
!
! print *, ">>>", nContacts
! allocate (contacts_array(2, nContacts))
!
! contacts_array = contactsTmp(:,1:nContacts)
!endsubroutine  read_contacts

endmodule
