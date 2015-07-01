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

	! Button connecting to both b-sheets
	type :: button_type
		integer  :: residue  
		real     :: energyBeta1, energyBeta2
		integer  :: side   !! Which side of the button was unfolded
	endtype

	! Seam formed by two b-sheets
    type :: seam_type
        integer  :: id          ! Number of elements (Contacts by 2 [AA1, AA2], [AA1, AA4] : 4 elements)
        integer  :: n           ! Number of elements (Contacts by 2 [AA1, AA2], [AA1, AA4] : 4 elements)
        integer  :: x (2, 500)  ! Array of contacts (AA1, AA2)
        real     :: energy      ! Total energy of the seeam
        integer  :: segments (4)! Delimiting segments seam (beta1: y1,y2, beta2: x1,x2)
		integer  :: nButtons
		type (button_type) ::  buttons (500)  ! Array of buttons
    endtype 

	type :: barrel_type
		integer :: nSeams  ! number of barrels_array
		type (seam_type) :: seams (100)
	endtype

  type (barrel_type), allocatable  :: barrels_array (:)

  !call test() !To check the module

CONTAINS
!----------------------------------------------------------------------------
! Only to test the main function of the module
!----------------------------------------------------------------------------
subroutine test (filename)
	character (len=*)                 :: filename
	type (barrel_type), allocatable	  :: barrels_array (:)
	integer, allocatable              :: allButtons (:)
 
	! Test read and write a barrel structure from .pdb and .sas files
	call geofold_seams_read (filename, barrels_array)
	!call writeBarrels (barrels_array)

	! Test get all buttons
	call getAllButtons (barrels_array, allButtons)
  	do i=1, size (allButtons)
  		print *, allButtons (i)
 	enddo
!	call getarg (2, contactsFilename)
!	call read_contacts (contactsFilename, contacts_array)
!
!	print *, contacts_array (1,:)
!	print *, contacts_array (2,:)
!
!	call geofold_seams_find_bridges (barrels_array, contacts_array, bridges_array)
!
!	print *, "Bridges"
!	print *, bridges_array 
!	print *, bridges_array (1,:)
!	print *, bridges_array (2,:)
endsubroutine
!----------------------------------------------------------------------------
! Return the points which delimit the region (A point: minX, minY, maxX, maxY)
! The file as the format " N : X1 Y1 X2 Y2 ...."
! N: Number of values in the list of coordinates
! X1, Y1: coordintates
!----------------------------------------------------------------------------
subroutine geofold_seams_read (seamsFilename, barrels_array)
	implicit none
	character (*), intent(in)			:: seamsFilename
	type (barrel_type), allocatable, intent (out):: barrels_array (:)
	type (seam_type)                             :: seam
	type (button_type)                           :: button
	character (len=500)                          :: aline, idStr*9
	integer                                      :: nBarrels, nSeams, nButtons, id, k, i, j, ios

	open(99,file=seamsFilename,status="old",form="formatted")
	do
		! Read the line
		read(99,'(a)',iostat=ios) aline
		if (ios /= 0) exit
		if (aline(1:1)=="#") cycle ! a comment
		if (aline(1:6)/="BARREL") cycle 

		read(aline, "(A, I5)") idStr, nBarrels ! read line ">>>>>>" nbb
		if (.not. allocated (barrels_array)) allocate (barrels_array (nBarrels))
		! For each barrel 
		do k=1, nBarrels
			read(99,'(a)',iostat=ios) aline
			read(aline, "(A,I5,I5)") idStr, id, nSeams
			barrels_array (k)%nSeams = nSeams
			do i=1, nSeams
				read(99,'(a)',iostat=ios) aline
				read(aline, "(A,I5,I5,f10.3,4I5)") idStr, seam%id, seam%nButtons, seam%energy, seam%segments 
				nButtons = seam%nButtons
				do j=1, nButtons
					read(99,'(a)',iostat=ios) aline
					read(aline, "(A,I5,f10.3,f10.3)") idStr, button%residue, button%energyBeta1, button%energyBeta2
					seam%buttons (j) = button
				enddo
				barrels_array (k)%seams (i) = seam
			enddo
		enddo
	enddo
	close(99)
endsubroutine 
!-----------------------------------------------------------------------
! Get all buttons from  barrels
!-----------------------------------------------------------------------
subroutine getAllButtons (barrels_array, allButtons)
	implicit none
	type (barrel_type), intent (in)	:: barrels_array (:)
	integer, allocatable, intent (out):: allButtons (:)
	type(barrel_type), target		  :: barrel 
	type(seam_type), target		      :: seam ! cycle data type record 
	type(button_type), target		  :: button 
	integer			                  :: nBarrels, nSeams, i, j, k, nButtons
	integer                           :: buttonsTmp (500), nButtonsTmp

	nBarrels = size (barrels_array)
  	buttonsTmp = 0
	nButtonsTmp = 0
	do i=1, nBarrels
		barrel = barrels_array (i)
		nSeams = barrel%nSeams
		do j=1, nSeams
			seam = barrel%seams (j)
			nButtons = barrel%seams(j)%nButtons
			do k=1, nButtons
				button = seam%buttons (k)
				if (any (buttonsTmp == button%residue)) cycle
				nButtonsTmp = nButtonsTmp + 1
				buttonsTmp (nButtonsTmp) = button%residue
				!write (*, "(A8, I5, f10.3, f10.3, I5)") &
				!	"BUTTON: ", button%residue, button%energyBeta1, button%energyBeta2, 0
			enddo
		enddo
	enddo
	allocate (allButtons (nButtonsTmp))
  	allButtons = buttonsTmp (1:nButtonsTmp)
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
function getEnergySeam (barrelsArray, contactsFilename, nResidues, iBarrel, iSeam, side, f) result(energySeam) 
	implicit none
		type (barrel_type), intent (in)  :: barrelsArray (:)
	  	character (*), intent (in)	     :: contactsFilename
		integer                          :: nResidues, iBarrel, iSeam
		integer, intent (out)            :: side

		!CHARACTER, dimension(1600), intent (in) :: iflags     ! flags that define the intermediate
		!character, intent (in), optional :: iflags (:)       ! Unfolded residue not taken into account
	  	type (intermediate), POINTER         :: f
		real                             :: energySeam, energyBeta1, energyBeta2
		real, allocatable  			     :: energyMatrix (:,:)
		integer                          :: nContacts, k, contact (2), c1, c2
		type(seam_type), target          :: seam
		integer                          :: nButtons, residue, x1,x2,y1,y2

	call loadContactsToEnergyMatrix (contactsFilename, energyMatrix, nResidues)

	seam = barrelsArray(iBarrel)%seams(iSeam)
	nContacts = seam%n
	energySeam = 0.0
	do k=1, nContacts
		contact = seam%x (:, k)
		c1 = contact (1)
		c2 = contact (2)
		
		if (f%iflag (c1) == '.' .or. f%iflag (c2) == '.') cycle
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
		if (f%iflag (residue) == '.') cycle !! broken contact
		energyBeta1 = energyBeta1 + getEnergyButtonBeta (residue, y1, y2, energyMatrix, f)
		energyBeta2 = energyBeta2 + getEnergyButtonBeta (residue, x1, x2, energyMatrix, f)
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
function getEnergyButtonBeta (residue, x1, x2, energyMatrix, f) result (energy)
	implicit none
	integer, intent (in)              :: residue, x1, x2 ! coordinates beta
	real, intent (in)			      :: energyMatrix (:,:)
	!character, intent (in), optional  :: iflags              ! Unfolded residue not taken into account
	type (intermediate), POINTER, optional         :: f
	integer, allocatable              :: allContacts (:)! residue forming contacts with input "residue" 
	real                              :: energy
	integer                           :: contact, i, j

	call getContactsResidue (residue, energyMatrix, allContacts)
	energy = 0
	do i=1, size (allContacts)
		contact = allContacts (i)   ! get contact 

		if (contact < x1 .or. contact > x2) cycle  ! not in beta
		if (present (f) .and. f%iflag (contact) == '.') cycle  ! Broken

		energy = energy + energyMatrix (residue, contact)
	enddo
endfunction

!------------------------------------------------------------------------------
! Get all contacts for an residue acid
! OPTIMIZE by reading or receiving directly the contact list
!------------------------------------------------------------------------------
subroutine getContactsResidue (residue, energyMatrix, allContacts)
	integer, intent (in)               :: residue
	real, intent (in)			       :: energyMatrix (:,:)
	integer, allocatable, intent (out) :: allContacts (:)
	integer                            :: allContactsTmp (500), nContacts

	n = size (energyMatrix, 1)
	!call createSeq (contactsSeq)
	nContacts=0
	do i=1, n
		if (energyMatrix (residue, i) /= 0 ) then
			nContacts = nContacts + 1
			allContactsTmp	(nContacts) = i
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
	  character (*), intent (in)	  :: contactsFilename
	  real, allocatable, intent (out)			  :: energyMatrix (:,:)
	  integer, intent (in)			  :: nResidues
	  character (len=50)			  :: aline
	  integer						  :: aa1, aa2, ios
	  real							  :: energy, x
	  
	allocate (energyMatrix (nResidues, nResidues))
	energyMatrix = 0
	open(99,file=contactsFilename,status="old",form="formatted")
	do
		read(99,'(a)',iostat=ios) aline
		if (ios /= 0) exit
		if (index (aline, "#") /= 0) cycle ! a comment

		read(aline, *), aa1, aa2, energy, x

		energyMatrix (aa1, aa2) = energy
		energyMatrix (aa2, aa1) = energy
	enddo
	close(99)
endsubroutine

!----------------------------------------------------------------------------
! Print the all the beta barrels from the barrels_array
! Only used to test visually the geofold_seams_read
!----------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Print the entire sequence to the screen (or a file)
!---------------------------------------------------------------------------
subroutine writeBarrels (barrelsArray)
	implicit none
	type (barrel_type), allocatable, intent (in)	:: barrelsArray (:)
	integer			   :: nBarrels, nSeams, i, j, k, nButtons
	type(seam_type), target		  :: seam ! cycle data type record 
	type(barrel_type), target		  :: barrel 
	type(button_type), target		  :: button 

	nBarrels = size (barrelsArray)

	write (*,"(A)") "# Format"
	write (*,"(A)") "# NBARRELS  <Number of Barrels>"
	write (*,"(A)") "# BARREL :  <id>  <Number of Seams>"
	write (*,"(A)") "# SEAM :    <id>  <Number of Buttons> <Total Energy> <x1-beta1> <x2-beta1> <x1-beta2> <x2-beta2>"
	write (*,"(A)") "# BUTTON :  <id>  <Energy-Beta1> <Energy-beta2>"

	write (*,"(A8,I5)") "NBARRELS ", nBarrels
	do i=1, nBarrels
		barrel = barrelsArray (i)
		nSeams = barrel%nSeams
		write (*, '("BARREL ",I5," ",I5)') i, nSeams
		do j=1, nSeams
			seam = barrel%seams (j)
			write (*, '("SEAM ",,I5,I5,f10.3,4I5)')  &
					   seam%id, seam%nButtons, seam%energy, seam%segments
			nButtons = barrel%seams(j)%nButtons
			do k=1, nButtons
				button = barrel%seams(j)%buttons (k)
				write (*, '("BUTTON ", I5, 2f10.3)') &
					button%residue, button%energyBeta1, button%energyBeta2
			enddo
		enddo
	enddo
endsubroutine


subroutine printBetaBarrelOld (barrels_array)
	implicit none
	type (barrel_type), intent (in)	:: barrels_array (:)
	type (barrel_type)				:: bb
	type (seam_type)				   :: bs
	integer							 :: nbb, nbs, k, i

	print *, ">>> ", "Printing BetaBarrel"
	nbb = size (barrels_array)
	print *, ">>> Num of Barrels: ", nbb
	do k=1, nbb
		bb = barrels_array (k)
		nbs = bb%nSeams
		print *, ">>> Num of seams: ", nbs
		do i=1, nbs
			bs = bb%seams (i)
			write (*,"(I5, A, f8.3, A, 4I5)"),bs%id, " ", bs%energy, " ", bs%segments
		enddo
	enddo
endsubroutine
!!====================================================================================
!---------------------------------------------------------------------------------
! Read a contact list file. Only for testing
!---------------------------------------------------------------------------------
subroutine read_contacts(contactsFilename, contacts_array)
	implicit none
	character (*), intent(in)			:: contactsFilename
	integer, allocatable, intent (out):: contacts_array (:,:)
	character (len=500)					:: aline
    integer                           :: aax,aay,reason
    real                              :: enr
    character (len=50)                :: strLine
	integer                                  :: contactsTmp (2,2000), nContacts

	! Read the file
	nContacts = 0
	open (10, file=contactsFilename)
	do
		read (10, "(A)", iostat=reason) strLine
		if (reason /=0) exit
		if (index (strLine, "#") > 0) cycle

		nContacts = nContacts + 1
		read (strLine,*, iostat=reason) aax, aay, enr

		contactsTmp (:, nContacts) = (/aax, aay/)
	enddo
	close(10)

	print *, ">>>", nContacts
	allocate (contacts_array (2, nContacts))

	contacts_array = contactsTmp (:,1:nContacts)
endsubroutine 

endmodule
