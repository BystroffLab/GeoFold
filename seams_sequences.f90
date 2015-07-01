module seams_sequences
    use seams_sequence_generic
    !-------------------------------------------------------------------------
    ! Structures used by the sequences module: 
    !               list of betasheets (region), list of contacts (int2d)
    !-------------------------------------------------------------------------
	
	! Button connecting to both b-sheets
	type :: button_type
		integer  :: residue  
		real     :: energyBeta1, energyBeta2
	endtype

	! Seam formed by two b-sheets
    type :: seam_data
        integer  :: id           ! Number of elements (Contacts by 2 [AA1, AA2], [AA1, AA4] : 4 elements)
        integer  :: n           ! Number of elements (Contacts by 2 [AA1, AA2], [AA1, AA4] : 4 elements)
        integer  :: x (2, 500)  ! Array of contacts (AA1, AA2)
        real     :: energy      ! Total energy of the bsheet
        integer  :: segments (4)  ! Delimiting segments bsheet (minX, minY, maxX, maxY)
		integer  :: nButtons
		type (button_type) ::  buttons (500)  ! Array of buttons
    endtype 

    ! A container for storing region pointers
    type :: seam_pointer
        type(seam_data), pointer :: p
    endtype 

    !-------------------------------------------------------------------------
    ! General names of sequence operations (overloaded for regions and int2D)
    !-------------------------------------------------------------------------
    ! Create a empty (null) sequence
    interface createSeq
        module procedure createSeqReg 
        module procedure createSeq1D
        module procedure createSeq2i
    endinterface 

    ! Append an element to the sequence
    interface appendSeq
        module procedure appendSeqReg
        module procedure appendSeq1D
        module procedure appendSeq2i
    endinterface
    
    interface putAtSeq
        module procedure putAtSeqReg
    endinterface
    
    ! Get the "k" element of the sequence
    interface getAtSeq
        module procedure getAtSeqReg
        module procedure getAtSeq1D
        module procedure getAtSeq2i
    endinterface

    ! Get the number of elements of the sequence
    interface lengthSeq
        module procedure lengthSeqReg
        module procedure lengthSeq1D
        module procedure lengthSeq2i
    endinterface

    ! Checks if the sequence is empty
    interface isEmptySeq
        module procedure isEmptySeq1D
        module procedure isEmptySeq2i
    endinterface

    ! Get the "k" element and remove it from the sequence
    interface popSeq
        module procedure popSeq1D
        module procedure popSeq2i
    endinterface
    
    ! Check (True/False) if an element exist in the sequence
    interface existsSeq
        module procedure existsSeq1D
        module procedure existsSeq2i
    endinterface
    
    ! Extend the sequence with the elements of another sequence
    interface extendsSeq
        module procedure extendsSeq1D
        module procedure extendsSeq2i
    endinterface

    ! Reverse the sequence
    interface reverseSeq
        module procedure reverseSeqReg
        module procedure reverseSeq1D
        module procedure reverseSeq2i
    endinterface

    ! Print to screen the elements of the sequence
    interface printSeq
        module procedure printSeqReg
        module procedure printSeq1D
        module procedure printSeq2i
    endinterface

    ! Write to file the elements of the sequence
    interface writeSeq
        module procedure writeSeqReg
    endinterface

CONTAINS
    !-------------------------------------------------------------------------
	! Subroutine only for testing parts of the module
    !-------------------------------------------------------------------------
	subroutine testArray1D 
		integer, allocatable :: arr1d (:), a2 (:)
		integer element

		print *, "Testing array 1D"
		print *, "Creating..."
		call createSeq (arr1d)
		call createSeq (a2)

		print *, "Appending.."
		call appendSeq (arr1d, 11)
		call appendSeq (arr1d, 22)

		call appendSeq (a2, 111)
		call appendSeq (a2, 222)

		print *, "Printing..."
		call printSeq (arr1d)

		print *, "Poping 3..."
		element = popSeq (arr1d,3)
		print *, "Element: ", element

		print *, "Printing..."
		call printSeq (arr1d)

		print *, "Expanding "
		call extendsSeq (arr1d, a2)
		call printSeq (arr1d)

		print *, "Reversing "
		call reverseSeq (arr1d)
		call printSeq (arr1d)

		print *, "Exists 44?"
		print*,  existsSeq (arr1d, 44)
	endsubroutine

!--Sequence regions of contacts---------------------------------------------
    !public :: createSeq, appendSeq, printSeq, getAtSeq, lengthSeq
    !public :: popSeq, existsSeq, isEmptySeq, extendsSeq

    !public :: seam_data
    !public :: seam_pointer

    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    ! A sequence derived type for storing region (contact maps) data
    ! A wraper for the generic sequence (in this module)
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Create a null sequence
    !------------------------------------------------------------------------------
    subroutine createSeqReg (seq)
        type(sequence_node), pointer :: seq 

        call createSeqG (seq)
    endsubroutine
        
    !------------------------------------------------------------------------------
    ! Append the region (array) to the sequence seq
    !------------------------------------------------------------------------------
    subroutine appendSeqReg (seq, array)
        type(sequence_node), pointer :: seq 
        integer, intent (in)   :: array (:,:)
        type(seam_pointer) :: ptr
        integer :: n, sh (2)

        sh = shape (array)
        n = sh (2)

        allocate (ptr%p)
        ptr%p%n = n
        ptr%p%x = 0
        ptr%p%x (:,1:n) =  array(:,1:n)
        ptr%p%energy = -1.0

        minX = minval (array (1, 1:n))
        minY = minval (array (2, 1:n))
        maxX = maxval (array (1, 1:n))
        maxY = maxval (array (2, 1:n))

        ptr%p%segments = (/minY, maxY, minX,  maxX/)

        call appendSeqG (seq, DATA=transfer(ptr, sequence_data))
    endsubroutine

    !------------------------------------------------------------------------------
    ! Append the region (array) to the sequence seq
    !------------------------------------------------------------------------------
    subroutine putAtSeqReg (seq, pos, reg)
        type(sequence_node), pointer    :: seq 
        integer, intent (in)            :: pos
        type(seam_data), target       :: reg
        type(seam_pointer)            :: rPtr
        integer :: n, sh (2)

        !regPtr%p = reg 
        allocate (rPtr%p)
        rPtr%p%n      = reg%n
        rPtr%p%x      = reg%x
        rPtr%p%energy = reg%energy
        rPtr%p%segments = reg%segments

        call putAtSeqG (seq, pos, DATA=transfer(rPtr, sequence_data))
    endsubroutine

    !---------------------------------------------------------------------------
    ! Get the REGION of the sequence at position "pos"
    !---------------------------------------------------------------------------
    function getAtSeqReg (seq, pos) result (reg)
        type(sequence_node), pointer     :: seq
        integer, intent (in)           :: pos
        type(seam_data), target           :: reg
        type(seam_pointer)                 :: ptr

        ptr = transfer(getAtSeqG(seq, pos), ptr)
        reg = ptr%p
     endfunction

    !---------------------------------------------------------------------------
    ! Reverse the sequence
    !---------------------------------------------------------------------------
    subroutine reverseSeqReg (seq) 
        type(sequence_node), pointer  :: seq

        call reverseSeqG (seq)
    endsubroutine
     
    !---------------------------------------------------------------------------
    ! Get the Length (number of elements) of the sequence
    !---------------------------------------------------------------------------
    function lengthSeqReg (seq) result (length)
        type(sequence_node), pointer  :: seq
        integer                    :: length

        length = lengthSeqG (seq)
    endfunction
     
    !---------------------------------------------------------------------------
    ! Return the DATA stored in the node SELF
    !---------------------------------------------------------------------------
    function getSeqData (seq) result (reg)
        type(sequence_node), pointer :: seq 
        type(seam_data), target :: reg
        type(seam_pointer) :: ptr
        
        ptr = transfer(getSeqG(seq), ptr)
        reg = ptr%p
    endfunction

    !---------------------------------------------------------------------------
    ! Print the entire sequence to the screen (or a file)
    ! Uses the Generic Print with the callback function printOneReg
    !---------------------------------------------------------------------------
   subroutine printOneReg (data)
       implicit none
         integer                     :: data (:)
         type(seam_pointer)        :: ptr
         type(seam_data), target   :: reg
       
       ptr = transfer(data, ptr)
       reg = ptr%p
       write (*, "(I5, A, f8.3, A, 4I5)", advance='no'), &
       reg%n, " : ", reg%energy, " : ", reg%segments
       write (*,*), ""
   endsubroutine 

   subroutine printSeqReg (seq)
        type(sequence_node), pointer :: seq

        call printSeqG (seq, printOneReg)
    endsubroutine

    !---------------------------------------------------------------------------
    ! Write the entire sequence to file
    !---------------------------------------------------------------------------
    subroutine writeSeqReg (seq, filename)
        implicit none
          type(sequence_node), pointer :: seq
          character (*), optional, intent (in)  :: filename
          type(sequence_node), pointer :: pSelf
          type(seam_data), target :: reg
          integer                   :: array (2,500), j, len, n, ios, dunit
          character (50)            :: fmt

        dunit = 6
        if (present (filename)) then
            dunit = 99
            open(dunit,file=filename,status="replace",form='formatted',iostat=ios)
            if (ios/=0) stop 'Error in printSeqReg: writen sequence of regions to a file'
        endif

        pSelf => seq
        
        write (dunit, *) "# List of 4 points delimiting regions (square) of beta sheets"
        write (dunit, *) "# FORMAT: 'number of contacts : Energy : minX, minY, maxX, maxY'"
        do while (associated (pSelf))
            reg = getSeqData (pSelf)
            len = reg%n
            array (:,1:len) = reg%x(:,1:len)

            n = len*2
            write (dunit, '(I5, A, f8.3, A)', advance='no'),n," : ", reg%energy, " : "
            fmt = "("//trim (catSS (strII (n),"(I5)"))//")"
            !write (dunit, fmt) (array(1:2,j), j=1,len)
            write (dunit, fmt) reg%segments

            pSelf => nextSeqG (pSelf)
        enddo

        if (present (filename)) close (dunit)
    endsubroutine
    
    !----------------------------------------------------------------------------
	! Maybe OBSOLETE and UNUSED
    ! Return the points which delimit the region (A point: minX, minY, maxX, maxY)
    ! The file as the format " N : X1 Y1 X2 Y2 ...."
    ! N: Number of values in the list of coordinates
    ! X1, Y1: coordintates
    !----------------------------------------------------------------------------
    subroutine loadPointsSeqReg (regFilename, pointsRegion)
        implicit none
        character (*), intent(in)           :: regFilename
        integer, allocatable, intent (out)  :: pointsRegion (:,:)
        integer                             :: n, minX, minY, maxX, maxY, ios
        character (len=80)                  :: aline
        character (len=3)                   :: sep
        integer, allocatable                :: minSeq (:,:), maxSeq (:,:)

        call createSeq (minSeq)  ! Sequence 2D for min Points
        call createSeq (maxSeq)  ! Sequence 2D for max Points

        open(99,file=regFilename,status="old",form="formatted")
        do
            read(99,'(a)',iostat=ios) aline
            if (ios /= 0) exit

            if (index (aline, "#") /= 0) cycle ! a comment

            read(aline, "(I5,A,4I5)"), n, sep, minX, minY, maxX, maxY 

            call appendSeq (minSeq, (/minX, minY/))
            call appendSeq (maxSeq, (/maxX, maxY/))
            
        enddo
        close(99)

        n = lengthSeq (minSeq)
        allocate (pointsRegion (4,n))
        pointsRegion (1:2,:) = minSeq
        pointsRegion (3:4,:) = maxSeq
    endsubroutine 

    !----------------------------------------------------------------------------
    ! Get the string from the integer "n"
    !----------------------------------------------------------------------------
    function strII (n) result (str)
        implicit none
        integer :: n
        character (len=100) :: s1
        character (len=100) :: str
        
        write (s1, "(i10)") n
        write (str, "(A)")trim(adjustl(s1))
    end function

    !----------------------------------------------------------------------------
    ! Paste the string "s1" and "s2"
    !----------------------------------------------------------------------------
    function catSS (s1, s2) result (str)
        integer :: n
        character (len=*):: s1, s2
        character (len=100) :: str
        
        write (str, *)trim (s1)//trim (s2)
    endfunction

!--Sequence Integer Array 1 D-----------------------------------------------
    !---------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Sequence for integer array 1 dimensional
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine createSeq1D (array)
        implicit none
        integer,  allocatable, intent (inout) :: array (:)

        if (allocated (array)) then
            deallocate (array)
        endif

        allocate (array (0))
    endsubroutine
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine appendSeq1D (array, element)
        implicit none
        integer,  allocatable, intent (inout) :: array (:)
        integer   :: element, n
        integer,  allocatable  :: tmp (:)

        if (.not. allocated (array)) Stop "Error in appendSeq1D: not allocated array"

        if (isEmptySeq1D (array)) then
            deallocate (array)
            allocate (array(1))
            array (1) = element
        else
            n = lengthSeq1D (array)
            allocate (tmp (n))
            tmp = array (1:n)

            deallocate (array)
            allocate (array (n+1))

            array (1:n) = tmp
            array (n+1) = element
            deallocate (tmp)
        endif
    endsubroutine

    ! Return true if the "array" is empty.
    logical function isEmptySeq1D  (array) result (value)
        implicit none
        integer, allocatable, intent (in) :: array (:)

        value = (lengthSeq1D (array) == 0)
    end function

    ! Returns the element "k" from the 2d array
    function getAtSeq1D (array, k) result (value)
        integer,  intent (in) :: array (:)
        integer               :: k, value 

        value = array (k)
    endfunction

    ! Returns and remove the element "k" from the 2d array
    function popSeq1D (array, k) result (value)
        integer,  allocatable, intent (inout) :: array (:)
        integer,  allocatable                :: tmp (:)
        integer     :: n, k, value 

        n = lengthSeq1D (array)

        value = array (k)

        allocate (tmp (n-1))
        tmp =  reshape ( (/array(1:k-1),array(k+1:n)/),  shape (tmp))

        deallocate (array)
        allocate (array (n-1))
        array = tmp
        deallocate (tmp)
    endfunction
    
    function lengthSeq1D (array) result (value)
        implicit none
        integer,  intent (in) :: array (:)
        integer  :: value
        
        value = size (array)
    endfunction

	subroutine reverseSeq1D (array)
        implicit none
        integer,  intent (inout) :: array (:)
        integer, allocatable   :: arrayRev (:)
		integer                :: n, i
		
		n = lengthSeq1D(array)
		allocate (arrayRev (n))

		do i=1, n
			arrayRev (i) =  array (n-i+1)	
		enddo

		array = arrayRev 

		deallocate (arrayRev)
	endsubroutine

    logical function existsSeq1D (array, element) result (value)
        implicit none
        integer,  allocatable, intent (in) :: array (:)
        integer   :: element, n, i

        n = lengthSeq1D (array)
        value = .false.
        do i=1, n
            if (any (array (:) == element)) then
                value = .true.
                return
            endif
        enddo
    endfunction

    ! Extends (appends) the current sequence with another lists
    subroutine extendsSeq1D (array, arrayExt)
        implicit none
        integer,  allocatable, intent (inout) :: array (:)
        integer,  allocatable, intent (in) :: arrayExt (:)
        integer   :: n,k

        n = lengthSeq1D (arrayExt)
        do k=1, n
            call appendSeq1D (array, getAtSeq1D (arrayExt, k))
        enddo
    endsubroutine

    ! Print  the 2d "array" to screen plainly (x1,y1,x2,y2...)
    subroutine printSeq1D (array)
        implicit none
        integer, intent (in)  :: array (:)
        integer i,j,n, rows
        character (len=20) :: fmt
        fmt = "(I5)"
        n = lengthSeq1D (array)
        do i=1, n
            write (*, fmt, advance='no') array (i)
        enddo
		write (*,*) ""
    endsubroutine

!--Sequence Integer Array 2 D-----------------------------------------------
    !---------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Sequence for integer array 2 dimiensionales
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine createSeq2i (array)
        implicit none
        integer,  allocatable, dimension (:,:), intent (inout) :: array

        if (allocated (array)) then
            deallocate (array)
        endif

        allocate (array (2,0))
    endsubroutine

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine appendSeq2i (array, element)
        implicit none
        integer,  allocatable, dimension (:,:), intent (inout) :: array
        integer   :: element(2), n, s(2)
        integer,  allocatable, dimension (:,:)  :: tmp

        if (.not. allocated (array)) Stop "Error in appendSeq2i: not allocated array"

        if (isEmptySeq2i (array)) then
            deallocate (array)
            allocate (array(2,1))
            array (:,1) = element
        else
            n = lengthSeq2i (array)
            allocate (tmp (2, n))
            tmp = array (:,1:n)

            deallocate (array)
            allocate (array (2, n+1))

            array (:, 1:n) = tmp
            array (:, n+1) = element
            deallocate (tmp)
        endif
    endsubroutine


    ! Return true if the "array" is empty.
    logical function isEmptySeq2i  (array) result (value)
        implicit none
        integer, allocatable,  dimension (:,:), intent (in) :: array

        value = (lengthSeq2i (array) == 0)
    end function

    ! Returns the element "k" from the 2d array
    function getAtSeq2i (array, k) result (value)
        integer,  dimension (:,:), intent (in) :: array
        integer     :: k, value (2)

        value = array (:, k)
    endfunction

    ! Returns and remove the element "k" from the 2d array
    function popSeq2i (array, k) result (value)
        integer,  allocatable, dimension (:,:), intent (inout) :: array
        integer,  allocatable, dimension (:,:)                 :: tmp
        integer     :: n, k, value (2)

        n = lengthSeq2i (array)

        value = array (:, k)

        allocate (tmp (2, n-1))
        tmp =  reshape ( (/array(:,1:k-1),array(:,k+1:n)/),  shape (tmp))

        deallocate (array)
        allocate (array (2, n-1))
        array = tmp
        deallocate (tmp)
    endfunction
    
    function lengthSeq2i (array) result (value)
        implicit none
        integer,  intent (in) :: array (:,:)
        integer  :: sh (2), value
        
        sh = shape (array)
        value = sh (2)
    endfunction

	subroutine reverseSeq2i (array)
        implicit none
        integer,  intent (inout) :: array (:,:)
        integer, allocatable   :: arrayRev (:,:)
		integer                :: n, i
		
		n = lengthSeq2i(array)
		allocate (arrayRev (2,n))

		do i=1, n
			arrayRev (:,i) =  array (:,n-i+1)	
		enddo

		array = arrayRev 

		deallocate (arrayRev)
	endsubroutine

    logical function existsSeq2i (array, element) result (value)
        implicit none
        integer,  allocatable, intent (in) :: array (:,:)
        integer   :: element(2), s(2), n, i

        n = lengthSeq2i (array)
        value = .false.
        do i=1, n
            if (all (array (:,i) == element)) then
                value = .true.
                return
            endif
        enddo
    endfunction

    ! Extends (appends) the current sequence with another lists
    subroutine extendsSeq2i (array, arrayExt)
        implicit none
        integer,  allocatable, intent (inout) :: array (:,:)
        integer,  allocatable, intent (in) :: arrayExt (:,:)
        integer   :: n,k, element(2)

        n = lengthSeq2i (arrayExt)
        do k=1, n
            call appendSeq2i (array, getAtSeq2i (arrayExt, k))
        enddo
    endsubroutine

    ! Print  the 2d "array" to screen plainly (x1,y1,x2,y2...)
    subroutine printSeq2i (array)
        implicit none
        integer, intent (in)  :: array (:,:)
        integer i,j,n, rows
        character (len=20) :: fmt
        
    !        n = size (array)
    !        rows = nrows (array)
    !
    !        fmt = "("//trim (catS (strI (n),"(I5)"))//")"
    !        write (*, fmt) (array(1:2,j), j=1,rows)
    !
        fmt = "(I5,I5)"
        n = lengthSeq2i (array)
        do i=1, n
            write (*, fmt) array (:,i)
        enddo
    endsubroutine
endmodule

