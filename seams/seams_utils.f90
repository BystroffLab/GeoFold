!---------------------------------------------------------------------------
! Utilities for strings, file, and matrix handling
!---------------------------------------------------------------------------

module seams_utils
  public :: nrows, writeMatrix, str2real
CONTAINS
!------------------------ S Y S T E M ---------------------------------------
! system utilities
!----------------------------------------------------------------------------
  ! Execute a external command
  subroutine execute (command, outputArray) 
    implicit none
    character (*), intent (in) :: command
    character (len=200), allocatable :: outputArray (:)
    character (len=1000) :: tmpFilename,LName
    character (len=1000) :: newCmm
    character (len=200) :: tmpdir
    integer  :: i, stat
    
    call get_environment_variable("TMPDIR",tmpdir)
    if(tmpdir=="") tmpdir = "./tmp"
    call get_environment_variable("LNAME",LName)
    if(LName=="")LName="_"
    tmpFilename = trim(adjustl(tmpdir)) // "/tmp" // trim(adjustl(LName)) // ".tmp"
    write(0,*)"tmpFilename"
    write(0,*) tmpFilename
    newCmm = trim (command) // " > " // trim(adjustl(tmpFilename))
    write(0,*) newCmm
    call system (newCmm, stat)
        if (stat/=0) then
          write(0,*) "seams_utils.f90:: error on system execution, command=", trim(newCmm)
          stop
        endif

    call readFile (tmpFilename, outputArray)

    newCmm = "rm " // tmpFilename
    call system (newCmm, stat)
        if (stat/=0) then
          write(0,*) "seams_utils,f90:: error on system execution, command=", trim(newCmm)
          stop
        endif
    
  endsubroutine

  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  subroutine readFile (inputFile, outputArray)
    implicit none
    character (len=*), intent(in) :: inputFile 
    character (len=200), allocatable :: outputArray (:)
    integer, parameter :: maxrecs = 10000
    integer :: J, NR, ios
    character(LEN=1) :: junk
    character (len=200)    :: line

    !Determine total number of lines in file
    NR = 0
    OPEN(UNIT=1,FILE=inputFile)
    do J=1,maxrecs
      read(1,'(a)',iostat=ios) line
      if (ios /= 0) exit

      if (J == maxrecs) &
            stop "seams_utils.f90:: Error (readFile): Max number of line exceeded..."

      NR = NR + 1
    enddo

    rewind(1)
    !Now we can allocate data variables
    allocate (outputArray(NR))
    !Now read data into mydata
        outputArray = " "
    do J=1,NR
      read (1, '(a)') line
      outputArray(J) = trim(line)
    enddo
    close(1)
  endsubroutine

  !----------------------------------------------------------------------------
  ! F I L E S
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Get a free unit to open a file
  !----------------------------------------------------------------------------
  integer function pickUnit(iunit)
    implicit none
    integer,intent(in) :: iunit
    logical :: alreadyused
    integer :: ounit
    ounit = iunit
    inquire(unit=ounit,opened=alreadyused)
    do while (alreadyused)
      ounit = ounit + 1
      inquire(unit=ounit,opened=alreadyused)
    enddo
    pickunit = ounit
  end function pickunit

  !---------------------------r-------------------------------------------------
  ! Return the unit to open a file 
  !----------------------------------------------------------------------------
  function openFile (filename)
    character (len=*), intent(in) :: filename 
    integer, save :: dunit

    ! dunit = pickUnit(10)
    open (newunit=dunit, file=filename, iostat=ierr)
    if (ierr > 0 ) then
      print *, "Error ", ierr, " opening file ", filename
      Stop
    endif

    openFile = dunit
  end function openFile

  !----------------------------------------------------------------------------
  ! Return the number of lines of a file
  !----------------------------------------------------------------------------
  integer function lengthFile (filename)
    character (len=*), intent(in) :: filename 
    character(len=1000) :: aline
    integer  :: n, dunit, ioFlag

    n=0
    dunit = openFile (filename) 
    if ((dunit) < 0) then
      print *, "Error opening ", filename
      Stop
    endif

    open (dunit, file=filename)
    do
      read (dunit, *, iostat=ioFlag) aline

      if (ioFlag /=0) exit
      n = n+1
    enddo
    close (dunit)

    lengthFile = n
  endfunction

  !----------------------------------------------------------------------------
  ! M A T R I X
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! writeMatrix: 
  ! Write a matrix (NxN) structure to a file
  !----------------------------------------------------------------------------
  subroutine writeMatrix (matrix, filename)
    integer, dimension (:,:)  :: matrix
    character (len=*)     :: filename
    integer           :: i, j, n, dim (2)
    integer           :: dunit

    dunit = openFile (filename)
    !open (dunit, file=filename)
    dim = shape (matrix)

    do i=1, dim (1)
      do j=1, dim (2)
        write (unit=dunit, FMT="(I1)", advance="no"), matrix (j,i)
      enddo
      write (unit=dunit, FMT=*) ""
    enddo
    close(dunit)
  endsubroutine

  !----------------------------------------------------------------------------
  ! Read a NxN binary matrix from a file
  !----------------------------------------------------------------------------
  subroutine readMatrix (filename, matrix)
    implicit none
    character(len=*), intent (in) :: filename
    character(len=1000) :: line   ! Max number of residues of a protein
    integer       :: n, i, j,  io, dunit, value
    integer, allocatable, dimension (:,:) :: matrix
    character (len=10)  :: strFormat

    n = lengthFile (filename)
    call getStringFormat (n, strFormat)

    dunit = openFile (Filename)
    if (dunit < 0) Stop "Error opening file matrix"
    
    allocate (matrix (n, n))
    do i=1, n
      read (dunit, *, iostat=io) line
      do j=1, n
        read (line (j:j), *) value
        matrix (j,i) = value
      enddo
    end do
    close (dunit)
  endsubroutine readMatrix

  !----------------------------------------------------------------------------
  ! Return the number of columns of the 2d array
  !----------------------------------------------------------------------------
  function ncols (matrix) result (value)
      integer, dimension (:,:), intent (in) ::  matrix
      integer        :: value, shp(2)
      
      shp   = shape (matrix)
      value = shp (1)
  endfunction
  !----------------------------------------------------------------------------
  ! Return the number of rows of the 2d array
  !----------------------------------------------------------------------------
  function nrows (matrix) result (value)
      integer, dimension (:,:), intent (in) ::  matrix
      integer        :: value, shp(2)
      
      shp   = shape (matrix)
      value = shp (2)
  endfunction

  !------------------------------ S T R I N G S -------------------------------
  ! string utilities 
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !Paste the integer "n" with the string "s1" 
  !----------------------------------------------------------------------------
  function catI (n, s) result (str)
    integer :: n
    character (len=*):: s
    character (len=100) :: str
    
    write (str, *)trim (catS (int2str(n),s))
  endfunction

  !----------------------------------------------------------------------------
  ! Paste the string "s1" and "s2"
  !----------------------------------------------------------------------------
  function catS (s1, s2) result (str)
    integer :: n
    character (len=*):: s1, s2
    character (len=100) :: str
    
    write (str, *)trim (s1)//trim (s2)
  endfunction

  !----------------------------------------------------------------------------
  ! Get the string from the integer "n"
  !----------------------------------------------------------------------------
  function int2str (n) result (str)
    implicit none
    integer :: n
    character (len=100) :: s1
    character (len=100) :: str
    
    write (s1, "(i10)") n
    write (str, "(A)")trim(adjustl(s1))
  end function

  !----------------------------------------------------------------------------
  ! Get the integer from the string
  !----------------------------------------------------------------------------
  function str2int (str) result (n)
    implicit none
    character (*), intent (in)  :: str
    integer           :: n, ios

    read (str, '(I5)',iostat=ios) n
        if (ios/=0) then
          write(0,*) "seams_utils.f90:: str2int: cant convert to int: ", trim(str)
          stop
        endif
  end function
  !----------------------------------------------------------------------------
  ! Get the real from the string
  !----------------------------------------------------------------------------
  function str2real (str) result (n)
    implicit none
    character (*), intent (in)  :: str
    real            :: n

    read (str, '(f8.3)'), n
  endfunction
  
  !----------------------------------------------------------------------------
  ! Get the string format to write a number "n"
  !----------------------------------------------------------------------------
  subroutine getStringFormat (n, s1) 
    integer, intent (in) :: n
    character (len=10), intent (out) :: s1
    character (len=10) :: s2
    integer :: m
    
    write (s2, "(I5)") n
    write (s1, *) "(a",trim(adjustl(s2)),")"

  endsubroutine 

  !----------------------------------------------------------------------------
  ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
  ! the delimiters contained in the string 'delims'. Preceding a delimiter in
  ! 'str' by a backslash (\) makes this particular instance not a delimiter.
  ! The integer output variable nargs contains the number of arguments found.
  !----------------------------------------------------------------------------
  subroutine parse(str,delims,args,nargs)
    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:) :: args

    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
      args(i)=' '
    end do  
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    k=0

    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       call split(str,delims,args(nargs))
       call removebksl(args(nargs))
    end do   
    str=strsav
  endsubroutine
  !----------------------------------------------------------------------------
  ! Converts multiple spaces and tabs to single spaces; deletes control characters;
  ! removes initial spaces.
  !----------------------------------------------------------------------------
  subroutine compact(str)

  character(len=*):: str
  character(len=1):: ch
  character(len=len_trim(str)):: outstr

  str=adjustl(str)
  lenstr=len_trim(str)
  outstr=' '
  isp=0
  k=0

  do i=1,lenstr
    ch=str(i:i)
    ich=iachar(ch)
    
    select case(ich)
    
    case(9,32)     ! space or tab character
      if(isp==0) then
      k=k+1
      outstr(k:k)=' '
      end if
      isp=1
      
    case(33:)    ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0
      
    end select
    
  end do

  str=adjustl(outstr)

  end subroutine compact

  !-----------------------------------------------------------------------
  ! Routine finds the first instance of a character from 'delims' in the
  ! the string 'str'. The characters before the found delimiter are
  ! output in 'before'. The characters after the found delimiter are
  ! output in 'str'. The optional output character 'sep' contains the 
  ! found delimiter. A delimiter in 'str' is treated like an ordinary 
  ! character if it is preceded by a backslash (\). If the backslash 
  ! character is desired in 'str', then precede it with another backslash.
  !-----------------------------------------------------------------------

  subroutine split(str,delims,before,sep)

    character(len=*) :: str,delims,before
    character,optional :: sep
    logical :: pres
    character :: ch,cha

    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if(lenstr == 0) return    ! string str is empty
    k=0
    ibsl=0            ! backslash initially inactive
    before=' '
    do i=1,lenstr
     ch=str(i:i)
     if(ibsl == 1) then      ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
     end if
     if(ch == '\') then      ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
     end if
     ipos=index(delims,ch)       
     if(ipos == 0) then      ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
     end if
     if(ch /= ' ') then      ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
     end if
     cha=str(i+1:i+1)      ! character is a space delimiter
     iposa=index(delims,cha)
     if(iposa > 0) then      ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
     else
      str=str(i+1:)
      if(pres) sep=ch
      exit
     end if
    end do
    if(i >= lenstr) str=''
    str=adjustl(str)        ! remove initial spaces
    return
  endsubroutine 
  !-----------------------------------------------------------------------
  ! Removes backslash (\) characters. Double backslashes (\\) are replaced
  ! by a single backslash.
  !-----------------------------------------------------------------------
  subroutine removebksl(str)
    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0
    ibsl=0            ! backslash initially inactive

    do i=1,lenstr
    ch=str(i:i)
    if(ibsl == 1) then      ! backslash active
     k=k+1
     outstr(k:k)=ch
     ibsl=0
     cycle
    end if
    if(ch == '\') then      ! backslash with backslash inactive
     ibsl=1
     cycle
    end if
    k=k+1
    outstr(k:k)=ch        ! non-backslash with backslash inactive
    end do

    str=adjustl(outstr)
  endsubroutine 

  !-----------------------------------------------------------------------
  ! Removes spaces, tabs, and control characters in string str
  !-----------------------------------------------------------------------
  subroutine removesp(str)
    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0

    do i=1,lenstr
    ch=str(i:i)
    ich=iachar(ch)
    select case(ich)    
      case(0:32)  ! space, tab, or control character
         cycle     
      case(33:)  
      k=k+1
      outstr(k:k)=ch
    end select
    end do

    str=adjustl(outstr)

  end subroutine removesp

endmodule
