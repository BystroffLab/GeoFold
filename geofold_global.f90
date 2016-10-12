!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! geofold_global.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE geofold_global
  private
  !!--------- global types and parameters for GEOFOLD ------------!!
  !! INTEGER, PARAMETER :: maxres = 80 !maximum number of residues
  INTEGER, PARAMETER :: MAXRES=1600
  integer,parameter :: MAXCHAIN=26  !! Maximum number of separate chains.
  integer,save :: maxsplit=8
  real,parameter ::  FULSPLITDEPTH=4.  !! recursion depth limit for maximum splitting.
  integer,parameter :: MAXBARREL=8
  real,parameter :: MAXHINGEANG=30.0 ! atoms must hinge at least this much
  character,dimension(MAXCHAIN) :: chainletters=(/"A","B","C","D","E","F","G","H","I","J","K","L","M","N",&
                                                  "O","P","Q","R","S","T","U","V","W","X","Y","Z"/)
  TYPE intermediate
     INTEGER :: idnum                          ! unique number identifier of the intermediate state
     CHARACTER, dimension(maxres) :: iflag     ! flags that define the intermediate
     TYPE (intermediate), POINTER :: next       ! next intermediate
     INTEGER :: state                           ! breakflag, pivotflag, hingeflag or meltflag
     integer :: sym                             ! copy number for symmetric proteins
     integer :: axis                            ! index of vectorball
     integer, dimension(MAXBARREL) :: barrel    ! Array of barrels (when protein has barrels)
  END TYPE intermediate
  TYPE tstate
     INTEGER :: id                              ! id of the tstate
     INTEGER :: tp                              ! type of the tstate (hinge, pivot, etc.)
     REAL :: entropy                            ! entropy value 0.-1.
     integer :: parent, child1, child2        ! idnum from type (intermediate)
     TYPE (tstate), POINTER :: next             ! next tstate
     REAL :: energy                             !
     integer :: axis                            ! index to vectorball, axis of + rotation
     integer :: seam                            ! If this is a seam move, this is the seam number.
  END TYPE tstate

  ! Button connecting to both b-sheets
  type :: button_type  !  OBSOLETE?
    integer  :: residue
    real     :: energyBeta1, energyBeta2
    integer  :: side   !! Which side of the button was unfolded
  endtype

  ! Seam formed by two b-strands in a barrel
    type :: seam_type
        integer  :: id            ! seam identifier
        CHARACTER, dimension(maxres) :: u1flag, u2flag     ! flags that define the residues that separate
        integer  :: n             ! Number of elements (Contacts by 2 [AA1, AA2], [AA1, AA4] : 4 elements)
        integer  :: x (2, 200)    ! Array of contacts (AA1, AA2)   OBSOLETE?
        real     :: energy        ! Total energy of the seam
        integer  :: segments(4)   ! Delimiting segments seam (beta1: y1,y2, beta2: x1,x2) OBSOLETE?
    integer  :: nButtons      ! OBSOLETE?
    type (button_type) ::  buttons(200)  ! Array of buttons   OBSOLETE?
    endtype


  type :: barrel_type
    integer :: nSeams  ! number of barrels_array
    type (seam_type) :: seams(100)
  endtype
  type (barrel_type), allocatable, target  :: barrels_array(:)

  TYPE (intermediate), POINTER :: ilistroot   ! declare global intermediate list
  TYPE (tstate), POINTER :: tlistroot         ! declare global tstate list
  REAL, dimension(3, maxres) :: allcoords      ! contain all the c-a coordinates
  REAL :: nativeenergy                         ! contains information about native state
  TYPE (intermediate), POINTER :: iptr         ! pointer for ilistroot (for later)
  CHARACTER, dimension(maxres) :: masterchains ! contains all the chain IDs
  integer,dimension(maxres) :: seq, resseq
  INTEGER,parameter ::  breakflag=1, pivotflag=2, hingeflag=3, seamflag=4, meltflag=5
  INTEGER ::  geofold_nres, geofold_split
  logical :: verbose=.true.
  real :: hcutoff=0.5, bcutoff=0.2, pcutoff=0.01, scutoff=0.50 !! max number of seam buttons
  integer       :: seamcut = 8          ! MAX number of barrels
  integer, dimension(:,:), allocatable :: geofold_hb
  integer, dimension(:), allocatable :: geofold_ss
  integer ::  geofold_nhbonds
  !!
  !---------------------- INTERFACES -------------------------------------------
    interface pickunit
        module procedure geofold_pickunit
    endinterface
    interface zerointermediate
        module procedure geofold_zerointermediates, geofold_zerointermediate_ptr
    endinterface
    interface geofold_readparameter
        module procedure read_parameter_str, read_parameter_real, &
                         read_parameter_bool, read_parameter_int
    end interface
!        subroutine read_parameter_str(dunit,keyword,str,required)
!          integer,intent(in) :: dunit
!          character(len=*),intent(in) :: keyword
!          character(len=*),intent(out) :: str
!          integer,intent(in),optional :: required
!        end subroutine read_parameter_str
!        subroutine read_parameter_real(dunit,keyword,val,low,high,default,required)
!          integer,intent(in) :: dunit
!          character(len=*),intent(in) :: keyword
!          real,intent(out) :: val
!          real,intent(in),optional :: low, high, default
!          integer,intent(in),optional :: required
!        end subroutine read_parameter_real
!        subroutine read_parameter_bool(dunit,keyword,bool,default,required)
!          integer,intent(in) :: dunit
!          character(len=*),intent(in) :: keyword
!          logical,intent(out) :: bool
!          logical,intent(in),optional :: default
!          integer,intent(in),optional :: required
!        endsubroutine read_parameter_bool
!        subroutine read_parameter_int (dunit,keyword,intg,low,high,default,required)
!          integer,intent(in) :: dunit
!          character(len=*),intent(in) :: keyword
!          integer,intent(out) :: intg,low,high
!          integer,intent(in),optional :: default
!          integer,intent(in),optional :: required
!        endsubroutine read_parameter_int
!    endinterface
  !!----------------------------------------------------------------------------------
!---------------------- PUBLIC -----------------------------------------------
  public :: button_type, seam_type, barrel_type, barrels_array
  public :: geofold_nres, geofold_split, intermediate, tstate, ilistroot
  public :: tlistroot, allcoords, masterchains, breakflag
  public :: pivotflag, hingeflag, seamflag,  meltflag, verbose, hcutoff, bcutoff, pcutoff, scutoff
  public :: MAXRES, MAXCHAIN, maxsplit, MAXHINGEANG, FULSPLITDEPTH, MAXBARREL
  public :: seq, resseq, chainletters
  public :: geofold_readpdb, geofold_writepdb, geofold_hb, geofold_ss, geofold_nhbonds, pickunit
  public :: zerointermediate, geofold_readparameter
  !-----------------------------------------------------------------------------
CONTAINS
  !-----------------------------------------------------------------------------
  integer function aa2num(aa)
    character(len=3),intent(in) :: aa
    character(len=80),parameter :: allaa="ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR "
    aa2num = index(allaa,aa)/4 + 1
  end function aa2num
!!====================================================================================
  integer function geofold_pickunit(iunit)
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
    geofold_pickunit = ounit
  end function geofold_pickunit
!!====================================================================================
  SUBROUTINE geofold_readpdb(dunit)
  implicit none
  integer,intent(in) :: dunit
  ! character(len=*),intent(in) :: filename
  integer ::  ierr, ires
  character(len=1000) :: aline
  character :: chainID
  character(len=3) :: aa
  !!
  !! masterchains is global
  masterchains = '.'
  ires = 0
  DO
     read(dunit,'(a)', iostat=ierr) aline !read in pdb file
     IF ( ierr /=0 ) EXIT
     IF (aline(1:5) /= "ATOM ".and.aline(1:7) /= "HETATM ") CYCLE
     IF (aline(13:16) /=" CA " ) CYCLE
     ires = ires + 1
     ! Get Coordinates and Chain ID (_ 'Char')
     ChainID = aline(22:22)
     IF ( ChainID == ' ') THEN
        masterchains(ires) = 'b'
     ELSE
        masterchains(ires) = ChainID
     END IF
     read( aline(31:54),'(3f8.3)' ) allcoords(1:3, ires)
     aa = aline(18:20)
     seq(ires) = aa2num(aa)
     read( aline(23:26),'(i4)' ) resseq(ires)
  END DO
  geofold_nres = ires
  END SUBROUTINE geofold_readpdb
  !----------------------------------------------------------------------------------
  SUBROUTINE geofold_writepdb(dunit)
    integer,intent(in) :: dunit
    character(len=80) :: aline
    integer :: ires, ios
    character(len=80) :: aa="ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR "
    aline = "ATOM      2  CA  MET A   1       2.666  13.148  30.331  1.00 55.69"
    !--------123456789012345678901234567890123456789012345678901234567890
    do ires=1,geofold_nres
      write(aline(8:11),'(i4)',iostat=ios) ires; if (ios/=0) stop 'BUG ires'
      aline(18:20) = aa((seq(ires)-1)*4+1:seq(ires)*4-1)
      if (masterchains(ires)=='b') then
        aline(22:22) = ' '
      else
        aline(22:22) = masterchains(ires)
      endif
      write(aline(23:26),'(i4)',iostat=ios) resseq(ires); if (ios/=0) stop 'BUG ires'
      write(aline(31:54),'(3f8.3)',iostat=ios) allcoords(1:3,ires)
      if (ios/=0) stop 'BUG allcoords'
      write(aline(55:66),'(2f6.2)' ) 0.0, 0.0 ; if (ios/=0) stop 'BUG B'
      write(dunit,'(a)',iostat=ios) trim(aline) ; if (ios/=0) stop 'BUG aline'
    enddo
  END SUBROUTINE geofold_writepdb
  !----------------------------------------------------------------------------------
  subroutine geofold_zerointermediates(f)
    type(intermediate),dimension(:),intent(inout) :: f
    integer :: i
    do i=1,size(f)
      f(i)%idnum = 0
      f(i)%iflag = " "
      f(i)%state = 0
      f(i)%sym = 0
      f(i)%axis = 0
      f(i)%barrel(:) = 0
      nullify(f(i)%next)
    enddo
  end subroutine geofold_zerointermediates
  !--
  subroutine geofold_zerointermediate_ptr(f)
    type(intermediate),pointer :: f
    if (.not.associated(f)) &
      stop 'geofold_global:: zerointermediate: attempting to initialize unassociated pointer.'
    f%idnum = 0
    f%iflag = " "
    f%state = 0
    f%sym = 0
    f%axis = 0
    f%barrel(:) = 0
    nullify(f%next)
  end subroutine geofold_zerointermediate_ptr
  !---
  !subroutine geofold_zerointermediate(f)
  !  type(intermediate),intent(inout) :: f
  !  f%idnum = 0
  !  f%iflag = " "
  !  f%state = 0
  !  f%sym = 0
  !  f%axis = 0
  !  f%barrel(:) = 0
  !  nullify(f%next)
  !end subroutine geofold_zerointermediate
  function allcaps(card) result (caps)
  character(len=*),intent(in) ::  card
  character(len=len(card)) ::  caps
  integer :: I,icap
  icap = ichar('A') - ichar('a')
  do I=1,len(card)
    if ('a'<=card(I:I).and.card(I:I)<='z') then
      caps(I:I) = char(ichar(card(I:I))+icap)
    else
      caps(I:I) = card(I:I)
    endif
  enddo
  end function allcaps
  !!----------------------------------------------------------------------------------
  subroutine read_parameter_str (dunit,keyword,str,low,high,default,required)
    implicit none
    integer,intent(in) :: dunit
    character(len=*),intent(in) :: keyword
    character(len=*),intent(out) :: str
    integer,intent(in),optional :: required,low,high,default
    character(len=len(keyword)) :: key
    character(len=1000) :: aline
    integer :: ios=0,i
    rewind(dunit)
    key = allcaps(keyword)
    str = " "
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      i=index(aline,' ')-1
      if (allcaps(aline(1:i))/=trim(key)) cycle
      str = trim(adjustl(aline(i+1:)))
    enddo
    if (str/=" ") return
    if (present(required)) then
      write(0,'("Required parameter is missing: ",a)') trim(key)
      stop 'geofold_global.f90 :: read_parameter_str, required field missing.'
    endif
  end subroutine read_parameter_str
  !!----------------------------------------------------------------------------------
  subroutine read_parameter_real(dunit,keyword,val,low,high,default,required)
    implicit none
    integer,intent(in) :: dunit
    character(len=*),intent(in) :: keyword
    real,intent(out) :: val
    real,intent(in),optional :: low, high, default
    integer,optional :: required
    character(len=len(keyword)) :: key
    character(len=1000) :: aline
    integer :: ios=0,i
    rewind(dunit)
    key = allcaps(keyword)
    write(0,*) "in read_parameter_real"
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      i=index(aline,' ')-1
      if(i == 0) cycle
      !write(0,*) len(key), i
      if (allcaps(aline(1:i))/=trim(key(1:len(key)))) cycle
      !stop 'works'
      read(aline(i+1:),*,iostat=ios) val
      if (ios/=0) then
        stop 'geofold_global.f90 :: read_parameter_real. Error parsing line.'
      endif
      if (present(low)) then
        if (val<low) then
          if (present(default)) then
            val = default
          else
            write(0,*) "geofold_global.f90 :: read_parameter_real. Value too low. must be >=", low
            stop      'geofold_global.f90 :: read_parameter_real. Value too low. '
          endif
        endif
      endif
      if (present(high)) then
        if (val>high) then
          write(0,*) "geofold_global.f90 :: read_parameter_real. Value too high. must be <=", high
          if (present(default)) then
            val = default
            write(0,*) "Setting the value for ",trim(key)," to ",val
          else
            stop      'geofold_global.f90 :: read_parameter_real. Value too high. '
          endif
        endif
      endif
      return
    enddo
    if (present(required)) then
      if (present(default)) then
        val = default
        write(0,'("Required parameter is missing: ",a," Setting to default",f9.3)') trim(key), val
      else
        write(0,'("Required parameter is missing: ",a)') trim(key)
        stop 'geofold_global.f90 :: read_parameter_str, required field missing.'
      endif
    endif
  end subroutine read_parameter_real
  !!----------------------------------------------------------------------------------
  subroutine read_parameter_bool(dunit,keyword,bool,low,high,default,required)
    implicit none
    integer,intent(in) :: dunit
    character(len=*),intent(in) :: keyword
    character(len=len(keyword)) :: key
    logical,intent(out) :: bool
    logical,intent(in),optional :: default
    integer,optional :: required, low, high
    integer :: ios=0,i,j
    character(len=100) :: str
    character(len=1000) :: aline
    rewind(dunit)
    key = allcaps(keyword)
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      i=index(aline,' ')-1
      if (allcaps(aline(1:i))/=trim(key)) cycle
      read(aline(i+1:),'(L)',iostat=ios) bool
      if (ios/=0) then
        read(aline(i+1:),*,iostat=ios) j
        if (ios/=0) then
          read(aline(i+1:),*,iostat=ios) str
          if (ios/=0) then
            if (present(default)) then
              bool = default
            else
              bool = .true.
            endif
          elseif ((str(1:1)=="T").or.(str(1:1)=="t")) then
            bool = .true.
          elseif ((str(1:1)=="F").or.(str(1:1)=="f")) then
            bool = .false.
          elseif (present(default)) then
            bool = default
          else
            bool = .true.
          endif
        elseif (j==0) then
          bool = .false.
        else
          bool = .true.
        endif
      endif
      return
    enddo
    if (present(required)) then
      write(0,'("Required parameter is missing: ",a)') trim(key)
      stop 'geofold_global.f90 :: read_parameter_bool, required field missing.'
    endif
  end subroutine read_parameter_bool
  !!----------------------------------------------------------------------------------
  subroutine read_parameter_int (dunit,keyword,intg,low,high,default,required)
    implicit none
    integer,intent(in) :: dunit
    integer,intent(out) :: intg
    integer,intent(in),optional :: low,high,default,required
    character(len=*),intent(in) :: keyword
    character(len=len(keyword)) :: key
    integer :: ios=0,i
    character(len=1000) :: aline
    rewind(dunit)
    key = allcaps(keyword)
    do
      read(dunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      i=index(aline,' ')-1
      if (allcaps(aline(1:i))/=trim(key)) cycle
      read(aline(i+1:),*,iostat=ios) intg
      if (ios/=0) stop 'geofold_global.f90 :: read_parameter_int: bad int.'
      return
    enddo
    if (present(required)) then
      write(0,'("Required parameter is missing: ",a)') trim(key)
      stop 'geofold_global.f90 :: read_parameter_str, required field missing.'
    elseif (present(default)) then
      intg = default
    endif
  end subroutine read_parameter_int
  !!----------------------------------------------------------------------------------
END MODULE geofold_global
