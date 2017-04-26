!! ========= PDB2HB ============
!! Find backbone H-bonds. Output
!! donor and acceptor id numbers :  idonor iacceptor value
!! Number is sequential as used in GeoFOLD programs.
!! Input PDB must be "clean". Use getchain.f90
!! value is not used, yet.
!! Modified from hbcontactmap.f90
!! --------- in more detail ----------
!! READ PDB file  (output from getchain.f90)
!!   Loop over all residues
!!     For backbone N
!!       Find nearest backbone O under 3.6A
!!       save residue number
!!   enddo
!! -----------------------------------
!! C.Bystroff Sat Nov  8 16:00:14 EST 2008
!!  hbcontactmap.f90 C.Bystroff Thu Apr 13 06:14:50 EDT 2006
!!  Side-chain H-bonds added B. Walcott Tue May 6 14:57:34 EDT 2014
!! -----------------------------------
program pdb2hb
implicit none
integer :: iargc, jarg, ios, n,m
character(len=80) :: pdbfile,paramfile,aline
integer :: nres=0, nsg=0, nd = 0, na = 0
!real,dimension(:,:),pointer :: xyz
real,dimension(:,:),pointer :: sg
real, parameter :: pi = 3.1415927
logical :: sidechains = .true.
logical :: reduced = .false. !decide whether or not to count disulfide bonds

!Declaring donor and and acceptor types	- B Walcott 02.05.2014 12:07:39
type hb_donor
  integer :: res_num = 0				!residue number containing this group
  character(len=3) :: res_name = 'XXX'			!3-letter abbreviation for this residue
  character(len=3) :: atom ='XXX'				!atom participating in the hydrogen bond (PDB format)
  real, dimension(3) :: xyz = 999., axyz = 999.		!coordinates of the atom and the antecedent one (used to calculate hydrogen bond angle)
  type(hb_donor),pointer :: next	=> null()	!pointer to the next donor in the list
end type hb_donor

type hb_acceptor						!variables here are basically the same as donor
  integer :: res_num = 0
  character(len=3) :: res_name = 'XXX'
  character(len=3) :: atom = 'XXX'
  real, dimension(3) :: xyz = 999.
  type(hb_acceptor),pointer :: next => null()
end type hb_acceptor
type(hb_acceptor),pointer :: a_head
type(hb_donor),pointer :: d_head

!chain = '?'
jarg = iargc()
if (jarg < 2) then
  write(0,'("Usage: xpdb2hb paramfile pdbfile > hblist ")')
  stop 'pdb2hb.f90 Sat Nov  8 16:01:44 EST 2008'
endif

call getarg(2,pdbfile)

!---- announce
write(*,'("!-----------------------------------------------------!")')
write(*,'("! pdb2hb.f90 : This program extracts Hbonding infor-  !")')
write(*,'("!            mation to stdout                         !")')
write(*,'("!            C.Bystroff 8-NOV-2008                    !")')
write(*,'("!-----------------------------------------------------!")')
write(*,'("! donor acceptor H-or-S ")')
write(*,'("! pdb file = ",a)') trim(pdbfile)

  allocate(d_head)
  allocate(a_head)

  !read SIDECHAINHB from parameters file and REDUCING
  call getarg(1,paramfile)
  open(8, file=paramfile,status='old',form='formatted',action='read',iostat=ios)
  if(ios/=0) stop 'Error reading parameter file'
  do
    read(8,'(a)',iostat=ios)aline
    if(ios/=0) exit
    if(aline(1:11)/='SIDECHAINHB' .and. aline(1:8)/='REDUCING') cycle
    if(aline(1:11)=='SIDECHAINHB' .and. aline(13:13)=='0') sidechains = .false.
    if(aline(1:8)=='REDUCING' .and. aline(10:10)=='1') reduced = .true.
  enddo


  !! read PDB file
!  write(*,*) "Start READPDB"
  call READPDB(pdbfile, d_head, a_head, nres, sg, nsg, nd, na, sidechains)
!for debugging
!  call readout(d_head, a_head)
!  call readoutss(sg,nres)
  !! find backbone H-bonds
  call GETHBONDS(d_head, a_head)
  call cleanup(d_head, a_head)
!  write(*,*) "Lists cleaned"
  if (nsg/=0 .and. .not.  reduced) call GETSSBONDS(sg,nres)

!if (associated(xyz)) deallocate(xyz)
if (associated(sg)) deallocate(sg)

CONTAINS

!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
subroutine GETSSBONDS(sg,nres)
  integer,intent(in) :: nres
  real,dimension(:,:),pointer :: sg
  integer :: ires,jres,ios,nss
  real,dimension(3) :: vec
  real :: dd
  real,parameter :: SSCUT=2.5
  if (.not.associated(sg)) stop 'Bug in GETSSBONDS 5'
  !! loop over all residues
  nss = 0
  do ires=1,nres
    if (sg(1,ires)==0.000) cycle
    do jres=ires+1,nres
      if (sg(1,jres)==0.000) cycle
      vec = sg(1:3,ires) - sg(1:3,jres)
      dd = sqrt(sum(vec*vec))
      if (dd > SSCUT) cycle
      nss = nss + 1
      write(*,'(i7," SG  ",i7," SG  "," S ")') ires, jres
      exit
    enddo
  enddo
  deallocate(sg)
end subroutine GETSSBONDS
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!!---------Subroutines added by B. Walcott 02.05.2014 12:08:07----------!
!-----------------------------------------------------------------------!
!-------------new READPDB B. Walcott 05.05.2014 17:06:26----------------!
!!READPDB rewritten by B. Walcott 02.05.2014 14:40:15!

!simplifying assumption made: amino groups behave
!only as hydrogen bond donor

subroutine READPDB(pdbfile, d_head, a_head, nres, sg, nsg, nd, na, schb)
  implicit none
  character(len=*),intent(in) :: pdbfile
  type(hb_donor),pointer :: d_head, i_d
  type(hb_acceptor),pointer :: a_head, i_a
  real,dimension(:,:),pointer :: sg
  integer,intent(out) :: nres, nsg, nd, na
  integer :: ios
  logical, intent(in) :: schb
  character(len=80) :: card,filename
  character(len=4) ::  atomname
  character(len=4) ::  acid
  character(len=6) :: last
!  character(len=1) :: ocode
  integer :: jres, ires
  real,dimension(3) :: coords2, prevO
  real,dimension(8,3) :: coordsC  !coordinates of previous C, CA, CB, CG, CD, CE, CZ, CH
  real :: b,p !,chi1,chi2,x
!  integer :: I,J,K,L,i1,i2,ist,startres
  character(len=80),parameter :: pdbinfmt="(BZ,6x,5x,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2)"

!  write(*,*) "Variables declared"

  !** initialization
  i_d=>d_head
  nullify(i_d%next)
  i_a=>a_head
  nullify(i_a%next)
  last = '      '
  nsg = 0
  open(1, file=pdbfile,status='old',form='formatted',action='read',iostat=ios)
  if(ios/=0) then
    write(0,*) 'Error opening PDB file: ', trim(pdbfile)
    stop 'Error opening PDB file'
  endif
!  ocode = ' '
  jres = 0
  ires = 0

  !Count the number of residues in file
  do
    read(1, '(a80)',iostat=ios) card
    if(ios/=0) exit
    if (jres > 0 .and. card(1:5) == 'ENDMD') exit
    if (jres > 0 .and. card(1:5) == 'TER  ') exit
    if (card(1:5) /=  'ATOM ') cycle
    if (card(22:27) /= last) then
      jres = jres + 1
      last = card(22:27)
    endif
  enddo
  nres = jres
  write(*, '("!",i9," residues in PDB file.")') nres

  !start from the beginning of file, more initialization
  rewind(1)
  allocate(sg(3,nres))
  sg = 0.0
  last = '      '
  jres = 0
  ios = 0
  coordsC(1:8,1:3) = 999.
  coords2 = 999.

  !Parse the PDB file into lists of acceptor and donor groups and an array of cys sulfides
  do
    read(1,'(a80)', iostat=ios) card
    if(ios/=0) exit
    if (jres > 0 .and. card(1:5) == 'ENDMD') exit
    if (jres > 0 .and. card(1:5) == 'TER  ') exit
    if (card(1:5) /= 'ATOM ') cycle
    if (atomname(2:2)=='C') then
	  select case(atomname)
    case(' C  ')
		  coordsC(1,1:3) = coords2
		case(' CA ')
		  coordsC(2,1:3) = coords2
		case(' CB ')
		  coordsC(3,1:3) = coords2
		case(' CG ')
		  coordsC(4,1:3) = coords2
		case(' CD ')
		  coordsC(5,1:3) = coords2
		case(' CE ')
		  coordsC(6,1:3) = coords2
		case(' CZ ')
		  coordsC(7,1:3) = coords2
		case(' CH ')
		  coordsC(8,1:3) = coords2
		case default
	  end select
	endif
    read(card,pdbinfmt,iostat=ios) atomname,acid,ires,coords2,b,p
    if(ios/=0) cycle
    if(card(22:27) /= last) then	!new residue starting
      jres = jres + 1
      if(jres>nres) exit
      last = card(22:27)
    endif

      !debugging
!      write(*,'(a3, ", residue number ",i7," of ",i7)') trim(acid),jres,nres

      !backbone hydrogen bond donors/acceptors
      !Pro is an imino acid and doesn't donate a hydrogen bond
      if(acid /= 'PRO ') then
!        write(*,*) "Not PRO"
        if(atomname == ' N  ') then
          if(i_d%res_num /= 0) then
            allocate(i_d%next)
            i_d=>i_d%next
          endif
          i_d%res_num = jres
          i_d%res_name = trim(acid)
          i_d%atom = atomname(2:4)
          i_d%xyz = coords2
          i_d%axyz = prevO
          nd=nd+1
          cycle
        endif
      endif
      if(atomname == ' O  ' .or. atomname == ' OXT') then
        if(i_a%res_num /= 0) then
          allocate(i_a%next)
          i_a=>i_a%next
        endif
        i_a%res_num = jres
        i_a%res_name = trim(acid)
        i_a%atom = atomname(2:4)
        i_a%xyz = coords2
        prevO = coords2
        na=na+1
        cycle
      endif

      if(schb .or. trim(acid)=='CYS') then
        !sidechain hydrogen bond donors/acceptors
        select case(trim(acid))

          !hydrogen bond accepting sidechains
          case('ASP')
  			if(atomname == ' OD1' .or. atomname == ' OD2') then
          if(i_a%res_num /= 0) then
            allocate(i_a%next)
            i_a=>i_a%next
          endif
  			  i_a%res_num = jres
  			  i_a%res_name = trim(acid)
  			  i_a%atom = atomname(2:4)
  			  i_a%xyz = coords2
  			  na=na+1
  			endif
          case('GLU')
  			if(atomname == ' OE1' .or. atomname == ' OE2') then
          if(i_a%res_num /= 0) then
            allocate(i_a%next)
            i_a=>i_a%next
          endif
  			  i_a%res_num = jres
  			  i_a%res_name = trim(acid)
  			  i_a%atom = atomname(2:4)
  			  i_a%xyz = coords2
  			  na=na+1
  			endif

          !hydrogen bond donating/accepting sidechains
          case('ASN', 'GLN')
            !Carbonyl Oxygen acceptor
            if(atomname(2:2) == 'O') then
              if(i_a%res_num /= 0) then
                allocate(i_a%next)
                i_a=>i_a%next
              endif
              i_a%res_num = jres
              i_a%res_name = trim(acid)
              i_a%atom = atomname(2:4)
              i_a%xyz = coords2
              na=na+1
            endif
            !amino nitrogen donor
            if (atomname(2:2) == 'N' ) then
              if(i_d%res_num /= 0) then
                allocate(i_d%next)
                i_d=>i_d%next
              endif
              i_d%res_num = jres
              i_d%res_name = trim(acid)
              i_d%atom = atomname(2:4)
              i_d%xyz = coords2
              if(atomname(3:3) == 'D') i_d%axyz = coordsC(4,1:3)
              if(atomname(3:3) == 'E') i_d%axyz = coordsC(5,1:3)
              nd=nd+1
            endif

          !sidechains with hydroxyl groups
          case('SER', 'THR', 'TYR')
            if(atomname(2:2) /= 'O') cycle
            !acceptor
            if(i_a%res_num /= 0) then
              allocate(i_a%next)
              i_a=>i_a%next
            endif
            i_a%res_num = jres
            i_a%res_name = trim(acid)
            i_a%atom = atomname(2:4)
            i_a%xyz = coords2
            na=na+1
            !donor
            if(i_d%res_num /= 0) then
              allocate(i_d%next)
              i_d=>i_d%next
            endif
            i_d%res_num = jres
            i_d%res_name = trim(acid)
            i_d%atom = atomname(2:4)
            i_d%xyz = coords2
            if(trim(acid)=='TYR') then
              i_d%axyz = coordsC(7,1:3)
            else
                i_d%axyz = coordsC(3,1:3)
            endif
                nd=nd+1

          !hydrogen bond donating sidechains
          case('ARG', 'LYS', 'TRP', 'HIS')
            if(atomname(2:2) /= 'N' ) cycle
            if(i_d%res_num /= 0) then
              allocate(i_d%next)
              i_d=>i_d%next
            endif
            i_d%res_num = jres
            i_d%res_name = trim(acid)
            i_d%atom = atomname(2:4)
            i_d%xyz = coords2
            if(atomname(3:3)=='D') then
              i_d%axyz = coordsC(4,1:3)
            else if(atomname(3:3)=='E') then
              i_d%axyz = coordsC(5,1:3)
            else if(atomname(3:3)=='Z') then
              i_d%axyz = coordsC(6,1:3)
            else
              i_d%axyz = coordsC(7,1:3)
            endif
                nd=nd+1

          !Cysteine add sulfur group
          case('CYS')
  		  if(atomname /= ' SG ') cycle
  		  nsg = nsg + 1
  		  sg(1:3,jres) = coords2

          !Do nothing
          case default
        end select
      endif

  enddo
  close(1)
end subroutine READPDB
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
subroutine GETHBONDS(d_head, a_head)
  implicit none
  type(hb_donor),pointer,intent(in) :: d_head
  type(hb_acceptor),pointer,intent(in) :: a_head
  type(hb_donor),pointer :: i_d
  type(hb_acceptor),pointer :: i_a
  integer :: ios
  real, parameter :: HBCUT=3.5 !, pi = 3.1415927
  real :: dd, angle
  real,dimension(3) :: vec

  !!initialize
  i_d=>d_head
  i_a=>a_head
  angle = 0.0

  !!loop over all donors
  do while(associated(i_d))
    if(i_d%xyz(1)>998.) then
      i_d=>i_d%next
      cycle
    endif
    i_a=>a_head
    do while(associated(i_a))
      if(i_d%res_num == i_a%res_num) then
        i_a=>i_a%next
        cycle
      endif
     if(i_a%xyz(1)>998.) then
        i_a=>i_a%next
       cycle
      endif
      vec = i_d%xyz - i_a%xyz
      dd = sqrt(sum(vec*vec))
     if(dd > HBCUT) then
        i_a=>i_a%next
        cycle
      endif
      !! for all close acceptors, calculate angle defined by it, the donor
      !  and the atom preceeding the donor
      call calc_angle(i_d%xyz, i_a%xyz, i_d%axyz, angle)
      if(angle <(pi/2) .or. angle > pi) then
        i_a=>i_a%next
        cycle
      endif
      !if angle is between pi/2 and pi (90 and 180 deg), H-bond formed

      !example output: 4 N 9 OE2 H
      write(*,'(i7," ",a3," ",i7," ",a3," "," H ")') i_d%res_num, i_d%atom, i_a%res_num, i_a%atom
      i_a=>i_a%next
      !if you only want 1 h-bond per acceptor
      !exit
    enddo
    i_d=>i_d%next
  enddo
end subroutine GETHBONDS
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
subroutine calc_angle(don, acc, ant, angle)			!calculates the angle (rad) between three points using law of cos
  real,dimension(3), intent(in):: don, acc, ant		!coordinates of the donor, acceptor, and antecedent
  real,dimension(3) :: temp = 0.0					!temporary variable
  real,intent(out) :: angle							!the angle between these three points
  real :: acc_ant, acc_don, ant_don					!sides of the imaginary triangle formed by these three points

  angle = 0.0
  !Calculate acc_ant
  temp = acc(1:3)-ant(1:3)
  acc_ant = sqrt(sum(temp*temp))

  !Calculate acc_don
  temp = acc(1:3)-don(1:3)
  acc_don = sqrt(sum(temp*temp))

  !Calculate ant_don
  temp = ant(1:3)-don(1:3)
  ant_don = sqrt(sum(temp*temp))

  !Calculate angle
  angle = acos(((ant_don*ant_don)+(acc_don*acc_don)-(acc_ant*acc_ant))&
               /(2*ant_don*acc_don))

end subroutine calc_angle
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
subroutine cleanup(d_head, a_head)			!clears the memory from the linked lists
  type(hb_donor),pointer :: d_head, i_d, j_d
  type(hb_acceptor),pointer :: a_head, i_a, j_a
  i_d=>d_head
  do while(associated(i_d%next))
    j_d=>i_d%next
    deallocate(i_d)
    i_d=>j_d
  enddo
  i_a=>a_head
  do while(associated(i_a%next))
    j_a=>i_a%next
    deallocate(i_a)
    i_a=>j_a
  enddo
  nullify(d_head)
  nullify(a_head)
  nullify(j_d)
  nullify(j_a)
  nullify(i_d)
  nullify(i_a)
end subroutine cleanup
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
subroutine readout(d_head, a_head)
  type(hb_donor),pointer,intent(in) :: d_head
  type(hb_acceptor),pointer,intent(in) :: a_head
  type(hb_donor),pointer :: i_d
  type(hb_acceptor),pointer :: i_a

  i_d => d_head
  i_a => a_head
  write(*,*) "Writing out Donors list"
  do
    if(.not. associated(i_d)) exit
    if(i_d%res_name=='XXX')exit
    write(*,'(a3,i7," ",a3, "("e10.3","e10.3","e10.3")")')i_d%res_name, i_d%res_num, i_d%atom, i_d%xyz
    i_d => i_d%next
  enddo
  write(*,*) "Writing out Acceptors list"
  do
    if(.not. associated(i_a)) exit
    if(i_a%res_name=='XXX')exit
    write(*,'(a3,i7," ",a3, "("e10.3","e10.3","e10.3")")')i_a%res_name, i_a%res_num, i_a%atom, i_a%xyz
    i_a => i_a%next
  enddo
end subroutine readout
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

!For debugging disulfide stuff, writes out all of the SG groups
subroutine readoutss(sg,nres)
  integer,intent(in)::nres
  real,dimension(:,:),pointer :: sg
  integer :: i

  do i=1,nres
    if(sg(1,i) == 0.0) cycle
    write(*,'("SG ",i3,"(",e10.3","e10.3","e10.3")")') i, sg(1:3,i)
  enddo

end subroutine readoutss
end program pdb2hb
