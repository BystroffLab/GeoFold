!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!          geofold_flory.f90 -B Walcott 26 September 2014          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module geofold_flory
  use geofold_global
  
  !private
  public
      
  type :: contact
    integer :: cost = 9999
    integer :: path = 0
    integer,dimension(:),allocatable :: contacts
  end type contact
  
  !!Reimplementing with dynamically allocated arrays

  public :: geofold_flory_calc_entropy,geofold_flory_get_subcontacts
  public :: contact, geofold_flory_all_contacts, Path
  
interface find_contacts
  module procedure find_contacts_noseam,find_contacts_seam
endinterface
  
contains

!Calculates the change in contact order for the transition state
!To be used in the Flory equation in place of a sequence separation
!contact order measures the avg seq separation of contacts within a protein
!So, the change in that should be reflective of the change in conformational
!entropy it confers as per the Flory equation
real function DCO(c_list,w,f,u1,u2)
  implicit none
  type(contact),dimension(:),allocatable,intent(in) :: c_list
  type(intermediate),pointer,intent(in) :: f,u1
  type(intermediate),pointer,intent(in),optional :: u2
  real,intent(in) :: w
  type(contact),dimension(:),allocatable :: f_contacts, broken
  real :: fco, uco,fn,un
  integer :: i,j,cj
  
  !initialize parent contact order
  fco = 0.
  fn = 0.
  !find the contacts that are broken by the transition and the parent subcontacts
  if(.not. present(u2)) then
    call find_contacts_seam(f,u1,c_list,broken,geofold_nres,f_contacts)
  else
    call find_contacts_noseam(f,u1,u2,broken, c_list,geofold_nres,f_contacts)
  endif
  !calculate the contact order of the parent best case: O(nres), worst case: O(nres^2), avg: O(nres*b)
  iloop: do i = 1, geofold_nres
    jloop: do j = 1, geofold_nres
      if(f_contacts(i)%contacts(j)==0) exit jloop
      cj = f_contacts(i)%contacts(j)
      fco = fco + abs(i-cj)
      fn = fn + 1.
    enddo jloop
  enddo iloop
  !initialize the transition contact order
  uco = fco
  un = fn  
  !subtract broken contacts from the transition contact order O(nres), worst case: O(nres^2), avg: O(nres*b)
  iloop2: do i = 1, geofold_nres
    jloop2: do j = 1, geofold_nres
      if(broken(i)%contacts(j)==0) exit jloop2
      cj = broken(i)%contacts(j)
      uco = uco - abs(i-cj)
      un = un - 1.
    enddo jloop2
  enddo iloop2
  !calculate the parent and transition contact orders
  fco = fco/fn
  uco = uco/un
  !calculate the difference between them
  DCO = abs(uco-fco)
end function DCO

recursive subroutine newLeff(a,b,subcontacts,wn,T,Ln,Ld, Leff,searched)
  implicit none
  integer,intent(in) :: a,b
  real, intent(in) :: wn,T
  type(contact),dimension(:),allocatable,intent(in) :: subcontacts
  real,parameter :: R = 8.314
  real :: Li
  real,intent(inout) :: Ln,Ld
  real,intent(inout) :: Leff
  integer :: i,cntct,j
  integer,save :: flory = 0.
  integer,dimension(10) :: cov
  integer, save :: start = 0
  integer :: diff
  logical,dimension(geofold_nres),intent(inout):: searched

  cov = 0
  j = 0
  !if(Leff /= Leff) stop "Leff NaN"
  if(start == 0) then
    flory = abs(a-b)
    start = a
    Ln = 0.
    Ld = 0.
  else
  endif
  searched(a) = .true.
  diff = abs(a-start)
  !base: a = b
  if(a == b) then
    Leff = Ln/Ld
    flory = 0
  else
    do i = 1, geofold_nres      
      Li = 0.
      cntct = subcontacts(a)%contacts(i)
      if(cntct == 0) exit
      if(start == a .and. cntct == b) cycle
      if(isCovalent(a,cntct)) then
        if(j == 10) cycle
        j = j + 1
        cov(j) = cntct
        cycle
      endif
      Li = real(abs(cntct-b)+wn)+diff
      if(Li > flory) cycle
      Ln = Ln + Li*exp(-Li/(R*T))
      Ld = Ld + exp(-Li/(R*T))
    enddo
    do i = 1, j
      if(searched(cov(i))) cycle
      call newLeff(cov(i),b,subcontacts,wn,T,Ln,Ld,Leff,searched)
      searched(cov(i)) = .true.
      if(cov(i) == b) exit
    enddo
  endif
  if(start == a) then
    start = 0
  endif
end subroutine newLeff

!depth-limited search to find all paths between residues a and b where depth
!limit is flory number.  Calculates Leff as it goes.
!Works #BDW 09.03.2015
real function Leff(a,b,subcontacts,wn, T)
  implicit none
  integer,intent(in) :: a,b
  real,intent(in) :: wn,T
  type(contact),dimension(:),allocatable,intent(in) :: subcontacts
  real, parameter :: R = 8.314 !8.314 J/molK using R instead of Boltzmann constant
  real :: Li,Ln,Ld 
  integer :: flory,i,j,prev,k
  integer,dimension(:,:),allocatable :: path
  logical,dimension( : ),allocatable :: searched
  real :: start,finish
  
  call cpu_time(start)
  
  !initialize searched
  if(allocated(searched)) deallocate(searched)
  allocate(searched(geofold_nres))
  searched = .false.
  searched(a) = .true.
  !calculate flory
  flory = abs(a-b)+1
  !initialize path
  if(allocated(path)) deallocate(path)
  allocate(path(2,flory))
  path = 0
  path(1,1) = a
  !path is a 2xflory matrix the first row are residues numbers within the path
  !second row are the index of those residues in the previous one's set of contacts
  !initialize Ln and Ld (the numerator and denominator of the Leff equation)
  Ln = 0.
  Ld = 0.
  i = 1
  loop1: do while(i /= 0)
    i = i + 1
    j = path(2,i)+1
    prev = path(1,i-1)
    if(i == 2 .and. subcontacts(prev)%contacts(j)==b) j = j + 1
    if(subcontacts(prev)%contacts(j)==0) then
      call backtrack(searched,path,i)
      cycle loop1
    endif
    do while(searched(subcontacts(prev)%contacts(j)))
      j = j + 1
      if(subcontacts(prev)%contacts(j)==0) then
        call backtrack(searched,path,i)
        cycle loop1
      endif
    enddo
    path(1,i) = subcontacts(prev)%contacts(j)
    path(2,i) = j
    !call fprintArray(path)
    searched(path(1,i)) = .true.
    if(searched(b)) then
      Li = 0.
      do k = 1, flory-1
        if(path(1,k+1) == 0) exit
        if(isCovalent(path(1,k),path(1,k+1))) then
          Li = Li + 1.
          !write(6,*) "+1"
        else
          Li = Li + wn
          !write(6,*) "+wn"
        endif
      enddo
      Ln = Ln + Li*exp(-Li/(R*T))
      Ld = Ld + exp(-Li/(R*T))
      call cpu_time(finish)
      !write(6,*) finish
      !write(6,*)"Li,Ln,Ld: ",Li,Ln,Ld
      !call fprintArray(path)
    endif
    if(i == flory .or. searched(b)) then  !searched(b) implies that path(1,i) = b
      call backtrack(searched,path,i)
      cycle loop1
    endif    
  enddo loop1
  Leff = Ln/Ld
  call cpu_time(finish)
  write(6,*) Leff,finish-start
end function Leff

subroutine backtrack(searched, path, i)
  implicit none
  logical,dimension(:),allocatable :: searched
  integer :: i
  integer,dimension(:,:),allocatable :: path
  
  if(.not. allocated(searched) .or. .not. allocated(path)) then
    stop "geofold_flory::backtrack: array not allocated"
  endif
  searched(path(1,i)) = .false.
  path(1,i) = 0
  path(2,i) = 0
  searched(path(1,i-1)) = .false.
  path(1,i-1) = 0
  i = i - 2
  !call fprintArray(path)
end subroutine backtrack

!removes paths from the paths array that contain only zeroes
!reallocates the array to be the correct size
!Works #BDW 26.02.2015
subroutine removeZeroes(paths)
  implicit none
  integer,dimension(:,:),allocatable,intent(inout) :: paths
  integer,dimension(:,:),allocatable :: tmp
  integer :: i, nzero,newsize,j
  integer, dimension(2) :: oldsize
  
  
  if(.not. allocated(paths)) stop "geofold_flory::removeZeroes: Error - &
    paths array not allocated."
  oldsize(1) = size(paths,1)
  oldsize(2) = size(paths,2)
  nzero = 0
  allocate(tmp(oldsize(1),oldsize(2)))
  j = 0
  do i = 1, oldsize(1)
    if(paths(i,1) == 0) then
      nzero = nzero + 1
    else
      j = j+1
      tmp(j,1:oldsize(2)) = paths(i,1:oldsize(2))
    endif    
  enddo
  newsize = oldsize(1) - nzero
  call reallocate(tmp,newsize)
  if(allocated(paths)) deallocate(paths)
  allocate(paths(newsize,oldsize(2)))
  paths = tmp
  if(allocated(tmp)) deallocate(tmp)
end subroutine removeZeroes
  
  

!reallocates a 2d array of integers, part of two subroutines to turn a
!2d array of integers into a vector
!Works #BDW 26.02.2015
subroutine reallocate(array,ns)
  implicit none
  integer,dimension(:,:),allocatable :: array,tmp
  integer :: oldsize,width,newsize,i
  integer, optional :: ns
  
  if(.not. allocated(array)) then
    stop "geofold_flory::reallocate: Error - Array not allocated"
  endif
  oldsize = size(array,1)
  width = size(array,2)
  if(present(ns)) then
    newsize = ns
  else 
    newsize = oldsize+10000
  endif
  if(allocated(tmp)) deallocate(tmp)
  allocate(tmp(oldsize,width))
  tmp = array
  if(allocated(array)) deallocate(array)
  allocate(array(newsize,width))
  array = 0
  if(newsize > oldsize) then
    do i = 1, oldsize
      array(i,1:size(array,2)) = tmp(i,1:size(array,2))
    enddo
  else
    do i = 1, newsize
      array(i,1:size(array,2)) = tmp(i,1:size(array,2))
    enddo
  endif
  if(allocated(tmp)) deallocate(tmp)
end subroutine reallocate

!adds value n to array, reallocates array if size is too small
!works #BDW 26.02.2015
subroutine append(array,n)
  implicit none
  integer,dimension(:,:),allocatable :: array
  integer,dimension(:),allocatable :: n
  integer :: i
  integer,save :: ind = 1
  
  if(size(array,2) /= size(n)) then
    stop "geofold_flory::append: Error - adding array of wrong size"
  endif
  do while(array(ind,1) /= 0)
    ind = ind + 1
    if(ind > size(array,1)) call reallocate(array)
  enddo    
  array(ind,1:size(array,2)) = n(1:size(array,2))
end subroutine append


!Works #BDW 26.02.2015
logical function isCovalent(a,b)
  implicit none
  integer :: a,b
  
  if(abs(a-b) == 1 .or. geofold_ss(a) == b) then
    isCovalent = .true.
  else
    isCovalent = .false.
  endif
end function isCovalent

  subroutine fprintArray(array)
    implicit none
    integer,dimension(:,:),allocatable :: array
    integer :: i,j
        
    write(*,'(i3,"x",i3," array")') size(array,1),size(array,2)
    do i = 1, size(array,1)
      write(*,*) array(i,1:size(array,2))
    enddo
  end subroutine fprintArray
 

!find all valid paths between contacts a and b
!a path is valid if it is <= default loop length
!and doesn't loop over itself
!needs to be rewritten for new implementation
subroutine findPaths(a,b,paths,subcontacts)
  implicit none
  integer :: a,b,flory,i,j,k,contct,point
  integer,dimension(:,:),allocatable,intent(inout) :: paths
  integer,dimension(:),allocatable :: tmp
  type(contact),dimension(:),allocatable,intent(in) :: subcontacts
  
  write(0,*) "finding paths...."
  
  !first, set flory
  flory = abs(a-b)+1
  !initialize paths
  if(allocated(paths)) deallocate(paths)
  allocate(paths(9000000,flory))
  allocate(tmp(flory))
  paths = 0
  paths(1,1) = a !start of first (and every) path
  do i = 1, flory-1
    j = 1
    do while(j<size(paths,1))
      if(paths(j,1)==0) then
        exit
      endif
      write(0,*) "j = ",j, " of ",size(paths,1)     
      write(0,*) "i ",i," out of ",flory-1
      if(.not. inPath(b,paths,j) .and. paths(j,i) /= 0) then
        if(i /= flory) then
          if(paths(j,i+1) /= 0)then
            j = j+1
            cycle
          endif
        endif
        point = paths(j,i)
        !check neighbors
        k = 1
        contct = subcontacts(point)%contacts(k)
        do while(contct /= 0)
          if(contct == b .and. i == 1) then
            k = k+1
            contct = subcontacts(point)%contacts(k)
            cycle
          endif
          !create new path for all that haven't been used
          if(.not. inPath(contct,paths,j)) then
            tmp(1:flory) = paths(j,1:flory)
            tmp(i+1) = contct
            call append(paths,tmp)
            !call fprintArray(paths)
            tmp = 0
          endif
          k = k+1
          contct = subcontacts(point)%contacts(k)
        enddo
      endif
      j = j+1      
      !if(mod(j,100000) == 0) call removeZeroes(paths)
    enddo
  enddo
  do i = 1, size(paths,1)
    if(.not. inPath(b,paths,i)) paths(i,1:flory) = 0
  enddo
  call removeZeroes(paths)
  call fprintArray(paths)
end subroutine findPaths

!determine the weighted average of the pathlength
!using Boltzmann weighting
!Works #BDW 27.02.2015
real function getLeff(wc,wn,paths,T) result(Leff)
  implicit none
  real,intent(in) :: wc, wn
  integer,dimension(:,:),allocatable,intent(in) :: paths
  real :: Ln, Ld,Li
  real, parameter :: R = 8.314
  real,intent(in) :: T
  integer :: i,j,a,b
  
  Ln = 0.
  Ld = 0.
  
  do i = 1, size(paths,1)    
    Li = 0.
    do j = 1, size(paths,2)-1
      a = paths(i,j)
      b = paths(i,j+1)
      if(isCovalent(a,b)) then
        Li = Li + wc
      else
        Li = Li + wn
      endif
    enddo
    Ln = Ln + Li*exp(-Li/(R*T))
    Ld = Ld + exp(-Li/(R*T))
  enddo
  Leff = Ln/Ld
end function getLeff

!test to see if path already contains a residue

logical function inPath(a,paths,ind)
  implicit none
  integer :: a,i,ind
  integer,dimension(:,:),allocatable :: paths

  inPath = .false.
  do i = 1, size(paths,2)
    if (paths(ind,i)==a) then
      inPath = .true.
      exit
    endif
  enddo
end function inPath
  
!get subcontacts contained within f

subroutine geofold_flory_get_subcontacts(c_list,subcontacts,f,nres)
  implicit none
  integer :: nres,i,j,k
  type(contact), dimension(:),allocatable, intent(in) :: c_list
  type(contact),dimension(:), allocatable,intent(out) :: subcontacts
  type(intermediate),pointer,intent(in) :: f
  
  if(.not.allocated(subcontacts)) then
  allocate(subcontacts(nres))
  do i = 1, nres
    allocate(subcontacts(i)%contacts(nres))
    subcontacts(i)%contacts = 0
  enddo
  endif
  do i = 1, nres
  if(f%iflag(i)==".") cycle
  do j = 1, nres
    if(c_list(i)%contacts(j)==0) exit
    k = c_list(i)%contacts(j)
    if(f%iflag(k)==".")cycle
    if(.not. inArray(k,subcontacts(i)%contacts,nres))&
    call addValue(k,subcontacts(i)%contacts,nres)
    if(.not. inArray(i,subcontacts(k)%contacts,nres))&
    call addValue(i,subcontacts(k)%contacts,nres)
    k = 0
  enddo
  enddo
end subroutine geofold_flory_get_subcontacts

!finds the contacts broken by a seam move, requires existence of barrels_array
subroutine find_contacts_seam(f,u1,c_list,broken,nres,subcontacts)
  implicit none
  type(contact), dimension(:),allocatable, intent(in) :: c_list
  type(contact),dimension(:), allocatable,intent(out) :: broken
  type(contact),dimension(:),allocatable,intent(out) :: subcontacts
  type(intermediate),pointer,intent(in) :: f,u1
  integer :: i,j,nres,seam,k
  character,dimension(maxres) :: u1flag, u2flag
  
  if(.not. allocated(barrels_array)) stop "geofold_flory::find_contacts&
  _seam: Error, barrels_array not allocated."
  call geofold_flory_get_subcontacts(c_list,subcontacts,f,nres)
  !determine which seam is broken
  do i = 1, maxbarrel 
  if(u1%barrel(i)==0) cycle
  seam = u1%barrel(i)
  exit
  enddo
  !set u1flag and u2flag
  u1flag = barrels_array(i)%seams(seam)%u1flag
  u2flag = barrels_array(i)%seams(seam)%u2flag
  !initializations
  if(allocated(broken)) deallocate(broken)
  allocate(broken(nres))
  do i = 1, nres
    allocate(broken(i)%contacts(nres))
    broken(i)%contacts = 0
  enddo
  !iterate through f and find the contacts it contains
  do i = 1, nres
  !check if contact exists in f
    do j = 1, nres
    if(subcontacts(i)%contacts(j)==0)exit
    k=subcontacts(i)%contacts(j)
    if(u1flag(i)/='.')then
    if(u1flag(k)/='.') cycle
    if(.not. inArray(k,broken(i)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
    call addValue(k,broken(i)%contacts,nres)
    endif
    if(.not. inArray(i,broken(k)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
    call addValue(i,broken(k)%contacts,nres)
    endif
    else if(u2flag(i)/='.') then
    if(u2flag(k)/='.') cycle
    if(.not. inArray(k,broken(i)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
    call addValue(k,broken(i)%contacts,nres)
    endif
    if(.not. inArray(i,broken(k)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
    call addValue(i,broken(k)%contacts,nres)
    endif
    endif
    enddo
  enddo
  
end subroutine find_contacts_seam

!checks if a node is in the visited array

logical function inArray(a,array,arraySize) result(answer)
  integer,intent(in) :: a,arraySize
  integer,dimension(:),intent(in) :: array
  integer :: i
  
  answer = .false.
  do i = 1, arraySize
    if(array(i)==a) then
      answer = .true.
      exit
    endif
    if (array(i)==0) exit
  enddo
end function inArray

!finds the contacts broken by a particular move (not seam)
!Operational 7 October 2014

subroutine find_contacts_noseam(f,u1,u2,broken, c_list,nres,subcontacts)
  implicit none
  type(contact),dimension(:),allocatable, intent(in) :: c_list
  type(contact),dimension(:),allocatable, intent(out) :: broken
  type(contact),dimension(:),allocatable,intent(out) :: subcontacts
  integer :: i,j,nres,k
  type(intermediate),pointer,intent(in) :: f,u1,u2
  
  !initializations
  if(allocated(broken)) deallocate(broken)
  allocate(broken(nres))
  do i = 1, nres
    allocate(broken(i)%contacts(nres))
    broken(i)%contacts = 0
  enddo
  call geofold_flory_get_subcontacts(c_list,subcontacts,f,nres)
  do i = 1, nres
  do j = 1, nres
    if(subcontacts(i)%contacts(j)==0)exit
    k=subcontacts(i)%contacts(j)   
    if(u1%iflag(i)/='.')then
      if(u1%iflag(k)/='.') cycle
      if(.not. inArray(k,broken(i)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
        call addValue(k,broken(i)%contacts,nres)
        endif
      if(.not. inArray(i,broken(k)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
      call addValue(i,broken(k)%contacts,nres)
      endif
    else if(u2%iflag(i)/='.') then
      if(u2%iflag(k)/='.') cycle
      if(.not. inArray(k,broken(i)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
      call addValue(k,broken(i)%contacts,nres)
      endif
      if(.not. inArray(i,broken(k)%contacts,nres))then
        if(k /= i-1 .and. k/=i+1)&
      call addValue(i,broken(k)%contacts,nres)
      endif
    endif
    enddo
  enddo
end subroutine find_contacts_noseam

!finds all contacts within the native state and creates a data structure
!to contain them
!Operational 6 October 2014

subroutine geofold_flory_all_contacts(mfile, cijfile, c_list, nres)
  implicit none
  type(contact),dimension(:),allocatable :: c_list
  character (len=*), intent(in) :: mfile, cijfile
  integer :: res1, res2, ierr, dunit = 31, nres,i,j
  real :: tmp
  character (len=3) :: datom, aatom
  character (len=1):: bond
  character (len=200) :: aline
  

  !determine nres
    !open mfile (hbonds file)
    dunit = pickunit(dunit)
    open(dunit,file=mfile,iostat=ierr,status='old',form='formatted')
    if(ierr /= 0) stop "geofold_flory::all_contacts: Error opening &
    hbonds file"
    !find line 8
    do i = 1, 8
      read(dunit, '(a)') aline
    enddo
    !read aline(2:10) as (i9) into nres
    read(aline(2:10),'(i9)',iostat=ierr) nres
    if(ierr/=0) stop 'geofold_flory::all_contacts: Error reading nres'
  !allocate nres contacts in c_list and allocate nres contacts in each 
  !contacts array
  if(allocated(c_list)) deallocate(c_list)
  allocate(c_list(nres))
  do i = 1, nres
    allocate(c_list(i)%contacts(nres))
    c_list(i)%contacts = 0
    !add previous and last i to the contacts array for each entry
    if(i > 1) then
      do j = 1, nres
        if(c_list(i)%contacts(j)/=0) cycle
        c_list(i)%contacts(j) = i-1
        exit
      enddo
    endif
    if(i<nres) then
      do j = 1, nres
        if(c_list(i)%contacts(j)/= 0) cycle
        c_list(i)%contacts(j) = i+1
        exit
      enddo
    endif
  enddo    
  !rewind mfile
  rewind(dunit)
  !read through mfile
  do
    read(dunit,'(a)',iostat=ierr) aline
    if(ierr /= 0) exit
    !ignore comment lines
    if(aline(1:1)=='!') cycle
    !ignore empty lines
    if(trim(aline)=='') cycle
    !read in res1 datom res2 aatom bond
    read(aline,*) res1, datom, res2, aatom, bond    
    !add res2 node to res1 contacts
    do i = 1, nres    
      if(inArray(res2,c_list(res1)%contacts,nres)) exit
      if(c_list(res1)%contacts(i) /= 0) cycle
      c_list(res1)%contacts(i) = res2
      exit
    enddo
    !add res1 node to res2 contacts
    do i = 1, nres
      if(inArray(res1,c_list(res2)%contacts,nres)) exit
      if(c_list(res2)%contacts(i) /= 0) cycle
      c_list(res2)%contacts(i) = res1
      exit
    enddo
  enddo
  !close(mfile)
  close(dunit)
  !open cijfile
  dunit = pickunit(dunit)
  open(dunit,file=cijfile,iostat=ierr,status='old',form='formatted')
  if(ierr /= 0) stop "geofold_flory::all_contacts: Error opening cij &
  file."
  !read through cijfile
  do
    read(dunit, '(a)',iostat=ierr) aline
    if(ierr /= 0) exit
    !ignore comment lines
    if(aline(1:1)=='!') cycle
    !ignore empty lines
    if(trim(aline)=='') cycle
    !read in res1 res2
    read(aline,'(2i6)') res1,res2    
    !add res2 to res1 contacts
    do i = 1, nres
      if(inArray(res2,c_list(res1)%contacts,nres)) exit
      if(c_list(res1)%contacts(i) /= 0) cycle
      c_list(res1)%contacts(i) = res2
      exit
    enddo
    !add res1 to res2 contacts
    do i = 1, nres
      if(inArray(res1,c_list(res2)%contacts,nres)) exit
      if(c_list(res2)%contacts(i) /= 0) cycle
      c_list(res2)%contacts(i) = res1
      exit
    enddo
  enddo
  !close(cijfile)
  close(dunit)
!end subroutine  
end subroutine geofold_flory_all_contacts



!adds a node to the visited array

subroutine addValue(a,visited,nres)
  implicit none
  integer :: a, i,nres
  integer, dimension(:),allocatable :: visited
  
  do i = 1, nres
    if(visited(i)/= 0) cycle
    visited(i) = a
    exit
  enddo
end subroutine addValue

!flips an integer array, excluding 0 values

subroutine flip(array, arraySize)
  implicit none
  integer, dimension(:), allocatable :: array,flop
  integer :: arraySize,i,j

  
  allocate(flop(arraySize))
  flop = array
  j = 1
  do i = arraySize, 1, -1
    if(flop(i)==0)cycle
    array(j) = flop(i)
    j = j+1
  enddo
  if(allocated(flop)) deallocate(flop)
end subroutine flip

logical function allVisited(j,c_list,visited,nres) result(answer)
  implicit none
  integer :: j,i,k,nres
  integer,dimension(:),allocatable,intent(in) :: visited
  type(contact),dimension(:),allocatable,intent(in) :: c_list
  logical :: notInArray = .false.
  
  answer = .false.
  do i = 1, nres
    k = c_list(j)%contacts(i)
    if(k==0) exit
    if(.not. inArray(k,visited,nres)) then
      notInArray = .true.
      exit
    endif
  enddo
  if(.not. notInArray) answer = .true.
  
end function allVisited

!finds the shortest path between two contacts
!uses Dijkstra algorithm with a maximum path constraint

!!Dijkstra's algorithm:
!!keep track of visited and unvisited nodes
!!start node cost is 0, all others are initialized to infinity (9999)
!!starting from node a, calculate cost of traveling to all unvisited 
!!neighbors
!!if this cost is lower than that neighbor's current cost, set it to the 
!!lower value
!!move node a from unvisited to visited.  Repeat this procedure from the
! node's lowest cost neighbor
!!When at goal node, stop
!!additionally, this algorithm doesn't search any node where the cost is 
!!greater than the sequence distance between the start and the finish as 
!!that is always possible and much shorter
!Operational 7 October 2014

subroutine shortest_path(a, b, c_list,cost,nres)
  implicit none
  integer :: a,b,cost,maxCost,i,j,current,k,l
  integer,intent(in) :: nres
  integer, dimension(:),allocatable :: visited, unvisited
  type(contact),dimension(:),allocatable :: c_list
  logical :: atgoal
  
  
  atgoal = .false.
  cost = 0
  do i = 1, nres
    c_list(i)%cost=9999
  enddo
  !start from a
  c_list(a)%cost=0
  !set max to |a-b|  
  maxCost = abs(a-b)
  !create data structures of visited and unvisited nodes
  allocate(visited(nres))
  visited = 0
  allocate(unvisited(nres))
  unvisited = 0
  do i = 1, nres
    unvisited(i) = i
  enddo
  !set current to start point
  current = a
  !do while not at goal
  do while(.not.atgoal)
    !if current is goal, atgoal is true
    if(current == b) then
      atgoal = .true.
      exit
    endif
    !exclude paths greater than maxCost |a-b|
    if(c_list(current)%cost>maxCost) then
      call addValue(current,visited,nres)
      unvisited(current)=0
      cycle
    endif
    !calculate distance from current residue to its unvisited neighbors 
    !through current path
    do i = 1, nres
      j = c_list(current)%contacts(i)
      !exit if j is 0, eg at the end of current's neighbors
      if(j == 0) exit
      if(inArray(j,visited,nres)) cycle
      !exclude goal if current is start
      if(current == a .and. j==b) cycle
      !determine cost of neighbors
      cost = c_list(current)%cost+1
      !save cost for neighbors if lower than previous value
      if(cost < c_list(j)%cost) then
        c_list(j)%cost = cost
        !also, adjust path to index of neighbor
        c_list(j)%path = current
      endif
    enddo
    !move this node to visited nodes list
    unvisited(current) = 0
    call addValue(current,visited,nres)
    !set current node to nearest neighbor
    k = maxCost
    l = 0
    do i = 1, nres
      j = c_list(current)%contacts(i)        
      !exclude goal if current is start
      if(current == a .and. j==b) cycle 
      if(j==0) then
        exit
      endif
      if(inArray(j,visited,nres)) then
        cycle
      endif
      if(allVisited(j,c_list,visited,nres)) then
        call addValue(j,visited,nres)
        cycle
      endif
      if(c_list(j)%cost < k) then
        if(j/=0) then
          k = c_list(j)%cost
          l = j
        endif
      endif
    enddo
    if(l/=0) current = l
    !no shorter path found, save straight walkalong
    if(l == 0) then
      cost = maxCost
      if(allocated(visited)) deallocate(visited)
      if(allocated(unvisited)) deallocate(unvisited)
      return
    endif
  !enddo
  enddo
  !reset costs and path pointers of every node
  do i = 1, nres
    c_list(i)%cost = 9999
    c_list(i)%path = 0
  enddo
  !deallocate visited and unvisited
  if(allocated(visited)) deallocate(visited)
  if(allocated(unvisited)) deallocate(unvisited)
end subroutine shortest_path


real function geofold_flory_calc_entropy(n,f,u1,u2,c_list,w,T) result(e)
  implicit none
  real, parameter :: R = 8.314 !J/molK
  type(intermediate),pointer,intent(in) :: f,u1
  type(intermediate),pointer,intent(in),optional :: u2
  integer :: nres,i,j,k
  integer, intent(in) :: n
  integer,dimension(:),allocatable :: path
  integer,dimension(:,:),allocatable :: paths
  real :: shortest,cost3,cost
  real,intent(in) :: w, T !weight for calibrating DCO, Temperature
  type(contact),dimension(:),allocatable,intent(in) :: c_list
  type(contact),dimension(:),allocatable :: broken,subcontacts
  real :: Ln,Ld
  logical,dimension(geofold_nres) :: searched
  
  nres = geofold_nres
  cost = 0.
  if(n == 3) then
    if(present(u2)) then
      cost3 = DCO(c_list,w,f,u1,u2)
    else
      cost3 = DCO(c_list,w,f,u1)
    endif
    if(.not. present(u2)) write(*,'("cost3: ",e10.3)') cost3
    e = 2.1+(1.5*R*log(w*cost3))
  else
    if(n==1) then 
      shortest = 99999.
    else if (n == 2) then
      shortest = 0.
    !else if (n == 3) then
     ! shortest = 99999.
    else
      e = 0
      return
    endif
    !call geofold_flory_get_subcontacts(c_list,subcontacts,f,nres)
    !seam move?
    if(.not. present(u2))then
      call find_contacts_seam(f,u1,c_list,broken,nres,subcontacts)
    else
      call find_contacts_noseam(f,u1,u2,broken,c_list,nres,subcontacts)
    endif
    do i =1, nres
      do j = 1, nres
        if(subcontacts(i)%contacts(j)==0) exit
      enddo
    enddo
    do i = 1, nres
      do j = 1, nres
        if(broken(i)%contacts(j)==0) exit
        k = broken(i)%contacts(j)
        if(n /= 3) then
          !call shortest_path(i,k,subcontacts, cost,nres)
          !newLeff(a,b,subcontacts,wn,T,Ln,Ld, Leff,d)
          searched = .false.
          call newLeff(i,k,subcontacts,w,T,Ln,Ld,cost,searched)
        else
          !Leff(a,b,subcontacts,wn, T)
          !cost3 = Leff(i,k,subcontacts,wn,T)
        endif
        if(n==2 .and. cost > shortest) shortest = cost
        if(n==1 .and. cost < shortest) shortest = cost
      !  if(n==3 .and. cost3 < shortest) shortest = cost3
      enddo
    enddo
    e = 2.1+(1.5*R*log(shortest))
  endif
end function geofold_flory_calc_entropy


end module geofold_flory
