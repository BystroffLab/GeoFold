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
  
  type :: pathNode
    integer :: sz = 0, maxsz = 0
    integer,dimension(:),pointer :: points => null()
    type(pathNode),pointer :: next => null(),prev => null()
  end type pathNode
  
  type :: Path
    type(pathNode),pointer :: root => null()
    integer :: sz = 0
  end type Path

!don't forget that geofold_nres is a thing
  
  
  


  public :: geofold_flory_calc_entropy,geofold_flory_get_subcontacts
  public :: contact, geofold_flory_all_contacts, Path
  
interface find_contacts
  module procedure find_contacts_noseam,find_contacts_seam
endinterface
  
contains

!cleanup subroutine for Path type
!Works #BDW 20/2/2015
subroutine cleanPath(paths)
  implicit none
  type(Path) :: paths
  type(pathNode),pointer :: itr,jtr
  
  itr => paths%root
  do while(associated(itr%next))
    jtr => itr%next
    deallocate(itr%points)
    deallocate(itr)
    itr => jtr
  enddo
  deallocate(itr%points)
  deallocate(itr)
  paths%sz = 0
  nullify(paths%root)
end subroutine cleanPath

!add new node to path given template (root), Path, and next step
!Works #BDW 19/2/2015
subroutine addStep(node, paths, step)
  implicit none
  integer :: step
  type(pathNode),pointer :: node,newNode
  type(Path) :: paths
  
  if(node%sz < node%maxsz) then
    !initialize newNode
    allocate(newNode)
    newNode%sz = node%sz
    newNode%maxsz = node%maxsz
    allocate(newNode%points(newNode%maxsz))
    !up until the new point newNode should be identical to node
    newNode%points = node%points

    !find the insertion point and insert step
    newNode%points(newNode%sz+1) = step
    newNode%sz = newNode%sz + 1

    !adjust pointers
    if(associated(node%next)) then
      newNode%next => node%next
    else
      nullify(newNode%next)
    endif
    newNode%prev => node
    node%next => newNode
    paths%sz = paths%sz + 1
  endif
end subroutine addStep

!remove a path from the list of paths if it does not meet the end goal
!Works #BDW 23/2/2015
subroutine removePath(paths, node)
  implicit none
  type(Path),intent(inout) :: paths
  type(pathNode),pointer,intent(inout) :: node
  type(pathNode),pointer :: itr,jtr


  if(associated(node,target=paths%root))then
    if(associated(node%next))then
      itr=>node%next
      deallocate(node%points)
      deallocate(node)
      nullify(node)
      paths%root=>itr
      nullify(itr%prev)
      paths%sz = paths%sz - 1
    else
      !only item in path
      deallocate(node%points)
      deallocate(node)
      nullify(node)
      nullify(paths%root)
      paths%sz = paths%sz - 1 !Should be zero
    endif
  else
  !The bug is with itr and node%prev....
    !itr => node%prev
    itr => paths%root
    do while(.not. associated(itr%next,target=node))
      itr=>itr%next
    enddo
    if(associated(node%next)) then
      jtr => node%next
    else
      nullify(jtr)
    endif
    deallocate(node%points)
    deallocate(node)
    nullify(node)
    itr%next => jtr
    if(associated(jtr)) then
      jtr%prev => itr
    else
      nullify(itr%next)
    endif
    paths%sz = paths%sz - 1
  endif
end subroutine removePath



!Works #BDW 19/2/2015
logical function isCovalent(a,b)
  implicit none
  integer :: a,b
  
  if(abs(a-b) == 1 .or. geofold_ss(a) == b) then
    isCovalent = .true.
  else
    isCovalent = .false.
  endif
end function isCovalent

subroutine fprintPaths(paths)
    implicit none
    type(Path) :: paths
    integer :: i,j1
    type(pathNode),pointer :: itr

    write(*,'( "Printing ",i9," paths...")')paths%sz
    itr => paths%root
    i = 0
    do while(associated(itr))! .and. i<paths%sz)
      i = i + 1
      write(6,'("Path ",i3," size = ",i9,":")')i,itr%sz
      write(6,*) itr%points
      flush(6)
      itr => itr%next
    enddo
  end subroutine fprintPaths

!find all valid paths between contacts a and b
!a path is valid if it is <= default loop length
!and doesn't loop over itself
!Works #BDW 25.02.2015
subroutine findPaths(a,b,paths,subcontacts)
  implicit none
  integer :: a,b,flory
  type(Path),intent(inout) :: paths
  type(contact),dimension(:),allocatable,intent(in) :: subcontacts
  !iterators
  integer :: i,j
  type(pathNode),pointer :: itr,jtr
  integer :: point,contct

  !initializations
  flory = abs(a-b)+1 !the loop length considering no other contacts, this is also the max path size
  !allocate and initialize first path
  if(associated(paths%root)) deallocate(paths%root)
  allocate(paths%root)
  allocate(paths%root%points(flory))
  paths%root%points = 0
  paths%root%points(1) = a
  paths%root%sz = 1
  paths%root%maxsz = flory
  paths%sz = 1

  !main loop of the subroutine
  do i = 1, flory
    itr => paths%root
    do while(associated(itr))
      !check if path is atGoal and is size i
      if(itr%sz == i .and. .not. atGoal(itr,b)) then
        !start at last non-zero value of each path not atGoal
        point = itr%points(i)
        !check all of its neighbors if they've been used in the path already
        j = 1
        contct = subcontacts(point)%contacts(j)
        do while(contct /= 0)
          if(contct == b .and. i == 1) then
            j = j + 1
            contct = subcontacts(point)%contacts(j)
            cycle
          endif
          !create new path for all that haven't been used
          if(.not. inPath(contct,itr)) then
            call addStep(itr,paths,contct)
          endif
          j = j + 1
          contct = subcontacts(point)%contacts(j)
        enddo
      endif
      !repeat
      itr=>itr%next
    enddo
  enddo
  itr => paths%root
  do while(associated(itr))
    if(.not. atGoal(itr,b)) then
      jtr => itr%next
      call removePath(paths,itr)
      if(.not. associated(jtr)) exit
      itr => jtr
    else
      itr => itr%next
    endif
  enddo
  !congratulations, you should have all the valid paths now!

end subroutine findPaths  

!determine the weighted average of the pathlength
!using Boltzmann weighting
real function getLeff(wc,wn,paths,T) result(Leff)
  implicit none
  real,intent(in) :: wc, wn
  type(Path),intent(in) :: paths
  real :: Ln, Ld
  real,dimension(:),pointer :: Li
  real, parameter :: R = 8.314
  real,intent(in) :: T
  type(pathNode),pointer :: itr
  integer :: i,j,a,b

  
  if(associated(Li)) deallocate(Li)
  allocate(Li(paths%sz))
  Li = 0.
  Ln = 0.
  Ld = 0.
  
  itr => paths%root
  do i = 1, paths%sz
    do j = 1, itr%sz-1
      a = itr%points(j)
      b = itr%points(j+1)
      write(*,'("a: ",i4," b: ",i4)') a,b
      if(isCovalent(a,b)) then
        Li(i) = Li(i) + wc
      else
        Li(i) = Li(i) + wn
      endif
    enddo
    Ln = Ln + Li(i)*exp(-(Li(i))/(R*T))
    Ld = Ld + exp(-(Li(i))/(R*T))
    itr => itr%next
  enddo
  Leff = Ln/Ld
end function getLeff
        
      
  
  

!determine if a path has reached its goal
!Works #BDW 19/2/2015
logical function atGoal(root,a)
  implicit none
  type(pathNode),pointer :: root
  integer :: a,i

  atGoal = .false.

  do i = 1, root%sz
    if(root%points(i)==a) then
      atGoal = .true.
      exit
    endif
  enddo
end function atGoal


!test to see if path already contains a residue
!Works #BDW 19/2/2015
logical function inPath(a,root)
  implicit none
  integer :: a,i
  type(pathNode),pointer :: root

  inPath = .false.
  do i = 1, root%sz
    if (root%points(i)==a) then
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

subroutine find_contacts_seam(f,u1,c_list,broken,nres)
  implicit none
  type(contact), dimension(:),allocatable, intent(in) :: c_list
  type(contact),dimension(:), allocatable,intent(out) :: broken
  type(contact),dimension(:),allocatable :: subcontacts
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
  integer,dimension(:),allocatable,intent(in) :: array
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

subroutine find_contacts_noseam(f,u1,u2,broken, c_list,nres)
  implicit none
  type(contact),dimension(:),allocatable, intent(in) :: c_list
  type(contact),dimension(:),allocatable, intent(out) :: broken
  type(contact),dimension(:),allocatable :: subcontacts
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
    ! dunit = pickunit(dunit)
    open(newunit=dunit,file=mfile,iostat=ierr,status='old',form='formatted')
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
  ! dunit = pickunit(dunit)
  open(newunit=dunit,file=cijfile,iostat=ierr,status='old',form='formatted')
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
  deallocate(flop)
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
      deallocate(visited)
      deallocate(unvisited)
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
  deallocate(visited)
  deallocate(unvisited)
end subroutine shortest_path


real function geofold_flory_calc_entropy(n,f,u1,u2,c_list,nres) result(e)
  implicit none
  real, parameter :: R = 8.314 !J/molK
  type(intermediate),pointer,intent(in) :: f,u1
  type(intermediate),pointer,intent(in),optional :: u2
  integer :: nres,i,j,cost=0,k
  integer, intent(in) :: n
  integer,dimension(:),allocatable :: path
  real :: shortest
  type(contact),dimension(:),allocatable,intent(in) :: c_list
  type(contact),dimension(:),allocatable :: broken,subcontacts
  
  if(n==1) then 
    shortest = 99999.
  else if (n == 2) then
    shortest = 0.
  else
    e = 0
    return
  endif
  call geofold_flory_get_subcontacts(c_list,subcontacts,f,nres)
  !seam move?
  if(.not. present(u2))then
  call find_contacts_seam(f,u1,c_list,broken,nres)
  else
  call find_contacts_noseam(f,u1,u2,broken,c_list,nres)
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
      call shortest_path(i,k,subcontacts, cost,nres)
      if(n==2 .and. cost > shortest) shortest = cost
      if(n==1 .and. cost < shortest) shortest = cost
    enddo
  enddo
  e = 2.1+(1.5*R*log(shortest))
end function geofold_flory_calc_entropy


end module geofold_flory
