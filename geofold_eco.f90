module geofold_eco
    use geofold_global
    ! private
    public !for testing purposes

    public :: geofold_eco_AvgECO, getECO !getECO only public for testing purposes

    type :: contact
        integer :: cost = 9999
        integer,dimension(:),allocatable :: neighbors
    end type contact

contains

    subroutine getbroken(f,u1,contacts,broken,nres,subcontacts)
        implicit none
        type(contact),dimension(:),allocatable,intent(in) :: contacts
        type(contact),dimension(:),allocatable,intent(out) :: broken,subcontacts
        type(intermediate),pointer,intent(in) :: f,u1
        integer,intent(in) :: nres
        integer :: i,j,seam,k
        character,dimension(maxres) :: u1flag,u2flag

        if(.not. allocated(barrels_array)) stop 'geofold_eco::getbroken: Error, &
        barrels_array not allocated'
        call getsubcontacts(contacts,subcontacts,f,nres)
        !determine which seam is broken
        do i = 1,maxbarrel
            if(u1%barrel(i) == 0) cycle
            seam = u1%barrel(i)
            exit
        enddo
        !set u1flag and u2flag
        u1flag = barrels_array(i)%seams(seam)%u1flag
        u2flag = barrels_array(i)%seams(seam)%u2flag
        !initializations
        if(allocated(broken)) deallocate(broken)
        allocate(broken(nres))
        do i = 1,nres
            allocate(broken(i)%neighbors(nres))
            broken(i)%neighbors = 0
        enddo
        !iterate through f and find the contacts it contains
        do i = 1,nres
            !check if contact exists in f
            do j = 1,nres
                if(subcontacts(i)%neighbors(j) == 0) exit
                k = subcontacts(i)%neighbors(j)
                if(u1flag(i) /= '.') then
                    if(u1flag(k) /= '.') cycle
                    if(.not. inArray(k,broken(i)%neighbors,nres)) then
                        if(k /= i-1 .and. k /= i+1)&
                            call addValue(k,broken(i)%neighbors,nres)
                    endif
                    if(.not. inArray(i,broken(k)%neighbors,nres)) then
                        if(k /= i-1 .and. k /= i+1)&
                            call addValue(i,broken(k)%neighbors,nres)
                    endif
                else if(u2flag(i) /= '.') then
                    if(u2flag(k) /= '.') cycle
                    if(.not. inArray(k,broken(i)%neighbors,nres)) then
                        if(k /= i-1 .and. k /= i+1)&
                            call addValue(k,broken(i)%neighbors,nres)
                    endif
                    if(.not. inArray(i,broken(k)%neighbors,nres)) then
                        if(k /= i-1 .and. k /= i+1)&
                            call addValue(i,broken(k)%neighbors,nres)
                    endif
                endif
            enddo
        enddo
    end subroutine getbroken

    subroutine getsubcontacts(contacts,subcontacts,f,nres)
        implicit none
        integer,intent(in) :: nres
        integer :: i,j,k
        type(contact),dimension(:),allocatable,intent(in) :: contacts
        type(contact),dimension(:),allocatable,intent(out) :: subcontacts
        type(intermediate),pointer,intent(in) :: f

        if(.not. allocated(subcontacts)) then
            allocate(subcontacts(nres))
            do i = 1,nres
                allocate(subcontacts(i)%neighbors(nres))
                subcontacts(i)%neighbors = 0
            enddo
        endif
        do i = 1,nres
            if(f%iflag(i) == '.') cycle
            do j = 1, nres
                if(contacts(i)%neighbors(j) == 0) exit
                k = contacts(i)%neighbors(j)
                if(f%iflag(k) == '.') cycle
                if(.not. inArray(k,subcontacts(i)%neighbors,nres))&
                    call addValue(k,subcontacts(i)%neighbors,nres)
                if(.not. inArray(i,subcontacts(k)%neighbors,nres))&
                    call addValue(i,subcontacts(k)%neighbors,nres)
                k = 0
            enddo
        enddo
    end subroutine getsubcontacts

    subroutine addValue(a, array, size)
        implicit none
        integer,intent(in) :: a,size
        integer :: i
        integer,dimension(:),allocatable,intent(inout) :: array

        do i = 1,size
            if(array(i) /= 0) cycle
            array(i) = a
            exit
        enddo
    end subroutine addValue

    subroutine getcontacts(mfile,cijfile,contacts,nres)
        implicit none
        type(contact),dimension(:),allocatable,intent(out) :: contacts
        character (len=*),intent(in) :: mfile,cijfile
        integer :: res1,res2,ierr,dunit,nres,i,j
        real :: tmp
        character (len=3) :: datom, aatom
        character (len=1) :: bond
        character (len=200) :: aline

        !determine nres
        !open mfile (hbonds file)
        open(newunit=dunit,file=mfile,iostat=ierr,status='old',form='formatted')
        if(ierr /= 0) stop "geofold_eco::getcontacts: Error opening hbonds file"
        !find line 8
        do i = 1,8
            read(dunit,'(a)') aline
        enddo
        !read aline(2:10) as (i9) into nres
        read(aline(2:10),'(i9)',iostat=ierr) nres
        if(ierr /= 0) stop 'geofold_eco_getcontacts: Error reading nres'
        !allocate nres contacts in contacts and allocate nres contacts in each
        !contacts array
        if(allocated(contacts)) deallocate(contacts)
        allocate(contacts(nres))
        do i = 1,nres
            allocate(contacts(i)%neighbors(nres))
            contacts(i)%neighbors = 0
            !add previous and last i to the contacts array for each entry
            if(i > 1)  call addValue(i-1,contacts(i)%neighbors,nres)
            if(i < nres) call addValue(i+1,contacts(i)%neighbors,nres)
        enddo
        !rewind mfile
        rewind(dunit)
        !read through mfile
        do
            read(dunit,'(a)',iostat=ierr) aline
            if(ierr /= 0) exit
            !ignore comment lines
            if(aline(1:1) == '!') cycle
            !ignore empy lines
            if(trim(aline) == '') cycle
            !read in res1 datom res2 aatom bond
            read(aline,*) res1,datom,res2,aatom,bond
            !add res2 node to res1 contacts
            if(.not. inArray(res2,contacts(res1)%neighbors,nres))&
                call addValue(res2,contacts(res1)%neighbors,nres)
            !add res1 node to res2 contacts
            if(.not. inArray(res1,contacts(res2)%neighbors,nres))&
                call addValue(res1,contacts(res2)%neighbors,nres)
        enddo
        !close(mfile)
        close(dunit)
        !open cijfile
        open(newunit=dunit,file=cijfile,iostat=ierr,status='old',form='formatted')
        if(ierr /= 0) stop 'geofold_eco::getcontacts: Error opeing cij file.'
        !read through cijfile
        do
            read(dunit,'(a)',iostat=ierr) aline
            if(ierr /= 0) exit
            !ignore comment lines
            if(aline(1:1)=='!') cycle
            !ignore empty lines
            if(trim(aline) == '') cycle
            !read in res1 res2
            read(aline,'(2i6)') res1,res2
            !add res2 to res1 contacts
            if(.not. inArray(res2,contacts(res1)%neighbors,nres))&
                call addValue(res2,contacts(res1)%neighbors,nres)
            !add res1 to res2 contacts
            if(.not. inArray(res1,contacts(res2)%neighbors,nres))&
                call addValue(res1,contacts(res2)%neighbors,nres)
        enddo
        !close(cijfile)
        close(dunit)
    end subroutine getcontacts

    integer function getECO(i,j,contacts,broken,nres) result(ECO)
        implicit none
        integer,intent(in) :: i,j,nres
        integer :: maxCost,cost,current,a,b,k,l
        integer,dimension(:),allocatable :: visited,unvisited
        type(contact),dimension(:),allocatable :: contacts,broken
        logical :: atGoal

        atGoal = .false.
        ECO = 0
        do a = 1,nres
            contacts(a)%cost = 9999
        enddo
        !start from i
        contacts(i)%cost=0
        !set max to |i-j|
        maxCost = abs(i-j)
        !create data structures for visited and unvisited nodes
        allocate(visited(nres))
        visited = 0
        allocate(unvisited(nres))
        unvisited = 0
        do a = 1,nres
            unvisited(i) = a
        enddo
        !set current to start point
        current = a
        !do while not at goal
        do while(.not. atGoal)
            !if current is goal, atGoal is true
            if(current == j) then
                atGoal = .true.
                exit
            endif
            !exclude paths greater than maxCost
            if(contacts(current)%cost > maxCost) then
                call addValue(current,visited,nres)
                unvisited(current) = 0
                cycle
            endif
            !calculate distance from current residue to its unvisited neighbors
            !through current path.  CAN'T USE BROKEN CONTACTS
            do a = 1, nres
                b = contacts(current)%neighbors(a)
                if(b == 0) exit
                !Exclude broken contacts
                if(inArray(b,broken(a)%neighbors,nres)) cycle
                !Exclude visited contacts
                if(inArray(b,visited,nres)) cycle
                !Exclude goal if current is start
                if(current == i .and. b == j) cycle
                !determine cost of neighbors
                cost = contacts(current)%cost + 1
                !save cost for neighbors if lower than previous value
                if(cost < contacts(b)%cost) contacts(b)%cost = cost
            enddo
            !move this node to visited nodes list
            unvisited(current) = 0
            call addValue(current,visited,nres)
            !set current node to nearest neighbor
            k = maxCost
            l = 0
            do a = 1,nres
                b = contacts(current)%neighbors(a)
                if(b == 0) exit
                !Exclude gola if current is start
                if(current == i .and. b == j) cycle
                if(inArray(b,visited,nres)) cycle
                if(inArray(b,broken(a)%neighbors,nres)) cycle
                if(allVisited(b,contacts,visited,nres)) then
                    call addValue(b,visited,nres)
                    cycle
                endif
                if(contacts(b)%cost < k) then
                    if(b /= 0) then
                        k = contacts(b)%cost
                        l = b
                    endif
                endif
            enddo
            if(l /= 0) current = l
            !no shorter path found, save straight walkalong
            if(l == 0) then
                ECO = maxCost
                if(allocated(visited)) deallocate(visited)
                if(allocated(unvisited)) deallocate(unvisited)
                return
            endif
        enddo
        !reset costs of every node
        do a = 1,nres
            contacts(i)%cost = 9999
        enddo
        if(allocated(visited)) deallocate(visited)
        if(allocated(unvisited)) deallocate(unvisited)
    end function getECO

    !works
    logical function isCovalent(a,b)
        implicit none
        integer, intent(in) :: a,b

        if(abs(a-b) == 1 .or. geofold_ss(a) == b) then
            isCovalent = .true.
        else
            isCovalent = .false.
        endif
    end function isCovalent

    !works
    subroutine printArray(array)
        implicit none
        integer,dimension(:,:),allocatable :: array
        integer :: i,j

        write(*,'(i3,"x",i3," array")') size(array,1),size(array,2)
        do i = 1, size(array,1)
            write(*,*) array(i,1:size(array,2))
        enddo
    end subroutine printArray

    logical function inArray(a,array,arraySize) result(answer)
        implicit none
        integer,intent(in) :: a,arraySize
        integer,dimension(:),intent(in) :: array
        integer :: i

        answer = .false.
        do i = 1, arraySize
            if(array(i) == a) then
                answer = .true.
                exit
            endif
            if(array(i) == 0) exit
        enddo
    end function inArray

    logical function allVisited(j,contacts,visited,nres) result(answer)
        implicit none
        integer :: j,i,k,nres
        integer,dimension(:),allocatable,intent(in) :: visited
        type(contact),dimension(:),allocatable,intent(in) :: contacts
        logical :: notInArray = .false.

        answer = .false.
        do i = 1,nres
            k = contacts(j)%neighbors(i)
            if(k == 0) exit
            if(.not. inArray(k,visited,nres)) then
                notInArray = .true.
                exit
            endif
        enddo
        if(.not. notInArray) answer = .true.
    end function allVisited

    real function getAvgECO(contacts,broken,nres) result(avgECO)
        implicit none
        type(contact),dimension(:),allocatable,intent(in) :: contacts,broken
        integer,intent(in) :: nres
        integer :: count, sumECO,i,j,k

        sumECO = 0
        count = 0

        do i = 1,nres
            do j = 1,nres
                k = broken(i)%neighbors(j)
                if(k == 0) cycle
                sumECO = sumECO + getECO(i,k,contacts,broken,nres)
                count = count + 1
            enddo
        enddo
        avgECO = real(sumECO)/real(count)
    end function getAvgECO

    real function geofold_eco_AvgECO(mfile,cijfile,f,u1,nres) result(avgECO)
        implicit none
        type(contact),dimension(:),allocatable :: contacts,broken,subcontacts
        integer,intent(in) :: nres
        character (len=*),intent(in) :: mfile,cijfile
        type(intermediate),pointer,intent(in) :: f,u1

        call getcontacts(mfile,cijfile,contacts,nres)
        call getbroken(f,u1,contacts,broken,nres,subcontacts)
        avgECO = getAvgECO(contacts,broken,nres)
    end function geofold_eco_AvgECO

end module geofold_eco
