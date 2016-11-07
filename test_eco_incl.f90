subroutine test_printArray()
    implicit none
    integer,dimension(:,:),allocatable :: array
    integer,parameter :: seed = 8311990
    integer :: i,j

    call srand(seed)

    write(*,*) "allocating a",10,"by",5,"array..."
    if(allocated(array)) deallocate(array)
    allocate(array(10,5))
    do i = 1,10
        do j = 1,5
            array(i,j) = int(rand()*10)
            ! write(*,*) array(i,j)
        enddo
        ! write(*,*) ''
    enddo
    call printArray(array)
    deallocate(array)
end subroutine test_printArray

subroutine test_isCovalent()
    implicit none
    integer :: a,b
    character, dimension(maxres) :: iflag

    iflag = 'A'
    iflag(7) = 'B'

    geofold_nres = 7
    allocate(geofold_ss(7))
    geofold_ss = 0
    geofold_ss(3) = 6
    geofold_ss(6) = 3
    if(.not. isCovalent(3,6,iflag)) stop "isCovalent(3,6) should be true"
    if(.not. isCovalent(6,3,iflag)) stop "isCovalent(6,3) should be true"
    if(.not. isCovalent(2,3,iflag)) stop "isCovalent(2,3) should be true"
    if(isCovalent(1,6,iflag)) stop "isCovalent(1,6) should be false"
    if(isCovalent(2,5,iflag)) stop "isCovalent(2,5) should be false"
    if(isCovalent(6,7,iflag)) stop "isCovalent(6,7) should be false"
    write(*,*) "test_isCovalent passed!"
    deallocate(geofold_ss)
end subroutine test_isCovalent

subroutine test_addValue()
    implicit none
    integer,dimension(:),allocatable :: array

    write(*,*) 'allocating a 10 by 1 array'
    allocate(array(10))
    array = 0
    array(1) = 5
    write(*,*) "adding 7 to array.  Should be added to array(2)"
    call addValue(7,array,10)
    write(*,*) "array(2)",array(2)
    deallocate(array)
end subroutine test_addValue

subroutine test_inArray()
    implicit none
    integer,dimension(:),allocatable :: array

    allocate(array(10))
    array = 0
    array(1:6) = 2
    array(7) = 5
    if(.not. inArray(5,array,10)) stop '5 should be in array'
    if(inArray(7,array,10)) stop '7 should not be in array'
    write(*,*) "test_inArray successful!"
    deallocate(array)
end subroutine test_inArray


subroutine test_getcontacts()
    implicit none
    type(contact),dimension(:),allocatable :: contacts
    integer :: dunit, i
    character (len=1000) :: pdbfile,hbfile,cijfile

    write (*,*) "======== STARTING TEST_GETCONTACTS ========"

    pdbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.pdb"
    hbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.hb"
    cijfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.cij"
    write (*,*) "cijfile:",cijfile
    write (*,*) "pdbfile:",pdbfile
    write (*,*) "hbfile:",hbfile

    ! stop "TEST STOP"

    open(newunit=dunit,file=pdbfile)
    call geofold_readpdb(dunit)
    close(dunit)
    if(geofold_nres /= 20) then
        write (*,*) "geofold_nres =",geofold_nres
        stop "geofold_nres should be 20"
    endif
    call geofold_hbonds_read(hbfile)
    if(size(geofold_hb,2) /= 20) then
        write (*,*) "size(geofold_hb,2):",size(geofold_hb,2)
        do i = 1, size(geofold_hb)
            write (*,*) i, geofold_hb(1,i),geofold_hb(2,i)
        enddo
        ! write(*,*) "geofold_hb =",geofold_hb
        stop "geofold_hb should be 20"
    endif
    do i = 1, geofold_nres
        if(geofold_ss(i) /= 0) then
            write(*,'("geofold_ss(",i2,") = ",i2)') i,geofold_ss(i)
            stop "There should be no disulfides"
        endif
    enddo
    call getcontacts(cijfile,contacts)
    !use inArray to test expected values
    !cij contact
    if(.not. inArray(19,contacts(17)%neighbors,geofold_nres)) then
        stop "Contact should exist between 19 and 17"
    endif
    !hb contact
    if(.not. inArray(6,contacts(17)%neighbors,geofold_nres)) then
        stop "Contact should exist between 6 and 17"
    endif
    !not a contact
    if (inArray(1,contacts(19)%neighbors,geofold_nres)) then
        stop "Contact should not exist between 1 and 19"
    endif

    write(*,*) "test_getcontacts completed!"

end subroutine test_getcontacts

subroutine test_getsubcontacts()
    implicit none
    type(contact),dimension(:),allocatable :: contacts,subcontacts
    type(intermediate), pointer :: f
    character (len=1000) :: cijfile,pdbfile,hbfile
    integer :: dunit,i

    write (*,*) "======= TEST_GETSUBCONTACTS ========"

    !SETUP
    pdbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.pdb"
    hbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.hb"
    cijfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.cij"
    open(newunit=dunit,file=pdbfile)
    call geofold_readpdb(dunit)
    close(dunit)
    call geofold_hbonds_read(hbfile)
    call getcontacts(cijfile,contacts)

    !Setup hypothetical F
    allocate(f)
    f%idnum = 6
    nullify(f%next)
    f%state = pivotflag
    f%sym = 0
    f%barrel = 0
    f%iflag = '.'
    do i = 7,15
        f%iflag(i) = 'A'
    enddo

    call getsubcontacts(contacts,subcontacts,f)
    !contact that should still exist
    if(.not. inArray(9,subcontacts(7)%neighbors,geofold_nres)) &
        stop "7 and 9 should still be in contact"
    !contact that should no longer exist
    if(inArray(5,subcontacts(7)%neighbors,geofold_nres)) &
        stop "contact between 5 and 7 should be broken"
    write (*,*) "test_subcontacts passed"
end subroutine test_getsubcontacts

subroutine test_getbroken()
    implicit none
    type(intermediate),pointer :: f,u1,u2
    type(tstate),pointer :: transition
    type(contact),dimension(:),allocatable :: contacts,broken,subcontacts
    character (len=1000) :: cijfile,pdbfile,hbfile
    integer :: dunit,i

    write (*,*) "======== TEST_GETBROKEN ========"

    !setup
    pdbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.pdb"
    hbfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.hb"
    cijfile = "/Users/walcob/GeoFold/test_1l2y/tmp/1L2Y_PC25.cij"
    open(newunit=dunit,file=pdbfile)
    call geofold_readpdb(dunit)
    close(dunit)
    call geofold_hbonds_read(hbfile)
    call getcontacts(cijfile,contacts)

    allocate(f)
    allocate(u1)
    allocate(u2)
    allocate(transition)
    allocate(barrels_array(0))

    !setup f
    f%idnum = 2
    f%iflag = '.'
    f%iflag(1:15) = 'A'
    f%state = pivotflag
    f%sym = 0
    f%axis = 0
    f% barrel = 0
    nullify(f%next)

    !u1
    u1%idnum = 3
    u1%iflag = '.'
    u1%iflag(1:6) = 'A'
    u1%state = pivotflag
    u1%sym = 0
    u1%axis = 0
    u1% barrel = 0
    nullify(u1%next)

    !u2
    u2%idnum = 6
    u2%iflag = '.'
    u2%iflag(7:15) = 'A'
    u2%state = pivotflag
    u2%sym = 0
    u2%axis = 0
    u2% barrel = 0
    nullify(u2%next)

    !transition
    transition%id = 3
    transition%tp = pivotflag
    transition%entropy = 0.41
    transition%parent = 2
    transition%child1 = 3
    transition%child2 = 6
    nullify(transition%next)
    transition%energy = 0.
    transition%axis = 0
    transition%seam = 0

    call getbroken(f,u1,u2,transition,contacts,broken,subcontacts)



    !contact should be broken
    if(.not. inArray(5,broken(7)%neighbors,geofold_nres)) &
        stop "contact between 5 and 7 should be broken"

    !contact should not be broken
    if(inArray(9,broken(7)%neighbors,geofold_nres)) &
        stop "7 and 9 should still be in contact"

    write (*,*) "test_getbroken pass"
end subroutine test_getbroken
