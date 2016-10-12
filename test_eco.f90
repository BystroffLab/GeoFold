program test_eco
    use geofold_global
    use geofold_eco
    implicit none

    call test_printArray()
    call test_isCovalent()
    call test_addValue()
    call test_inArray()


contains

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
    integer :: a,b,nres

    nres = 7
    allocate(geofold_ss(7))
    geofold_ss = 0
    geofold_ss(3) = 6
    geofold_ss(6) = 3
    if(.not. isCovalent(3,6)) stop "isCovalent(3,6) should be true"
    if(.not. isCovalent(6,3)) stop "isCovalent(6,3) should be true"
    if(.not. isCovalent(2,3)) stop "isCovalent(2,3) should be true"
    if(isCovalent(1,6)) stop "isCovalent(1,6) should be false"
    if(isCovalent(2,5)) stop "isCovalent(2,5) should be false"
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

end program test_eco
