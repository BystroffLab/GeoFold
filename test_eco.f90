program test_eco
    use geofold_global
    use vectormath
    use geofold_seams
    use geofold_pivots
    use geofold_hbonds
    use geofold_eco
    implicit none

    call test_printArray()
    call test_isCovalent()
    call test_addValue()
    call test_inArray()
    call test_getcontacts()
    call test_getsubcontacts()
    call test_getbroken()

end program test_eco
