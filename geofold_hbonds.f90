!!====================================================================================
!! geofold_hbonds
!!====================================================================================
MODULE geofold_hbonds
 use geofold_global
 use geofold_pivots
 private
 public ::  geofold_hbonds_get, geofold_hbonds_read, geofold_hbonds_cleanup
 public ::  geofold_hbonds_disulfidepresent, geofold_hbonds_ssinseg
 real,public :: geofold_hbonds_eperbond

 interface geofold_hbonds_get
   module procedure geofold_hbonds_getbetween
   module procedure geofold_hbonds_getwithin
 endinterface


 CONTAINS
!!----------------------------------------------------------------------
!! Get hbonds between two sides of a seam
!!
  subroutine geofold_hbonds_getbetween(aseam,hbonds,seamchar)
    implicit none
    type(seam_type),pointer :: aseam
    character,dimension(MAXRES) :: f1,f2
    integer,intent(out) :: hbonds
    character,intent(in),optional :: seamchar
    integer :: n,ires,jres,da,h,kres
    n = 0
    if (present(seamchar)) then
      do ires=1,geofold_nres
        if (aseam%u1flag(ires)/=seamchar) cycle
        do da=1,2    !da = donor/acceptor
          do h = 1, size(geofold_hb, 2)
            jres = geofold_hb(da,h)
            if (jres/=ires) cycle
            kres = geofold_hb((2/da),h) ! note 2/da flips 1 to 2, 2 to 1
            if (aseam%u2flag(kres)/=seamchar) cycle
            n = n + 1
          enddo
        enddo
      enddo
    else
      do ires=1,geofold_nres
        if (aseam%u1flag(ires)==".") cycle
        do da=1,2
          do h = 1, size(geofold_hb, 2)
            jres = geofold_hb(da,h)
            if (jres/=ires) cycle
            kres = geofold_hb((2/da),h) ! note 2/da flips 1 to 2, 2 to 1
            if (aseam%u2flag(kres)==".") cycle
            n = n + 1
          enddo
        enddo
      enddo
    endif
    hbonds = n   ! every H-bond is counted once only
  end subroutine geofold_hbonds_getbetween
!!====================================================================================
! get_hbonds
! count the number of hbonds within intermediate f
!
  subroutine geofold_hbonds_getwithin(f, hbonds,seamchar)
    implicit none
    type(intermediate), POINTER :: f
    integer,intent(out) :: hbonds
    character,intent(in),optional :: seamchar
    integer :: n,ires,jres,da,h,kres,bar_itr,seam_itr
    character, dimension(maxres) :: u1flag, u2flag
    n = 0
    ILOOP: do ires=1,geofold_nres
      if (f%iflag(ires)==".") cycle ILOOP
      do da=1,2
        do h = 1, size(geofold_hb, 2)
          jres = geofold_hb(da,h)
          if (jres/=ires) cycle
          kres = geofold_hb((2/da),h) !kres is the acceptor/donor residue where ires is the donor/acceptor
          if (f%iflag(kres)==".") cycle
          if (geofold_pivots_queryinseam(f, ires,kres)) cycle
          n = n + 1
        enddo
      enddo
    enddo ILOOP
!    
!    !Account for H-bonds broken by a seam move
!    !for each barrel in f%barrel
!    do bar_itr = 1, MAXBARREL
!      !check if barrel /= 0 --> open seam
!      if(f%barrel(bar_itr) == 0) cycle      
!      !if /= 0, find and remove bonds broken by seam
!      !iterate through seams to find the one with the proper id
!      do seam_itr = 1, barrels_array(bar_itr)%nSeams
!        if(barrels_array(bar_itr)%seams(seam_itr)%id /= f%barrel(bar_itr)) cycle
!        !set u1flag and u2flag
!        u1flag = barrels_array(bar_itr)%seams(seam_itr)%u1flag
!        u2flag = barrels_array(bar_itr)%seams(seam_itr)%u2flag
!        !iterate through every residue
!        SLOOP: do ires = 1, geofold_nres
!          if(f%iflag(ires)==".") cycle SLOOP
!          !a contact is broken if both residues are in opposite flags
!          do da = 1,2
!            do h = 1, size(geofold_hb, 2)
!              if(geofold_hb(da,h) /= ires) cycle
!              kres = geofold_hb((2/da),h)
!              if (f%iflag(kres)==".") cycle
!              if (geofold_pivots_queryinseam(f,ires,kres)) cycle
!              if (u1flag(ires) /= ".") then
!                if(u2flag(kres) /= ".") n = n-1
!                write (*,'("Bond removed", 2i5)') ires, kres
!              endif
!              if(u1flag(kres) /= ".") then
!                if(u2flag(ires) /= ".") n = n-1
!                write (*,'("Bond removed", 2i5)') ires, kres
!              endif
!            enddo
!          enddo
!        enddo SLOOP
!      enddo
!    enddo
    hbonds = n/2   ! because every H-bond is counted twice !
  end subroutine geofold_hbonds_getwithin

 !!----------------------
 !! read Hbonds file
 SUBROUTINE geofold_hbonds_read(mfile)
  implicit none
  ! integer,intent(in) :: harg
  !CHARACTER (len=*) :: mfile
  character (len=1000) :: mfile
  INTEGER :: res1, res2, ierr, dunit, n, ss, nhbonds = 0
  REAL :: engy, x
  CHARACTER (len=200) :: aline
  character :: bond !! "H" or "S" for Hbond or disulfide
  character(len=3) :: datom,aatom
  integer :: minres, maxres
  dunit = 34
  minres = 999
  maxres = -999
  if (allocated(geofold_hb)) deallocate(geofold_hb)
  allocate(geofold_ss(geofold_nres),stat=ierr)
  if (ierr/=0) stop 'read_hbonds:: error allocating geofold_ss'
  dunit = pickunit(dunit)
  write(0,'("dunit ",i3,", mfile ",a,", ierr ",i3)') dunit,mfile,ierr
  open(dunit, file = mfile, iostat=ierr, status="old", form="formatted")
  IF (ierr /= 0 ) STOP "geofold:: geofold_hbonds_read:  error opening file!"
  do
    read(dunit, '(a)', iostat=ierr) aline
    write (*,*) aline
    if(ierr/=0) exit
    if (aline(1:1)=='!') cycle
    if (trim(aline)=="") cycle
    read(aline,*, iostat=ierr) res1, datom, res2, aatom, bond
    if(ierr/=0) stop 'parsing error 1'
    if (bond == "H" .or. bond == "h") then
      nhbonds = nhbonds+1
    endif
  enddo
  allocate(geofold_hb(2,nhbonds),stat=ierr)
  if(ierr/=0) stop 'read_hbonds:: error allocating geofold_hb'
  rewind(dunit)

  geofold_hb = 0.
  geofold_ss = 0.
  n = 0
  ss = 0
  DO
     !! res1 is donor, res2 is acceptor
     !! geofold_hbond_hb(1,i) is donor to i
     !! geofold_hbond_hb(2,i) is acceptor to i
     !! i is the hydrogen bond number
     read(dunit,'(a)',iostat=ierr) aline
     write (*,*) aline
     IF (ierr /=0) then
      exit
     endif
     if (aline(1:1)=='!') cycle
     read(aline,*, iostat=ierr) res1, datom, res2, aatom, bond
     IF (ierr /=0) then
      exit
     endif
     select case (bond)
     case ("H","h")
       n = n + 1
       geofold_hb(1, n) = res1
       geofold_hb(2, n) = res2
       n = n + 1
     case ("S","s")
       geofold_ss(res1) = res2
       geofold_ss(res2) = res1
       ss = ss + 1
     end select
     if (res1 < minres) minres = res1
     if (res2 < minres) minres = res2
     if (res1 > maxres) maxres = res1
     if (res2 > maxres) maxres = res2
  END DO
  !!! segfault on following line?  It apparently seems to be caused by the trim() function
  write(0,'("Number of HBONDs read from ",a,i8)') mfile(1:20), n
  write(0,'("Number of SS-bonds read from ",a,i8)') mfile(1:20), ss
  write(0,'("HBOND contact data found for residue range ",2i8)') minres, maxres
  write(0,'("Last Hbond  ",2i8)') res1,res2
  close(dunit) ! close hbonds file
END SUBROUTINE geofold_hbonds_read

!!====================================================================================
  subroutine geofold_hbonds_cleanup()
    if (allocated(geofold_hb)) deallocate(geofold_hb)
    if (allocated(geofold_ss)) deallocate(geofold_ss)
  end subroutine geofold_hbonds_cleanup
!!====================================================================================
end MODULE geofold_hbonds
