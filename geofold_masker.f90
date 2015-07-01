!!====================================================================================
!! geofold_masker
!!====================================================================================
MODULE geofold_masker
 use geofold_global
 use geofold_hbonds
 use geofold_pivots
 use geofold_seams
 private
 REAL, dimension(:,:), allocatable :: sasnrg,sasvalues !contains energies of all the c-a residues
 real,dimension(20)::lambda=(/ &
      0.00, 0.55, 1.25, 1.81, 0.58, 0.00, 0.96, 0.89, 1.94, 0.78, &
      1.61, 1.57, 0.00, 2.11, 2.03, 1.71, 1.63, 0.51, 0.97, 0.98/)
 real,dimension(20)::maxsas=(/ &
      107.05,136.83,150.35,170.56,211.48,90.30,168.59,154.89,176.46,156.33, &
      168.88,149.62,127.45,203.74,200.29,120.02,137.17,136.02,205.44,191.31/)
 REAL, dimension(3, maxres) :: voidcoords      
 real,public :: geofold_masker_omega=0.0   !! solvation free energy
 real,public :: geofold_masker_lambdaweight=0.0   !! sidechain entropy weight
 real,public :: geofold_masker_spervoid=0.0   !! void entropy per void.
 type voidtype
     INTEGER :: nn
     integer,dimension(:),pointer :: wall
 end type voidtype
 type(voidtype),dimension(maxres) :: avoid
 integer :: geofold_masker_nvoid
 interface geofold_masker_getscenergy
   module procedure getscenergy
   module procedure getscenergy_seam
 endinterface
 public :: geofold_masker_energy , geofold_masker_read, geofold_masker_intermediates
 public :: geofold_masker_setvoids, geofold_masker_readvoids, geofold_masker_seamenergy
 public :: geofold_masker_getscenergy, geofold_masker_seamsas
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE geofold_masker_readvoids(dunit)
    implicit none
    integer,intent(in) :: dunit
    character(len=1000) :: aline
    integer :: ierr, ivoid 
    !!
    !! Get voids from PDB file, if present
    rewind(dunit)
    ivoid = 0
    DO
      read(dunit,'(a)', iostat=ierr) aline !read in pdb file
      IF (ierr /=0 ) EXIT 
      IF (aline(1:7) /= "HETATM ") CYCLE
      IF (aline(18:20) /="HOH" ) CYCLE
      IF (aline(22:22) /="V" ) CYCLE
      ivoid = ivoid + 1
      ! Get Coordinates and Chain ID (_ 'Char')
      read( aline(31:54),'(3f8.3)' ), voidcoords(1:3, ivoid)   
    END DO 
    geofold_masker_nvoid = ivoid
  END SUBROUTINE geofold_masker_readvoids
  !-----------------------------------------------------------------------------
  ! getvoids
  ! count the number of voids that have all of their neighbors
   subroutine getvoids(f, voids)
     implicit none
     type(intermediate), POINTER :: f
     integer,intent(out) :: voids
     integer :: i,j,k,n
     n = 0
     ILOOP: do i=1,geofold_masker_nvoid
       if (avoid(i)%nn == 0) cycle
       do j=1,avoid(i)%nn
         if (avoid(i)%wall(j)<=0) stop "BUG in getvoids"
         if (f%iflag(avoid(i)%wall(j))=='.') cycle ILOOP
         if (any(f%barrel/=0)) then
           do k=j+1,avoid(i)%nn
             if (geofold_pivots_queryinseam(f, k,j)) cycle ILOOP
           enddo
         endif
       enddo
       n = n + 1
     enddo ILOOP
     voids = n
   end subroutine getvoids
  !-----------------------------------------------------------------------------
  SUBROUTINE geofold_masker_read(mfile)
    ! integer,intent(in) :: marg
    INTEGER :: res1, res2, ierr, dunit
    REAL :: engy, x
    CHARACTER (len=*) :: mfile
    integer :: minres, maxres
    ! WRITE(0,*) 'in read_maskerfile',"  geofold_nres=",geofold_nres
    dunit = 33
    minres = 999
    maxres = -999
    ! call getarg(marg, mfile)
    allocate(sasnrg(geofold_nres, geofold_nres),stat=ierr)
    if (ierr/=0) stop 'read_maskerfile:: error allocating sasnrg'
    allocate(sasvalues(geofold_nres,geofold_nres),stat=ierr)
    if (ierr/=0) stop 'read_maskerfile:: error allocating sasvalues'
    dunit = pickunit(dunit)
    open(dunit, file = mfile, iostat=ierr, status="old", form="formatted")
    IF (ierr /= 0 ) then
      write(0,*) "geofold_masker:: error opening masker file! ", trim(mfile)
      STOP "geofold_masker:: error opening masker file!"
    endif
    sasnrg = 0.
    sasvalues = 0.
    DO
       !read(dunit,*, iostat=ierr) res1, res2, engy, x !read in masker file
       read(dunit,'(2i6,2f9.3)',iostat=ierr) res1, res2, x, engy
       IF (ierr /=0) EXIT
       sasnrg(res1, res2) = engy !store information in data structure
       sasvalues(res1,res2) = x
       if (res1 < minres) minres = res1
       if (res2 < minres) minres = res2
       if (res1 > maxres) maxres = res1
       if (res2 > maxres) maxres = res2
    END DO
    close(dunit) ! close masker file
    write(*,'("MASKER SAS contact data found for residue range ",2i8)') minres, maxres
    write(*,'("Last data point  ",2i8,f9.4)') res1,res2,engy
  END SUBROUTINE geofold_masker_read
  !!--------------------------------------------------------------
  SUBROUTINE geofold_masker_seamenergy(seam, energy,nc,seamchar)
    ! get the amount of buried surface area exposed by opening a seam
    type(seam_type),pointer :: seam
    real,intent(out) :: energy
    integer,intent(out),optional :: nc
    character,intent(in),optional :: seamchar
    integer :: i,j, k
    energy = 0
    k = 0  !! contact count
    if (present(seamchar)) then
      DO i = 1, geofold_nres
         IF ( seam%u1flag(i) /= seamchar ) cycle
         DO j = 1, geofold_nres ! calculate the energy value
           IF ( seam%u2flag(j) /= seamchar ) cycle
           energy = energy + sasnrg(i, j)
           k = k + 1
         enddo
      enddo
    else
      DO i = 1, geofold_nres
         IF ( seam%u1flag(i) == '.' ) cycle
         DO j = 1, geofold_nres ! calculate the energy value
           IF ( seam%u2flag(j) == '.' ) cycle
           energy = energy + sasnrg(i, j)
           k = k + 1
         enddo
      enddo
    endif
    if (present(nc)) then
      nc = k
    endif
  ENDSUBROUTINE geofold_masker_seamenergy

  subroutine geofold_masker_seamsas(seam,sas,nc,seamchar)
    ! get the amount of buried surface area exposed by opening a seam
    type(seam_type),pointer :: seam
    real,intent(out) :: sas
    integer,intent(out),optional :: nc
    character,intent(in),optional :: seamchar
    integer :: i,j, k
    energy = 0
    k = 0  !! contact count
    if (present(seamchar)) then
      DO i = 1, geofold_nres
         IF ( seam%u1flag(i) /= seamchar ) cycle
         DO j = 1, geofold_nres ! calculate the energy value
           IF ( seam%u2flag(j) /= seamchar ) cycle
           !energy = energy + sasnrg(i, j)
           energy = energy + sasvalues(i,j)
           k = k + 1
         enddo
      enddo
    else
      DO i = 1, geofold_nres
         IF ( seam%u1flag(i) == '.' ) cycle
         DO j = 1, geofold_nres ! calculate the energy value
           IF ( seam%u2flag(j) == '.' ) cycle
           !energy = energy + sasnrg(i, j)
           energy = energy + sasvalues(i,j)
           k = k + 1
         enddo
      enddo
    endif
    if (present(nc)) then
      nc = k
    endif
  ENDSUBROUTINE geofold_masker_seamsas
  !!--------------------------------------------------------------
  SUBROUTINE geofold_masker_energy(f, energy, sas,sumsas)
    ! this involves looking at all pairs of residues in that intermediate
    ! Total Energy = Summation over all i and j
    ! multiplied by omega, the solvation free energy
    implicit none
    type(intermediate), POINTER :: f
    REAL,intent(out) :: energy
    real, optional :: sas,sumsas
    REAL :: sumenergy
    INTEGER :: i, j, k
    energy = 0
    sumenergy = 0
    k = 0
    if(present(sumsas)) sumsas = 0
    DO i = 1, geofold_nres
       IF ( f%iflag(i) /= '.' ) THEN ! look at intermediate f: f%iflag
          k = k + 1
          DO j = i, geofold_nres ! calculate the energy value
             IF ( f%iflag(j) /= '.' ) THEN 
               !! skip contacts that are broken by a barrel opening move
               if (geofold_pivots_queryinseam (f, i,j)) cycle
               sumenergy = sumenergy + sasnrg(i, j)
               if(present(sumsas)) sumsas = sumsas + sasvalues(i,j)
             ENDIF
          END DO
       END IF
    END DO
    if (present(sas)) sas = sumenergy
    energy = sumenergy*geofold_masker_omega
  END SUBROUTINE geofold_masker_energy
!!====================================================================================
!
  subroutine geofold_masker_setvoids()
    implicit none
    integer :: ivoid,i,j,k,ires, ios, nn
    real :: dd
    real,parameter :: dcut=7.0
    !!
    do ivoid=1,geofold_masker_nvoid
      nn = 0
      do ires=1,geofold_nres
        dd = sqrt(sum((allcoords(1:3,ires)-voidcoords(1:3,ivoid))**2))
        if (dd < dcut) then
          nn = nn + 1
        endif
      enddo
      avoid(ivoid)%nn = nn
      if (associated(avoid(ivoid)%wall)) deallocate(avoid(ivoid)%wall)
      allocate(avoid(ivoid)%wall(nn),stat=ios)
      if (ios/=0) stop 'geofold_masker.f90:: setvoids : ERROR aloocating wall.'
    enddo
    do ivoid=1,geofold_masker_nvoid
      nn = 0
      do ires=1,geofold_nres
        dd = sqrt(sum((allcoords(1:3,ires)-voidcoords(1:3,ivoid))**2))
        if (dd < dcut) then
          nn = nn + 1
          avoid(ivoid)%wall(nn) = ires
        endif
      enddo
    enddo
  end subroutine geofold_masker_setvoids

  !!------------------------------------------
  SUBROUTINE geofold_masker_intermediates(dunit)
    implicit none
    integer,intent(in) :: dunit
    REAL :: returnedenergy, scenergy, sas,sumsas
    integer ::  voids=0
    type(intermediate), POINTER :: itemp,f,u1,u2
    integer :: nres, ierr, hbonds, seams(MAXBARREL)
    nres = geofold_nres
    !!
    ! call geofold_masker_read(2)
    ! call geofold_hbonds_read(3)
    !!
    itemp => ilistroot
    WRITE(dunit, '(a)') 'HEADER ISEGMT isegmt_no recursion_depth symops Gsolv SAS '//&
                        'sc_entropy Voids Hbonds concentration seamsX8'
    DO while ( ASSOCIATED(itemp%next) ) 
       itemp => itemp%next
       call geofold_masker_energy(itemp, returnedenergy, sas=sas,sumsas=sumsas)
       call getscenergy(itemp, scenergy)
       call getvoids(itemp, voids)
       call geofold_hbonds_get(itemp, hbonds)
       call getseams(itemp, seams)
       !---- print out intermediate information and associated energy
       !---- 0.00 is a placeholder for concentration, to be added by unfoldsim.f90
       WRITE(dunit, '(a,i7,2i5,3f12.2,i8,i6,f8.5,10i4)') 'ISEGMT ', itemp%idnum, itemp%state, &
           itemp%sym, sas,sumsas, scenergy, voids, hbonds, 0.00, seams
       WRITE(dunit, '(2000(a1))') itemp%iflag(1:geofold_nres)
    END DO
    
    !!
    deallocate(sasnrg,stat=ierr)
    if (ierr/=0) stop 'Error deallocating!!'
    call geofold_hbonds_cleanup()
    !!
  CONTAINS
!!====================================================================================
! getseams
! retrieve the seam numbers for this intermediate
     subroutine getseams(f, seams)
       implicit none
       type(intermediate), POINTER :: f
       integer,dimension(:) :: seams
       integer :: i
       i = size(seams)
       seams = f%barrel(1:i)
     end subroutine getseams

  END SUBROUTINE geofold_masker_intermediates 
  
  !!-----------------------------------------------------
  subroutine getscenergy(f, scentropy, energy, seamchar)   
    implicit none
    type(intermediate), POINTER :: f
    real,intent(out) :: scentropy
    real,intent(out),optional :: energy
    character,intent(in),optional :: seamchar
    real :: sumsas, sumi, saszero
    integer :: i, j
    logical :: countit
    scentropy = 0.
    DO i = 1, geofold_nres
      sumsas = 0.
      sumi = 0.
      DO j = 1, geofold_nres ! calculate the energy value
        IF (f%iflag(i)/='.'.and.f%iflag(j)/='.') then 
          countit = .true.
          if (any(f%barrel/=0)) then
             countit = .not.geofold_pivots_queryinseam(f, i,j) 
          endif
          if (countit) sumi = sumi + sasvalues(i, j)
        END IF
        sumsas = sumsas + sasvalues(i, j)
      END DO
      saszero = maxsas(seq(i)) - maxsas(6)  !! subtract a Gly to get sc SAS
      if (saszero > 0.00000001) then
         scentropy = scentropy + lambda(seq(i))*(sumi/saszero)
      endif
    END DO
    if (present(energy)) energy=scentropy*geofold_masker_lambdaweight
  end subroutine getscenergy
  !!-----------------------------------------------------
  subroutine getscenergy_seam(aseam, scentropy, energy, seamchar)   
    !! returns side chain entropy gain with sem unfolding,
    !! based on buried surface ara exposed.
    implicit none
    type(seam_type),pointer :: aseam
    real,intent(out) :: scentropy
    real,intent(out),optional :: energy
    character,intent(in),optional :: seamchar
    character :: seamch
    real :: sumsas, sumi, saszero
    integer :: i, j
    logical :: countit
    scentropy = 0.
    seamch = "A"
    if (present(seamchar)) seamch=seamchar
    DO i = 1, geofold_nres
      !! get % exposure of residue i side chain
      sumsas = 0.
      sumi = 0.
      DO j = 1, geofold_nres 
        !! if j is on the other split of the seam, count the exposed sas.
        IF     (aseam%u1flag(i)==seamch.and.aseam%u2flag(j)==seamch) then 
          sumi = sumi + sasvalues(i, j)
        elseIF (aseam%u2flag(i)==seamch.and.aseam%u1flag(j)==seamch) then 
          sumi = sumi + sasvalues(i, j)
        END IF
        sumsas = sumsas + sasvalues(i, j)
      END DO
      saszero = maxsas(seq(i)) - maxsas(6)  !! subtract a Gly to get sc SAS
      if (saszero > 0.00000001) then
         scentropy = scentropy + lambda(seq(i))*(sumi/saszero)
      endif
    END DO
    if (present(energy)) energy=scentropy*geofold_masker_lambdaweight
  end subroutine getscenergy_seam
!!====================================================================================
END MODULE geofold_masker
