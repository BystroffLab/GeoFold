!put geofold_flory_calc_entropy in the get subroutines for each cut
!or check subroutine?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FORTRAN90 Geofold 
! by Chris Bystroff, Suzanne Matthews, Luis Garreta
! Latest version: 
!  Tue Jul  1 15:02:32 EDT 2014
!--------------
! Based on GEOFOLD by Saeed Salem, 
! Vibin Ramakrishnan, Chris Bystroff and Mohammed Zaki
! 2007
!-------------
! Based on UNFOLD, by Mohammed Zaki, Vinay Nadimpaly, 
! Deb Bardham, and Chris Bystroff
! 2005
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!====================================================================================
!! MODIFICATIONS:
!!  2-JUL-2008  Forced melting of left-over intermediates  C.B.
!!
!!  8-NOV-2008  Read and save H-bond donor acceptor list
!!              Modified to keep more strict '_'-based object naming.
!!  
!!  Wed Mar 25 17:04:12 EDT 2009 Calculate masker energy as a
!!              sorting element for selecting pivots, hinges, breaks,
!!              after passing the entropy cutoff.  C.B.
!!  Thu May 21 21:36:26 PDT 2009
!!              Command line changes. Parameter file is used to read
!!              file names and entropy cutoff values.  C.B.
!!  Mon May  3 11:06:29 EDT 2010
!!              Outputs SEQUENCE line to DAG file. C.B.
!!  Tue Jun 4 2013  Adding SEAMS. Additions by Luis Garreta.
!!  Fri Jun 5 2013: Adding SEAMS module, and read_seams. L.G.
!!  Thu Aug 1 2013 Bug fixes in seam energies. C.B.
!!  Tue Jul 1 2014 Speed/Bug fixed in seam-associated routines. C.B.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM geofold
  USE geofold_global    !! geofold_global.f90
  USE geofold_pivots    !! geofold_pivots.f90
  USE geofold_masker    !! geofold_global.f90
  USE geofold_hbonds    !! geofold_global.f90
  USE geofold_seams     !! geofold_seams.f90
  USE geofold_flory     !! geofold_flory.f90
  ! mpi added by SAN
  use mpi
  implicit none
  ! omp added by SAN
  include 'omp_lib.h'
  

  INTEGER :: n_procs, mpierr, proc_id
  
  INTEGER :: nres, ierr !number of residues
  CHARACTER, dimension(:), allocatable:: chainid !gets passed in 
  CHARACTER(len=200) :: aline
  ! CHARACTER :: chainID
  INTEGER :: ires, i, j, jarg, ios, ivoid, dunit, flory
  CHARACTER(len=1000) :: filename=" ",dagfile=" ",parfile=" ",cmfile=" ",hbfile=" ",seamfile=" "
  type(intermediate), POINTER :: gptr          ! points to the global intermediate
  type(intermediate), TARGET :: Native
  type(contact),dimension(:),allocatable :: contacts
  real :: w,T
  logical :: all_seams = .true.
  ! variables added by san
  CHARACTER(len=1000) :: timerfile=" ", msg = " "
  DOUBLE PRECISION :: timer, start_timer, end_timer
  ! character(len=3) :: aa
  !------------------------------ COMMAND LINE ----------------
  ! mpi added by SAN
  !  Initialize MPI.
  call MPI_INIT ( mpierr )
  !  Get the number of processors.
  call MPI_COMM_SIZE ( MPI_COMM_WORLD, n_procs, mpierr )
  !  Get the rank of this processor.
  call MPI_COMM_RANK ( MPI_COMM_WORLD, proc_id, mpierr )
  
  ires = 0
  jarg = iargc() ! get the number of command line arguments
  IF ( jarg < 3 ) THEN
  	IF ( proc_id == 0 ) THEN
        write(*,*) "Usage: xgeofold <input PDB> <output DAG> <parameter file>"
        ! write(*,*) "Usage: xgeofold <inputpdb> <mas file> <hbond file> <outputdag> [<minhinge> <minpivot> <minbreak>]" ! write usage line
        write(*,*) "<inputpdb> Pre-processed input coordinates in PDB format (xGetChain), including Void atoms(Voidmask)."
        write(*,*) "<outputdag> Output directed acyclic graph, for input to xUnfoldSim"
        write(*,*) "<parameter file> Keyworded file containing the following keywords (others keywords ignored)"
        write(*,*) "  CONTACTS <filename> -- Buried surface areas from Masker program ContactMask.f90"
        write(*,*) "  HBONDS <Hbond file> H-bond file from the utility pdb2hb.f90"
        write(*,*) "  HINGECUT <hcutoff> Minimum possible hinge rotation as a fraction of", 2*MAXHINGEANG, &
                    "   degrees. default=",hcutoff," [optional]"
        write(*,*) "  PIVOTCUT <pcutoff> Minimum number of possible pivot axes as a fraction of", NVB,&
                    "   default=",pcutoff," [optional]"
        write(*,*) "  BREAKCUT <bcutoff> Minimum number of possible translation vectors as a fraction of", NVB,&
                    " default=",bcutoff," [optional]"
        write(*,*) "  SEAMCUT <scutoff> Minimum number of possible translation vectors as a fraction of", NVB,&
                    " default=",scutoff," [optional]"
        write(*,*) "  MAXSPLIT <n> n={2,4,6,8} is the maximum number of children for each parent."
        write(*,*) "  default=",maxsplit
        write(*,*) "  MINSEG <n> n={2..20} is the Minimum size of a terminal fragment to unfold."
        write(*,*) "  default=",pivottail
        stop 'geofold.f90 v.  Wed Sep 30 12:00:57 EDT 2009'
    END IF
  END IF
  call getarg(1, filename) ! get the name of input pdb from command line
  call getarg(2, dagfile) ! get the name of output file
  call getarg(3, parfile) ! get the name of parameters file
  !!=============== READ PARAMETERS FILE ================
  dunit = pickunit(20)
  ! MPI I/O added by San
 
  IF ( proc_id == 0 ) THEN
	  !call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
	  !call MPI_FILE_OPEN(MPI_COMM_WORLD, parfile, MPI_MODE_RDONLY, MPI_INFO_NULL, dunit, mpierr)
	  open(dunit, file=parfile, status='old', form='formatted', iostat=ios)
	  IF (ios/=0) STOP 'geofold.f90:: parameters file not found.'
	!  write(0,*) 'Reading parameters'
	  call geofold_readparameter(dunit,"HINGECUT",hcutoff,low=0.0,high=1.0,default=hcutoff)
	  call geofold_readparameter(dunit,"PIVOTCUT",pcutoff,low=0.0,high=1.0,default=pcutoff)
	  call geofold_readparameter(dunit,"BREAKCUT",bcutoff,low=0.0,high=1.0,default=bcutoff)
	  call geofold_readparameter(dunit,"CONTACTS",cmfile,required=1)
	  call geofold_readparameter(dunit,"HBONDS",hbfile,required=1)
	  call geofold_readparameter(dunit,"SEAMS",seamfile,required=1)
	  call geofold_readparameter(dunit,"MINSEG",pivottail,default=pivottail)
	  call geofold_readparameter(dunit,"VERBOSE",verbose,default=.false.)
	  call geofold_readparameter(dunit,"MAXSPLIT",maxsplit,default=maxsplit)
	  call geofold_readparameter(dunit,"HBONDENERGY",geofold_hbonds_eperbond,low=0.00,default=100.)
	  call geofold_readparameter(dunit,"OMEGA",geofold_masker_omega,low=0.00,default=1.)
	  call geofold_readparameter(dunit,"SIDECHAINENTROPY",geofold_masker_lambdaweight,low=0.00,default=1.)
	  call geofold_readparameter(dunit,"VOIDENTROPY",geofold_masker_spervoid,low=0.0,default=0.0)
	  call geofold_readparameter(dunit,"FLORY",flory,default=0)
	  call geofold_readparameter(dunit,"FLORYW",w,default=1.)
	  call geofold_readparameter(dunit,"TEMPERATURE",T,default=300.)
	  call geofold_readparameter(dunit,"SEAMCUT",scutoff,low=0.0,high=1.0,default=scutoff)
	  call geofold_readparameter(dunit,"ALLSEAMS",all_seams,default=.true.)
	  !call MPI_FILE_CLOSE(dunit, mpierr)
	  close(dunit)
	!  write(0,*) 'Read parameters'
	  !------------------------------ READ INPUT PDB and other FILEs ---------
	  dunit = pickunit(10)
	  ! MPI I/O added by San
	  !call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
	  !call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, dunit, mpierr)
	  open(dunit, file=filename, form="formatted", status="old", iostat=ierr)
	  IF (ierr > 0 ) STOP "geofold:: Error! File not found!"
	  !write(0,*) 'geofold_readpdb'
	  call geofold_readpdb(dunit)
	  !write(0,*) 'geofold_masker_readvoids'
	  call geofold_masker_readvoids(dunit)   !! read from same PDB file
	  !call MPI_FILE_CLOSE(dunit, mpierr)
	  close(dunit)
	  !write(0,*) 'geofold_masker_read'
	  call geofold_masker_read(cmfile)
	  write(0,'("hbfile: ",a)') hbfile

	!  write(0,*) "TESTING............."
	  !write(0,*) 'geofold_hbonds_read'
	  call geofold_hbonds_read(hbfile)
	  !write(0,*) 'geofold_seams_read'
	  call geofold_seams_read(seamfile)
	  !------------------------------ INITIALIZE   ----------------
	  nres = geofold_nres 
	  Native%iflag = masterchains  
	  Native%idnum = 1
	  Native%state = 1 
	  Native%axis = 0  
	  Native%barrel = 0
	  NULLIFY(Native%next)
	  gptr => Native
	  nullify(ilistroot)
	  !write(0,*) 'initpivot'
  
    CALL geofold_initpivot(allcoords,nres,chainid=masterchains)
    !write(0,*) 'geofold_masker_setvoids()'
    CALL geofold_masker_setvoids()
    !write(0,*) 'geofold_flory_all_contacts'
    call geofold_flory_all_contacts(hbfile,cmfile,contacts,geofold_nres)
    !------------------------------ WORK   ----------------
    !write(0,*) 'getcutpoints'
    !timer_write added by SAN
    timerfile = "/bach1/home/sanw/GeoFold/tmp/timer.txt"
    msg = "time used by RECURSIVE SUBROUTINE getcutpoints"
    start_timer = omp_get_wtime()
    call getcutpoints(gptr,contacts,flory,T)
    end_timer = omp_get_wtime()
    timer = end_timer-start_timer
  
    call timer_write(timerfile, msg, timer)
    !------------------------------ FINISH UP   ----------------
    !write(0,*) 'dag_write'
    call dag_write(dagfile,ounit=dunit) 
    !  write(0,*) 'geofold_seams_write'
    call geofold_seams_write(ounit=dunit)
    close(dunit)
    !  write(0,*) 'cleanuplists'
    call cleanuplists()
    !  close(45)
  END IF
  call MPI_FINALIZE(mpierr)
CONTAINS

!!====================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine getcutpoints
!takes a folded portion of protein
!recurses until unfolded
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE getcutpoints( f ,contacts,flory,T)
  !! ------- cut point algorithm ---------
  !! 
  !! If the intermediate is multichain, break if possible.
  !! If not, then if the chain(s) is long enough, pivot if possible.
  !! If not, then if the chain is long enough, hinge if possible.
  !! If not, then if barrels are present, then seam if possible.
  !! If not, if the chain is greater than length 1, then melt
  !! If not, then it is length 1, return.
  !!
  implicit none
  type(intermediate), POINTER :: f
  type(intermediate), target :: child1, child2
  type(intermediate), pointer :: u1, u2
  type(contact),dimension(:),allocatable,intent(in) :: contacts
  INTEGER :: lastpos,lastbreak,lasthinge
  INTEGER :: interval
  INTEGER :: flag,ios=0
  integer :: flory
  INTEGER :: cuttype, nchain, nres, mres, nseam
  integer :: nbreak,ibreak,npivot,ipivot,nhinge,ihinge, iseam, nmove
  REAL :: entropy
  real,intent(in) :: T
  ! logical,optional,intent(out) :: done
  ! logical :: u1done,u2done
  character,dimension(MAXCHAIN) :: uniqchains
  type(intermediate),dimension(:),allocatable :: allu1
  real,dimension(:),allocatable :: allentropy
  type (seammove_type), pointer :: seammove (:)
  !-----------
  nullify(u1, u2)
  nullify(seammove)
  if (allocated(allentropy)) deallocate(allentropy)
  allocate(allentropy(maxsplit+1),stat=ios)
  if (ios/=0) stop 'geofold:: getcutpoints: error allocating allentropy'
  if (allocated(allu1)) deallocate(allu1)
  allocate(allu1(maxsplit+1),stat=ios)
  if (ios/=0) stop 'geofold:: getcutpoints: error allocating allu1'
  nres = geofold_nres
  if (.not. associated(f)) then
    write(0,'("Warning: called getcutpoints with an unassociated pointer!!")')
    stop 'Error: called getcutpoints with an unassociated pointer'
  endif
  interval = oldintermediate(f)
  !if (verbose) then
  !   if (interval==0) then
  !     WRITE (0,*) 'in getcutpoints:: NO oldintermediate found for ',f%idnum
  !   else
  !     WRITE (0,*) 'in getcutpoints:: oldintermediate number = ', interval
  !   endif
  !endif
  f%idnum = interval
  IF (interval /= 0 )  RETURN
  IF ( count(f%iflag(1:geofold_nres) /= '.') == 0 ) RETURN ! empty
  !!---- CYCLE FROM HERE as long as no cutpoints can be found
  ! CUTOFFLOOP: DO
  !!---- if intermediate is too small, melt it.
!  write(0,*) 'getallchains'
  call getallchains(f%iflag,uniqchains(1:nres),nchain,nres)  !! returns chain chars in uniqchains
  mres = count(f%iflag(1:geofold_nres) /= '.')
  if (nchain==1.and.mres<=(2*pivottail+3)) then
!    write(0,*) 'getmelting'
    call getmelting(f)
    return
  endif
!  write(0,*) 'saveintermediate'
  CALL saveintermediate(f)
  IF ( count(f%iflag(1:geofold_nres) /= '.') == 1 ) RETURN     ! leaf node
  !! ---- two new children are created for each recursion depth ----
  u1 => child1
  u2 => child2
  u1%idnum = 0
  u2%idnum = 0
  u1%state = f%state + 1    ! recursion depth
  u2%state = f%state + 1    ; if (f%state > 5000) stop 'geofold:: getcutpoints: RECURSION DEPTH LIMIT EXCEEDED!!!'
  u1%iflag = '.'
  u2%iflag = '.'
  u1%axis = 0
  u2%axis = 0
  u1%barrel=f%barrel
  u2%barrel=f%barrel
  NULLIFY(u1%next)
  NULLIFY(u2%next)
  !! Split the pathway 
  geofold_split = nint(real(maxsplit)/(f%state / FULSPLITDEPTH ) )
  if (geofold_split==0) geofold_split = 1
  if (geofold_split>maxsplit) geofold_split = maxsplit
  !!------------- BREAKS ----------------
!  write(0,*) 'getbreaks'
  call getbreaks(f,allu1,allentropy,nbreak,contacts,flory)
  if (verbose.and.nbreak > 0) write(*,*) '============>>> found ',nbreak,' BREAKs'

  breakloop: DO ibreak=1,nbreak
     entropy = allentropy(ibreak)
     u1%iflag = allu1(ibreak)%iflag
     f%axis = allu1(ibreak)%axis
!     write(0,*) 'calling getcutpoints from break u1'

     CALL getcutpoints(u1,contacts,flory,T)   !!  recusrively find cutpoints
     
     u2%iflag = f%iflag
     where (u1%iflag/='.') u2%iflag = '.'
!     write(0,*) 'calling getcutpoints from break u2'

     CALL getcutpoints(u2,contacts,flory,T)

     cuttype = breakflag
!     write(0,*) 'calling savetstate from break'

     CALL savetstate(f,u1,u2=u2,t=cuttype, ent=entropy)
  enddo breakloop

  if (nbreak > 0) return
  !!------------- PIVOTS ----------------
!  write(0,*) 'getpivots'
  call getpivots(f,allu1,allentropy,npivot,contacts,flory)
!  write(0,*) 'got pivots'
  if (verbose.and.npivot > 0) write(*,*) '============>>> found ',npivot,' PIVOTS'
  pivotloop: DO ipivot=1,npivot 
!     write(0,*) 'in pivotloop'
     entropy = allentropy(ipivot)
     u1%iflag = allu1(ipivot)%iflag
     f%axis = allu1(ipivot)%axis
!     write(0,*) 'pu1 getcutpoints'
     CALL getcutpoints(u1,contacts,flory,T)   !!  recusrively find cutpoints
     u2%iflag = f%iflag
     where (u1%iflag/='.') u2%iflag = '.'
!     write(0,*) 'pu2 getcutpoints'
     CALL getcutpoints(u2,contacts,flory,T)
     cuttype = pivotflag
!     write(0,*) 'p savetstate'
     CALL savetstate(f,u1,u2=u2,t=cuttype, ent=entropy)
  enddo pivotloop
!  write(0,*) 'out of pivot loop'
!  write(0,'(i3," pivots found.")') npivot
  if (npivot > 0) return
  !!------------- SEAMS  ----------------
  !!-1: left, 1: right  --> allside
  !This is apparently where the segfault is happening.... What do?
  write(0,*) 'deallocating seammove'
  write(0,*) associated(seammove)
  if (associated(seammove)) deallocate(seammove)
  write(0,*) 'seammove deallocated'
  nseam = maxsplit
  !!!All seams
  if(all_seams) then
    nseam = 0
    do i = 1, size(barrels_array)
      nseam = nseam + barrels_array(i)%nSeams
    enddo
  endif
  !!!
  allocate(seammove(nseam),stat=ios); if (ios/=0) stop 'geofold:: getcutpoints: error allocating seammove.'
  seammove(:)%barrel=0;seammove(:)%seam=0;seammove(:)%energy=0;seammove(:)%side=0;
!  write(0,*) 'getseams'
  call getseams(f, seammove, nseam,contacts,flory,w,T) 
  if (verbose.and.nseam > 0) then
     write(*,*) '============>>> found ',nseam,' SEAMS'
     do iseam=1,nseam
       write(*,*) seammove(iseam)%barrel, seammove(iseam)%seam, seammove(iseam)%energy
     enddo
  elseif (verbose) then
     write(*,*) '============>>> found no SEAMS'
  endif
  seamloop: DO iseam=1,nseam
     if (seammove(iseam)%barrel==0) cycle
     if (seammove(iseam)%barrel<0) then
       write(0,'("ERROR: seammove(",i3,")%barrel=",i3)') iseam,seammove(iseam)%barrel
       stop "geofold.f90:: error in seamloop."
     endif
     entropy = seammove(iseam)%energy  !! contact energy
     u1%iflag = f%iflag
     u1%barrel(seammove(iseam)%barrel) = seammove(iseam)%seam !! set flag
     write(*,*) "seamloop:: u1%barrel(seammove(",iseam,")%barrel)=",u1%barrel(seammove(iseam)%barrel)
!     write(0,*) 'su1 getcutpoints'
     CALL getcutpoints(u1,contacts,flory,T)   !!  recursively find cutpoints
     cuttype = seamflag
     !Flory entropy
     !entropy = entropy + geofold_flory_calc_entropy(flory,f,u1,c_list=contacts,w=w,T=T)
!     write(0,*) 's savetstate'
     CALL savetstate(f,u1,t=cuttype, ent=entropy, iseam=u1%barrel(seammove(iseam)%barrel)) !! add to list of intermediates
  enddo seamloop
  if (associated(seammove)) deallocate(seammove)
  if (nseam > 0) return
  !!------------- HINGES  ----------------
!  write(0,*) 'gethinges'
  call gethinges(f,allu1,allentropy,nhinge,contacts,flory)
  if (verbose.and.nhinge > 0) then
    write(*,*) '============>>> found ',nhinge,' HINGES'
  elseif (verbose) then
    write(*,*) '============>>> found no HINGES'
  endif
  hingeloop: DO ihinge=1,nhinge
     entropy = allentropy(ihinge)
     u1%iflag = allu1(ihinge)%iflag
     f%axis = allu1(ihinge)%axis
!     write(0,*) 'hu1 getcutpoints'
     CALL getcutpoints(u1,contacts,flory,T)   !!  recusrively find cutpoints
     u2%iflag = '.'
     where (f%iflag/='.'.and.u1%iflag=='.') u2%iflag = f%iflag
!     write(0,*) 'hu2 getcutpoints'
     CALL getcutpoints(u2,contacts,flory,T)
     cuttype = hingeflag
!     write(0,*) 'h savetstate'
     CALL savetstate(f,u1,u2=u2,t=cuttype, ent=entropy)
  enddo hingeloop
  if (nhinge > 0) return
  !!------------- MELT    ----------------
  ! WRITE (*,*) 'ERROR: No cutpoints found for ',f%iflag(1:geofold_nres)
  ! stop 'geofold:: BUG or ERROR. No cutpoints found. '
  WRITE(*,*) 'WARNING: No cutpoints found for ',f%iflag(1:geofold_nres)
  if (any(f%barrel(:)/=0)) write(*,*) "NOTE: It is an open barrel."
  WRITE(*,*) 'MELTING it.'
!  write(0,*) 'getmelting'
  call getmelting(f,force=.true.)
  return
END SUBROUTINE getcutpoints

!!====================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function getbreaks
!Attempts to find chains tht can separate
!by simple translation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE getbreaks(f,allu1,allentropy,nbreak,contacts,flory)
  implicit none
  integer,intent(out) :: nbreak
  type(intermediate), POINTER :: f
  type(contact),dimension(:),allocatable,intent(in) :: contacts
  type(intermediate),dimension(:),intent(out) :: allu1   !! changed to assumed-shape 5/3/10
  real,dimension(:),intent(out) :: allentropy   !! changed to assumed-shape 5/3/10
  type(intermediate),pointer :: u1, u2
  character,dimension(MAXCHAIN) :: uniqchains
  real,dimension(:),allocatable :: allenergy
  real :: entropy, energy, fenergy, u1energy, u2energy
  integer :: ires, nchain, nres, bpoint, ibreak, mbreak, nsplit, bvec, i, j, ios
  integer :: flory
  if (allocated(allenergy)) deallocate(allenergy)
  allocate(allenergy(geofold_split+1),stat=ios); if (ios/=0) stop 'geofold:: getbreaks: BUG allocating allenergy'
  call geofold_masker_energy(f,fenergy)
  allocate(u1,u2,stat=ios); if (ios/=0) stop 'geofold:: error allocating u1,u2 in getpivots'
  call zerointermediate(u1)
  call zerointermediate(u2)
  call zerointermediate(allu1)
  entropy = 0.
  nbreak = 0
  bvec = 0
  nres = geofold_nres
  call getallchains(f%iflag,uniqchains,nchain,nres)  !! returns chain chars in uniqchains
  if (nchain <= 1) then
    if(associated(u1)) deallocate(u1)
    if(associated(u2)) deallocate(u2)
    return
  endif
  mbreak = 2**(nchain-2)
  bpoint = 0
  nullify(u1%next, u2%next)
  nsplit = geofold_split
  allenergy = -99999999.
  allentropy = 0.
  do while (bpoint<mbreak)
    call geofold_getnextbreak(f=f,calpha=allcoords,chainid=f%iflag,u1=u1%iflag,u2=u2%iflag, &
                            nres=nres,breakpoint=bpoint,entropy=entropy,bvec=bvec)
    if(entropy /=-1) &
      entropy = entropy + geofold_flory_calc_entropy(flory,f,u1,u2,contacts,w,T)
    if (entropy>bcutoff) then
      call geofold_masker_energy(u1,u1energy)
      call geofold_masker_energy(u2,u2energy)
      energy = u1energy + u2energy - fenergy
      if (nbreak<nsplit) nbreak = nbreak + 1
      ibreak = nbreak + 1
      !!do while (entropy>allentropy(ibreak-1))
      do while (energy>allenergy(ibreak-1))
        allenergy(ibreak) = allenergy(ibreak-1)
        allentropy(ibreak) = allentropy(ibreak-1)
        allu1(ibreak)%iflag = allu1(ibreak-1)%iflag
        allu1(ibreak)%axis = allu1(ibreak-1)%axis
        ibreak = ibreak - 1
        if (ibreak==1) exit
      enddo
      allenergy(ibreak) = energy
      allentropy(ibreak) = entropy
      allu1(ibreak)%iflag = u1%iflag
      allu1(ibreak)%axis = bvec
    endif
  enddo
  if(associated(u1)) deallocate(u1)
  if(associated(u2)) deallocate(u2)
  if(allocated(allenergy)) deallocate(allenergy)
END SUBROUTINE getbreaks

!!====================================================================================
! get pivots finds positions within the chain that permit
! free rotation of the N-term versus the C-term
! Separate chains are distributed to N or C-term all possible ways.
!
SUBROUTINE getpivots(f,allu1,allentropy,npivot,contacts,flory)
  implicit none
  integer,intent(out) :: npivot
  type(intermediate), POINTER :: f
  type(contact),dimension(:),allocatable,intent(in) :: contacts
  type(intermediate),dimension(:),intent(out) :: allu1   !! assumed shape 5/3/10
  real,dimension(:),intent(out) :: allentropy   !! assumed shape 5/3/10
  type(intermediate),pointer :: u1, u2
  ! character,dimension(MAXCHAIN) :: uniqchains
  real :: entropy, energy, fenergy, u1energy, u2energy
  real,dimension(:),allocatable :: allenergy
  integer :: ires, nchain, nres, bpoint, ipivot,nsplit, bvec, ios,flory
  !!
  allocate(allenergy(geofold_split+1),stat=ios); if (ios/=0) stop 'geofold:: getpivots: BUG allocating allenergy'
  call geofold_masker_energy(f, fenergy)
  !diagnostic
  !write(0,*) fenergy
  !write(0,*) "getpivots!!!"
  entropy = 0.
  npivot = 0
  nres = geofold_nres
  bpoint = 0
  allocate(u1,u2,stat=ios); if (ios/=0) stop 'geofold:: error allocating u1,u2 in getpivots'
  call zerointermediate(allu1)
  call zerointermediate(u1)
  call zerointermediate(u2)
  nsplit = geofold_split
  allentropy = 0.
  allenergy = -99999999.
  bvec = 0
  write(0,*) 'getpivots stuff initialized'
  do while (bpoint<nres)
    write(0,*) bpoint, " ", nres
    write(0,*) 'geofold_getnextpivot'
    call geofold_getnextpivot(f=f,calpha=allcoords,chainid=f%iflag,u1=u1%iflag,u2=u2%iflag, &
                             nres=nres,pivotpoint=bpoint,entropy=entropy,bvec=bvec)
    if(entropy /= -1)&
      entropy = entropy + geofold_flory_calc_entropy(flory,f,u1,u2,contacts,w,T)
    if (entropy<pcutoff) then
      write(0,*) 'entropy too low...'
      if(associated(u1)) deallocate(u1)
      if(associated(u2)) deallocate(u2)
      write(0,*) 'returning...'
      return
    endif
    write(0,*) 'geofold_masker_energy'
    call geofold_masker_energy(u1,u1energy)
    call geofold_masker_energy(u2,u2energy)
    energy = u1energy + u2energy - fenergy
    if (entropy>pcutoff) then
      write(0,*) 'energy > pcutoff'
      if (npivot<nsplit) npivot = npivot + 1
      ipivot = npivot + 1
      !!do while (entropy>allentropy(ipivot-1))
      do while (energy>allenergy(ipivot-1))
        allentropy(ipivot) = allentropy(ipivot-1)
        allenergy(ipivot) = allenergy(ipivot-1)
        allu1(ipivot)%iflag = allu1(ipivot-1)%iflag
        allu1(ipivot)%axis = allu1(ipivot-1)%axis
        ipivot = ipivot - 1
        if (ipivot==1) exit
      enddo
      allenergy(ipivot) = energy
      allentropy(ipivot) = entropy
      allu1(ipivot)%iflag = u1%iflag 
      allu1(ipivot)%axis = bvec
    endif
  enddo
  if(associated(u1)) deallocate(u1)
  if(associated(u2)) deallocate(u2)
  if(allocated(allenergy)) deallocate(allenergy)
  write(0,*) 'end of getpivots'
END SUBROUTINE getpivots
!!====================================================================================
! gethinges attempts to find pairs of positions in the
! same chain that permit the rotation around the
! hinge axis of the segment between the hine points.
!
SUBROUTINE gethinges(f,allu1,allentropy,nhinge,contacts,flory)
  implicit none
  type(contact),dimension(:),allocatable,intent(in) :: contacts
  integer,intent(out) :: nhinge
  type(intermediate), POINTER :: f
  type(intermediate),dimension(:),intent(out) :: allu1   !! assumed shape 5/3/10
  real,dimension(:),intent(out) :: allentropy   !! assumed shape 5/3/10
  real :: energy, fenergy, u1energy, u2energy
  real,dimension(:),allocatable :: allenergy
  type(intermediate),pointer :: u1, u2
  character,dimension(MAXCHAIN) :: uniqchains
  real :: entropy, mat(3)
  integer :: ires, nchain, nres, bpoint, ihinge, nsplit,flory
  integer :: lasti, lastj
  if (allocated(allenergy)) deallocate(allenergy)
  allocate(allenergy(geofold_split+1),stat=ios); if (ios/=0) stop 'geofold:: gethinges: BUG allocating allenergy'
  call geofold_masker_energy(f, fenergy)
  entropy = 0.
  nhinge = 0
  nres = geofold_nres
  bpoint = 0
  allocate(u1,u2,stat=ios); if (ios/=0) stop 'geofold:: error allocating u1,u2 in gethinges'
  call zerointermediate(allu1)
  call zerointermediate(u1)
  call zerointermediate(u2)
  nsplit = geofold_split
  lasti = 1
  lastj = 0
  allentropy = 0.
  allenergy = -99999999.
  allu1(:)%axis = 0
  do while (lasti<nres)
   call geofold_getnexthinge(f=f,calpha=allcoords,chainid=f%iflag,u1=u1%iflag,u2=u2%iflag, &
                             nres=nres,hinge1=lasti,hinge2=lastj,entropy=entropy)
    if(entropy /= -1)&
      entropy = entropy + geofold_flory_calc_entropy(flory,f,u1,u2,contacts,w,T)
    if (entropy<hcutoff) then
      if(associated(u1)) deallocate(u1)
      if(associated(u2)) deallocate(u2)
      return
    endif
    if (entropy>hcutoff) then
      call geofold_masker_energy(u1,u1energy)
      call geofold_masker_energy(u2,u2energy)
      energy = u1energy + u2energy - fenergy
      if (nhinge<nsplit) nhinge = nhinge + 1
      ihinge = nhinge + 1
      !!do while (entropy>allentropy(ihinge-1))
      do while (energy>allenergy(ihinge-1))
        allenergy(ihinge) = allenergy(ihinge-1)
        allentropy(ihinge) = allentropy(ihinge-1)
        allu1(ihinge)%iflag = allu1(ihinge-1)%iflag
        allu1(ihinge)%axis = allu1(ihinge-1)%axis
        ihinge = ihinge - 1
        if (ihinge==1) exit
      enddo
      allenergy(ihinge) = energy
      allentropy(ihinge) = entropy
      !! NOTE: we use u2 from geofold_getnexthinge because it has new flags, 
      !! while u1 retains flags from f. New flags are needed in calling
      !! routine, getcutpoints.
      allu1(ihinge)%iflag = u2%iflag
      !! watch out. tricky business: store residue numbers for hinge in "axis" element.
      allu1(ihinge)%axis = 10000*lasti + lastj
    endif
  enddo
  if(associated(u1)) deallocate(u1)
  if(associated(u2)) deallocate(u2)
  if(allocated(allenergy)) deallocate(allenergy)

END SUBROUTINE gethinges
!---------------------------------------------------------------------------------
! getseams finds seams in beta barrels
!---------------------------------------------------------------------------------

!!  Need to modify this to play nice with Flory
subroutine getseams (f, seammove, nMove,contacts,flory,w,T)
  implicit none
      type(contact),dimension(:),allocatable,intent(in) :: contacts
      type (intermediate), POINTER         :: f,u1
      type (seammove_type), pointer        :: seammove (:)
      integer,intent (inout)               :: nMove
      integer                              :: nBarrels, iBarrel, nseams, iSeam, i, side, nseam
      type (seammove_type)                 :: tmpMove
      type (seam_type),pointer             :: aseam
      real                                 :: energy
      integer, intent(in)                  :: flory
      real, intent(in)                     :: w,T
    CHARACTER, dimension(1600) :: flags     ! flags that define the intermediate
!testing something out here...  This will make it more consistent with everything else
  seammove(:)%energy = HUGE(0.0)
  !seammove(:)%energy = -HUGE(0.0)
  nseam = nMove
  nMove = 0
    side = 0

  nBarrels = size (barrels_array)
    if (nBarrels==0) return
  do iBarrel=1, nBarrels
    if (f%barrel(iBarrel) /= 0) cycle  !! open already. No more seam moves on this barrel.
    do iseam=1, barrels_array(iBarrel)%nseams
            aseam => barrels_array(iBarrel)%seams(iseam)
      energy = getEnergySeam(aseam,T=T)
      tmpMove%barrel = iBarrel
      tmpMove%seam   = iSeam
      tmpMove%energy = energy

      if(flory /= 0) then
        allocate(u1)
        u1%idnum = 0
        u1%iflag = f%iflag
        nullify(u1%next)
        u1%state = seamflag
        u1%sym = 0
        u1%axis = iSeam
        u1%barrel = f%barrel
        u1%barrel(iBarrel)=iSeam      
        tmpMove%energy = tmpMove%energy - &
        T*geofold_flory_calc_entropy(flory,f,u1,c_list=contacts,w=w,T=T)
        if(associated(u1)) deallocate(u1)
        energy = tmpMove%energy
      endif
      if (energy < maxval(seammove(1:nseam)%energy,dim=1)) then
        i = maxloc(seammove(1:nseam)%energy,dim=1)
!      if(energy > minval(seammove(1:nseam)%energy,dim=1)) then
!        i = minloc(seammove(1:nseam)%energy,dim=1)
        seammove(i) = tmpMove
        nMove = nMove + 1
      endif
    enddo
  enddo
    !if (nMove > nseam) nMove = nseam
endsubroutine getseams
!!====================================================================================
! getmelting divides a short segment into
! smaller pieces until it consists of single residues.
!
recursive subroutine getmelting(f,force)
  type(intermediate), POINTER :: f
  logical,optional :: force
  type(intermediate), target :: child1, child2
  type(intermediate), pointer :: u1, u2
  type(intermediate), POINTER :: ihead
  integer :: found, i, id, n, z1,z2, check
  logical,parameter :: veryverbose=.true.
  !!
  if (.not.associated(ilistroot)) stop 'geofold.f90:: getmelting BUG 1'
  IF ( count(f%iflag(1:geofold_nres) /= '.') == 0 ) RETURN ! empty
  id =  oldintermediate(f)
  if (id/=0) then
    f%idnum = id
    if (verbose.and.veryverbose) write(0,*) "getmelting:: Found an old intermediate ",id
    if (.not.present(force)) return
  endif
  if (id==0) call saveintermediate(f)
  IF ( count(f%iflag(1:geofold_nres) /= '.') == 1 ) RETURN !checks to see if it is a single aa
  ihead => ilistroot
  found = 0
  n = 0
  outloop: do while ( associated(ihead%next) )
     ihead => ihead%next
     n = n + 1
     if (f%idnum==ihead%idnum) cycle
     if (all((ihead%iflag(1:geofold_nres)=='.').or.(f%iflag(1:geofold_nres)/='.'))) then
       !! check to see that the first or last flag is set
       inloop1: do i=1,geofold_nres
         if (f%iflag(i)/='.') then
           if (ihead%iflag(i)/='.') then
             found = 1
             exit outloop
           else
             exit inloop1
           endif
         endif
       enddo inloop1
       inloop2: do i=geofold_nres,1,-1
         if (f%iflag(i)/='.') then
           if (ihead%iflag(i)/='.') then
             found = 1
             exit outloop
           else
             exit inloop2
           endif
         endif
       enddo inloop2
     endif
  enddo outloop
  u1 => child1
  u2 => child2
  u1%iflag = '.'
  u2%iflag = '.'
  u1%axis = 0
  u2%axis = 0
  u1%state = f%state + 1    ! recursion depth
  u2%state = f%state + 1
  u1%barrel = f%barrel
  u2%barrel = f%barrel
  if (found==1) then
    do i=1,geofold_nres
      u1%iflag(i) = ihead%iflag(i)
      if (f%iflag(i)/='.' .and. ihead%iflag(i)=='.') then
        u2%iflag(i) = f%iflag(i)
      endif
    enddo
    ! call getmelting(u1)  !! not necessary since u1 = ihead
  else  
    !  write(0,*) "Melting by one character  n= ", n
    !! no existing subsets in the intermediate list. peel one residue.
    u2%iflag = f%iflag
    do i=1,geofold_nres
       if (f%iflag(i)/='.') then
         u1%iflag(i) = f%iflag(i)
         u2%iflag(i) = '.'
         exit
       endif
    enddo
  endif
  !set idnums for u1 and u2
  u1%idnum = oldintermediate(u1)
  u2%idnum = oldintermediate(u2)
  z1 = count(u1%iflag(1:geofold_nres)/='.')
  z2 = count(u2%iflag(1:geofold_nres)/='.')
  if (z1==0) call getmelting(u1)
  if (z2==0) call getmelting(u2)
  if (z1/=0.and.z2/=0) then
    !if u1 or u2 are intermediate #0, saveintermediate
    if(u1%idnum == 0) call saveintermediate(u1)
    if(u2%idnum == 0) call saveintermediate(u2)
    call savetstate(f,u1,u2=u2,t=meltflag, ent=0.1)
  endif
  !!  
end subroutine getmelting

!!====================================================================================

!!====================================================================================
! Changed u2 to optional. For seams 
!!====================================================================================
! savetstate adds  one elemental subsystem, composed
! of segment f and segments u1 and u2, to a lnked list of
! saved transitions states.
!
SUBROUTINE savetstate(f, u1, u2, t, ent, iseam)
  implicit none
  type(intermediate), POINTER :: f, u1
  type(intermediate), POINTER, optional ::  u2
  INTEGER, intent(in) :: t
  REAL, intent(in) :: ent
  integer, optional, intent(in) :: iseam
  INTEGER, save :: tstatecounter=1
  TYPE (tstate), POINTER  :: ttail   ! declare tail of the tstate list
  character(len=5),parameter :: cuttypes="bphsm"
  if (.not. associated(f)) STOP 'save tstate: error! f not associated!'
  if (.not. associated(u1)) STOP 'save tstate: error! u1 not associated!'
  if (present(u2)) then
    if (.not. associated(u2)) STOP 'save tstate: error! u2 not associated!'
  endif
  if (f%idnum==0) then
    write(0,*) 'savetstate: WARNING! attempting to save f%idnum=0'
    return
  endif
  if (u1%idnum==0) then
    write(0,*) 'savetstate: WARNING! attempting to save u1%idnum=0'
    return
  endif
  if (present(u2)) then
    if (u2%idnum==0) then
      write(0,*) 'savetstate: WARNING! attempting to save u2%idnum=0, cuttype=',t
      return
    endif
  endif
  if (u1%idnum==f%idnum) then
    write(0,*) 'savetstate: WARNING!  attempting to save u1%idnum=f%idnum cuttype=',t
    !debug return
    !reassign u1%idnum
    call saveintermediate(u1)
    !if it fails, return
    if (u1%idnum==f%idnum) return
  elseif(present(u2)) then
    if (u2%idnum==f%idnum) then
      write(0,*) 'savetstate: WARNING! attempting to save u2%idnum=f%idnum cuttype=',t
      !debug return
      !reassign u2%idnum
      call saveintermediate(u2)
      !if it fails, return
      if (u2%idnum==f%idnum) return
    endif
  endif
  IF (.NOT. ASSOCIATED(tlistroot)) THEN !check to see if list is empty
     ALLOCATE(tlistroot)
     NULLIFY(tlistroot%next)
  ENDIF
  ttail => tlistroot
  do while (associated(ttail%next))
     ttail => ttail%next
     if (ttail%parent/=f%idnum) cycle
     if (t==seamflag) then
       if (ttail%child1/=u1%idnum) cycle
     else
       if (ttail%child1/=u1%idnum.and.ttail%child2/=u1%idnum) cycle
       if (ttail%child1/=u2%idnum.and.ttail%child2/=u2%idnum) cycle
     endif
     !! this tstate is identical to a previous one
     if (verbose) then
       write(0,*) "savetstate:: tstate identical to previous one"
     endif
     return
  enddo
  ALLOCATE(ttail%next)
  NULLIFY(ttail%next%next)
  ttail => ttail%next
  ttail%tp = t
  ttail%id = tstatecounter
  ttail%entropy = ent
  ttail%parent = f%idnum
  ttail%child1 = u1%idnum
  ttail%child2 = 0
  if (present (u2)) ttail%child2 = u2%idnum
  ttail%axis = f%axis
  ttail%seam = 0
  if (present(iseam)) ttail%seam = iseam
  tstatecounter = tstatecounter + 1   !need tstate number
  if (verbose) then
     write(0,*) "TSTATE ", tstatecounter, ttail%parent, ttail%child1, ttail%child2, cuttypes(t:t)
  endif
END SUBROUTINE savetstate
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine saveintermediate
!saves intermediate f to a linked list
!of intermediates.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE saveintermediate(f)
  ! save a copy of f to a linked list
  ! We NOW assume that this is a new intermediate!!  Tue Aug 20 19:45:26 EDT 2013
  type(intermediate), POINTER :: f
  TYPE (intermediate), POINTER, save :: itail ! declare tail of the intermediate list
  INTEGER, save :: intercounter=1, id=0
  logical,parameter :: semiverbose=.true.
  character(len=1500) :: tmpstr
  integer :: i
  logical :: debugging=.true.
  !!--
  if (.not.associated(f)) then
     write(0,*) "saveintermediate:: BUG!  nothing to save!! "
     return
  endif
  if (debugging) then
    id = oldintermediate(f)
    if (id /= 0) then
      f%idnum = id
      write(0,*) "WARNING: FOUND AN UNEXPECTED OLD INTERMEDIATE while in saveintermediates() ..."
      return
    endif
  endif
  IF (.NOT. ASSOCIATED(ilistroot)) THEN
     ALLOCATE (ilistroot)
     NULLIFY(ilistroot%next)
     itail => ilistroot   !! empty root 
  END IF
  ALLOCATE(itail%next)
  NULLIFY(itail%next%next)
  itail => itail%next
  f%idnum = intercounter
  f%sym = countsyms(f%iflag,geofold_nres)
  itail = f
  NULLIFY(itail%next)
  if (verbose.or.semiverbose) then
     tmpstr = " "
     do i=1,geofold_nres
       tmpstr = trim(tmpstr)//f%iflag(i)
     enddo
     write(0,'(a,i4,a,a,$)') "Intermediate ",f%idnum," ",trim(tmpstr) 
     write(0,'(a,i4,$)') " depth",f%state
     if (f%sym/=0) write(0,'(a,i2,$)') " sym=",f%sym
     if (any(f%barrel(:)/=0)) write(0,'(a,$)') " seams: "
     do i=1,MAXBARREL
       if (any(f%barrel(i:)/=0)) write(0,'(i2,$)') f%barrel(i)
     enddo
     write(0,*)
  endif
  intercounter = intercounter + 1
END SUBROUTINE saveintermediate
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print intermediates
!to stdout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE printintermediates
  type(intermediate), POINTER :: pptr
  pptr => ilistroot
  output: DO while (associated(pptr%next) )
     pptr => pptr%next
     WRITE(*, *) 'INTERMEDIATE ', pptr%idnum, ' ', pptr%state
  END DO output
END SUBROUTINE printintermediates
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check to see whether f is already in 
! the intermediates list. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION oldintermediate(f)
  type(intermediate), POINTER :: f
  type (intermediate), POINTER,save :: ihead
  integer :: ios=0, n = 0, i, iseam, y1,y2,x1,x2
  logical :: thesame=.false.
  if (.not.associated(f)) then
    write(0,'("Warning:: called oldintermediate(f) with an unassociated pointer!!")')
    stop "Warning:: called oldintermediate(f) with an unassociated pointer!!"
    oldintermediate = 0
    return
  endif
  if (associated(ilistroot)) then
    ihead => ilistroot
  else
    write(0,'("Warning:: oldintermediate :: head of linked list (ilistroot) is missing!!")')
    oldintermediate = 0
    return
  endif
  n = 0
  DO while ( associated(ihead%next) ) 
     ihead => ihead%next
     thesame = .false.
     IF ( checkflags(ihead%iflag,f%iflag) ) THEN
       thesame = .true.
       do i=1,MAXBARREL  !! = 8
         if ( ihead%barrel(i)/=f%barrel(i) ) then ! not the same seam
           thesame = .false.  !!! added to fix bug resulting in several seam moves
           iseam = abs(f%barrel(i))
           if (iseam /= 0) then
             y1=barrels_array(i)%seams(iseam)%segments(1)
             y2=barrels_array(i)%seams(iseam)%segments(2)
             x1=barrels_array(i)%seams(iseam)%segments(3)
             x2=barrels_array(i)%seams(iseam)%segments(4)
             if (any(f%iflag(y1:y2)/=".").and.any(f%iflag(x1:x2)/=".")) then
               thesame = .false.
             endif
           endif
           iseam = abs(ihead%barrel(i))
           if (iseam /= 0) then
             y1=barrels_array(i)%seams(iseam)%segments(1)
             y2=barrels_array(i)%seams(iseam)%segments(2)
             x1=barrels_array(i)%seams(iseam)%segments(3)
             x2=barrels_array(i)%seams(iseam)%segments(4)
             if (any(f%iflag(y1:y2)/=".").and.any(f%iflag(x1:x2)/=".")) then
               thesame = .false.
             endif
           endif
         endif
       enddo
     ENDIF
     if (thesame) then
       oldintermediate = ihead%idnum !! found old intermediate
       return
     endif
     n = n + 1
  ENDDO
  oldintermediate = 0
  !if (verbose) then
  !   write(0,*) "New intermediate added. Number of intermediate so far =",n
  !endif
  return
END FUNCTION oldintermediate
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!determine if newf, oldf are equivalent
!return true or false based on result
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
logical FUNCTION checkflags(newf, oldf)
  CHARACTER, dimension(maxres) :: newf, oldf 
  CHARACTER, dimension(maxres) :: symf
  integer :: nsym, isym, nres
  nres = geofold_nres
  nsym = geofold_nsym
  do isym=0,nsym
    call applysymmetry(isym,newf,symf,nres)
    !  if ( all(symf(1:geofold_nres) == oldf(1:geofold_nres))) then
    if ( all( ((oldf(1:geofold_nres)/='.').and.(symf(1:geofold_nres)/='.')).or. &
              ((oldf(1:geofold_nres)=='.').and.(symf(1:geofold_nres)=='.')) ) ) then
      checkflags = .true.
      return
    endif
  enddo
  checkflags = .false.
END FUNCTION checkflags
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! countsyms()
! count how many syms there are of a given intermediate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
integer FUNCTION countsyms(flags,nres)
  integer,intent(in) :: nres
  CHARACTER, dimension(nres),intent(in) :: flags
  CHARACTER, dimension(nres) :: symf
  integer :: isym, nn, nsym
  nn = 0
  nsym = geofold_nsym
  do isym=1,nsym
    call applysymmetry(isym,flags,symf,nres)
    if ( .not. any((flags/='.').and.(symf/='.'))) then
      nn = nn + 1
    endif
  enddo
  countsyms = nn
END FUNCTION countsyms
!!====================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Apply symmetry operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine applysymmetry(isym,f,symf,nres)
  ! use geofold_pivots, only : geofold_symop
  implicit none
  integer,intent(in) :: isym, nres
  character, dimension(nres), intent(in) :: f
  character, dimension(nres) , intent(out):: symf
  integer :: i
  !!
  if (isym==0) then
    symf(1:nres) = f(1:nres)
    return
  endif
  symf = '.'
  do i=1,geofold_nres
    symf(i) = f(geofold_symop(isym,i))
  enddo
end subroutine applysymmetry

!!====================================================================================

SUBROUTINE cleanuplists()
  TYPE (tstate), POINTER :: tptr, thead ! pointer to tlistroot
  TYPE (intermediate), POINTER ::iptr, ihead
  !!
  if (associated(tlistroot)) then
    tptr => tlistroot
    !!
    do while (associated (tptr%next) )
       thead => tptr%next
       if(associated(tptr)) deallocate(tptr)
       tptr => thead
    end do
    if (associated(tptr)) deallocate(tptr)
  endif
  !!
  if (associated(ilistroot)) then
    iptr => ilistroot
    !!
    do while (associated (iptr%next) )
       ihead => iptr%next
       if(associated(iptr)) deallocate(iptr)
       iptr => ihead
    end do
    if (associated(iptr)) deallocate(iptr)
  endif
  !if (associated(seq)) deallocate(seq)
  !if (associated(resseq)) deallocate(resseq)
  call geofold_cleanup()

END SUBROUTINE cleanuplists

!!====================================================================================

!!Timer added by SAN
SUBROUTINE timer_write(ofile, msg, timer)
  implicit none
  CHARACTER(len=*) :: ofile, msg 
  DOUBLE PRECISION :: timer 
  integer :: out_unit = 8 
  open (unit = out_unit, file = ofile, action = "write", status = "replace")
  write (out_unit, '(A, d15.7)') msg, timer
  close (out_unit)
END SUBROUTINE timer_write


SUBROUTINE dag_write(ofile,ounit)
  implicit none
  TYPE (tstate), POINTER :: tptr !pointer for tlistroot (for later)
  character(len=*),intent(in) :: ofile
  integer,optional,intent(out) :: ounit
  integer :: dunit=8,ios=0,ires,i
  character :: cutchar
  character(len=200) :: pdbcode, lname
  real,dimension(3,3) :: mat
  real,dimension(3) :: vec
  real :: ntrp
  character(len=20),parameter :: aa1="ACDEFGHIKLMNPQRSTVWY"
  !!
  dunit = pickunit(dunit)
  if (present(ounit)) ounit=dunit
  open(dunit,file=ofile,form='formatted',status='replace',iostat=ios)
  if (ios/=0) stop 'geofold:: error opening file for writing. Permissions?'
  !!
  pdbcode = " "
  !! No big deal if these environment varables are not set...
  call getenv("TARG",pdbcode)
  if (pdbcode/=" ") then
     WRITE(dunit, '(a,a)') "HEADER PDBcode= ",trim(pdbcode)
  endif
  lname = " "
  call getenv("LNAME",lname)
  if (lname/=" ") then
     WRITE(dunit, '(a,a)') "HEADER jobname= ",trim(lname)
  endif
  !!--- 1-letter code amino acid sequence
  write(dunit,'("SEQUENCE ",i4," ",$)') geofold_nres
  do i=1,geofold_nres
    write(dunit,'(a1,$)') aa1(seq(i):seq(i))
  enddo
  write(dunit,*)
  !!--- Write PDB format alpha carbon coordinates
  WRITE(dunit, '(a)') "HEADER C-alpha coordinates with original numbering"
  call geofold_writepdb(dunit)
  WRITE(dunit, '(a)') "HEADER tstate tstate_no f u1 u2 entropy cuttype axis"
  if (.not.associated(tlistroot)) then
    write(0,'("Warning: tlistroot is an unassociated pointer!!")')
  else
    tptr => tlistroot ! point to head of tstate list
    outloop: DO while (ASSOCIATED(tptr%next))
       tptr => tptr%next
       select case (tptr%tp)
       case (breakflag)
         cutchar='b'
         ntrp = tptr%entropy!*(1-bcutoff)
         call geofold_getballvec(bvec=tptr%axis,vec=vec)
         WRITE(dunit, '(a,4i7,f12.2,3x,a1,i5)') "TSTATE ",tptr%id, tptr%parent, tptr%child1, &
                                       tptr%child2, ntrp, cutchar, 0
         !! Note: MATRIX lines were intended to be used to generate a movie directly
         !! from the DAG file, but this has not been done. Commenting out for now.
         ! WRITE(dunit, '(a,i7,1x,a1,9f7.3)') "MATRIX ",tptr%id, cutchar, vec
       case (pivotflag)
         cutchar='p'
         ntrp = tptr%entropy!*(1-pcutoff)
         call geofold_getballmat(bvec=tptr%axis,mat=mat)
         WRITE(dunit, '(a,4i7,f12.2,3x,a1,i5)') "TSTATE ",tptr%id, tptr%parent, tptr%child1, &
                                       tptr%child2, ntrp, cutchar, 0
         ! WRITE(dunit, '(a,i7,1x,a1,9f7.3)') "MATRIX ",tptr%id, cutchar, mat
       case (hingeflag)
         cutchar='h'
         ntrp = tptr%entropy!*(1-hcutoff)
         call geofold_gethingemat(ij=tptr%axis,mat=mat,vec=vec)
         WRITE(dunit, '(a,4i7,f12.2,3x,a1,i5)') "TSTATE ",tptr%id, tptr%parent, tptr%child1, &
                                       tptr%child2, ntrp, cutchar, 0
         ! WRITE(dunit, '(a,i7,1x,a1,12f7.3)') "MATRIX ",tptr%id, cutchar, mat, vec
       case (seamflag)
         cutchar='s'
         ntrp = tptr%entropy!*(1-scutoff)
         WRITE(dunit, '(a,4i7,f12.2,3x,a1,i5)') "TSTATE ",tptr%id, tptr%parent, tptr%child1, &
                                       0, ntrp, cutchar, tptr%seam
       case (meltflag)
         cutchar='m'
         WRITE(dunit, '(a,4i7,f12.2,3x,a1,i5)') "TSTATE ",tptr%id, tptr%parent, tptr%child1, &
                                       tptr%child2, tptr%entropy, cutchar, 0
       case default
         cutchar ='u'
       end select 
    END DO outloop
  endif
  !!
  call geofold_masker_intermediates(dunit) ! calculates ISEGMT energies
  !!
  if (.not.present(ounit)) close(dunit)
END SUBROUTINE dag_write
!-----------------------------------------------------------------------------------------------
  function getEnergySeam (aseam,sas,nhb,sce,T) result(energy)
    !! Seam energy includes buried surface area, H-bonds and side chain entropy terms only.
    use geofold_masker
    use geofold_hbonds
    implicit none
    type(seam_type),pointer :: aseam
    real,optional,intent(out) :: sas, sce
    integer,intent(out),optional :: nhb
    integer :: nhbonds
    real,intent(in) :: T
    real :: sasa, energy, hbenergy, scentropy
    !!----------------
    if (.not.associated(aseam)) stop 'geofold.f90:: getEnergySeam BUG,aseam not associated !'
    !! get amount of buried surface area exposed by opening seam
    call geofold_masker_seamenergy(aseam,sasa)
    if (present(sas)) sas = sasa
    !! get number of hydrogen bonds broken by opening seam
    call geofold_hbonds_get(aseam,nhbonds,"A")
    if (present(nhb)) nhb = nhbonds
    !! get sidechain entropy gained on opening seam
    call geofold_masker_getscenergy(aseam,scentropy,seamchar="A")
    if (present(sce)) sce = scentropy
    energy = sasa*geofold_masker_omega + nhbonds*geofold_hbonds_eperbond -  &
             T*(scentropy*geofold_masker_lambdaweight)
  end function getEnergySeam 
!-----------------------------------------------------------------------------------------------
END PROGRAM geofold

