
!! ===================== SPLITSEAMS.F90 ==================================
!! This is a utlity program to process the seams file (from pdb2seams.f90)
!! to a new format for speedier and more accurate modeling of barrels.
!! The old format uses BUTTONs, which are residues that bridge
!! two beta strands in a SEAM. The new format splits the buttons
!! into two sets of flags in the seam_type, u1flag and u2flag,
!! so that some go to beta1 and some go to beta2. The split is
!! heuristic, maximizing the contiguity of the button residues
!! and minimizing the seam unfolding energy.
!!
!! This routine calculates the energy of a seam, however, it uses only surface area
!! and Hbonds. No voids, no sidechain entropy.
!!
!!------------------------------------------------------------------------
!! C. Bystroff  Thu Jun 26 10:46:41 EDT 2014
!!------------------------------------------------------------------------
program splitseams
  use geofold_global !  , only : geofold_nres, barrels_array, seam_type
  use geofold_seams
  use geofold_masker
  use geofold_hbonds
  !----
  implicit none
  integer,parameter :: MAXSEG=20
  character :: flag(MAXRES,MAXSEG)
  integer,parameter :: segment_tolerance=2
  type(seam_type),pointer :: seam, aseam, hems
  integer :: mergeit(2)
  integer :: iseg, jseg, mseg, i,j,k,m,n,t, ires, jres, kres, nh1,nh2
  real :: x, energy, u1x,u2x, eperbond, omega
  integer :: iargc, jarg, dunit, ounit, ios, ib, ierr, nc1=1, nc2=0
  character(len=1000) :: seamsfile, sasfile, newfile, filename, hbfile, parfile, pdbfile
  nullify(hems,seam,aseam)
  hbfile = " "
  jarg = iargc()
  if (jarg < 5) then
    write(*,'("usage: xsplitseams <pdbfile> <seamsfile> <sasfile> '//&
       '<hbondsfile> <parametersfile> <splitseamsfile>")')
    stop 'splitseams.f90 v.  Wed Jun 25 23:02:45 EDT 2014'
  endif
  call getarg(1,pdbfile)
  call getarg(2,seamsfile)
  call getarg(3,sasfile)
  call getarg(4,hbfile)
  call getarg(5,parfile)
  call getarg(6,newfile)
  allocate(aseam)
  aseam%u1flag="."
  aseam%u2flag="."
  dunit = pickunit(10)
  open(dunit, file=pdbfile, form="formatted", status="old", iostat=ierr)
  IF (ierr > 0 ) STOP "splitseams.f90:: Error! File pdbfile not found!"
  write(*,*) "Reading input files: pdb, sas, seams..."
  call geofold_readpdb(dunit)     !! this sets geofold_nres global.
  call geofold_masker_readvoids(dunit)
  close(dunit)
  !! Read seams file, result in barrels_array()
  call geofold_seams_read(seamsfile)
  !! Read hbonds file. Private array used in geofold_hbonds_get
  !!!SEGFAULT HERE, CAUSE UNKNOWN.  CHECKING geofold_hbonds_read
  call geofold_hbonds_read(hbfile)
  !!  Read SAS file, result in  private variable sasnrg.
  !!  Accessed through geofold_masker_seamenergy
  call geofold_masker_read(sasfile)
  dunit = pickunit(10)
  open(dunit, file=parfile, form="formatted", status="old", iostat=ierr)
  IF (ierr > 0 ) STOP "splitseams.f90:: Error! File parfile not found!"
  !!! SEGFAULT HERE TOO
  call geofold_readparameter(dunit,"HBONDENERGY",geofold_hbonds_eperbond)
  eperbond = geofold_hbonds_eperbond
  call geofold_readparameter(dunit,"OMEGA",geofold_masker_omega)
  omega = geofold_masker_omega
  close(dunit)
  !!  Generate u1flag and u2flag arrays
  write(*,*) "Dividing buttons into contiguous segments... tolerance = ",segment_tolerance
  t = segment_tolerance
  !SEGFAULT HERE barrels_array not initialized?
  do k=1,size(barrels_array)
    !! For each barrel, each seam, split buttons between beta1 and beta2 sides.
    write(*,*) "Barrel number ",k," ..."
    do i=1,barrels_array(k)%nSeams
      !! assign all non-seam, non-button residues that contact to u1hem or u2hem
      !! divide all residues into segments that are maximally contiguous
      !! first, initialize two segments for beta1 and beta2
      write(*,*) "   Seam number ",i, barrels_array(k)%seams(i)%segments(1:4), &
                 " has ",barrels_array(k)%seams(i)%nButtons," buttons."
      flag(:,:) = "."
      flag(barrels_array(k)%seams(i)%segments(1):barrels_array(k)%seams(i)%segments(2),1) = "A"
      flag(barrels_array(k)%seams(i)%segments(3):barrels_array(k)%seams(i)%segments(4),2) = "A"
      call sethems(hems,flag(:,1),flag(:,2))
      iseg = 2
      mseg = 2
      BLOOP: do ib=1,barrels_array(k)%seams(i)%nButtons
        !! place each button into a segment
        m = barrels_array(k)%seams(i)%buttons(ib)%residue
        !! remove buttons from hems
        hems%u1flag(m) = '.'
        hems%u2flag(m) = '.'
        !! diagnostic
        write(*,'("                  checking button ",i5,$)') m
        do t=1,segment_tolerance
          !! first try to add buttons to seam betas
          if     (any(flag(m-t:m+t,1)=="A")) then
            flag(m,1) = "A"
            write(*,'("     segment.... ",i5)') 1
            cycle BLOOP
          elseif (any(flag(m-t:m+t,2)=="A")) then
            flag(m,2) = "A"
            write(*,'("     segment.... ",i5)') 2
            cycle BLOOP
          else
            !! then try to add buttons to existing disjoint segments
            iseg = 3
            do while (iseg<=mseg)
              if (any(flag(m-t:m+t,iseg)=="A")) then
                flag(m,iseg) = "A"
                write(*,'("     segment.... ",i5)') iseg
                cycle BLOOP
              endif
              iseg = iseg + 1
            enddo
          endif
          !! if we did not exit loop, increment tolerance and try again
        enddo
        !! if here, then no place was found for this button. Start a new segment
        mseg = mseg + 1
        flag(m,mseg) = "A"
        !! diagnostic
        write(*,'(" new segment.... ",i5)') mseg
      enddo BLOOP !! buttons
      write(*,*) "    segmented into  ",mseg," segments."
      !! merge segments that are within segment_tolerance
      TLOOP: do t=1,segment_tolerance
        iseg = 2
        do while (iseg < mseg)
          iseg = iseg + 1
          do m=1,geofold_nres
            if (flag(m,iseg)==".") cycle
            do jseg=1,iseg-1
              !! see if jseg has a residue within tolerance of iseg. Merge.
              if (any(flag(m-t:m+t,jseg)=="A")) then
                write(*,*) "    merging         ",iseg," with ",jseg
                where (flag(:,iseg)=="A") flag(:,jseg) = "A"
                !! jseg now merged with iseg. Copy last segment to iseg's place.
                flag(:,iseg) = flag(:,mseg)
                !! decrementing mseg, since it has now been copied to iseg
                mseg = mseg - 1
              endif
            enddo
          enddo
        enddo
      enddo TLOOP
      write(*,*) "    merged    into  ",mseg," segments by tolerance."
      !! add hems for energy calculations. hems include buttons.
      where (flag(:,1)=='.') flag(:,1)=hems%u1flag(:)
      where (flag(:,2)=='.') flag(:,2)=hems%u2flag(:)
      !! Get all break energies and select the smallest.
      !! mseg is the current number of segments. When this is 2, we are done.
      do while (mseg > 2)
        !! merge segments in order of decreasing energy
        energy = -999.
        mergeit = (/0,0/)
        do jseg=3,mseg
          aseam%u2flag(1:geofold_nres)=flag(1:geofold_nres,jseg)
          !! try merging with u1
          aseam%u1flag(1:geofold_nres)=flag(1:geofold_nres,1)
          call geofold_masker_seamenergy(aseam,u1x,nc=nc1)
          call geofold_hbonds_get(aseam,nh1)
          !! try merging with u2
          aseam%u1flag(1:geofold_nres)=flag(1:geofold_nres,2)
          call geofold_masker_seamenergy(aseam,u2x,nc=nc2)
          call geofold_hbonds_get(aseam,nh2)
          u1x = u1x + eperbond*nh1
          u2x = u2x + eperbond*nh2
          x = abs(u1x-u2x)
          if (x > energy) then
            energy = x
            iseg = maxloc((/u1x,u2x/),dim=1)
            mergeit = (/iseg,jseg/)
            if (energy<0.001) then
              write(*,'("WARNING: near zero energy for ",i3,i3," nc=",i3)') &
              iseg, jseg, nc1+nc2
            endif
          endif
        enddo
        !! merge mergeit(1:2), note mergeit(1) < mergeit(2)
        !! diagnostic
        iseg=mergeit(1); jseg=mergeit(2)
        if ((iseg==0).or.(jseg==0)) stop 'BUG in splitseams.f90 :: mergeit'
        where (flag(:,jseg)=="A") flag(:,iseg) = "A"
        if (jseg /= mseg) flag(:,jseg) = flag(:,mseg)
        mseg = mseg - 1
        write(*,*) "    merging         ",jseg," with ",iseg," mseg=",mseg," energydiff=",energy
      enddo
      !! fill in u1, u2 segments if gaps are within tolerance.
      t = segment_tolerance
      ires = 0
      do while (ires < geofold_nres)
        ires = ires + 1
        do iseg=1,2
          if (count(flag(ires:ires+t,iseg)=="A")>=2) then
            jres = ires
            do while (flag(jres,iseg)/="A")
              jres = jres + 1
              if (jres > ires+t) stop "BUG in splitseams.f90"
            enddo
            kres = ires + t
            do while (flag(kres,iseg)/="A")
              kres = kres - 1
              if (kres<jres+1) stop 'BUG in splitseams 2'
            enddo
            if (jres+1<kres-1) flag(jres+1:kres-1,iseg) = "A"
          endif
        enddo
      enddo
      !! Save u1/u2 division in u1flag and u2flag
      seam => barrels_array(k)%seams(i)
      seam%u1flag(1:geofold_nres) = flag(1:geofold_nres,1)
      seam%u2flag(1:geofold_nres) = flag(1:geofold_nres,2)
      !! erase hems for energy calculation
      where (flag(1:geofold_nres,1)=='b') flag(1:geofold_nres,1) = '.'
      where (flag(1:geofold_nres,2)=='b') flag(1:geofold_nres,2) = '.'
      aseam%u1flag(1:geofold_nres) = flag(1:geofold_nres,1)
      aseam%u2flag(1:geofold_nres) = flag(1:geofold_nres,2)
      !! get new seam energy, uses private variable sasnrg
      call geofold_masker_seamenergy(aseam, energy) !! energy of unfolding.
      call geofold_hbonds_get(aseam,nh1)
      energy = energy*omega + eperbond*nh1
      seam%energy = energy
      write(*,'(a,$)') "      U1 =  "
      do ires=1,geofold_nres;if (seam%u1flag(ires)=="A") write(*,'(i5,$)') ires; enddo; write(*,*)
      write(*,'(a,$)') "      U2 =  "
      do ires=1,geofold_nres;if (seam%u2flag(ires)=="A") write(*,'(i5,$)') ires; enddo; write(*,*)
      write(*,*) "      Energy =  ", energy
    enddo
  enddo
  dunit = pickunit(10)
  !stop 'works'
  write (0,*) 'newfile ',trim(newfile),', dunit ',dunit,', ios ',ios
  !SEGFAULT HERE?
  open(dunit,file=newfile,status='replace',iostat=ios)
  !stop 'works'
  if (ios/=0) stop 'splitseams.f90 :: error opening new file. Permissions?'
  write(dunit,'("# Processed by splitseams.f90 ")')
  call geofold_seams_write(dunit)
  close(dunit)
  deallocate(aseam)
  dunit = pickunit(10)
  open(dunit, file=parfile, form="formatted", status='old', access="append", iostat=ierr)
  IF (ierr > 0 ) STOP "splitseams.f90:: Error opening parfile to append."
  write(dunit,'(a)') "SEAMS "//trim(newfile)
  close(dunit)
CONTAINS
  subroutine sethems(hems,f1,f2)
    !! assign all non-seam residues to u1 or u2 side of seam, based on energy.
    !! f1 and f2 are the sides of the seam. hems is returned with 'b' in hem residues.
    type(seam_type),pointer :: hems, hseam
    character,dimension(MAXRES),intent(in) :: f1, f2
    real :: u1x, u2x, energy, x
    real,parameter :: tinysas=0.1
    if (.not.associated(hems)) then
      allocate(hems)
    endif
    hems%u1flag = "."
    hems%u2flag = "."
    allocate(hseam)
    do ires=1,geofold_nres
      if (f1(ires)/='.') cycle
      if (f2(ires)/='.') cycle
      hseam%u1flag = f1(1:geofold_nres)
      hseam%u2flag = '.'
      hseam%u2flag(ires) = 'A'
      call geofold_masker_seamenergy(hseam,u1x)
      hseam%u1flag = f2(1:geofold_nres)
      hseam%u2flag = '.'
      hseam%u2flag(ires) = 'A'
      call geofold_masker_seamenergy(hseam,u2x)
      if (u1x<tinysas.and.u2x<tinysas) cycle
      if (u1x>u2x) then
        hems%u1flag(ires) ="b"
      else
        hems%u2flag(ires) ="b"
      endif
    enddo
    deallocate(hseam)
  endsubroutine sethems
  !! --------- read parameter file
  subroutine readparameters(filename,keyword,x)
    implicit none
    character(len=*),intent(in) :: filename,keyword
    real,intent(out) :: x
    integer :: iunit,iw,ios
    character(len=1000) :: aline
    iunit = pickunit(11)
    open(iunit,file=filename,status='old',iostat=ios)
    if (ios/=0) stop 'splitseams.f90:: input file missing.'
    do
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      iw = index(aline,' ')
      if (aline(1:iw-1)/=trim(keyword)) cycle
      read(aline(iw:),*,iostat=ios) x
      if (ios/=0) STOP 'splitseams.f90:: bad float'
      close(iunit)
      return
    enddo
    write(*,'("KEYWORD not found: ",a)') trim(keyword)
    stop
  end subroutine readparameters
end program splitseams
