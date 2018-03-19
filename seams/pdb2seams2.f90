program pdb2seams2
  implicit none
  !!----------------------------------------------------------
  !! REPLACES pdb2seams.f90
  !! C.Bystroff
  !! Wed Oct 11 16:35:08 EDT 2017
  !!----------------------------------------------------------
  !! OUTLINE
  !! Read a Hbond file and a PDB file.
  !! Find seams using anti and parallel beta Hbond operators
  !! Connect pairs of seams that share a residue (edges)
  !! Find barrels = cycles of edges.
  !! Return keyworded SEAM information 
  !! ------------------------------------------------
  !! Hbonds file format described in pdb2hb.f90
  !!----------------------------------------------------------
  type hbtype
    integer:: donornum, acceptnum
    character(len=4) :: donorat, acceptat
    character(len=1) :: bondtype
  end type
  type seamtype
    integer :: idx, start(2), end(2)
    integer :: nbbhb, nschb, nbulge,b1num, b2num
    character(len=12) :: orient
    type(hbtype),dimension(:),pointer :: bbhb
    type(hbtype),dimension(:),pointer :: schb
    type(seamtype),pointer :: next
  end type
  type(hbtype),dimension(:),pointer :: hb
  type(seamtype),dimension(:),pointer :: seams
  type(seamtype),pointer :: seamsroot, aseam
  integer,dimension(:,:),allocatable :: overlap
  integer :: nhb,nseams, this,ios=0,i,j,jarg,nb,nbarrel
  integer,parameter :: MINSEAM=2 !! minimum number of beta H-bonds to count as a seam
  integer,parameter :: MAXBARREL=22 !! Maximum number of seams in a barrel
  integer,dimension(MAXBARREL) :: barrel
  logical :: isthere
  character(len=1000) :: tmpfile,hbfile
  character(len=1000) :: aline
  !!----------------------
  jarg = command_argument_count()
  if (jarg < 1) then
    write(*,*) 'Usage: xpdb2seams hbond_file'
    stop 'pdb2seams2.f90 v.  Tue Oct 24 15:53:57 EDT 2017'
  endif
  call get_command_argument(1,hbfile)
  inquire(file=hbfile,exist=isthere)
  if (.not.isthere) stop 'Hbond_file is missing'
  !!----------------------
  allocate(seamsroot)
  seamsroot%idx = 0
  nullify(seamsroot%next)
  nullify(seamsroot%bbhb)
  nullify(seamsroot%schb)
  call gethbonds(hb,hbfile,nhb)
  call findseams(hb,nhb,seamsroot,nseams)
  call mergeseams(seamsroot,nseams)
  allocate(overlap(nseams,nseams))
  allocate(seams(nseams))
  overlap = 0
  aseam => seamsroot
  do i=1,nseams
    seams(i) = aseam
    if (.not.associated(aseam%next)) exit
    aseam => aseam%next
  enddo
  call connectseams(seams,overlap,nseams)
  this = 1
  nb = 0
  barrel = 0
  nbarrel = 0
  ! It is bad to have a static tmp file
  ! If multiple jobs are running at once, they create a race condition.
  ! Things would end poorly. Let's salt it with a random integer
  call random_seed()
  write(tmpfile,'(a,i10.10,a)') "/tmp/pdb2seams2_",irand(),".tmp"
  open(11,file=trim(tmpfile),status='replace',iostat=ios)
  if (ios/=0) stop 'pdb2seams2.f90 :: permission problems. Do not open a new file.'
  do while (any(overlap/=0)) 
    this = 1
    !! find the first seam that has any overlap with another seam
    do while (all(overlap(this,:)==0))
      this = this + 1
    enddo
    !! diagnostic
    !! do i=1,nseams
    !!   do j=1,nseams
    !!     write(*,'(i1,$)') overlap(i,j)
    !!   enddo
    !!   write(*,*)
    !! enddo
    nb = 1
    barrel(nb) = this
    j = minloc(overlap(this,:),mask=(overlap(this,:)/=0),dim=1)
    if (findbarrels(nseams,overlap,barrel,nb,j)) then
      nbarrel = nbarrel + 1
      call outputDAGlines(11,nseams,seams,barrel,nb,nbarrel)
    endif
  enddo
  close(11)
  write (*,"(A)") "# ========== Results from pdb2seams2.f90 =========="
  write (*,"(A)") "# NBARRELS  <Number of Barrels>"
  write (*,"(A)") "# BARREL :  <id>  <Number of Seams>"
  write (*,"(A)") "# SEAM :    <id>  <Number of Buttons> <Total Energy> <x1-beta1> <x2-beta1> <x1-beta2> <x2-beta2>"
  write (*,"(A)") "# BUTTON :  <id>  <Energy-Beta1> <Energy-beta2>"
  write (*,'("NBARRELS  ",I5)') nbarrel
  open(11,file=tmpfile,status='old',iostat=ios)
  if (ios/=0) stop 'pdb2seams2.f90:: failed to open temporary file /tmp/pdb2seams.tmp'
  do 
    read(11,'(a)',iostat=ios) aline
    if (ios/=0) exit
    write(*,'(a)') trim(aline)
  enddo
  close(11,status='DELETE')
  deallocate(overlap)
  deallocate(seams)
  aseam => seamsroot
  do while (associated(aseam%next))
    seamsroot => aseam%next
    deallocate(aseam)
    aseam => seamsroot
  enddo
  deallocate(aseam)
  stop
  !!----------------------
CONTAINS
  recursive logical function findbarrels(nseams,overlap,barrel,nb,this) result(isabarrel)
    integer,intent(in) :: nseams, this
    integer,intent(inout) :: nb
    integer,dimension(nseams,nseams),intent(inout) :: overlap
    integer,dimension(nseams),intent(inout) :: barrel
    integer :: i,j,last,nbin
    !! ---------------
    !! Trace the connected seams to look for a cycle.
    !! If the connections in overlap() loop back
    !! then we found a barrel.
    !! overlap(i,j) is non-zero if seam i shares any residues with seam j.
    !! barrel() returns the seams composing the barrel.
    !! nb returns the number of seams in the barrel.
    !! ---------------
    nbin = nb
    isabarrel = .false.
    if (this==0) return
    last = barrel(nb)
    !! diagnostic
    !! write(*,*) "LOOKING for barrel adding to ",nb," -th seam =",last
    if (overlap(last,this)==0) stop "BUG BUG BUG calling findbarrels without seam (last,this)"
    overlap(last,this) = 0
    if (overlap(this,last)==0) stop "BUG BUG BUG calling findbarrels without seam (last,this)"
    overlap(this,last) = 0
    !! add this to barrel
    nb = nb + 1
    if (nb <= nseams) barrel(nb) = this  !! guards against outofbounds when adding last seam for 2nd time.
    !! look for a cycle
    if (nb > 2) then
      if (any(barrel(1:nb-2)==this)) then
        isabarrel = .true.
        return
      endif
    endif
    !! no cycle, go to the next edge of this
    j = 1
    do while (j/=0)
      j = minloc(overlap(this,:),mask=(overlap(this,:)/=0),dim=1)
      isabarrel = findbarrels(nseams,overlap,barrel,nb,j)
      if (isabarrel) return
    enddo
    nb = nbin
    return   
  end function findbarrels
  !!---------------------
  !! Read H-bonds from a file. H-bonds are generated
  !! as an edge list, from pdb2hb.f90
  !!
  subroutine gethbonds(hb,hbfile,nhb)
    implicit none
    character(len=*),intent(in) :: hbfile
    integer,intent(out) :: nhb
    type (hbtype),dimension(:),pointer :: hb
    integer :: iunit,ios
    character(len=200) :: aline
    open(newunit=iunit,file=hbfile,status='old',iostat=ios)
    if (ios/=0) stop 'pdb2seams2.f90:: gethbonds: error opening hbfile'
    nhb = 0
    do
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:1)/='!') nhb = nhb + 1
    enddo
    rewind(iunit)
    allocate(hb(nhb),stat=ios)
    if (ios/=0) stop 'pdb2seams2.f90:: gethbonds: error allocating hb'
    nhb = 0
    do
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      !     11 N        99 OH   H
      if (aline(1:1)=='!') cycle
      nhb = nhb + 1
      read(aline,'(i7,1x,a4,i7,1x,a4,1x,a1)') hb(nhb)%donornum,hb(nhb)%donorat, &
        hb(nhb)%acceptnum, hb(nhb)%acceptat, hb(nhb)%bondtype
    enddo
    close(iunit)
  end subroutine gethbonds
  
  subroutine printStack2(stack,hb,nhb)
      implicit none
      integer,intent(in) :: nhb
      type(hbtype),intent(in),dimension(:),pointer :: hb
      integer,intent(in),dimension(:),allocatable :: stack
      integer :: i
      
      write(0,*) "PRINTING STACK"
      do i = 1, nhb
          if(stack(i) == 2) write(0,*) hb(i)%donornum, hb(i)%acceptnum
      enddo
      
  end subroutine printStack2
  !!---------------------
  subroutine findseams(hb,nhb,seamsroot,nseams)
    implicit none
    !! =======================================================
    !! Identify donor-acceptor patterns for beta sheets
    !!       Antiparallel: (i,j) , (j,i), (i-2,j+2), (j+2,i-2)
    !!       Parallel :  (i+1,j),(j+2,i+1),(j,i-1),(i-1,j-2)
    !! =======================================================
    integer,intent(in) :: nhb
    type (hbtype),dimension(:),pointer :: hb
    type(seamtype),pointer :: seamsroot, myseam
    integer,intent(out) :: nseams
    !! =======================================================
    !! Find seams pseudocode
    !!   stack1 = all H-bonds
    !!   until stack1 is empty
    !!     stack2 = {}
    !!     pop top hb from stack1 => stack2
    !!     while (hb=parallel(stack2,stack1) /= 0)
    !!       push hb to stack2
    !!       pop hb from stack1
    !!     end
    !!     if (size(stack2) >= minimum beta sheet) then
    !!       write parallel seam
    !!       cycle
    !!     else
    !!       pop all hb from stack2 back to stack1
    !!       pop top hb from stack1 => stack2
    !!     endif
    !!     while (hb=antiparallel(stack2,stack1) /= 0)
    !!       push hb to stack2
    !!       pop hb from stack1
    !!     end
    !!     if (size(stack2) >= minimum beta sheet) then
    !!       write antiparallel seam
    !!     else
    !!       pop all hb from stack2 back to stack1
    !!       pop top hb from stack1 => stack2
    !!     endif
    !!   end
    !! =======================================================
    integer,dimension(:),allocatable :: stack
    integer :: ihb,jhb,khb,icyc,nbulge
    character(len=12) :: ostr, bstr
    !! ---------------
    allocate(stack(nhb),stat=ios)
    if (ios/=0) stop 'pdb2seams2.f90:: findseams: error allocating stack'
    stack = 1  !! unused hbonds set to 1, used hbonds set to 0, sidechain hbonds set to 3, 
               !! If hbond is part of a seam, temporarily change from stack=1 to  stack=2.
               !! All stack=2 set to 0 after seam is completed.
    !! mark all sidechain Hbonds for later use and to exclude them from seam tests in nexthbond
    do ihb=1,nhb
      if ((hb(ihb)%donorat/="N   ").or.(hb(ihb)%acceptat/="O   ")) stack(ihb) = 3   !! save these sidechain Hbonds in stack=3 for later use.
    enddo
    nseams = 0
    do while (any(stack==1))
      !! ------------------------------------------------------------
      !! starting from the next unused Hbond, search upseam and downseam
      !! Parallel:
      !!  upseam: parallelN --> <repeat>
      !!  downseam: parallelC --> <repeat>
      !! Antiparallel:
      !!   upseam: [bulgeC] --> antiN --> [bulgeN] --> antiC --> <repeat>
      !!   downseam: [bulgeN] --> antiC --> [bulgeC] --> antiN --> <repeat>
      !!   [optional]
      !! ------------------------------------------------------------
      where (stack==2) 
        stack = 0  !! used or discarded as not part of a seam
      endwhere
      !! ---- find first unused backbone Hbond. All unused bb hbonds have stack=1
      ihb = 1
      do while (stack(ihb)/=1); ihb=ihb+1; if (ihb>nhb) exit; enddo
      if (ihb>nhb) exit
      stack(ihb) = 2  !! mark this hbond as part of a potential seam.
      !! !! if (hb(ihb)%donornum==216) write(*,*) "====216 to ",hb(ihb)%acceptnum
      !!----- Parallel seam? Look for next Hbond upseam
      jhb = nexthbond(ihb,nhb,stack,hb,orient="parallelN")
      if (jhb/=0) then  
        !! upseam parallel --------- keep going upseam parallel 
!!         if (hb(ihb)%donornum==216) write(*,*) "====216 is parallelN"
        do while (jhb /= 0) 
          do while (jhb /= 0) 
            khb = jhb
            stack(jhb) = 2
            jhb = nexthbond(jhb,nhb,stack,hb,orient="parallelN")
          enddo
          !! end of parallelN series, look for a single bulge in parallelN sequence
          !! if a bulge exists, see if parallelN continues from there.
          jhb = nexthbond(khb,nhb,stack,hb,orient="bulgeN")
        enddo
      endif
      !! downseam parallel --------- go back and look downseam parallel 
      jhb = nexthbond(ihb,nhb,stack,hb,orient="parallelC")
      if (jhb/=0) then  
        !! upseam parallel --------- keep going upseam parallel 
        !! if (hb(ihb)%donornum==216) write(*,*) "====216 is parallelC"
        do while (jhb /= 0) 
          do while (jhb /= 0) 
            khb = jhb
            stack(jhb) = 2
            jhb = nexthbond(jhb,nhb,stack,hb,orient="parallelC")
          enddo
          !! end of parallelC series, look for a single bulge in parallelC sequence
          !! if a bulge exists, see if parallelC continues from there.
          jhb = nexthbond(khb,nhb,stack,hb,orient="bulgeC")
        enddo
      endif
      if (count(stack==2)>=MINSEAM) then
        ! call printStack2(stack,hb,nhb)
        nseams = nseams+1
        call saveseam(seamsroot,nseams,orient="parallel",stack=stack,nhb=nhb,hb=hb)
        !! diagnostic
        !! write(*,'("Next hbond: ",i3,i5,i5,a,a,a,a,a,a,i8)') ihb,hb(ihb)%donornum, hb(ihb)%acceptnum,&
        !!   ">>",hb(ihb)%donorat,"<<",">>",hb(ihb)%acceptat,"<<  SEAM ",nseams
        cycle  !! no need to look for antiparallel seam if we found a parallel one.
      endif
      !!----- not parallel. Is it antiparallel?
      !! search upseam antiparallel ---------
      !!   upseam: [bulgeC] --> antiN --> [bulgeN] --> antiC --> <repeat>
      nbulge = 0
      jhb = nexthbond(ihb,nhb,stack,hb,orient="bulgeC")
      if (jhb==0) then
        !! no bulge, no worries, start again
        jhb = ihb
      else
        !! found a bulge, count it.
        stack(jhb) = 2
        nbulge = nbulge + 1
        !! if (hb(ihb)%donornum==216) write(*,*) "====216 is bulgeC"
      endif
      ostr = "antiN"
      jhb = nexthbond(jhb,nhb,stack,hb,orient=ostr)
      icyc = 0
      do while (jhb /= 0) 
        stack(jhb) = 2
        !! if (hb(ihb)%donornum==216) write(*,*) "====216 is ",trim(ostr)
        !! antiparallel beta is the same upseam and down,
        !! alternate antiN (or bulgeC), then antiC (or bulgeN)
        icyc = mod(icyc+1,2)
        if (icyc==0) then
          bstr = "bulgeC"
          ostr = "antiN"
        else
          bstr = "bulgeN"
          ostr = "antiC"
        endif
        khb = nexthbond(jhb,nhb,stack,hb,orient=bstr)
        if (khb==0) then
          khb = jhb
        else
          !! found a bulge
          stack(khb) = 2
          nbulge = nbulge + 1
          !! if (hb(ihb)%donornum==216) write(*,*) "====216 is ",trim(bstr)
        endif
        jhb = nexthbond(khb,nhb,stack,hb,orient=ostr)
      enddo
      !! downseam antiparallel ---------
      !!   downseam: [bulgeN] --> antiC --> [bulgeC] --> antiN --> <repeat>
      jhb = nexthbond(ihb,nhb,stack,hb,orient="bulgeN")
      if (jhb==0) then
        !! no bulge, no worries, start again
        jhb = ihb
      else
        !! found a bulge
        stack(jhb) = 2
        nbulge = nbulge + 1
        !! if (hb(ihb)%donornum==216) write(*,*) "====216 is bulgeN"
      endif
      ostr = "antiC"
      jhb = nexthbond(jhb,nhb,stack,hb,orient=ostr)
      icyc = 0
      do while (jhb /= 0) 
        stack(jhb) = 2
        !! if (hb(ihb)%donornum==216) write(*,*) "====216 is ",trim(bstr)
        icyc = mod(icyc+1,2)
        if (icyc==0) then
          bstr = "bulgeN"
          ostr = "antiC"
        else
          bstr = "bulgeC"
          ostr = "antiN"
        endif
        khb = nexthbond(jhb,nhb,stack,hb,orient=bstr)
        if (khb==0) then
          khb = jhb
        else
          stack(khb) = 2
          nbulge = nbulge + 1
          !! if (hb(ihb)%donornum==216) write(*,*) "====216 is ",trim(bstr)
        endif
        jhb = nexthbond(khb,nhb,stack,hb,orient=ostr)
      enddo
      !! done looking for a antiparallel seam. Did we find it?
      if (count(stack==2)>=MINSEAM) then
        ! call printStack2(stack,hb,nhb)
        nseams = nseams+1
        call saveseam(seamsroot,nseams,orient="antiparallel",stack=stack,nhb=nhb,nbulge=nbulge,hb=hb)
        !! diagnostic
        !! write(*,'("Next hbond: ",i3,i5,i5,a,a,a,a,a,a,i8)') ihb,hb(ihb)%donornum, hb(ihb)%acceptnum,&
        !!   ">>",hb(ihb)%donorat,"<<",">>",hb(ihb)%acceptat,"<<  SEAM ",nseams
      else
        !! write(*,'("Next hbond: ",i3,i5,i5,a,a,a,a,a,a,i8)') ihb,hb(ihb)%donornum, hb(ihb)%acceptnum,&
        !!   ">>",hb(ihb)%donorat,"<<",">>",hb(ihb)%acceptat,"<<  NOT A SEAM "
      endif
    enddo
    !! Done finding seams. Saved in seamsroot
    !!    
    !! ------------ get start and end of each beta strand
    !! ------------ get sidechain hbonds
    myseam => seamsroot
    if(nseams > 0) then
        ! write(0,*) "431"
        ! write(0,*) nseams
        call getstartend(myseam)
        call getschbonds(myseam,hb,nhb,stack)
        do while (associated(myseam%next)) 
          myseam => myseam%next
          ! write(0,*) "436"
          call getstartend(myseam)
          call getschbonds(myseam,hb,nhb,stack)
        enddo
    endif
    !! diagnostic
    !! write(*,*) "Found ",nseams," seams."
    !!       
    deallocate(stack)
  end subroutine findseams
  !!-----------------------------------------------------------------------------
  !!-----------------------------------------------------------------------------
  subroutine getschbonds(myseam,hb,nhb,stack)
        implicit none
        integer,dimension(nhb),intent(in) :: stack
        type (seamtype),pointer :: myseam
        integer,intent(in) :: nhb
        type (hbtype),dimension(:),pointer :: hb
        integer :: nsc,i
        !! use global stack(), in which sc hbonds are set to 3
        nsc = 0
        do i=1,nhb
          if (stack(i)/=3) cycle
          if (hb(i)%donornum>=myseam%start(1).and.hb(i)%donornum<=myseam%end(1)) then
            if (hb(i)%acceptnum>=myseam%start(2).and.hb(i)%acceptnum<=myseam%end(2)) then
              nsc = nsc + 1
            endif
          endif
          if (hb(i)%donornum>=myseam%start(2).and.hb(i)%donornum<=myseam%end(2)) then
            if (hb(i)%acceptnum>=myseam%start(1).and.hb(i)%acceptnum<=myseam%end(1)) then
              nsc = nsc + 1
            endif
          endif
        enddo
        myseam%nschb = nsc
        allocate(myseam%schb(nsc))
        nsc = 0
        do i=1,nhb
          if (stack(i)/=3) cycle
          if (hb(i)%donornum>=myseam%start(1).and.hb(i)%donornum<=myseam%end(1)) then
            if (hb(i)%acceptnum>=myseam%start(2).and.hb(i)%acceptnum<=myseam%end(2)) then
              nsc = nsc + 1
              myseam%schb(nsc) = hb(i)
            endif
          endif
          if (hb(i)%donornum>=myseam%start(2).and.hb(i)%donornum<=myseam%end(2)) then
            if (hb(i)%acceptnum>=myseam%start(1).and.hb(i)%acceptnum<=myseam%end(1)) then
              nsc = nsc + 1
              myseam%schb(nsc) = hb(i)
            endif
          endif
        enddo
  end subroutine getschbonds
  !!-----------------------------------------------------------------------------
  !!-----------------------------------------------------------------------------
  subroutine getstartend(myseam)
        implicit none
        type (seamtype),pointer :: myseam
        integer :: i,j,k,l,m
        !!pseudocode------
        !! ---- find the min -> start(1)
        !! ---- find the max less that div -> end(1)
        !! ---- find the min greater than div -> start(2)
        !! ---- find the max -> end(2)
        j = 999; k=-999; l=999; m=-999
        do i=1,myseam%nbbhb
           j = min(j,    myseam%bbhb(i)%donornum,myseam%bbhb(i)%acceptnum)
           k = max(k,min(myseam%bbhb(i)%donornum,myseam%bbhb(i)%acceptnum))
           l = min(l,max(myseam%bbhb(i)%donornum,myseam%bbhb(i)%acceptnum))
           m = max(m,    myseam%bbhb(i)%donornum,myseam%bbhb(i)%acceptnum)
           !! diagnostic
           !! if (myseam%idx==7) then
           !!   write(*,*) "idx=7 ",i, myseam%bbhb(i)%donornum, myseam%bbhb(i)%acceptnum,&
           !!                          myseam%bbhb(i)%donorat, myseam%bbhb(i)%acceptat
           !! endif
        enddo
        myseam%start(1) = j
        myseam%end(1) = k
        myseam%start(2) = l
        myseam%end(2) = m
        !! diagnostic
        !! write(*,*) "In getstartend nbbhb=",myseam%idx,myseam%nbbhb,j,k,l,m
      end subroutine getstartend
      !!-----------------------------------------------------------------------------
      !!-----------------------------------------------------------------------------
      subroutine saveseam(seamsroot,nseams,orient,stack,nhb,nbulge,hb)
        implicit none
        integer,intent(in) :: nseams
        character(len=*),intent(in) :: orient
        integer,optional,intent(in) :: nbulge
        type(seamtype),pointer :: seamsroot
        type(seamtype),pointer :: myseam
        integer,intent(in) :: nhb
        type (hbtype),dimension(:),pointer :: hb
        integer,dimension(nhb),intent(in) :: stack
        integer :: ihb,mhb
        !! uses globals seams(), stack(), nhb
        myseam => seamsroot
        do while (associated(myseam%next)) 
          myseam => myseam%next
        enddo
        if (myseam%idx==0) then  
              !! seamsroot. Populate the root node.
        else  !! Last seam saved. Make next.
          allocate(myseam%next)
          myseam => myseam%next
          nullify(myseam%next)
          nullify(myseam%bbhb)
          nullify(myseam%schb)
        endif
        myseam%start(1) = 0
        myseam%start(2) = 0
        myseam%end(1) = 0
        myseam%end(2) = 0
        myseam%idx = nseams
        mhb = count(stack==2)
        myseam%nbbhb = mhb
        myseam%nschb = 0  !! add later
        myseam%nbulge = 0  !! add later
        if (present(nbulge)) myseam%nbulge = nbulge
        allocate(myseam%bbhb(mhb))
        ihb = 0
        do i=1,size(stack)
          if (stack(i)/=2) cycle
          ihb = ihb + 1
          if (ihb>nhb) stop 'BUG. ihb > nhb'
          myseam%bbhb(ihb) = hb(i)  
          !! diagnostic
          !! write(*,*) "In saveseam: ",i,stack(i),hb(i)%donornum,hb(i)%acceptnum,&
          !!            hb(i)%donorat,hb(i)%acceptat
        enddo
        myseam%orient = trim(orient)
        !! seam saved.
        !! diagnostic
        !! write(*,*) " idx  start1  end1  start2 end2 / nbbhb  orient"
        !! write(*,*) myseam%idx,myseam%start(1),myseam%end(1),myseam%start(2),myseam%end(2)
        !! write(*,*) myseam%nbbhb, myseam%orient
  end subroutine saveseam
  !!-----------------------------------------------------------------------------
  !!-----------------------------------------------------------------------------
  integer function nexthbond(jhb,nhb,stack,hb,orient) result(idx)
        implicit none
        integer,intent(in) :: nhb, jhb
        integer,dimension(nhb),intent(in) :: stack
        type (hbtype),dimension(:),pointer :: hb
        character(len=*),intent(in) :: orient
        !!----------------------------
        !! ask whether hb contains a donor acceptor
        !! pair that is related to hb(jhb) by one of the
        !! following operators and not already
        !! marked.
        !! Parallel :  (i+1,j),(j+2,i+1),(j,i-1),(i-1,j-2)
        !!  (i,j) => (j,i-2)
        !!  (i+1,j) => (j,i-1)
        !!  (j+2,i+1) => (i+1, j)
        !!  (j,i-1) => (i-1,j-2)
        !! ParallelN operator
        !!  (0  1  |  0 ) i
        !!  (1  0  | -2 ) j  = ( j, i-2 )
        !! ParallelC operator
        !!  (0  1  |  2 ) i
        !!  (1  0  |  0 ) j    = ( j+2, i)
        !! antiN operator
        !!  (0  1  |  0 ) i
        !!  (1  0  |  0 ) j    = ( j, i)
        !! antiC operator
        !!  (0  1  | +2 ) i
        !!  (1  0  | -2 ) j    = ( j+2, i-2)
        !! bulgeN operator
        !!  (1  0  | -1 ) i
        !!  (0  1  |  0 ) j    = ( i-1, j)
        !! bulgeC operator
        !!  (1  0  | +1 ) i
        !!  (0  1  |  0 ) j    = ( i+1, j)
        !!----------------------------
        integer :: ihb
        if (jhb==0) stop 'BUG BUG BUG in nexthbond. jhb==0'
        idx = 0
        do ihb=1,nhb
          if (ihb==jhb) cycle
          if (stack(ihb)/=1) cycle
          select case (trim(orient))
          case ("parallelN")
            if (hb(ihb)%donornum == hb(jhb)%acceptnum .and. hb(ihb)%acceptnum == hb(jhb)%donornum-2) then
              idx = ihb
              return
            endif
          case ("parallelC")
            if (hb(ihb)%donornum == hb(jhb)%acceptnum+2 .and. hb(ihb)%acceptnum == hb(jhb)%donornum) then
              idx = ihb
              return
            endif
          case ("antiN")
            if (hb(ihb)%donornum == hb(jhb)%acceptnum .and. hb(ihb)%acceptnum == hb(jhb)%donornum) then
              idx = ihb
              return
            endif
          case ("antiC")
            if (hb(ihb)%donornum == hb(jhb)%acceptnum+2 .and. hb(ihb)%acceptnum == hb(jhb)%donornum-2) then
              idx = ihb
              return
            endif
          case ("bulgeN")
            if (hb(ihb)%donornum == hb(jhb)%donornum-1 .and. hb(ihb)%acceptnum == hb(jhb)%acceptnum) then
              idx = ihb
              return
            endif
          case ("bulgeC")
            if (hb(ihb)%donornum == hb(jhb)%donornum+1 .and. hb(ihb)%acceptnum == hb(jhb)%acceptnum) then
              idx = ihb
              return
            endif
          end select
        enddo
  end function nexthbond
  !!--------------------
  subroutine connectseams(seams,overlap,nseams)
      implicit none
      !!------------------------------------
      !! Input seams(nseams) is an array of pointers of seamtype.
      !! Reads seams start and end residues for strands 1 and 2
      !! Draws a connection between seams if either strand has
      !!   overlapping residue ranges.
      !! Returns a matrix of 0,1,2 indicating no overlap, strand 1 or strand 2
      !!   overlap between i (row) and j (column).
      !! PARAMETER: novr is the additional number of residues of overlapping
      !!   range required. If one residue is considered too weak,
      !!   then set novr=1 to require at least 2 residues overlap.
      !!-----------------------------------
      integer,intent(in) :: nseams
      type(seamtype),dimension(:),pointer :: seams
      integer,dimension(nseams,nseams),intent(inout) :: overlap
      integer,parameter :: novr=1  !! additional required overlap in sequence
      integer :: i,j
      do i=1,nseams
        overlap(i,i) = 0
        do j=1,i-1
          if     (seams(i)%start(1)+novr<=seams(j)%end(1)  .and. &
                  seams(i)%end(1)-novr>=seams(j)%start(1)) then
            overlap(i,j) = 1  !! seam i connects to j, strand 1 to strand 1
            overlap(j,i) = 1
            !! diagnostic
            !! write(*,*) "1__",seams(i)%idx, seams(i)%start(1), seams(i)%end(1), &
             !!                 seams(j)%idx, seams(j)%start(1), seams(j)%end(1)
         elseif (seams(i)%start(2)+novr<=seams(j)%end(1)  .and. &
                  seams(i)%end(2)-novr>=seams(j)%start(1)) then
            overlap(i,j) = 2  !! seam i connects to j, strand 2 to strand 1
            overlap(j,i) = 1
            !! write(*,*) "2__",seams(i)%idx, seams(i)%start(2), seams(i)%end(2), &
            !!                  seams(j)%idx, seams(j)%start(1), seams(j)%end(1)
        elseif (seams(i)%start(1)+novr<=seams(j)%end(2)  .and. &
                  seams(i)%end(1)-novr>=seams(j)%start(2)) then
            overlap(i,j) = 1  !! seam i connects to j, strand 1 to strand 2
            overlap(j,i) = 2
            !! write(*,*) "3__",seams(i)%idx, seams(i)%start(1), seams(i)%end(1), &
            !!                  seams(j)%idx, seams(j)%start(2), seams(j)%end(2)
        elseif (seams(i)%start(2)+novr<=seams(j)%end(2)  .and. &
                  seams(i)%end(2)-novr>=seams(j)%start(2)) then
            overlap(i,j) = 2  !! seam i connects to j, strand 2 to strand 2
            overlap(j,i) = 2
            !! write(*,*) "4__",seams(i)%idx, seams(i)%start(2), seams(i)%end(2), &
            !!                  seams(j)%idx, seams(j)%start(2), seams(j)%end(2)
          else
            overlap(i,j) = 0
            overlap(j,i) = 0
          endif
        enddo
      enddo
      do i = 1, nseams
          ! write(0,*) overlap(i,:)
      enddo
  end subroutine connectseams
  !!--------------------
  subroutine outputDAGlines(iunit,nseams,seams,barrel,nb,nbarrel)
      implicit none
      integer,intent(in) :: iunit
      integer,intent(in) :: nseams, nb,nbarrel
      type (seamtype),dimension(:),pointer :: seams
      integer,dimension(nseams),intent(in) :: barrel
      integer :: ib
      !----------------------
      write(iunit,'("BARREL ",I5," ",I5)') nbarrel, nb-1
      do ib=1,nb-1
        aseam => seams(barrel(ib))
        write (iunit, '("SEAM    ",I5,I5,f10.3,4I5)')  &
             aseam%idx, 0, 0.0, aseam%start(1), aseam%end(1), &
             aseam%start(2), aseam%end(2)
      enddo
  end subroutine outputDAGlines
  
  subroutine mergeseams(root,nseams)
      implicit none
      integer,intent(inout)::nseams
      type(seamtype),intent(inout),pointer :: root
      type(seamtype),pointer :: itr,jtr
      integer :: i,j,k
      integer,parameter :: over = 1
      integer :: newstart(2),newend(2)
      
      ! mess with the linked list
      itr => root
      loop1: do while(associated(itr%next))
          jtr => itr%next
          do while(associated(jtr))
              ! write(0,'("itr%idx: ",i4," jtr%idx: ",i4)')itr%idx,jtr%idx
              !1-3, 2-4
              if(itr%end(1)-jtr%start(1) <= over .and. itr%end(1)-jtr%start(1) >= 0&
              .and. itr%start(2)-jtr%start(2) <= over .and. itr%start(2)-jtr%start(2) >= 0) then
                
                ! start
                newstart(1) = itr%start(1)
                newstart(2) = itr%start(2)
                ! end
                newend(1) = jtr%end(1)
                newend(2) = jtr%end(2)
                call merge(itr,jtr,root,newstart,newend)
                if(.not. associated(jtr)) exit loop1
              ! 1-4, 2-3
          elseif(itr%end(1)-jtr%start(2) <= over .and. itr%end(1)-jtr%start(2) >= 0&
              .and. itr%end(2)-jtr%start(1) <= over .and. itr%end(2)-jtr%start(1) >= 0) then
                ! start
                newstart(1) = itr%start(1)
                newstart(2) = itr%start(2)
                ! end
                newend(1) = jtr%end(2)
                newend(2) = jtr%end(1)
                call merge(itr,jtr,root,newstart,newend)
                if(.not. associated(jtr)) exit loop1
              ! 3-1, 4-2
          elseif(jtr%end(1)-itr%start(1) <= over .and. jtr%end(1)-itr%start(1) >= 0&
              .and. jtr%end(2)-itr%start(2) <= over .and. jtr%end(2)-itr%start(2) >= 0) then
                ! start
                newstart(1) = jtr%start(1)
                newstart(2) = jtr%start(2)
                ! end
                newend(1) = itr%end(1)
                newend(2) = itr%end(2)
                call merge(itr,jtr,root,newstart,newend)
                if(.not. associated(jtr)) exit loop1
              ! 4-1, 3-2
          elseif(jtr%end(2)-itr%start(1) <= over .and. jtr%end(2)-itr%start(1) >= 0&
              .and. jtr%end(1)-itr%start(2) <= over .and. jtr%end(1)-itr%start(2) >= 0) then
                ! start
                newstart(1) = jtr%start(2)
                newstart(2) = jtr%start(1)
                ! end
                newend(1) = itr%end(1)
                newend(2) = itr%end(2)
                call merge(itr,jtr,root,newstart,newend)
                if(.not. associated(jtr)) exit loop1
              else  
                jtr => jtr%next
              endif
          enddo
          itr => itr%next
      enddo loop1
      nullify(itr)
      nullify(jtr)
      ! write(0,*) "Mergeseams complete"
  end subroutine mergeseams
  
  subroutine merge(itr,jtr,root,newstart,newend)
      !Do the common operations for merging seams itr and jtr into seam
      implicit none
      integer,intent(in) :: newstart(2),newend(2)
      type(seamtype),intent(inout),pointer :: itr,jtr
      type(seamtype),intent(inout),pointer :: root
      type(seamtype),pointer :: ktr
      integer :: i,j,k
      type(hbtype),dimension(:),pointer :: bbhb,schb
      ! write(0,'("Merging seams: ",2i4)')itr%idx,jtr%idx
      
      
      ! nbbhb, nschb, nbulge
      itr%nbulge = itr%nbulge + jtr%nbulge
      allocate(bbhb(itr%nbbhb+jtr%nbbhb))
      i = 1
      do while(i <= itr%nbbhb)
          bbhb(i) = itr%bbhb(i)
          i = i + 1
      enddo
      j = 0
      do while(j < jtr%nbbhb)
          bbhb(i+j) = jtr%bbhb(j+1)
          j = j + 1
      enddo
      itr%nbbhb = itr%nbbhb+jtr%nbbhb
      if(associated(itr%bbhb)) deallocate(itr%bbhb)
      itr%bbhb => bbhb
      ! schb
      allocate(schb(itr%nschb+jtr%nschb))
      i = 1
      do while(i <= itr%nschb)
          schb(i) = itr%schb(i)
          i = i + 1
      enddo
      j = 0
      do while(j < jtr%nschb)
          schb(i+j) = jtr%schb(j+1)
          j = j + 1
      enddo
      itr%nschb = itr%nschb+jtr%nschb
      if(associated(itr%schb)) deallocate(itr%schb)
      itr%schb=>schb
      
      ! upstream reindexing
      if(associated(itr%next)) then
          ktr => itr
          k = itr%idx
          do while(associated(ktr%next))
              if(associated(ktr%next,jtr)) then
                  if(associated(jtr%next)) then
                      ktr%next => jtr%next
                  else
                      nullify(ktr%next)
                      exit
                  endif
              endif
              ktr => ktr%next
              k = k + 1
              ktr%idx = k
          enddo
          if(associated(ktr)) nullify(ktr)
      endif
      ktr => jtr%next
      if(associated(jtr)) deallocate(jtr)
      jtr => ktr
      if(associated(ktr)) nullify(ktr)
  end subroutine merge
!!------------------------------------------------------------------
!! NOTE: outputDAGseams should be synched with geofold_seams.f90
!! writeBarrels()
!!------------------------------------------------------------------
!=====================================================================
end program pdb2seams2
