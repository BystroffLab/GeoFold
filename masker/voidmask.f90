program VoidMask
use masker
!! Find the 1.4A and 1.3A surfaces of the input atoms
!! Using only the centers of the re-entrant surfaces, substract
!! the centers of the 1.4A true-solvent surface from the centers of the
!! 1.3A surface to get the centers that correspond to maximal inaccessible spaces.
!! INPUT: coordinates in PDB format:
!!  Each atom field must match one of the lines in ATOMLIB
!!  ATOM lines and HETATM lines are accepted.
!! C.Bystroff
!!Sun Jul  2 11:36:42 EDT 2006
!! 
!! Modifications:---------
!! 23 June 2008
!! Added two routines, escapevoid() to check for voids on the surface,
!! and vcluster() to reduce multiple voides to a single one, with radius of gyration.
!! C.B.
!!
character(len=80) :: xyzfile, aline,outfile
real,dimension(:,:),allocatable :: xyz
integer :: i,j,k,icycle,ncycle,naccept,nseed,jarg,jj,n,nat
real :: tt, sig=0.3,x,y,z,nrg,scl,xmin,xmax,ymin,ymax,zmin,zmax,rwi
real,dimension(3) :: vec
real,dimension(:,:),pointer :: wats, watv, maskxyz
real :: vcutoff, ff
real,dimension(200) :: rg
integer :: iargc, time, nwatv, nwats, nv
logical :: drawingatall = .false., docking=.false.
xmin = 999.; ymin = 999.; zmin = 999.
xmax = -999.; ymax = -999.; zmax = -999.

!masker_rendering = .false.
!masker_smoothing = .false.
!masker_drawing = drawingatall
masker_saying = .false.  !! masker public variable
!! Readcommand line
tt = 100.
ncycle = 1
nseed = time()
jarg = iargc()
lbox = 0.00
rw = 1.4  !! global
rwi = 1.3
vcutoff = 0.5
if (jarg < 2) then
  write(*,'("Usage: xvoidmask coordfile outfile [water-proberadius inner-proberadius vcut]")')
  write(*,'(a,f7.2)') 'water proberadius rw=',rw
  write(*,'(a,f7.2)') 'inner proberadius rwi=',rwi
  write(*,'(a,f7.2)') 'cutoff for centers   =',vcutoff
  write(*,'("CURRENT MASK SIZE: ",i5)') MAXATOM
  write(*,*) 'voidmask.f90 Sun Jul  2 12:06:18 EDT 2006'
  stop
endif
call getarg(1,xyzfile)
call getarg(2,outfile)
if (jarg >=3 ) then
  call getarg(3,aline)
  read(aline,*,iostat=ios) rw
  if (ios/=0) stop 'Bad argument 3: rw'
endif
if (jarg >=4 ) then
  call getarg(4,aline)
  read(aline,*,iostat=ios) rwi
  if (ios/=0) stop 'Bad argument 4: rwi'
endif
if (jarg >=5 ) then
  call getarg(5,aline)
  read(aline,*,iostat=ios) vcutoff
  if (ios/=0) stop 'Bad argument 5: vcutoff'
endif
lbox = 0.
tt = 1.
!! Initialize masks
call masker_initmasks
call masker_render(.false.,smooth=.false.)
!! Read coordinates
open(12,file=xyzfile,status='old',form='formatted')
i = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(22:22)=="V") cycle
  if (aline(1:5)=='ATOM ') i = i + 1
  if (aline(1:6)=='HETATM') i = i + 1
enddo
nat = i
!! diagnostic
write(*,*) 'Number of atoms:', nat
rewind(12)
allocate(xyz(3,nat),stat=ios)
if (ios/=0) stop 'Error allocating xyz'
call masker_allocatoms(nat)
i = 0
!! vec = 0.
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    if (aline(22:22)=="V") then
      write(*,'("WARNING:: ignoring chain V in input file!!!!!!!!!!!")')
      cycle
    endif
    i = i + 1
    read(aline(31:54),'(3f8.3)',iostat=ios) xyz(1:3,i)
    if (ios/=0) stop 'Format error xyzfile'
    jj = 0
    if (aline(17:17)/=" ") then
      write(*,'("Warning!! alt loc indicator: ",a)') trim(aline)
      aline(17:17) = " "
    endif
    if (aline(14:17)=="OXT ") aline(14:17) = "OE1 "
    !!!  WZS HERE
    if(.not. wzs) then
      do j=1,nattype
        if (aline(14:17)==atomlib(j)%name(1:4)) then
          jj = j
          exit
        endif
      enddo
    else
      jj = getatype(aline)
    endif
    if (jj==0) then
      write(*,'("Unknown atom type ignored: ",a)') trim(aline)
      atype(i) = -1   !! dont draw the surface
    else
      atype(i) = jj
    endif
    !! diagnostic
    !! write(*,*) 'Atom ',i,jj,xyz(1:3,i)
    !! vec = vec + xyz(1:3,i)
  endif
enddo
write(*,*) i,' Atoms read.'
nat = i
!!

write(*,*) "Getting solvent accessible surface...."
call masker_getms(xyz,nat,x)
write(*,'(a,$)') "Finding reentrant surfaces...."
call masker_getwaters(wats,nwats)
write(*,'(i9,a)') nwats," reentrant points found."
rw = rwi
write(*,*) "Projecting surface inward to probe radius=",rwi
call masker_getms(xyz,nat,x)
write(*,'(a,$)') "Finding deep reentrant surfaces...."
call masker_getwaters(watv,nwatv)
write(*,'(i9,a)') nwatv," reentrant points found."

rewind(12)
open(13,file=outfile,status='replace',form='formatted',iostat=ios)
if (ios/=0) stop 'VoidMask:: error opening file for output. Permissions??'
i = 0
do
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(22:22)=="V") cycle
  if ((aline(1:5)=='ATOM ').or.(aline(1:6)=='HETATM')) then
    write(13,'(a)') trim(aline)
    i = i + 1
  else
    write(13,'(a)') trim(aline)
  endif
enddo
close(12)
write(*,'(i9," model atoms output.")') i

write(*,*) "Finding void spaces.....", nwatv
n = 0
do i=1,nwatv
  if (.not.any(sqrt( (watv(1,i)-wats(1,1:nwats))**2 + &
                     (watv(2,i)-wats(2,1:nwats))**2 + &
                     (watv(3,i)-wats(3,1:nwats))**2 ) < vcutoff)) then
    n = n + 1
    watv(1:3,n) = watv(1:3,i)
  endif
enddo
write(*,*) "Clustering void spaces.....", n
call vcluster(watv,n,nv,rg)
call masker_getmxyz(maskxyz)
n = 0
do i=1,nv
    !! See if this void is buried. Translate it. If it is free to float away, it's not buried.
    call escapevoid(watv(1:3,i), xyz, nat, maskxyz, MAXATOM, ff)   
    if (ff > 0.2) cycle
    write(13,'("HETATM",i5,2x,a1,3x,a3,1x,a1,i4,4x,3f8.3,2f6.2)') &
                9000+i,"O", "HOH","V",i,watv(1:3,i), ff, rg(i)*10.
    n = n + 1
enddo
close(13)
  
write(*,'(i9," void atoms output as HOH, chain V.")') n
write(*,*) "To see just the protein and the void volumes in RasMol, use:"
write(*,*) "restrict protein "
write(*,*) "wireframe"
write(*,*) "select :V"
write(*,*) "spacefill"
if (associated(wats)) deallocate(wats)
if (associated(watv)) deallocate(watv)
if (associated(maskxyz)) deallocate(maskxyz)
call masker_deallo()

CONTAINS
!----------------------------------------------------!
    subroutine escapevoid(v, xyz, N, maskxyz, M, fout)   
       !! translate the void (v) in all (N) directions (maskxyz)
       !! check distance to all (N) atoms (xyz). If the void
       !! can translate 10A without a collision, it has escaped
       !! to solvent and is not a buried void.
       !! CB 17-JEN-2008
       implicit none
       integer,intent(in) :: N,M
       real,dimension(3),intent(in) :: v
       real,dimension(3,N),intent(in) :: xyz
       real,dimension(3,M),intent(in) :: maskxyz
       real,intent(out) :: fout
       real,dimension(3) :: x, y
       real :: dd
       integer :: i,j,k, nout
       nout = 0
       MLOOP: do i=1,M
         do j=2,10,2
           x = v + maskxyz(1:3,i)*j
           do k=1,N
             y = x - xyz(1:3,k)
             dd = sqrt(sum(y*y))
             if (dd < 2.0) cycle MLOOP
           enddo
         enddo
         nout = nout + 1
       enddo MLOOP
       fout = real(nout)/real(M)
       !! diagnostic
       !! write(0,*) "escapes: ",nout,M,fout
      
    end subroutine escapevoid
!----------------------------------------------------!
    subroutine vcluster(watv,n,nv,rg)
      integer,intent(in) :: n
      integer,intent(out) :: nv
      real,dimension(3,n),intent(inout) :: watv
      real,dimension(3) :: v, wtmp
      real,dimension(n),intent(out) :: rg
      real,dimension(n,n) :: cntcts
      integer :: h,k,i,j,cluster(n),flag(n),nc
      real :: x, dd
      wtmp = 0.
      cntcts = 999.
      do i=1,n 
        do j=i+1,n 
          v = watv(1:3,i) - watv(1:3,j)
          dd = sqrt(sum(v*v))
          cntcts(i,j) = dd
          cntcts(j,i) = dd
        enddo
      enddo
      cluster = 0
      flag = 0
      nv = 0
      nc = 0
      do i=1,n
        if (flag(i)/=0) cycle
        k = 1
        cluster(k) = i
        flag(i) = k
        nc = 1
        wtmp(1:3) = watv(1:3,i)
        do j=1,n
          if (i==j) cycle
          if (flag(j)/=0) cycle
          if (cntcts(j,i)<2*rwi) then
            k=k+1
            cluster(k) = j
            flag(j) = k
            nc = nc + 1
            wtmp(1:3) = wtmp(1:3) + watv(1:3,j)
          endif
        enddo
        x = 0
        do j=1,k-1
          do h=j+1,k
            v = watv(1:3,cluster(j)) - watv(1:3,cluster(h))
            x = x + sum(v*v)
          enddo
        enddo
        if (k>1) then
          x = sqrt(x*2/real(k*(k-1)))
        else
          x = sqrt(x)
        endif
        nv = nv + 1
        wtmp = wtmp/real(nc)
        watv(1:3,nv) = wtmp
        rg(nv) = x
      enddo 
      !! returns nv values for void cluster radius of gyration, rg().
    end subroutine vcluster
!----------------------------------------------------!
!For WZS modeling
integer function getatype(aline)
  implicit none
  character(len=*),intent(in) :: aline
  character(len=3) :: res
  character(len=4) :: atom
  integer,dimension(20,14) :: types
  character(len=3),dimension(20) :: residues
  character(len=4),dimension(20,14) :: atomnames
  integer :: i,j

  !initialize residues
  residues = (/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',&
               'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'/)

  !initialize types
  types = 0
  types(1,1:5) = (/21,8,11,27,13/) 
  types(2,1:6) = (/21,8,11,27,10,29/) 
  types(3,1:8) = (/21,8,11,27,13,1,23,23/)
  types(4,1:9) = (/21,8,11,27,13,13,1,23,23/)
  types(5,1:11) = (/21,8,11,27,13,15,15,15,15,15,15/)
  types(6,1:4) = (/21,8,11,27/)
  types(7,1:10) = (/21,8,11,27,13,2,17,17,2,17/)
  types(8,1:8) = (/21,8,11,27,13,13,13,13/)
  types(9,1:9) = (/21,8,11,27,13,13,13,3,18/)
  types(10,1:8) = (/21,8,11,27,13,13,13,13/)
  types(11,1:8) = (/21,8,11,27,13,10,28,10/)
  types(12,1:8) = (/21,8,11,27,13,6,25,20/)
  types(13,1:7) = (/21,8,11,27,7,7,7/)
  types(14,1:9) = (/21,8,11,27,13,13,6,25,20/)
  types(15,1:11) = (/21,8,11,27,13,13,4,19,4,19,19/)
  types(16,1:6) = (/21,8,11,27,5,24/)
  types(17,1:7) = (/21,8,11,27,5,24,13/)
  types(18,1:7) = (/21,8,11,27,13,13,13/)
  types(19,1:14) = (/21,8,11,27,13,12,12,14,22,14,14,14,14,14/)
  types(20,1:12) = (/21,8,11,27,13,30,30,30,30,30,9,26/)

  !initialize atomnames
  atomnames = "    "
  atomnames(1,1:5) = (/'N   ','CA  ','C   ','O   ','CB  '/)
  atomnames(2,1:6) = (/'N   ','CA  ','C   ','O   ', 'CB  ','SG  '/)
  atomnames(3,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','OD1 ','OD2 '/)
  atomnames(4,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                       'OE1 ','OE2 '/)
  atomnames(5,1:11) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                        'CD2 ','CE1 ','CE2 ','CZ  '/)
  atomnames(6,1:4) = (/'N   ','CA  ','C   ','O   '/)
  atomnames(7,1:10) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','ND1 ',&
                        'CD2 ','CE1 ','NE2 '/)
  atomnames(8,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG1 ','CG2 ','CD1 '/)
  atomnames(9,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                       'CE  ','NZ  '/)
  atomnames(10,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ','CD2 '/)
  atomnames(11,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','SD  ','CE  '/)
  atomnames(12,1:8) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','OD1 ','ND2 '/)
  atomnames(13,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  '/)
  atomnames(14,1:9) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                        'OE1 ','NE2 '/)
  atomnames(15,1:11) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD  ',&
                         'NE  ','CZ  ','NH1 ','NH2 '/)
  atomnames(16,1:6) = (/'N   ','CA  ','C   ','O   ', 'CB  ','OG  '/)
  atomnames(17,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','OG1 ','CG2 '/)
  atomnames(18,1:7) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG1 ','CG2 '/)
  atomnames(19,1:14) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                         'CD2 ','NE1 ','CE2 ','CE3 ','CZ2 ','CZ3 ','CH2 '/)
  atomnames(20,1:12) = (/'N   ','CA  ','C   ','O   ', 'CB  ','CG  ','CD1 ',&
                         'CD2 ','CE1 ','CE2 ','CZ  ','OH  '/)



  read(aline(13:16),*) atom
  read(aline(18:20),*) res
  do i = 1, 21
    if(i == 21) then
      getatype = 0
      exit
    endif
    if(residues(i) == res) exit
  enddo
  do j = 1, 15
    if(j == 15 .or. i == 21) then
      getatype = 0
      exit
    endif
    if(atomnames(i,j) == atom) exit
  enddo
  if(i /= 21 .and. j /= 15) getatype = types(i,j)
end function getatype
  

end program VoidMask
