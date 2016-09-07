program renumber_one
!! Modified. Remove residues that don't have all backbone atoms.
!! Renumber a PDB file starting from 1
implicit none
character(len=1000) ::aline,nline,caline,cline,oline
integer :: I,J,ios,NNN
integer :: K,L,M,N,ires,kres,lres,nback
integer,parameter :: MAXLINE=1000
character(len=1),dimension(MAXLINE) :: aaseq
character(len=80) :: seqfile,linkfile
character(len=20) :: res1='ACDEFGHIKLMNPQRSTVWY'
character(len=60) :: &
res3='ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYR'
integer,dimension(1000) :: seq
character(len=1) :: achar,chainID
character(len=3) :: ares
character(len=80) :: infile,outfile
real :: x,y
integer :: iatm,firstres,jarg
real,parameter :: cfcut=0.50
logical :: complete

integer :: iargc
!! NNN is the number of the first residue
!! seqfile contains the AA sequence, unadulterated

jarg = iargc()
if (jarg < 2) then
  write(*,*) 'Usage xrenumber_one inputfile outputfile '
  write(*,*) 'WARNING: residues with incomplete backbone will be removed.'
  stop 'v. 30-JUL-10'
endif

call getarg(1,infile)
call getarg(2,outfile)

open(11,file=infile,status='old',iostat=ios)
if (ios/=0) stop 'file not found'
!! get true sequence
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5).ne.'ATOM ') cycle
  exit
enddo
write(*,*) trim(aline)
NNN = 1
!read(aline(23:26),*,iostat=ios) firstres
!write(*,'(a,i4,a,$)') "The first residue is ",firstres," New starting number? > "
!read(*,*) NNN
!write(*,'(a,$)') "New chain ID? > "
!read(*,'(a1)') chainID
rewind(11)
iatm = 0
open(12,file=outfile,status='replace',form='formatted',iostat=ios)
if (ios/=0) stop 'Error opening output file. Permissions?'
ires = 0
kres = 0
lres = 0
complete=.false.
nback = 0
do
  read(11,'(a)',iostat=ios) aline
  if (ios/=0) exit
  if (aline(1:5).ne.'ATOM ') then
    write(12,'(a)') trim(aline)
  else
    read(aline(23:26),*,iostat=ios) ires
    iatm = iatm + 1
    write(aline(7:11),'(i5)') iatm
    if (ires/=lres) then
      nback = 0
      complete = .false.
      lres = ires
    endif
    if (.not.complete) then
      select case (aline(13:16))
      case (" N  ")
        nback = nback + 1
        nline = trim(aline)
      case (" CA ")
        nback = nback + 2
        caline = trim(aline)
      case (" C  ")
        nback = nback + 4
        cline = trim(aline)
      case (" O  ")
        nback = nback + 8
        oline = trim(aline)
      end select 
    endif
    if (nback==15) then
      complete = .true.
      nback = 0
      kres = kres + 1
      write(nline(23:26),'(i4)') kres
      write(caline(23:26),'(i4)') kres
      write(cline(23:26),'(i4)') kres
      write(oline(23:26),'(i4)') kres
      write(12,'(a)') trim(nline)
      write(12,'(a)') trim(caline)
      write(12,'(a)') trim(cline)
      write(12,'(a)') trim(oline)
      nback = 0
    elseif (complete) then
      write(aline(23:26),'(i4)') kres
      write(12,'(a)') trim(aline)
    endif
    ! aline(22:22) = chainID
  endif
enddo
write(*,*) 'all done. atoms out=', iatm
close(11)
    
close(12)
end program renumber_one
