! Modified based on xeijcij2ps.f90,cij2ps.f90,eij2ps.f90 
! to read delta-G, convert to probability
! and output in colorscale.
!  X. Yuan  Wed Jul  7 13:26:33 EDT 2004 
!=========================================================
! PATHWAY2PS was generated from EIJCIJ2PS
! This program generates an "age plot" for a folding/unfolding
! pathway. An Age Plot is a square symmetric matrix
! where each value A(i,j) is the relative age of the i,j contact.
! Ages are converted to colors for output as an image.
! Output is a postscript file 
!=========================================================
! LAST MODIFIED: Tue Mar  3 18:28:36 EST 2009
!=========================================================

program pathway2ps
implicit none
integer :: I,J,t,q,N1,N2,NQ,ounit=12,L,k,M,ii,jj,nc,nc2
real,dimension(:,:),allocatable :: eij
integer,dimension(:,:),allocatable :: cij
real,dimension(:),allocatable :: dummy
character(len=200) :: cijfile,seqfile,psfile,str,eijfile
character(len=800) :: seq1,seq2
real :: x,idscore=0.1,xmax,xmin,scale,e,rt,y
integer :: jarg,ios=0,N,nres
integer :: A(1000),B(1000)
character(len=2850) :: aline
integer :: getfields, iargc, istart
integer :: ilocal=4,im,jm
logical :: incolor=.true.
real :: maxx,minx,prev

seq1=' '
seq2=' '
e = exp(1.)
rt = 0.1
istart = 1

jarg=iargc()
if (jarg < 4) then
  write(*,'("Usage: xpathway2ps seqfile agefile cijfile psfile [ilocal RT start B/W] ")')
  write(*,*) 'seqfile(in) has one-let aa'
  write(*,*) 'agefile(in) contains i j time(1/2) (from unfoldsim.f90)'
  write(*,*) 'cijfile(in) contains i j 1.0'
  write(*,*) 'psfile(out) '
  write(*,*) 'Default ilocal =',ilocal
  write(*,*) 'Default rt =',rt
  write(*,*) 'Default start =',istart
  stop 'pathway2ps.f90 v.  Tue Mar  3 18:30:37 EST 2009'
endif
!write(*,*) "===================================================="
!write(*,*) 'eijcij2ps.f90 v.1  Modified by X. Yuan. July-05-2004'
!write(*,*) '              v.2  C.Bystroff Tue Aug 10 10:42:21 EDT 2004'
!write(*,*) "===================================================="

call getarg(1,seqfile)
call getarg(2,eijfile)
call getarg(3,cijfile)
call getarg(4,psfile)
if (jarg >= 5) then
  call getarg(5,str)
  read(str,*,iostat=ios) ilocal
  if (ios/=0) stop 'Bad value for ilocal'
endif
if (jarg >= 6) then
  call getarg(6,str)
  read(str,*,iostat=ios) rt
  if (ios/=0) stop 'Bad value for rt'
endif
if (jarg >= 7) then
  call getarg(9,str)
  read(str,*,iostat=ios) istart
  if (ios/=0) stop 'Bad value for start'
endif
 if (jarg >= 8) then
   incolor = .false.
 endif

!!-------------------- read sequence ----------------------
open(12,file=seqfile,status='old',form='formatted',iostat=ios)
nres = 0
if (ios==0) then
  read(12,'(a)',iostat=ios) aline
  if (ios/=0) stop "Error reading first line of sequence file"
  N = 0
  do 
    if (aline(1:1)/=">".and.aline(1:1)/='#') then
      do i=1,index(aline,'       ')
        if (aline(i:i) >= 'A' .and. aline(i:i) <= 'Y') then
          N = N + 1
          seq1(N:N) = aline(i:i)
        endif
      enddo
    endif
    read(12,'(a)',iostat=ios) aline
    if (ios/=0) exit
  enddo
  nres = N
  write(*,'("Sequence length = ",i6)') nres
  close(12)
  seq2 = seq1
else
  write(0,'("Sequence file is missing: ",a)') trim(seqfile)
  write(0,'("Getting length from Age file: ",a)') trim(eijfile)
endif
         

!!-------------------- read Eij if necessary to get length----------------------
write(*,'(a)') "Opening Age file: ",trim(eijfile)
open(13,file=eijfile,status='old',form='formatted',iostat=ios)
if (ios/=0) stop 'pathway2ps.f90:: ERROR agefile missing!'
if (nres == 0) then
  do 
    read(13,'(a)',iostat=ios) aline
    if (ios/=0) exit
    read(aline,*,iostat=ios) i,j
    if (ios/=0) then
      write(*,'("Error reading line:",/,a)') trim(aline)
      stop 'Error reading line arg 1  0001'
    endif
    if (i > nres) then
      nres = i
    endif
    if (j > nres) then
      nres = j
    endif
  enddo
  rewind(13)
endif

!!-------------------- read Eij ----------------------
allocate(eij(nres,nres),stat=ios)
if (ios/=0) stop 'Error allocating eij'

eij = 0.0
maxx = -999
minx = 999
do 
  read(13,'(a)',iostat=ios) aline
  if (ios/=0) exit
  read(aline,*,iostat=ios) i,j,x
  if (ios/=0) then
    write(*,'("Error reading line:",/,a)') trim(aline)
    stop 'Error reading line arg 1 0002'
  endif
  if ((i > nres).or.(j > nres)) then
    write(*,*) "Out of range",i,j
    stop 'residue number too large 0003'
  endif
  if (abs(i-j) < ilocal) cycle
  eij(i,j)=x
  eij(j,i)=x
  if (x<0.) stop 'Negative concentration found. Is this a pathway (.path) file???'
  if (x>maxx) maxx = x
  if (x<minx) minx = x
enddo
close(13)  
write(*,*) "Done reading pathway."
!!----------------- read true contacts -----------------------
!!---------------- NOTE: for casp6 the .cij file is in RR format ---------

allocate(cij(nres,nres),stat=ios)
if (ios/=0) stop 'Error allocating cij'
open(13,file=cijfile,status='old',form='formatted',iostat=ios)
if (ios/=0) stop 'pathway2ps:: ERROR. cijfile missing!'

nc = 0
cij = 0
do 
  read(13,'(a)',iostat=ios) aline
  if (ios/=0) exit
  !! L = getfields(10,aline,A,B)
  read(aline,*,iostat=ios) i,j
  if (ios/=0) then
    write(*,'("Error reading line:",/,a)') trim(aline)
    stop 'Error reading line arg 1 0004'
  endif
  if ((i > nres).or.(j > nres)) then
    write(*,*) i,j
    stop 'Cij file residue number too large'
  endif
  if (abs(i-j)<ilocal) cycle
  !!=== Plot all listed i,j
  cij(i,j) = 1
  cij(j,i) = 1
  nc = nc + 1
enddo
close(13)
write(*,*) "Done reading true contacts. nc=",nc

!! Select just the contact pairs, then sort them by concentration (eij).
!! Low concentration will be colored red, high concentration blue.
!! Values between 0 and 1 map to HSB values (0.67,1,1) to (0., 1.,1.), 
!! blue to red through green.
do i=1,nres
  do j=i,min(i+ilocal-1,nres)
    eij(i,j) = 0.00
    eij(j,i) = 0.00
  enddo
enddo
!! assign eij to temporary values x = (hue - 1)/0.67
x = -100.   !! don't change this.
y = 1.
M = 0
prev = 0.
nc2 = getnc2(eij,nres)
do while (y > 0.000000001)
  y = 0.0
  !! get highest value of eij
  do i=1,nres
    do j=i+ilocal,nres
      if (cij(i,j)==0) cycle
      if (eij(i,j) > y) then
        y = eij(i,j)
        im = i
        jm = j
      endif
    enddo
  enddo
  !! replace the eij value with a rank
  !!
  !write(0,*) eij(im,jm)
  !write(0,*) prev 
  if(eij(im,jm) /= prev) x = x + (1/real(nc2))
  !x = x + (1/real(nc))
  prev = eij(im,jm)
  eij(im,jm) = x
  !write (0,*) x
  eij(jm,im) = x
  M = M + 1
  ! write(0,*) M,y,x
  if (M > nc) exit
enddo
eij = eij + 100.0
write(*,*) "Done assigning hues.", M
!! now contacts have sorted values from 0 (highest conc) to 1 (lowest conc)
!! These map to blue (highest conc) and red (lowest conc)

!! resign the lower right triangal eij with cij
!! which can be another way to draw the plot with whole eij
!do jj=1,nres
!  do ii=jj+ilocal,nres
!    if (cij(ii,jj)==-1) then !!never goes here
!      eij(ii,jj) = eij(jj,ii)
!    elseif (cij(ii,jj)==1) then !!cij=1
!      eij(ii,jj) = -1.
!    else !!cij=0
!      eij(ii,jj) = 0.
!    endif
!  enddo
!enddo
!! write(*,*) "Done setting colors."

call WRITEPS(psfile,eij,nres,nres)

CONTAINS
!---------------------------------------------------------------!
subroutine WRITEPS(filename,eij,N,M)
! Max values set to black, min set to white
implicit none
character(len=80),intent(in) :: filename
integer,intent(in) :: N,M
real,dimension(N,M),intent(in) :: eij
integer :: I,J,K
real,parameter :: cm=30 ! points/cm
real,parameter :: page=550 ! points/page
real,parameter :: margin=1. ! cm margin
character(len=10) :: phobic="ACFILMPVWY"
character(len=6) :: polar="GHNQST"
real :: boxsize,x,y,g,ticksize,xx,yy,boxline,boxlinegray
character(len=5) :: ticklabel

open(22,file=filename,status="replace",form="formatted",iostat=ios)
if (ios/=0) stop "Error opening psfile"
ticksize = 0.3
boxsize = (((page-ticksize*cm)/max(N,M))/cm)
boxline = boxsize/5.0
boxlinegray = 0.0
write(22,'("%!PS")')
write(22,'("% -- Created by eijcij2ps.f90  C.Bystroff, X.Yuan  Aug, 2004  ---- %")')
write(22,'("% ---- define box procedure")')
write(22,'("/cm { ",i4," mul } def")') int(cm)
!!
write(22,'("/uptick")')
write(22,'("{ newpath moveto 0 ",f8.3," cm rlineto closepath stroke } def ")') ticksize
!!
write(22,'("/rtick")')
write(22,'("{ newpath moveto ",f8.3," cm 0 rlineto closepath stroke } def ")') ticksize
!!
write(22,'("% ---- outline procedures")')
write(22,'("/boxbottomline")')
write(22,'("{ newpath setgray moveto ",f8.3," cm 0 rlineto ")') boxsize
write(22,'("  ",f8.3," cm setlinewidth closepath stroke } def")') boxline
!!
write(22,'("/boxleftline")')
write(22,'("{ newpath setgray moveto 0 ",f8.3," cm rlineto ")') boxsize
write(22,'("  ",f8.3," cm setlinewidth closepath stroke } def")') boxline
!!
write(22,'("/boxrightline")')
write(22,'("{ newpath setgray moveto ",f8.3," cm 0 rmoveto ")') boxsize
write(22,'("  0 ",f8.3," cm rlineto ")') boxsize
write(22,'("  ",f8.3," cm setlinewidth closepath stroke } def")') boxline
!!
write(22,'("/boxtopline")')
write(22,'("{ newpath setgray moveto 0 ",f8.3," cm rmoveto ")') boxsize
write(22,'("  ",f8.3," cm 0 rlineto ")') boxsize
write(22,'("  ",f8.3," cm setlinewidth closepath stroke } def")') boxline
!!
write(22,'("/box")')
write(22,'("{ ",f8.3," cm 0 rlineto")') boxsize
write(22,'("  0 ",f8.3," cm rlineto")') boxsize
write(22,'("  ",f8.3," cm 0 rlineto")') -boxsize
write(22,'("  closepath } def")')
!!
write(22,'("/graybox")')
write(22,'("{ newpath setgray moveto box fill } def")')
!!
write(22,'("/hsbbox")')
write(22,'("{ newpath moveto box sethsbcolor fill } def")')
!!
write(22,'("/Times_Roman findfont")')
write(22,'(f8.3," cm scalefont setfont")') ticksize
write(22,'("% -- begin eijcij2ps program ---- %")')
!  draw sequence 1 horizontal
! at the top
y = margin + N*boxsize +  boxsize
yy = margin + N*boxsize +  boxsize + ticksize
!! top margin ticks
do I=1,N
  J = I + istart - 1
  if (mod(J,10)==0) then
    write(ticklabel,'(i5)') J
    x = margin + (I-1)*boxsize 
    write(22,'(f8.2," cm ",f8.2," cm  moveto ",$)') x,yy
    write(22,'("(",a5,") show")') adjustl(ticklabel)
    x = margin + (I-1)*boxsize
    write(22,'(f8.2," cm ",f8.2," cm  uptick")') x,y
  endif
enddo
x = margin + N*boxsize +  boxsize
xx = margin + N*boxsize + boxsize + ticksize
do I=1,N
  J = I + istart - 1
  if (mod(J,10)==0) then
    write(ticklabel,'(i5)') J
    y = margin + (I-1)*boxsize 
    write(22,'(f8.2," cm ",f8.2," cm  moveto ",$)') xx,y
    write(22,'("(",a5,") show")') adjustl(ticklabel)
    y = margin + (I-1)*boxsize
    write(22,'(f8.2," cm ",f8.2," cm  rtick")') x,y
  endif
enddo
! draw upper triangle
  x = margin
  y = margin
  write(22,'("  newpath ")')
  write(22,'(f8.2," cm ",f8.2," cm  moveto")') x,y
  x = margin + N*boxsize
  y = margin + N*boxsize
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  x = margin 
  y = margin + N*boxsize
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  x = margin 
  y = margin 
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  write(22,'("  closepath stroke ")')
!
! draw lower triangle
  x = margin
  y = margin
  write(22,'("  newpath ")')
  write(22,'(f8.2," cm ",f8.2," cm  moveto")') x,y
  x = margin + N*boxsize
  y = margin + N*boxsize
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  x = margin + N*boxsize
  y = margin 
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  x = margin
  y = margin
  write(22,'(f8.2," cm ",f8.2," cm  lineto")') x,y
  write(22,'("  closepath stroke ")')
!!
!! draw the input .eij in the upper left triangle
do J=1,M
  ! y in units of cm
  y = margin + (J-1)*boxsize
  do I=1,J !! one way to draw only the upper  lefttriangle
    if (abs(I-J) < ilocal) cycle
    x = margin + (I-1)*boxsize
    !if (I < J) then !! another way to draw upper triangle
    if (cij(I,J)==0) cycle
      write(0,*) eij(I,J)
      if (eij(I,J) < 0.) then  !! draw a black box
        write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.
      elseif (eij(I,J) <= 1.000 ) then
        g = 0.670*(1.00 - eij(I,J))   !! blue to red. Low eij is blue.
        if (incolor) then
          write(22,'(3f6.3,f8.2," cm ",f8.2," cm  hsbbox")') g,1.0,1.0,x,y
        else
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,(1.-eij(I,J))
        endif
      endif
    !else
      !g = 0.670*(1.00 - eij(I,J))   !! blue to red.
      !if (eij(I,J) < 0.998) then
      !  if (incolor) then
      !    write(22,'(3f6.3,f8.2," cm ",f8.2," cm  hsbbox")') g,1.0,1.0,x,y
      !  else
      !    write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,(1.-eij(I,J))
      !  endif
      !endif
    !endif
  enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!
!! draw the input .cij in the lower right triangle
do J=1,M
  ! y in units of cm
  y = margin + (J-1)*boxsize
  do I=1,M
    x = margin + (I-1)*boxsize
    if (I > J) then !! draw lower right triangle
      g = 1.000 - cij(I,J)
      ! experimental: Add score for identities
      ! if (seq1(I:I)==seq2(J:J).and.g>=idscore) g = g - idscore
      if (abs(i-j) < ilocal) then
        !if (g < 0.997) then
        !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.4
        !else
        !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.8
        !endif
      elseif (g < 0.997) then
        if ((I-J) == 3) then
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.5
        else
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,g
        endif
      endif
    else  !! draw lines around contacts in upper right triangle
      if (cij(I,J)==1) then
        if (I==1) then !! leftline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxleftline")') x,y,boxlinegray
        elseif (cij(I-1,J)==0) then !! leftline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxleftline")') x,y,boxlinegray
        endif
        if (J==1) then !! bottomline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxbottomline")') x,y,boxlinegray
        elseif (cij(I,J-1)==0) then !! bottomline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxbottomline")') x,y,boxlinegray
        endif
        if (I==M) then !! rightline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxrightline")') x,y,boxlinegray
        elseif (cij(I+1,J)==0) then !! rightline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxrightline")') x,y,boxlinegray
        endif
        if (J==M) then !! topline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxtopline")') x,y,boxlinegray
        elseif (cij(I,J+1)==0) then !! topline
          write(22,'(f8.2," cm ",f8.2," cm ",f6.3," boxtopline")') x,y,boxlinegray
        endif
      endif
    endif
  enddo
enddo

!! Draw sequence along diagonal, if present
!! This chan be the chain ID characters
write(22,'("/Times_Roman findfont")')
write(22,'(f8.3," cm scalefont setfont")') boxsize
if (seq1 /= " ") then
  do I=1,M
    x = margin + (I-1)*boxsize
    y = x + boxsize
    !y = x - boxsize
    write(22,'(f8.2," cm ",f8.2," cm  moveto")') x,y
    write(22,'("(",a1,") show")') seq1(I:I)
    !! don't draw bars !!------
    !if (seq1(I:I) == "G") then
    !  y = x + boxsize
    !  write(22,'(3f6.3,f8.2," cm ",f8.2," cm  hsbbox")') 0.05,0.9,1.0,x,y
    !  y = x 
    !  write(22,'(3f6.3,f8.2," cm ",f8.2," cm  hsbbox")') 0.05,0.9,1.0,x,y
    !  y = x - boxsize
    !  write(22,'(3f6.3,f8.2," cm ",f8.2," cm  hsbbox")') 0.05,0.9,1.0,x,y
    !elseif (index(phobic,seq1(I:I)) /= 0) then
    !  y = x + boxsize
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.0
    !  y = x 
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.0
    !  y = x - boxsize
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.0
    !elseif (index(polar,seq1(I:I)) /= 0) then
    !  y = x + boxsize
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.65
    !  y = x 
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.65
    !  y = x - boxsize
    !  write(22,'(f8.2," cm ",f8.2," cm ",f6.3," graybox")') x,y,0.65
    !endif
  enddo
endif

!!----- draw legend

write(22,'("showpage")')
close(22)
end subroutine WRITEPS
!---------------------------------------------------------------!

integer function getnc2(eij,n)
  implicit none
  real,dimension(n,n),intent(in) :: eij
  integer,intent(in) :: n
  real,dimension(n) :: values
  integer :: i,j,k
  
  values = 0.
  getnc2 = 0
  do  i = 1, n
    loop: do j = i+ilocal, n
      do k = 1, n
        if(values(k) == 0.) then
          values(k) = eij(i,j)
          cycle loop
        endif
        if(eij(i,j) == values(k)) cycle loop
      enddo
    enddo loop
  enddo
  do k = 1, n
    if(values(k) == 0.) exit
    getnc2 = getnc2 + 1
  enddo
        
  write(0,'("getnc2 = ",i4)')getnc2
end function getnc2

end program pathway2ps
