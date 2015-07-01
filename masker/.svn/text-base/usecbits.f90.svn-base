program countbits
integer(kind=4) :: word4
integer,parameter :: KND=2
integer(kind=KND) :: ib,jb
integer :: i,j,k
integer,parameter :: N=100000
integer,parameter :: M=2**(KND*8-1)
integer,parameter :: ilo=(-M), ihi=(M-1)
integer(kind=1),dimension(ilo:ihi) :: cbits
!!
write(*,*) "************** reading cbits ***************"
open(13,file="cbits32.bin",form="unformatted",status="old")
read(13) cbits(ilo:ihi)
close(13)
do
  write(*,'("INT= ",$)')
  read(*,*) ib
  write(*,'("bits=",i5)') cbits(ib)
enddo
end program countbits
