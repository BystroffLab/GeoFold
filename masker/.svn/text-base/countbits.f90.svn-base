program countbits
integer(kind=4) :: word4
integer(kind=2) :: word2
integer,parameter :: KND=2
integer(kind=KND) :: ib,jb
integer :: i,j,k,ios=0
integer,parameter :: N=100000
integer,parameter :: M=2**(KND*8-1)
integer,parameter :: ilo=(-M), ihi=(M-1)
integer(kind=1),dimension(ilo:ihi) :: cbits
!!
do i=ilo,ihi
  ib = int(i,KND)
!   write(*,*)
!   write(*,'(i8,"  ",$)') ib
  k = 0
  do j=0,(KND*8-1)
   if (btest(ib,j)) then
!      write(*,'(i1,$)') 1
      k = k + 1
!   else
!      write(*,'(i1,$)') 0
   endif
  enddo
!  write(*,'(i3,$)') k
  cbits(i) = int(k,1)
enddo
!do i=1,N
!  j = btst(word4)
!enddo
!write(*,*) j
!write(*,*) "Starting using cbits, *********************"
!do i=1,N
!  j = usecbits(word4)
!enddo
!write(*,*) j
!write(*,*) "************** writing cbits ***************"
open(13,file="cbits.bin",form="unformatted",status="replace")
write(13) cbits
close(13)
write(*,*) "ALL DONE. Testing cbits. Enter an integer between -32768 and +32767 "
do
  write(*,'("I>> ",$)')
  read(*,*,iostat=ios) word2
  if (ios/=0) exit
  write(*,*) usecbits(word2)
enddo
write(*,*) "Bye."
CONTAINS
  !integer function btst(w)
  !integer(kind=4),intent(in) :: w
  !integer :: i,j
  !j = 0
  !do i=0,31
  !  if (btest(w,i)) j = j + 1
  !enddo
  !!btst = j
  !end function btst
  !! ====================
  integer function usecbits(w)
  integer(kind=KND),intent(in) :: w
  integer :: j
  j = 0
  !do i=1,4/KND
    j = j + cbits(w)
  !enddo
  usecbits = j
  end function usecbits
  !! ====================
end program countbits
