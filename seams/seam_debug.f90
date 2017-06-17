module seam_debug
  use seams_pdbtools
  use seams_sequences
  use seams_graph
  use seams_utils

CONTAINS
  !----------------------------------------------------------------------------
  ! writeMatrix
  subroutine dwriteMatrix(contactMatrix,nres,filename)
    implicit none
    character(len=*),intent(in) :: filename
    integer :: dunit
    integer, intent(in) :: nres
    integer,intent(in) :: contactMatrix(:,:)
    integer :: i,j,ios
    character(len=1),parameter :: tab = char(9)

    open(newunit=dunit,file=filename,form="formatted",status="replace",iostat=ios)
    if(ios/=0) stop "Couldn't open file to writeMatrix"

    do i = 1, nres
      do j = 1, nres
        if(contactMatrix(i,j)==0) cycle
        write(dunit,*) i,tab,j,tab,"1.000"
      enddo
    enddo
    close(dunit)
  endsubroutine dwriteMatrix
  !----------------------------------------------------------------------------
  subroutine itoa(num,ch)
    implicit none
    character(len=:),allocatable,intent(out) :: ch
    character(len=6000) :: tmpch
    integer,intent(in) :: num

    if(allocated(ch)) deallocate(ch)

    write(tmpch,'(i4.4)') num
    ch = trim(adjustl(tmpch))
  end subroutine itoa
end module seam_debug
