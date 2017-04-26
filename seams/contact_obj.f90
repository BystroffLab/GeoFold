module contact_obj
implicit none

private

public :: contact_list,print!,contact_node

type contact_node
  integer :: res1_,res2_,type_flag_
  type(contact_node),pointer :: next_,prev_
  end type contact_node

  type contact_list
    type(contact_node),pointer :: head_,tail_
    integer :: size_
    contains
      ! procedure :: cl_ctor
      procedure :: cl_dtor
      procedure :: get_indx
      procedure :: get_indx_values
      procedure :: append
      procedure :: pop
    end type contact_list
  interface contact_list
    module procedure init_contact_list
  end interface

CONTAINS
!
! type(contact_node) function init_contact_node(res1,res2,type_flag) result(self)
!   integer,intent(in),optional :: res1,res2,type_flag
!   !deallocate and nullify next_ and prev_
!   allocate(self%next_)
!   ! if(associated(self%next_)) deallocate(self%next_)
!   nullify(self%next_)
!   allocate(self%prev_)
!   ! if(associated(self%prev_)) deallocate(self%prev_)
!   nullify(self%prev_)
!   if(present(res1) .and. present(res2) .and. present(type_flag)) then
!     call self%cn_ctor(res1,res2,type_flag)
!   else
!     !set values to 0
!     call self%cn_ctor(0,0,0)
!   endif
! end function init_contact_node

type(contact_list) function init_contact_list() result(self)
  !deallocate and nullify head_ and tail_
  ! if(associated(self%head_)) deallocate(self%head_)
  nullify(self%head_)
  ! if(associated(self%tail_)) deallocate(self%tail_)
  nullify(self%tail_)
  !set size to 0
  self%size_ = 0
end function init_contact_list

!------------------------------------------------------------------------------!
!---------------------------- Contact Node Methods ----------------------------!
!----------------------------------- print ------------------------------------!
subroutine print(self)
  type(contact_node),intent(in) :: self

  write(0,'("res1: ",i4," res2: ",i4," flag: ",a)') self%res1_,self%res2_,get_type_flag_desc(self)
endsubroutine print
! !---------------------------------- set_next ----------------------------------!
! subroutine set_next(self,next)
!   type(contact_node),target,intent(inout) :: self
!   type(contact_node),target,intent(inout) :: next
!
!   self%next_ => next
!   next%prev_ => self
! endsubroutine set_next
! !---------------------------------- set_prev ----------------------------------!
! subroutine set_prev(self,prev)
!   type(contact_node),target,intent(inout) :: self
!   type(contact_node),target,intent(inout) :: prev
!
!   self%prev_ => prev
!   prev%next_ => self
! endsubroutine set_prev
! !---------------------------------- cn_dtor -----------------------------------!
! subroutine cn_dtor(self)
!   type(contact_node),target,intent(inout) :: self
!   type(contact_node),pointer :: self_ptr
!
!   self_ptr => self
!
!   if(associated(self%next_)) then
!     if(associated(self%prev_)) then
!       self%prev_%next_ => self%next_
!       self%next_%prev_ => self%prev_
!     else
!       nullify(self%next_%prev_)
!     endif
!   elseif(associated(self%prev_)) then
!     nullify(self%prev_%next_)
!   endif
!
!   nullify(self%next_)
!   nullify(self%prev_)
!   deallocate(self_ptr)
! endsubroutine cn_dtor
! !---------------------------------- cn_ctor -----------------------------------!
! subroutine cn_ctor(self,res1,res2,type_flag,next,prev)
!   type(contact_node),target,intent(inout) :: self
!   type(contact_node),target,intent(inout),optional :: next,prev
!   integer,intent(in) :: res1,res2,type_flag
!
!   ! self = contact_node()
!   self%res1_ = res1
!   self%res2_ = res2
!   self%type_flag_ = type_flag
!
!   if(present(next)) then
!     self%next_ => next
!     next%prev_ => self
!   endif
!   if(present(prev)) then
!     self%prev_ => prev
!     prev%next_ => self
!   endif
! endsubroutine cn_ctor
! !---------------------------------- get_next ----------------------------------!
! function get_next(self) result(next)
!   type(contact_node),intent(in) :: self
!   type(contact_node),pointer :: next
!
!   if(associated(self%next_)) then
!     next => self%next_
!   else
!     nullify(next)
!   endif
! end function get_next
! !---------------------------------- get_prev ----------------------------------!
! function get_prev(self) result(prev)
!   type(contact_node),intent(in) :: self
!   type(contact_node),pointer :: prev
!
!   if(associated(self%prev_)) then
!     prev => self%prev_
!   else
!     nullify(prev)
!   endif
! end function get_prev
! !---------------------------------- set_res1 ----------------------------------!
! subroutine set_res1(self,res1)
!   type(contact_node),intent(inout) :: self
!   integer,intent(in) :: res1
!
!   self%res1_ = res1
! endsubroutine set_res1
! !---------------------------------- set_res2 ----------------------------------!
! subroutine set_res2(self,res2)
!   type(contact_node),intent(inout) :: self
!   integer,intent(in) :: res2
!
!   self%res2_ = res2
! endsubroutine set_res2
! !------------------------------- set_type_flag --------------------------------!
! subroutine set_type_flag(self,type_flag)
!   type(contact_node),intent(inout) :: self
!   integer,intent(in) :: type_flag
!
!   self%type_flag_ = type_flag
! endsubroutine set_type_flag
!----------------------------- get_type_flag_desc -----------------------------!
function get_type_flag_desc(self) result(flag_desc)
  character(len=:),allocatable :: flag_desc
  type(contact_node),intent(in) :: self

  select case(self%type_flag_)
  case(0)
    flag_desc = "No Flag Set"
  case(1)
    flag_desc = "Disulfide Contact"
  case(2)
    flag_desc = "Backbone Hydrogen Bond"
  case(3)
    flag_desc = "Sidechain Hydrogen Bond"
  case(4)
    flag_desc = "CB-CB Contact"
  case default
    flag_desc = "Invalid Flag Set"
  end select
end function get_type_flag_desc
! !------------------------------- get_type_flag --------------------------------!
! integer function get_type_flag(self) result(flag)
!   type(contact_node),intent(in) :: self
!
!   flag=self%type_flag_
! end function get_type_flag
!------------------------------------------------------------------------------!
!---------------------------- Contact List Methods ----------------------------!
!---------------------------------- cl_dtor -----------------------------------!
subroutine cl_dtor(self)
  class(contact_list),intent(inout),target :: self
  type(contact_node),pointer :: current,next
  class(contact_list),pointer :: self_ptr

  self_ptr => self
  current => self%head_
  do while(associated(current%next_))
    next => current%next_
    if(associated(current)) deallocate(current)
    current => next
  enddo
  if(associated(current)) deallocate(current)
  if(associated(next)) deallocate(next)
  self%size_ = 0

  if(associated(self%head_)) deallocate(self%head_)
  if(associated(self%tail_)) deallocate(self%tail_)
  if(associated(self_ptr)) deallocate(self_ptr)
endsubroutine cl_dtor
!---------------------------------- get_indx ----------------------------------!
function get_indx(self,indx) result(node_ptr)
  class(contact_list),intent(in) :: self
  type(contact_node),pointer :: node_ptr
  integer,intent(in) :: indx
  integer :: i,half

  ! if(associated(node_ptr)) deallocate(node_ptr)

  if(indx > self%size_ .or. indx < 1) then
    write(0,*) "IndexError: Attempted to access contact out of bounds: ",indx
  else
    half = self%size_/2
    if(indx < half) then
      node_ptr => self%head_
      !Let's retain Fortran's 1 indexing for the sake of consistency
      do i = 1, indx
        if(i /= indx) node_ptr => node_ptr%next_
      enddo
    else
      node_ptr => self%tail_
      do i = self%size_,indx,-1
        if(i /= indx) node_ptr => node_ptr%prev_
      enddo
    endif
  endif
end function get_indx
!------------------------------ get_indx_values -------------------------------!
subroutine get_indx_values(self,indx,res1,res2,tflag)
  class(contact_list),intent(in) :: self
  integer,intent(in) :: indx
  integer,intent(out) :: res1,res2,tflag
  type(contact_node),pointer :: c_ptr

  c_ptr => self%get_indx(indx)
  res1 = c_ptr%res1_
  res2 = c_ptr%res2_
  tflag = c_ptr%type_flag_
endsubroutine get_indx_values
!----------------------------------- append -----------------------------------!
subroutine append(self,res1,res2,flag)
  class(contact_list),intent(inout) :: self
  integer,intent(in) :: res1,res2,flag
  type(contact_node),pointer :: node

  if(self%size_ == 0) then
    allocate(self%head_)
    node => self%head_
    nullify(node%prev_)
  else
    allocate(self%tail_%next_)
    node => self%tail_%next_
    node%prev_ => self%tail_
  endif
  node%res1_ = res1
  node%res2_ = res2
  node%type_flag_ = flag
  nullify(node%next_)
  self%tail_ => node
  self%size_ = self%size_ + 1

endsubroutine append
!------------------------------------ pop -------------------------------------!
subroutine pop(self,indx)
  class(contact_list),intent(inout) :: self
  integer,optional,intent(in) :: indx
  type(contact_node),pointer :: to_pop

  if(present(indx)) then
    to_pop => self%get_indx(indx)
  else
    to_pop => self%tail_
    if(associated(to_pop%prev_)) then
      self%tail_ => to_pop%prev_
    else !removing the only item in the list
      nullify(self%tail_)
      nullify(self%head_)
    endif
  endif
  if(.not. associated(to_pop)) then
    write(0,*) "IndexError: Attempting to pop a contact that does not exist: ",indx
  else
    if(indx == 1) then
      if(associated(to_pop%next_)) then
        self%head_ => to_pop%next_
      else !removing only item from list
        nullify(self%head_)
        nullify(self%tail_)
      endif
    elseif(indx == self%size_) then
      if(associated(to_pop%prev_)) then
        self%head_ => to_pop%prev_
      else !removing only item from list
        nullify(self%head_)
        nullify(self%tail_)
      endif
    endif
    deallocate(to_pop)
    self%size_ = self%size_ - 1
  endif
endsubroutine pop
!------------------------------------------------------------------------------!

end module
