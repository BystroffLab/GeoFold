module contact_obj
implicit none

private

public :: contact_list,print!,contact_node

type contact_node
  integer :: res1_,res2_,type_flag_
  type(contact_node),pointer :: next_,prev_
  end type contact_node

  type contact_list
    private
    type(contact_node),pointer :: head_,tail_
    integer :: size_
    contains
      ! procedure :: cl_ctor
      procedure,public :: cl_dtor
      procedure,public :: get_indx
      procedure,public :: get_indx_values
      procedure,public :: append
      procedure,public :: pop
      procedure,public :: size
    end type contact_list
  interface contact_list
    module procedure init_contact_list
  end interface

CONTAINS

type(contact_list) function init_contact_list() result(self)
  nullify(self%head_)
  nullify(self%tail_)
  self%size_ = 0
end function init_contact_list

!------------------------------------------------------------------------------!
!---------------------------- Contact Node Methods ----------------------------!
!----------------------------------- print ------------------------------------!
subroutine print(self)
  type(contact_node),intent(in) :: self

  write(0,'("res1: ",i4," res2: ",i4," flag: ",a)') self%res1_,self%res2_,get_type_flag_desc(self)
endsubroutine print
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
!---------------------------- Contact List Methods ----------------------------!
!------------------------------------ size ------------------------------------!
integer function size(self)
  class(contact_list),intent(in) :: self
  size = self%size_
end function size
!---------------------------------- cl_dtor -----------------------------------!
subroutine cl_dtor(self)
  class(contact_list),intent(inout),target :: self
  type(contact_node),pointer :: current,next
  class(contact_list),pointer :: self_ptr
  integer :: i

  current => self%head_
  do i = 1,self%size_
    call print(current)
    next => current%next_
    deallocate(current)
    current=>next
  enddo
  self%size_ = 0
  nullify(self%head_)
  nullify(self%tail_)
  nullify(current)
  nullify(next)
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
