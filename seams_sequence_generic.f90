!--Generic Secuence ------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    ! Generic Secuence 
    ! seasm_sequence_generic.f90 -- A Generic Linked SeqGuence Implementation in Fortran 95
    ! Copyright (C) 2013 Luis Garreta
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
module seams_sequence_generic
    ! A public variable used as a MOLD for transfer()
    integer, dimension(:), allocatable :: sequence_data

    ! Linked sequence node
    type :: sequence_node
        private
        integer, dimension(:), pointer :: data => null()
        type(sequence_node), pointer :: next => null()
    end type sequence_node

contains
    !---------------------------------------------------------------------------
    ! Check if the sequence is empty|
    !---------------------------------------------------------------------------
    logical function isEmptySeqG (self) result (value)
        type(sequence_node), pointer :: self

        value = (.not. associated (self))
    endfunction

    !---------------------------------------------------------------------------
    ! Create the sequence
    !---------------------------------------------------------------------------
    subroutine createSeqG (self)
        type(sequence_node), pointer :: self
        nullify  (self)
    endsubroutine

    !---------------------------------------------------------------------------
    ! Append a "data" element to the final
    !---------------------------------------------------------------------------
    subroutine appendSeqG (self, data)
        type(sequence_node), pointer :: self
        integer, intent(in), optional :: data (:)
        type(sequence_node), pointer :: ptr
        type(sequence_node), pointer :: next

        allocate(next)
        nullify(next%next)
        allocate(next%data(size(data)))
        next%data = data

        if (.not. associated (self)) then
            self => next
        else
            ptr => self
            do while (associated (ptr%next))
                ptr => ptr%next
            enddo
            ptr%next => next
         endif
    endsubroutine

    !---------------------------------------------------------------------------
    ! Return the length of the sequence
    !---------------------------------------------------------------------------
    function lengthSeqG (self) result (n)
        type(sequence_node), pointer  :: self
        type(sequence_node), pointer  :: ptr
        integer                     :: n

        n = 0
        ptr => self
        do while (associated (ptr))
            n = n + 1
            ptr => ptr%next
        enddo
    endfunction

    !---------------------------------------------------------------------------
    ! Reverse the sequence
    !---------------------------------------------------------------------------
    subroutine reverseSeqG (self) 
        type(sequence_node), pointer  :: self
        type(sequence_node), pointer  :: selfReversed
        integer, allocatable                 :: data (:)
        integer                     :: n, i

		call createSeqG (selfReversed)
		n = lengthSeqG (self)
		do i=n, 1
			data = getAtSeqG (self, i)
			call appendSeqG (selfReversed, data)
		enddo

		!call freeSeqG (self)
		self = selfReversed
    endsubroutine

    !---------------------------------------------------------------------------
    ! Get the DATA of the sequence at position "pos"
    !---------------------------------------------------------------------------
    function getAtSeqG (self, pos) result (data)
        type(sequence_node), pointer     :: self
        integer, intent (in)           :: pos
        integer, dimension(:), pointer :: data
        type(sequence_node), pointer     :: ptr
        integer                        :: i

        ptr => self

        i=1
        do while (i < pos)
            if (.not. associated (ptr)) STOP "Error: getSeqGAt out of bounds"
            ptr => ptr%next
            i = i + 1
        enddo

        data => ptr%data
    endfunction
     
    !---------------------------------------------------------------------------
    ! Return the DATA stored in the node SELF
    !---------------------------------------------------------------------------
    function getSeqG(self) result(data)
        type(sequence_node), pointer :: self
        integer, dimension(:), pointer :: data
        data => self%data
    end function

    !---------------------------------------------------------------------------
    ! Insert a sequence node after SELF containing DATA (optional)
    !---------------------------------------------------------------------------
    subroutine insertSeqG(self, data)
        type(sequence_node), pointer :: self
        integer, dimension(:), intent(in), optional :: data
        type(sequence_node), pointer :: next

        allocate(next)

        if (present(data)) then
        allocate(next%data(size(data)))
        next%data = data
        else
        nullify(next%data)
        end if

        next%next => self%next
        self%next => next
    end subroutine insertSeqG

    !---------------------------------------------------------------------------
    ! Store the encoded DATA in sequence node SELF
    ! Removed deallocate part as data is not allocatable (Luis Garreta)
    !---------------------------------------------------------------------------
    subroutine putSeqG(self, data)
        type(sequence_node), pointer :: self
        integer, dimension(:), intent(in) :: data

        if (associated(self%data)) then
            deallocate(self%data)
            nullify(self%data)
        endif
        allocate(self%data(size(data)))
        self%data = data
    end subroutine putSeqG

    !---------------------------------------------------------------------------
    ! Get the DATA of the sequence at position "pos"
    !---------------------------------------------------------------------------
    subroutine putAtSeqG (self, pos, data)
        type(sequence_node), pointer        :: self
        integer, intent (in)                :: pos
        integer, intent(in)                 :: data (:)
        integer                             :: i
        type(sequence_node), pointer        :: ptr

        ptr => self

        i=1
        do while (i < pos)
            if (.not. associated (ptr)) STOP "Error: putAtSeqG out of bounds"
            ptr => ptr%next
            i = i + 1
        enddo

        call putSeqG (ptr, data)
    endsubroutine
 
    !---------------------------------------------------------------------------
    ! Return the next node after SELF
    !---------------------------------------------------------------------------
    function nextSeqG(self)
        type(sequence_node), pointer :: self
        type(sequence_node), pointer :: nextSeqG
        nextSeqG => self%next
    end function nextSeqG

    !---------------------------------------------------------------------------
    ! Free the entire sequence and all data, beginning at SELF
    !---------------------------------------------------------------------------
    subroutine printSeqG (self, printFunction)
        implicit none
          type(sequence_node), pointer :: self
          type(sequence_node), pointer :: ptr
          integer, allocatable  :: data (:)
          integer :: lengthData

        interface
            subroutine printFunction (data)
                integer           :: data (:)
            endsubroutine
        endinterface

        if (.not. associated (self)) return

        ptr => self
        do while (associated (ptr))
            lengthData = size(transfer(ptr%data, data))
            allocate(data(lengthData))
            data = transfer(ptr%data, data)
            call printFunction (data)
            deallocate(data)

            ptr => ptr%next
        enddo
    endsubroutine
        
    !---------------------------------------------------------------------------
    ! Free the entire sequence and all data, beginning at SELF
    !---------------------------------------------------------------------------
    subroutine freeSeqG(self)
        type(sequence_node), pointer :: self
        type(sequence_node), pointer :: current
        type(sequence_node), pointer :: next

        current => self
        do while (associated(current))
           next => current%next
           if (associated(current%data)) then
              deallocate(current%data)
              nullify(current%data)
           end if
           deallocate(current)
           nullify(current)
           current => next
        end do
     endsubroutine
endmodule
