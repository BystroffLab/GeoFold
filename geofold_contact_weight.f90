module geofold_contact_weight
    use geofold_global
    use geofold_masker
    use geofold_hbonds
    
    private
    
    !Input: parent and child intermediates, logical array of contacts formed
    !Output: nxn array of reals where each entry is the weighted strength of a 
    !pair of residues in contact
    !Functional form: Aexp(\Delta{H})+1
    
    
    
    !NOTE: arrays to use from geofold_global
    ! geofold_hb : Hbonds
    ! geofold_ss : disulfides
    ! geofold_nres : number of residues
    !!geofold_masker
    ! sasnrg : solvation energy
    
    real, parameter :: HUGE = 9999999.0
    real, parameter :: A = 200. !Amplitude for energy function
    public :: add_contact, get_contact_weights
contains
    
    logical contacts_made(i,j,f,u1,u2,s) result(made)
        implicit none
        logical,dimension(geofold_nres,geofold_nres),intent(inout) :: contacts_made
        type(intermediate),pointer,intent(in) :: f,u1
        type(intermediate),pointer,intent(in),optional :: u2
        type(seam_type),pointer,intent(in),optional :: s
        character,dimension(MAXRES) :: fflag,u1flag,u2flag
        integer,intent(in) :: i,j
        
        
        fflag = f%iflag
        !not a seam move
        if(present(u2)) then
            u1flag = u1%iflag
            u2flag = u2%iflag
        !seam move
        else
            u1flag = s%u1flag
            u2flag = s%u2flag
                if((u1flag(i:i) /= "." .and. u1flag(j:j) /= ".") .or. &
                   (u2flag(i:i) /= "." .and. u2flag(j:j) /= ".") .or.
                   geofold_ss(i,j) == 1 .or. geofold_ss(j,i) == 1) then
                    made = .true.
                endif
            enddo
        enddo
        if(abs(i-j)==1) made = .true.
    end function contacts_made
        
    
    
    subroutine add_contact(i,j,contacts_made,contacts)
        !add contact between res i and j to contacts_made, should be called after
        !performing shortest path algorithm and adding contact.  Will adjust contact
        !weights as well
        implicit none
        integer,intent(in) :: i,j
        real,dimension(geofold_nres,geofold_nres),intent(inout) :: contacts
        real :: energy
        
        energy = geofold_hbonds_eperbond*geofold_hb(i,j) + sasnrg(i,j)
        energy = A*exp(energy)+1
        contacts(i,j) = energy
        contacts(j,i) = energy
    end subroutine add_contact
    
    !A*exp(energy)+1
    subroutine get_contact_weights(contacts,f,u1,u2,s,contacts_made)
        implicit none
        real, dimension(geofold_nres,geofold_nres), intent(inout) :: contacts 
        type(intermediate),pointer,intent(in) :: f,u1
        type(intermediate),pointer,intent(in),optional :: u2
        integer :: i,j
        real :: energy
        type(seam_type),pointer,intent(in),optional :: s
        
        contacts = HUGE
        do i = 1, geofold_nres-3
            do j = i+3, geofold_nres
                !check if these contacts have been added
                if(contacts_made(i,j,f,u1,u2,s) == .false.) cycle
                if(geofold_ss(i,j) == 1 or geofold_ss(j,i) == 1) then
                    contacts(i,j) = 1.
                    contacts(j,i) = 1.
                    cycle
                endif
                energy = geofold_hbonds_eperbond*(geofold_hb(i,j)) + sasnrg(i,j)
                energy = A*exp(energy)+1
            enddo
        enddo
        do i = 1,geofold_nres-1
            contacts(i,i+1) = 1.
            contacts(i+1,i) = 1.
        enddo
    end subroutine get_contact_weights
end module geofold_contact_weights