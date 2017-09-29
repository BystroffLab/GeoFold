module geofold_contact_weight
    use geofold_global
    use geofold_masker
    use geofold_hbonds
    
    private
    
    !Input: parent intermediate, solvation energies, h-bond info
    !Output: nxn array of reals where each entry is the weighted strength of a 
    !pair of residues in contact
    !Functional form: Aexp(\Delta{H})+1
    
    !TODO
    !calculate Delta_H
    !Determine A parameter
    !Do we want to be given f, u1, u2 or just do the whole protein????
    !Only F
    !What terms do we need to include in our Delta_H?
    !H-bond, Solvation
    !Is that it???
    
    !TODO
    !Create initializer subroutines that read hbond info and solvation into arrays
    !that only needs to be called once in GeoFold
    !the rest of this module will use those arrays coupled with a given f from getcutpoints
    !to calculate contact weights
    !NOTE: Hbond info is read into array geofold_hb in geofold_global --> WOOT
    
    !NOTE: arrays to use from geofold_global
    ! geofold_hb : Hbonds
    ! geofold_ss : disulfides
    ! geofold_nres : number of residues
    !!geofold_masker
    ! sasnrg : solvation energy
    
    real, parameter :: HUGE = 9999999.0
    real, parameter :: A = 200. !Amplitude for energy function
    public :: set_contacts_made, add_contact
contains
    subroutine set_contacts_made(contacts_made,f,u1,u2,s)
        implicit none
        logical,dimension(geofold_nres,geofold_nres),intent(inout) :: contacts_made
        type(intermediate),pointer,intent(in) :: f,u1
        type(intermediate),pointer,intent(in),optional :: u2
        type(seam_type),pointer,intent(in),optional :: s
        character,dimension(MAXRES) :: fflag,u1flag,u2flag
        integer :: i,j
        
        
        fflag = f%iflag
        !not a seam move
        if(present(u2)) then
            u1flag = u1%iflag
            u2flag = u2%iflag
        !seam move
        else
            u1flag = s%u1flag
            u2flag = s%u2flag
        endif
        contacts_made = .false.
        
        do i = 1, geofold_nres-3
            do j = i+3, geofold_nres
                if((u1flag(i:i) /= "." .and. u1flag(j:j) /= ".") .or. &
                   (u2flag(i:i) /= "." .and. u2flag(j:j) /= ".") .or.
                   geofold_ss(i,j) == 1 .or. geofold_ss(j,i) == 1) then
                    contacts_made(i,j) = .true.
                    contacts_made(j,i) = .true.
                endif
            enddo
        enddo
        do i = 1, geofold_nres -1
            contacts_made(i,i+1) = .true.
            contacts_made(i+1,i) = .true.
        enddo        
    end subroutine set_contacts_made
    
    
    subroutine add_contact(i,j,contacts_made)
        !add contact between res i and j to contacts_made
        implicit none
        integer,intent(in) :: i,j
        logical,dimension(geofold_nres,geofold_nres),intent(inout) :: contacts_made
        
        contacts_made(i,j) = .true.
        contacts_made(j,i) = .true.
    end subroutine add_contact
    
    !A*exp(energy)+1
    subroutine get_contact_weights(contacts,f,u1,u2,contacts_made)
        implicit none
        real, dimension(geofold_nres,geofold_nres), intent(out) :: contacts
        logical, dimension(geofold_nres,geofold_nres), intent(in) :: contacts_made 
        type(intermediate),pointer,intent(in) :: f,u1
        type(intermediate),pointer,intent(in),optional :: u2
        integer :: i,j
        real :: energy
        
        if(allocated(contacts)) deallocate(contacts)
        allocate(contacts(geofold_nres,geofold_nres))
        contacts = HUGE
        do i = 1, geofold_nres-3
            do j = i+3, geofold_nres
                !check if these contacts have been added
                if(contacts_made(i,j) == .false.) cycle
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