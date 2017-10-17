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
    
    real, parameter :: BIG = 9999999.0
    real, parameter :: A = 200. !Amplitude for energy function
    integer, public, dimension(:,:),pointer :: gcw_contacts
    integer,public,dimension(:),pointer :: gcw_nc
    real, public, dimension(:,:),pointer :: gcw_contact_weights
    public :: add_contact, get_contact_weights, gcw_init
contains
    
        
    subroutine gcw_init()
        implicit none
        if(.not.associated(gcw_contacts)) &
            allocate(gcw_contacts(geofold_nres,geofold_nres))
        if(.not.associated(gcw_nc)) & 
            allocate(gcw_nc(geofold_nres))
        if(.not.associated(gcw_contact_weights)) &
            allocate(gcw_contact_weights(geofold_nres,geofold_nres))
    end subroutine gcw_init


    subroutine append_contact(i,j)
        implicit none
        integer,intent(in) :: i,j
        
        gcw_nc(i) = gcw_nc(i) + 1
        gcw_nc(j) = gcw_nc(j) + 1
        gcw_contacts(i,gcw_nc(i)) = j
        gcw_contacts(j,gcw_nc(j)) = i
    end subroutine append_contact
            
    
    subroutine gcw_add_contact(i,j,T)
        !add contact between res i and j to contacts_made, should be called after
        !performing shortest path algorithm and adding contact.  Will adjust contact
        !weights as well
        implicit none
        integer,intent(in) :: i,j
        real,intent(in) :: T !Temperature
        real :: energy,scentropy
        
        call geofold_masker_getscenergy(i,j,scentropy,gcw_contacts)
        energy = geofold_hbonds_eperbond*geofold_hb(i,j) + &
            geofold_masker_getsasnrg(i,j) - T*scentropy
        energy = A*exp(energy)+1
        gcw_contact_weights(i,j) = energy
        gcw_contact_weights(j,i) = energy
        call append_contact(i,j)
    end subroutine gcw_add_contact
    
    
    
end module geofold_contact_weight