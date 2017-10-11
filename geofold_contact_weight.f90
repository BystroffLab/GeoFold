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
    integer, public, dimension(geofold_nres:geofold_nres) :: gcw_contacts = 0
    integer,public,dimension(geofold_nres) :: gcw_nc = 0
    real, public, dimension(geofold_nres:geofold_nres) :: gcw_contact_weights = BIG
    public :: add_contact, get_contact_weights
contains
    
        
    
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
        
        call geofold_masker_getscenergy(i,j,scentropy,contacts)
        energy = geofold_hbonds_eperbond*geofold_hb(i,j) + sasnrg(i,j) &
            - T*scentropy
        energy = A*exp(energy)+1
        gcw_contact_weights(i,j) = energy
        gcw_contact_weights(j,i) = energy
        call append_contact(i,j)
    end subroutine add_contact
    
    
    
end module geofold_contact_weights