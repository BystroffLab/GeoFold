module geofold_contact_weight
    
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
contains