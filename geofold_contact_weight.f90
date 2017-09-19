module geofold_contact_weight
    
    !Input: 
    !Output: nxn array of reals where each entry is the weighted strength of a 
    !pair of residues in contact
    !Functional form: Aexp(\Delta{H})+1
    
contains