#!/bin/csh
setenv N $1
setenv INCL vb$N.incl

   xmakemask $N vb.pdb
   echo 'integer,parameter :: NVB='$N' ' > $INCL   ## 200="N" in comment above
   echo 'real,dimension(3,NVB) :: vectorball=reshape((/ &' >> $INCL
   awk '{printf "%7.3f,%7.3f,%7.3f,",$6,$7,$8}' vb.pdb | fold | sed -e "s/\(.*\)/   \1 \&/"  >> $INCL
   echo '/),(/3,NVB/))' >> $INCL
   sed -e "s/, &\/)/\/)/" $INCL > junk
   mv junk $INCL
