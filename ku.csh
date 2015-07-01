#!/bin/csh
setenv PDBCODE $1
setenv CHAIN $2
setenv TARG $PDBCODE$CHAIN
setenv PARFILE $TARG.par

foreach OMEGA ( 5.  8.  11.  14.  )
     sed -e "s/OMEGA.*/OMEGA  "$OMEGA"/" parameters > $PARFILE
     ./xunfoldsim $TARG.dag $PARFILE > $TARG.tmp
     echo -n "$OMEGA   "
     grep ^"TIMECOURSE" $TARG.tmp | tail -1 
end


