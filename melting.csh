#!/bin/csh 

setenv DAG $1
setenv LAMBDA 0.2
setenv VOID 0.3
setenv SBREAK 100.
setenv TEMP  300.
setenv TIME 0.000000001
## NOTE: previously THETAM, now SOlidity
setenv SOLID 100.

setenv OUTDIR $DAG
if !(-e $OUTDIR) then
  mkdir $OUTDIR
endif

if ($#argv > 1) then
  @ OMEGA = $2
else
  @ OMEGA = 10
endif

if ($#argv > 2) then
  @ OMAX = $3
else
  @ OMAX = 50
endif

while ($OMEGA < $OMAX )

  ###### UNFOLDING  #########
  setenv LOGFILE $OUTDIR/$DAG.$OMEGA.u.log

  ./xunfoldsim $DAG.dag $SBREAK $TEMP $OMEGA $TIME 0 1. 0 $VOID 1000. $SOLID $LAMBDA > $LOGFILE

  echo -n "UNFOLDING $OMEGA     "
  grep ^"TIMECOURSE" $LOGFILE | tail -1 

  ###### FOLDING  #########
  setenv LOGFILE $OUTDIR/$DAG.$OMEGA.f.log

  ./xunfoldsim $DAG.dag $SBREAK $TEMP $OMEGA $TIME 0 1. 1 $VOID 1000. $SOLID $LAMBDA > $LOGFILE

  echo -n "FOLDING   $OMEGA     "
  grep ^"TIMECOURSE" $LOGFILE | tail -1 

  @ OMEGA ++

end
