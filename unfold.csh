#!/bin/csh 

setenv DAG $1
setenv OMEGA $2
setenv LAMBDA $3
setenv VOID $4
setenv SBREAK 100.
setenv TEMP  300.
setenv TIME 0.000000001
setenv THETAM 100.

  ./xunfoldsim $DAG.dag $SBREAK $TEMP $OMEGA $TIME 0 1. 0 $VOID 1000. $THETAM $LAMBDA > $DAG.unfold.log 
