#!/bin/csh

setenv NRES 1560
echo -n "Number of residue [$NRES]: "
setenv NNN $<
if ($NNN == "") then
else
  setenv NRES $NNN
endif
echo -n "Verbose? [y]/n : "
setenv NNN $<
if ($NNN == "n") then
  setenv TF false
else
  setenv TF true
endif

sed -e "s/MAXRES=.*/MAXRES=$NRES/" geofold_global.f90.bck | sed -e "s/verbose=.*/verbose=.$TF./" > geofold_global.f90

setenv NVB `grep "NVB=" vb.incl | head -1 | sed -e "s/.*NVB=\(....\).*/\1/"`
echo -n "Change the pivot angle resolution? [ $NVB ] ? "
setenv NNVB $<
if ($NNVB == "") exit
if ($NNVB == "n") exit
if ($NNVB == "N") exit
./vectorball.csh $NNVB
mv vb$NNVB.incl vb.incl

