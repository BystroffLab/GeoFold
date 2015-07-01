#!/bin/csh
## Modified: Wed Jan 11 10:33:16 EST 2012
setenv CONVERT /usr/bin/convert
setenv MOLSCRIPT /bach1/home/bystrc/bin/molscript
setenv MOLAUTO /bach1/home/bystrc/bin/molauto
setenv MOL $1
setenv MV /bin/mv
echo "Working directory = "
pwd
if !(-e $MOL.pdb ) then
  echo "$MOL.pdb not found."
  exit
endif
if ($#argv > 1) then
  setenv TODIR $2
else
  setenv TODIR " "
endif

( $MOLAUTO -nice $MOL.pdb > $MOL.mol ) >>& molscript.errors
@ GG = `grep -n "transform atom" $MOL.mol | awk -F : '{print $1}'`
@ HH = $GG + 1
@ SW = 25
@ TLT = 0

@ NN = 0
while ($NN < 360)

  @ MM = 10000 + $NN
  head -$GG $MOL.mol > $MOL.$MM
  echo "  transform atom * by rotation z ${TLT}.;" >> $MOL.$MM
  echo "  transform atom * by rotation y ${SW}.;" >> $MOL.$MM
  echo "  transform atom * by rotation z -${TLT}.;" >> $MOL.$MM
  echo "  transform atom * by rotation y ${NN}.;" >> $MOL.$MM
  tail -n +$HH $MOL.mol | head -n -1 >> $MOL.$MM
  ## echo 'set linewidth 20;' >> $MOL.$MM
  echo 'set bonddistance 2.2;' >> $MOL.$MM
  echo "ball-and-stick require in type CYS and either atom CA, atom CB or atom SG ;" >> $MOL.$MM
  tail -1 $MOL.mol >> $MOL.$MM
  ( $MOLSCRIPT < $MOL.$MM > $MOL.$MM.ps ) >>& molscript.errors
  @ NN += 5
  @ TLT += 5

end
setenv GIF $MOL.gif
$CONVERT -trim -delay 20 -geometry 300x300 $MOL.1????.ps $GIF
## CLEANUP
#rm -f $MOL.1????.ps $MOL.1????
if ($#argv > 1) then
  $MV $GIF $TODIR
endif
