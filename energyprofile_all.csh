#!/bin/csh
setenv GNUHELP /bach1/usr/local/share/gnuplot/4.2/gnuplot.gih
setenv GNUPLOT_PS_DIR /bach1/usr/local/share/gnuplot/4.2/PostScript
setenv GP /bach1/usr/local/bin/gnuplot
setenv GNAME $1
setenv HNAME $2
setenv GNUFILE ${GNAME}.all.gnu
echo "energyprofile_all.csh: Energy landscapes from ${GNAME}*.plot to gnuplot script $GNUFILE"
ls ${GNAME}*.plot
echo 'set terminal postscript portrait' > $GNUFILE
echo 'set size 1.4,0.7' >> $GNUFILE
echo 'set output "'${HNAME}.ps'"' >> $GNUFILE
echo 'set xlabel "unfolded <--------> folded\\nBuried surface area (1000*A^2)" font "TimesBold,18"' >> $GNUFILE
echo 'set ylabel "energy, kJ/mol" font "TimesBold,18" ' >> $GNUFILE
echo 'set pointsize 2' >> $GNUFILE
echo 'set nokey' >> $GNUFILE
@ NN = 0
foreach TARG (`ls ${GNAME}*.plot`)
  @ NN++
  echo 'NRG'$NN' = "'${TARG}'"' >> $GNUFILE
end
@ NN = 0
echo -n 'p ' >> $GNUFILE
foreach TARG (`ls ${GNAME}*.plot`)
  if ($NN > 0) echo -n ', '     >> $GNUFILE
  @ NN++
  echo -n ' NRG'$NN' w l lw 9 ' >> $GNUFILE
end
cat $GNUFILE
if ($NN == 0) then
  echo "WARNING....... plot files missing. NOT RUNNING gnuplot...."
else
  $GP < $GNUFILE
endif

