#!/bin/csh
setenv GNUHELP /bach1/usr/local/share/gnuplot/4.2/gnuplot.gih
setenv GNUPLOT_PS_DIR /bach1/usr/local/share/gnuplot/4.2/PostScript
setenv GP /bach1/usr/local/bin/gnuplot
setenv GNAME $1
setenv HNAME $2
echo 'set terminal postscript portrait' > ${GNAME}.nrg.gnu
echo 'set size 1.4,0.7' >> ${GNAME}.nrg.gnu
echo 'set output "'${HNAME}.ps'"' >> ${GNAME}.nrg.gnu
echo 'set xlabel "unfolded <--------> folded\\nBuried surface area (1000*A^2)" font "TimesBold,18"' >> ${GNAME}.nrg.gnu
echo 'set ylabel "energy, kJ/mol" font "TimesBold,18" ' >> ${GNAME}.nrg.gnu
echo 'set pointsize 2' >> ${GNAME}.nrg.gnu
echo 'set nokey' >> ${GNAME}.nrg.gnu
grep ^"ISTATE" ${GNAME}.dag.nrg | \
awk '{printf "%9.1f %9.1f\n", $3, $4}' > ${HNAME}.plot
echo 'NRG = "'${HNAME}.plot'"' >> ${GNAME}.nrg.gnu
setenv RNG `cat ${HNAME}.plot | awk 'BEGIN{x=999999.;y=-999999.;}{if ($2<x) {x=$2;}; if ($2>y) {y=$2;}}END{printf "%f",y-x}'`
grep ^"ISTATE" ${GNAME}.dag.nrg | \
   awk -v rg=$RNG 'BEGIN{tic=0.05*rg;}{printf "set label \"%s\" at %f,%f center\n",$5,$3,$4-tic}' >> ${GNAME}.nrg.gnu
grep ^"ISTATE" ${GNAME}.dag.nrg | \
   awk -v rg=$RNG 'BEGIN{tic=0.1*rg;}{printf "set label \"%s\" at %f,%f center\n",$6,$3,$4-tic}' >> ${GNAME}.nrg.gnu
##grep ^"TSTATE" ${GNAME}.dag.nrg | \
##   awk -v rg=$RNG 'BEGIN{tic=0.05*rg;}{printf "set label \"%s\" at %f,%f center\n",$6,$3,$4+tic}' >> ${GNAME}.nrg.gnu
echo 'p NRG w l lw 9 ' >> ${GNAME}.nrg.gnu
$GP < ${GNAME}.nrg.gnu

