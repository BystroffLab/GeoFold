#!/bin/csh
if ($#argv < 2) then
  echo "creategnuplot.csh  <uniquename> <intercept>"
  echo "The script expects <uniquename>.fit and writes <uniquename>.gp.ps, a postscript image."
  exit
endif
setenv GNUHELP /bach1/usr/local/share/gnuplot/4.2/gnuplot.gih
setenv GNUPLOT_PS_DIR /bach1/usr/local/share/gnuplot/4.2/PostScript
setenv GP /bach1/usr/local/bin/gnuplot
setenv LNAME $1
setenv WATER $2
setenv BB `grep -A 3 "Least-squares coefficients" ${LNAME}.fit | tail -2 | awk '{print $2}' | head -1`
setenv MM `grep -A 3 "Least-squares coefficients" ${LNAME}.fit | tail -2 | awk '{print $2}' | tail -1`
echo 'set terminal postscript portrait' > ${LNAME}.gnu
echo 'set size 1.4,0.7' >> ${LNAME}.gnu
echo 'set output "'${LNAME}.gp.ps'"' >> ${LNAME}.gnu
echo 'set xlabel "desolvation energy, omega (kJ/mol/A^2)"' >> ${LNAME}.gnu
echo 'set ylabel "unfolding rate, ln(k_u) s^-1"' >> ${LNAME}.gnu
echo 'set pointsize 2' >> ${LNAME}.gnu
echo 'LOG_KU = "'${LNAME}.plot'"' >> ${LNAME}.gnu
## echo 'p [0:'$WATER'] "'${LNAME}.plot'" w p pt 7 , '$BB' + ('$MM')*x' >> ${LNAME}.gnu
echo 'p [0:'$WATER'] LOG_KU w p pt 7 , '$BB' + ('$MM')*x' >> ${LNAME}.gnu
$GP < ${LNAME}.gnu
echo "Data from ${LNAME}.plot written to gnuplot script ${LNAME}.gnu"

