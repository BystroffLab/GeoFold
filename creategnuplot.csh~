#!/bin/csh
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
echo 'p [0:'$WATER'] "'${LNAME}.plot'" w p pt 7 , '$BB' + ('$MM')*x' >> ${LNAME}.gnu
gnuplot < ${LNAME}.gnu

