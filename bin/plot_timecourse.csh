#!/bin/csh
setenv LNAME $1
setenv TAG $2
setenv GNUFILE "${BASE}_${TAG}.gnu"
grep ^"TIMECOURSE" ${LNAME} | awk '{print $2, $3, $4, $5}' 
echo 'set terminal postscript portrait' > $GNUFILE
echo 'set size 1.4,0.7' >> $GNUFILE
echo 'set output "'${LNAME}.tc.ps'"' >> $GNUFILE
echo 'set xlabel "time (sec)"' >> $GNUFILE
echo 'set ylabel "% of total concentration"' >> $GNUFILE
echo 'set pointsize 2' >> $GNUFILE
echo 'p [0:'$WATER'] "'${LNAME}.plot'" w p pt 7 , '$BB' + ('$MM')*x' >> $GNUFILE
gnuplot < $GNUFILE

