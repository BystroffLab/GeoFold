#!/bin/csh
setenv GNUHELP /bach1/usr/local/share/gnuplot/4.2/gnuplot.gih
setenv GNUPLOT_PS_DIR /bach1/usr/local/share/gnuplot/4.2/PostScript
setenv GP /bach1/usr/local/bin/gnuplot
setenv INTM 0
setenv LOGF $1
if ($#argv == 0) then
  echo "Usage: createtrace.csh <unfoldsim log file w/TIMECOURSE lines> [<intermediate>]"
  exit
endif
if ($#argv > 1) then
  setenv INTM $2
endif
setenv LNAME `basename $LOGF .log`
setenv DNAME `dirname $LOGF`
setenv PSFILE $DNAME/${LNAME}.ps
setenv GNUFILE $DNAME/${LNAME}.gnu
echo '' > ${GNUFILE}
echo 'set terminal postscript portrait' >>  ${GNUFILE}
echo 'set size 1.4,0.7' >>  ${GNUFILE}
echo 'set output "'${PSFILE}'"' >>  ${GNUFILE}
echo 'set xlabel "time (sec)"' >>  ${GNUFILE}
echo 'set ylabel "rel. conc. (% total)"' >>  ${GNUFILE}
#echo 'set style line 3 lt 1 lw 4' >>  ${GNUFILE}
#echo 'set style line 2 lt 2 lw 4' >>  ${GNUFILE}
#echo 'set style line 1 lt 3 lw 4' >>  ${GNUFILE}
echo 'F = "< grep ^TIME '$LOGF' | awk '"'"'{print $2, $3}'"'"' " '  >>  ${GNUFILE}
echo 'U = "< grep ^TIME '$LOGF' | awk '"'"'{print $2, $4}'"'"' " '  >>  ${GNUFILE}
echo 'I = "< grep ^TIME '$LOGF' | awk '"'"'{print $2, $5}'"'"' " '  >>  ${GNUFILE}
if ($INTM > 0) then
  echo 'L = "< grep ^TIME '$LOGF' | awk '"'"'{print $2, $8}'"'"' " '  >>  ${GNUFILE}
  echo 'p [][0:100] F w line lt 1 lw 16, U w line lt 2 lw 16, I w line lt 3 lw 10, L w line lt 4 lw 16 '  >>  ${GNUFILE}
else
  echo 'p [][0:100] F w line lt 1 lw 16, U w line lt 2 lw 16, I w line lt 3 lw 12 '  >>  ${GNUFILE}
endif
$GP <  ${GNUFILE}

