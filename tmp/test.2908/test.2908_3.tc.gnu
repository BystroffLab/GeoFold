set terminal postscript enhanced portrait
set size 1.4,0.7
set output "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.tc.ps"
set encoding iso_8859_1
set xlabel "time (sec)"
set logscale x
set ylabel "% of total concentration"
set pointsize 2
set style line 80 lt rgb "#808080"
set style line 81 lt 0
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80
set xtics nomirror
set ytics nomirror
set key outside right
p "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.ftc.plot" w l lt rgb "red" lw 4 title "folded", "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.utc.plot" w l lt rgb "blue" lw 4 title "unfolded", "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.itc.plot" w l lt rgb "purple" lw 4 title "intermediate"
