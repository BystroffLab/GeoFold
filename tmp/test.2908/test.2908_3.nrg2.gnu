set term postscript enhanced portrait
set size 1.4,0.7
set title "test.2908_3"
set output "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.nrg.ps"
set xlabel "unfolded <---> folded\nBuried surface area (x1000A ^{2})"
set ylabel "energy (kJ/mol)"
set pointsize 2
unset key
set datafile separator ","
set style line 80 lt rgb "#808080"
set style line 81 lt 0
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80
set xtics nomirror
set ytics nomirror
p "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "orange" pt 0