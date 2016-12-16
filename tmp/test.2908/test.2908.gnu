set terminal postscript enhanced portrait
set size 1.4,0.7
set output "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908.gp.ps"
set encoding iso_8859_1
set xlabel "solvation weight, {/Symbol w}"
set ylabel "unfolding rate, ln(k_u) s^{-1}"
set pointsize 2
set style line 80 lt rgb "#808080"
set style line 81 lt 0
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80
set xtics nomirror
set ytics nomirror
set key right
p "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908.plot" w p pt 7, 15.5124870928 + 0.00000000000*x
