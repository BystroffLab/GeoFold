set terminal postscript enhanced portrait
set size 1.4,0.7
set output "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_all.nrg.ps"
set encoding iso_8859_1
set xlabel "unfolded {/Symbol «} folded\nBuried surface area ({/Symbol ´}1000{Å}^{2})" font "TimesBold,18"
set ylabel "energy (kJ/mol)" font "TimesBold,18"
set pointsize 2
set key outside right
set style line 80 lt rgb "#808080"
set style line 81 lt 0
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80
set xtics nomirror
set ytics nomirror
p "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_1.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "black" pt 0 title "{/Symbol w} = .04", "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_2.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "red" pt 0 title "{/Symbol w} = .05", "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_3.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "orange" pt 0 title "{/Symbol w} = .06"