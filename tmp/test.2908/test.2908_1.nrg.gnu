set term canvas size 1000,500 standalone mousing solid butt lw 1 jsdir "http://www.bioinfo.rpi.edu/geofold/output/gnuplot/"
set output '/bach1/home/sanw/public_html/geofold/output/test.2908/test.2908_1.nrg.html'
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
p "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_1.nrg.plot" u 1:2:3 w labels hypertext point pt 0 lw 2 lt rgb "black", "/bach1/home/sanw/GeoFold/tmp/test.2908/test.2908_1.nrg.plot" u 1:2 w lines lw 2 lt rgb "black"