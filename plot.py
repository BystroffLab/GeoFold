import commands
#make .plot in proper format
LName = "1VH0"
gnuplot = "/bach1/home/walcob/usr/bin/gnuplot"
nn = 5
#list of colors to use
colors = ['black','red','orange','yellow','green','blue','violet','cyan','magenta','pink','gold']
plot = [[]]
#open all the things
try:
  inFile = open("%s_%s.dag.nrg"%(LName,nn),'r')
except IOError:
  sys.exit("Failed to open file: %s_%s.dag.nrg"%(LName,nn))
try:
  gnuOut = open("%s_%s.nrg.gnu"%(LName,nn),'w+')
except IOError:
  inFile.close()
  sys.exit("Failed to open file: %s_%s.nrg.gnu"%(LName,nn))
try:
  plotOut = open("%s_%s.nrg.plot"%(LName,nn),'w+')
except IOError:
  inFile.close()
  gnuOut.close()
  sys.exit("Failed to open file: %s_%s.nrg.plot"%(LName,nn))

#Search .dag.nrg for ISTATE and corresponding TSTATE lines
#write out 3rd and 4th ISTATE lines and 4th and 5th of TSTATE
plotline=[]
for line in inFile:
  print(line)
  line = line.split()
  if line[0]=="ISTATE":
    plotline.append(line[2])
    plotline.append(line[3])
  if line[0]=="TSTATE" and len(plotline)==2:
    tmp = "%s %s"%(line[4],line[5])
    plotline.append(tmp)
  if len(plotline)==3:
    print("plotline")
    print(plotline)
    plot.append(plotline)
    plotline=[]
inFile.close()

plot.pop(0)
#write out to plotOut
for row in plot:
  print("row")
  print(row)
  plotOut.write("%s, %s, %s\n"%(row[0],row[1],row[2]))
plotOut.close()

#gnuplot script output
gnuOut.write('set term canvas size 1000,500 standalone mousing solid butt lw 1 jsdir "http://bach1.bio.rpi.edu/walcob/gnuplot"\n')
gnuOut.write("set output '%s_%s.html'\n"%(LName,nn))
#{/Symbol \253} = <----> {/Symbol \264} = multiplication symbol
gnuOut.write('set xlabel "unfolded <---> folded\\nBuried surface area (x1000A ^{2})"\n')
gnuOut.write('set ylabel "energy (kJ/mol)"\n')
gnuOut.write('set pointsize 2\n')
gnuOut.write('unset key\n')
gnuOut.write('set datafile separator ","\n')
#Cosmetics
gnuOut.write('set style line 80 lt rgb "#808080"\n')
gnuOut.write('set style line 81 lt 0\n')
gnuOut.write('set style line 81 lt rgb "#808080"\n')
gnuOut.write('set grid back linestyle 81\n')
gnuOut.write('set border 3 back linestyle 80\n')
gnuOut.write('set xtics nomirror\n')
gnuOut.write('set ytics nomirror\n')
#determine color and line type to use
kk = nn
count = 0
while kk > 11:
  count += 1
  kk = nn - 11
gnuOut.write('p "%s_%s.nrg.plot" u 1:2:3 w labels hypertext point pt %s lw 2 lt rgb "%s", "%s_%s.nrg.plot" u 1:2 w lines lw 2 lt rgb "%s"'%(LName,nn,count,colors[kk-1],LName,nn,colors[kk-1]))
gnuOut.close()
status,output = commands.getstatusoutput("%s<%s_%s.nrg.gnu"%(gnuplot,LName,nn))
print(output)
if status != 0:
  sys.exit(status)