
#!/usr/bin/env python
""" RUNGEOFOLD.py
This script runs the GeoFOLD suite of programs
It is a copy of the orginal C Shell server written by Dr. Chris Bystroff
The input files are (1) a parameters file and (2) a PDB file
All you need to do is modify the generic "parameters"
file so that it has the basename for the PDB file (i.e. usually the
4-letter PDB code), the chain letter(s) or "_" for
blank, and any modifications you might want.
Name the new file appropriately and then run this script:
./rungeofold.py myfile.par
output files will appear in the subdirectories output, tmp, log.
If these directories do not already exist, they are created.
The current directory must contain the geofold package,
installed by unpacking geofold.tgz
... then compiled... ????????????????
cd geofold
make clean
make all
... that's it.

-----------------------------------------------------
By modifying the directories below, you can
run this script and generate HTML output.
E. Gilbert 03-February-2014
B. Walcott 11 June 2014 """

# NOTE: If running this using python version 3,
# replace all occurrences of the commands module with subprocess

#### COMMAND LINE PARAMETERS ####
import os
import os.path
import sys
import commands
import time
import math
from georansac import fit

def makeZip(directory,LName):
  '''creates a zip archive of the directory and stores it within itself'''
  os.chdir("%s/.."%(directory))
  status,output = commands.getstatusoutput('zip -rv %s/%s.zip %s'%(LName,LName,LName))
  if status != 0:
    writeOut('%s: %s'%(status,output))
    runProgram('error')

def readConf(confFile):
  output = {}
  conf = open(confFile,'r')
  for line in conf:
    if line[0] != "#":
      line = line.split()
      if len(line) > 2:
        line[1] = " ".join(line[1:])
      if len(line) == 2:
        output[line[0]] = line[1]
  conf.close()
  return output


def writeOut(x):
  global htmlOut
  try:
    wout = open(htmlOut,'a')
  except IOError:
    sys.exit("Couldn't open file: %s"%(htmlOut))
  wout.write(x)
  wout.close()

def createForm(out,parfile):
  """write the output for the final form for do-over submission"""
  global parameters
  print parameters
  out.write('<FORM METHOD="POST" ACTION="../../geocgi.cgi" >\n')
  out.write('<input type="hidden" name="script" value=3>\n')
  ##New name for job
  out.write('<br><input type="text" name="keyword" value="" placeholder="Enter a new unique id for this job (avoid the words \'error\' and \'bug\')" size=80>\n')
  ##hidden parameters
  #Keep same pdb file
  #link to old parameters file to parsed by script
  out.write('<input type="hidden" name="oldParameters" value="%s">\n'%(parfile))
  #out.write('<input type="hidden" name="pdbcode" value="%s">\n'%(parameters["PDBCODE"]))
  #Keep same email
  #out.write('<input type="hidden" name="email_address" value="%s">\n'%(parameters["EMAIL"]))
  #submit button
  out.write('<br><input type="submit" name = "submit" value="submit"></form>')


#new energy profile function
def writeEnergyProfile(tmpDir,htmlDir,LName,nn):
  global gnuplot
  #list of colors to use
  colors = ['black','red','orange','yellow','green','blue','violet','cyan','magenta','pink','gold']
  plot = [[]]
  #open all the things
  try:
    inFile = open("%s/%s_%s.dag.nrg"%(tmpDir,LName,nn),'r')
  except IOError:
    return [1,"Failed to open file: %s_%s.dag.nrg"%(LName,nn)]
  try:
    gnuOut = open("%s/%s_%s.nrg.gnu"%(tmpDir,LName,nn),'w+')
  except IOError:
    inFile.close()
    return [1,"Failed to open file: %s_%s.nrg.gnu"%(LName,nn)]
  try:
    plotOut = open("%s/%s_%s.nrg.plot"%(tmpDir,LName,nn),'w+')
  except IOError:
    inFile.close()
    gnuOut.close()
    return [1,"Failed to open file: %s_%s.nrg.plot"%(LName,nn)]

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
      # print("plotline")
      # print(plotline)
      plot.append(plotline)
      plotline=[]
  inFile.close()

  plot.pop(0)
  #write out to plotOut
  for row in plot:
    # print("row")
    # print(row)
    plotOut.write("%s, %s, %s\n"%(row[0],row[1],row[2]))
  plotOut.close()

  #gnuplot script output
  gnuOut.write('set term canvas size 1000,500 standalone mousing solid butt lw 1 jsdir "http://www.bioinfo.rpi.edu/geofold/output/gnuplot/"\n')
  gnuOut.write("set output '%s/%s_%s.nrg.html'\n"%(htmlDir,LName,nn))
  print("set output '%s/%s_%s.nrg.html'\n"%(htmlDir,LName,nn))
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
  gnuOut.write('p "%s/%s_%s.nrg.plot" u 1:2:3 w labels hypertext point pt %s lw 2 lt rgb "%s", "%s/%s_%s.nrg.plot" u 1:2 w lines lw 2 lt rgb "%s"'%(tmpDir,LName,nn,count,colors[kk-1],tmpDir,LName,nn,colors[kk-1]))
  gnuOut.close()
  print("%s < %s/%s_%s.nrg.gnu"%(gnuplot,tmpDir,LName,nn))
  status,output = commands.getstatusoutput("%s < %s/%s_%s.nrg.gnu"%(gnuplot,tmpDir,LName,nn))
  print(status,output)
  if status != 0:
    return [status,output]
  else:
    try:
      gnuOut = open("%s/%s_%s.nrg2.gnu"%(tmpDir,LName,nn),'w+')
    except IOError:
      return [1,"Failed to open file %s_%s.nrg2.gnu"%(LName,nn)]
    gnuOut.write('set term postscript enhanced portrait\n')
    gnuOut.write('set size 1.4,0.7\n')
    gnuOut.write('set title "%s_%s"\n'%(LName,nn))
    gnuOut.write('set output "%s/%s_%s.nrg.ps"\n'%(tmpDir,LName,nn))
      # {/Symbol \253} = <----> {/Symbol \264} = multiplication symbol
    gnuOut.write('set xlabel "unfolded <---> folded\\nBuried surface area (x1000A ^{2})"\n')
    gnuOut.write('set ylabel "energy (kJ/mol)"\n')
    gnuOut.write('set pointsize 2\n')
    gnuOut.write('unset key\n')
    gnuOut.write('set datafile separator ","\n')
    # Cosmetics
    gnuOut.write('set style line 80 lt rgb "#808080"\n')
    gnuOut.write('set style line 81 lt 0\n')
    gnuOut.write('set style line 81 lt rgb "#808080"\n')
    gnuOut.write('set grid back linestyle 81\n')
    gnuOut.write('set border 3 back linestyle 80\n')
    gnuOut.write('set xtics nomirror\n')
    gnuOut.write('set ytics nomirror\n')
    print(colors[kk-1])
    gnuOut.write('p "%s/%s_%s.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "%s" pt %s'%(tmpDir,LName,nn,colors[kk-1],count))
    gnuOut.close()
    status,output = commands.getstatusoutput("%s < %s/%s_%s.nrg2.gnu"%(gnuplot,tmpDir,LName,nn))
    return [status,output]

#generate the molscript movie
def makeMolScript(mol, toDir = ' '):
  #programs used
  convert = 'convert'
  molScript = '/bach1/home/bystrc/bin/molscript'
  molAuto = '/bach1/home/bystrc/bin/molauto'

  print("makeMolScript")


  #run molAuto
  status,output = commands.getstatusoutput('%s -nice %s.pdb > %s.mol'%(molAuto,mol,mol))
  #error check
  if status != 0:
    return [status,output]
  #open the file created by molAuto and set gg = the line number of the "transform atom" line
  try:
    molOut = open('%s.mol'%(mol),'r')
  except IOError:
    return [1,"Couldn't open file %s.mol"%(mol)]
  gg = 0
  count = 0
  for line in molOut:
    count += 1
    if "transform atom" in line:
      gg = count
  hh = gg+1 ##hh is the line where the secondary structure data starts being stored
  ##nn, sw, and tlt are used to move the molecule around in the gif created
  sw = 25
  tlt = 0
  nn = 0
  while nn < 360:
    mm = 10000 + nn
    #write the first gg lines of mol.mol to mol.mm
    mmFile = open('%s.%s'%(mol,mm),'w+')
    molOut.seek(0,0)
    for i in range(0,gg):
      mmFile.write(molOut.readline())
    mmFile.write(" transform atom * by rotation z %s.;\n"%(tlt))
    mmFile.write(" transform atom * by rotation y %s.;\n"%(sw))
    mmFile.write(" transform atom * by rotation z -%s.;\n"%(tlt))
    mmFile.write(" transform atom * by rotation y %s.;\n"%(nn))

    #use hh to output structure data from molOut to mmFile
    molOut.seek(0,0)
    lines = []
    for line in molOut:
      lines.append(line)
    for i in range(hh,len(lines)):
      while 'end_plot' not in lines[i]:
        mmFile.write(lines[i])
        i += 1
    mmFile.write('set bonddistance 2.2;\n')
    mmFile.write('ball-and-stick require in type CYS and either atom CA, atom CB or atom SG ;\n')
    #write last line of mol.mol to mol.mm
    mmFile.write(lines[len(lines)-1])
    mmFile.close()
    #run molscript on mol.mm
    status,output = commands.getstatusoutput('%s < %s.%s > %s.%s.ps'%(molScript,mol,mm,mol,mm))
    if status != 0:
      return [status,output]
    nn += 5
    tlt += 5
  molOut.close()
  gif = '%s.gif'%(mol)
  runConvert = '%s -trim -delay 20 -geometry 300x300 %s.1????.ps %s'%(convert,mol,gif)
  status,output = commands.getstatusoutput(runConvert)
  if status != 0:
    return [status,output]
  if toDir != ' ':
    commands.getstatusoutput('mv %s %s'%(gif,toDir))

  #cleanup
  commands.getstatusoutput('rm -f %s.1????.ps %s.1????'%(mol,mol))
  return [0,'']

#gets the number of hydrogen and disulfide bonds from the hbfile
#returns a tuple of the # of h-bonds and disulfide bonds
#Used for the final HTML output
def getHS(LName):
  try:
    hbFile = open("%s.hb"%(LName),'r')
  except IOError:
    print("Couldn't open %s.hb"%(LName))
    return [-1,-1]
  h = 0
  s = 0
  for line in hbFile:
    if line[0] != '!':
      line = line.split()
      if line[4] == 'H':
        h += 1
      if line[4] == 'S':
        s += 1
  return [h,s]


#plots the lnku
def createGnuplot(LName, wat):
  #initial gnuplot setup
  #os.environ['GNUHELP']='/bach1/usr/local/share/gnuplot/4.2/gnuplot.gih'
  #os.environ['GNUPLOT_PS_DIR']='/bach1/usr/local/share/gnuplot/4.2/PostScript'
  global thermal
  global gnuplot

  print("createGnuplot")
  try:
    inFile = open("%s.fit"%(LName),'r')
  except IOError:
    return [1,"File not found: %s.fit"%(LName)]
  ##assign bb and mm
  line = inFile.readline()
  while "Least-squares coefficients" not in line:
    line = inFile.readline()
  line = inFile.readline()
  bb = inFile.readline()
  mm = inFile.readline()
  bb = bb.split()
  mm = mm.split()
  bb = bb[1]
  mm = mm[1]

  #open output file
  try:
    gnuOut = open("%s.gnu"%(LName),'w+')
  except IOError:
    return [1,"Couldn't open file %s.gnu"%(LName)]
  #initial gnuplot settings
  gnuOut.write('set terminal postscript enhanced portrait\n')
  gnuOut.write('set size 1.4,0.7\n')
  gnuOut.write('set output "%s.gp.ps"\n'%(LName))
  gnuOut.write('set encoding iso_8859_1\n') #allows for the use of special characters in gnuplot
  #{\305} = angstrom symbol {/Symbol \167} = omega
  #gnuOut.write('set xlabel "desolvation energy, {/Symbol \167} (J/mol/{\305}^2)"\n')
  if not thermal:
    gnuOut.write('set xlabel "solvation weight, {/Symbol \167}"\n')
  else:
    gnuOut.write('set xlabel "1/T, K^{-1}"\n')
  gnuOut.write('set ylabel "unfolding rate, ln(k_u) s^{-1}"\n')
  gnuOut.write('set pointsize 2\n')
  #Cosmetics
  gnuOut.write('set style line 80 lt rgb "#808080"\n')
  gnuOut.write('set style line 81 lt 0\n')
  gnuOut.write('set style line 81 lt rgb "#808080"\n')
  gnuOut.write('set grid back linestyle 81\n')
  gnuOut.write('set border 3 back linestyle 80\n')
  gnuOut.write('set xtics nomirror\n')
  gnuOut.write('set ytics nomirror\n')
  gnuOut.write('set key right\n')
  #The line that actually tells gnuplot what to plot
  gnuOut.write('p "%s.plot" w p pt 7, %s + %s*x\n'%(LName,bb,mm))
  gnuOut.close()
  status,output = commands.getstatusoutput("%s < %s.gnu"%(gnuplot,LName))
  return [status,output]



#plots the unfolding timecourse
def plotTimeCourse(LName,nn):
  global gnuplot
  print("plotTimeCourse")
  #initial setup
  try:
    gnuOut = open("%s_%s.tc.gnu"%(LName,nn),'w+')
  except IOError:
    return [1,"Couldn't open file %s_%s.tc.gnu"%(LName,nn)]
  # All info for this plot is read from the LName_nn.log file
  try:
    inFile = open("%s_%s.log"%(LName,nn),'r')
  except IOError:
    gnuOut.close()
    return [1,"Couldn't open file %s_%s.log"%(LName,nn)]
  #gnuplot environmental variables
  #os.environ['GNUPLOT_PS_DIR'] = '/bach1/usr/local/share/gnuplot/4.2/PostScript'
  #os.environ['GNUHELP'] = '/bach1/usr/local/share/gnuplot/4.2/gnuplot.gih'
  #initial gnuplot settings
  gnuOut.write('set terminal postscript enhanced portrait\n')
  gnuOut.write('set size 1.4,0.7\n')
  gnuOut.write('set output "%s_%s.tc.ps"\n'%(LName,nn))
  gnuOut.write('set encoding iso_8859_1\n')
  gnuOut.write('set xlabel "time (sec)"\n')
  gnuOut.write('set logscale x\n')
  gnuOut.write('set ylabel "% of total concentration"\n')
  gnuOut.write('set pointsize 2\n')
  #Cosmetics
  gnuOut.write('set style line 80 lt rgb "#808080"\n')
  gnuOut.write('set style line 81 lt 0\n')
  gnuOut.write('set style line 81 lt rgb "#808080"\n')
  gnuOut.write('set grid back linestyle 81\n')
  gnuOut.write('set border 3 back linestyle 80\n')
  gnuOut.write('set xtics nomirror\n')
  gnuOut.write('set ytics nomirror\n')
  gnuOut.write('set key outside right\n')
  #fOut = folded plot, uOut = unfolded plot, iOut = intermediate plot
  try:
    fOut = open("%s_%s.ftc.plot"%(LName,nn),'w+')
  except IOError:
    gnuOut.close()
    inFile.close()
    return [1,"Couldn't open file %s_%s.ftc.plot"%(LName,nn)]
  try:
    uOut = open("%s_%s.utc.plot"%(LName,nn),'w+')
  except IOError:
    gnuOut.close()
    inFile.close()
    fOut.close()
    return [1,"Couldn't open file %s_%s.utc.plot"%(LName,nn)]
  try:
    iOut = open("%s_%s.itc.plot"%(LName,nn),'w+')
  except IOError:
    gnuOut.close()
    inFile.close()
    fOut.close()
    uOut.close()
    return [1,"Couldn't open file %s_%s.itc.plot"%(LName,nn)]
  # Reads in data from TIMECOURSE lines
  for line in inFile:
    line = line.split()
    if len(line)!=0:
      if line[0]=='TIMECOURSE':
        fOut.write("%s %s\n"%(line[1],line[2]))
        uOut.write("%s %s\n"%(line[1],line[3]))
        iOut.write("%s %s\n"%(line[1],line[4]))
  fOut.close()
  uOut.close()
  iOut.close()
  inFile.close()
  #These next three lines tell gnuplot to plot the three curves together, each in a different color (red, blue, purple)
  #Note that it outputs as one line in gnuOut
  gnuOut.write('p "%s_%s.ftc.plot" w l lt rgb "red" lw 4 title "folded", '%(LName,nn))
  gnuOut.write('"%s_%s.utc.plot" w l lt rgb "blue" lw 4 title "unfolded", '%(LName,nn))
  gnuOut.write('"%s_%s.itc.plot" w l lt rgb "purple" lw 4 title "intermediate"\n'%(LName,nn))
  gnuOut.close()
  status,output = commands.getstatusoutput('%s < %s_%s.tc.gnu'%(gnuplot,LName,nn))
  return [status,output]

#Like writeEnergyProfile, but it superimposes the profile for all omega values onto one graph
def energyProfileAll(LName,omegaRange):
  print("energyProfileAll")
 #list of colors to use
  colors = ['black','red','orange','yellow','green','blue','violet','cyan','magenta','pink','gold']
  #gnuplot environmental variables
  # os.environ['GNUPLOT_PS_DIR'] = '/bach1/usr/local/share/gnuplot/4.2/PostScript'
  # os.environ['GNUHELP'] = '/bach1/usr/local/share/gnuplot/4.2/gnuplot.gih'
  global gnuplot
 #from geofold this program takes the arguments LName, omegaRange

  try:
    gnuOut = open("%s_all.nrg.gnu"%(LName),'w+')
  except IOError:
    return [1,"Couldn't open file %s_all.nrg.gnu"%(LName)]
 #initial gnuplot settings for plot
  gnuOut.write('set terminal postscript enhanced portrait\n')
  gnuOut.write('set size 1.4,0.7\n')
  gnuOut.write('set output "%s_all.nrg.ps"\n'%(LName))
  gnuOut.write('set encoding iso_8859_1\n')
  #{/Symbol \253} = <----> {/Symbol \264} = multiplication symbol
  gnuOut.write('set xlabel "unfolded {/Symbol \253} folded\\nBuried surface area ({/Symbol \264}1000{\305}^{2})" font "TimesBold,18"\n')
  gnuOut.write('set ylabel "energy (kJ/mol)" font "TimesBold,18"\n')
  gnuOut.write('set pointsize 2\n')
  gnuOut.write('set key outside right\n')
  #Cosmetics
  gnuOut.write('set style line 80 lt rgb "#808080"\n')
  gnuOut.write('set style line 81 lt 0\n')
  gnuOut.write('set style line 81 lt rgb "#808080"\n')
  gnuOut.write('set grid back linestyle 81\n')
  gnuOut.write('set border 3 back linestyle 80\n')
  gnuOut.write('set xtics nomirror\n')
  gnuOut.write('set ytics nomirror\n')
  nn = 0
  count = 0
  gnuOut.write('p ')
  for value in omegaRange:
    nn+=1
    kk = nn
    while kk > 11:
      count += 1
      kk = nn-11
      #{/Symbol \167} = lowercase omega
    gnuOut.write('"%s_%s.nrg.plot" u 1:2 w linespoints lw 4 lt rgb "%s" pt %s title "{/Symbol \167} = %s"'%(LName,nn,colors[kk-1],count,value))
    if nn != len(omegaRange):
      gnuOut.write(', ')
  gnuOut.close()
  status,output = commands.getstatusoutput('%s < %s_all.nrg.gnu'%(gnuplot,LName))
  return [status,output]

#adds link to non-canvas graph in the energy profile html file
def picLink(htmlDir,LName,nn):
  link = '<br><a href= "%s_%s.nrg.png">Click here if plot does not load.</a></body>'%(LName,nn)
  try:
    file = open("%s/%s_%s.nrg.html"%(htmlDir,LName,nn),'r')
  except IOError:
    return [1,"Failed to open file %s/%s_%s"%(htmlDir,LName,nn)]
  text =''
  for line in file:
    text += line
  file.close()
  text = text.split('</body>')
  text = text[0]+link+text[1]
  file = open("%s/%s_%s.nrg.html"%(htmlDir,LName,nn),'w+')
  for line in text:
    file.write(line)
  file.close()
  return [0,'']


#runs a command and exits the program if an error is found
#if the command is "error" then it is an internal error and
#runProgram is being called to handle it
def runProgram(command):
  global tmpWrite
  global outWrite
  global htmlOut
  global htmlTmp
  global debug
  print(command)
  if command != "error":
    status,output = commands.getstatusoutput(command)
  else:
    status = 1
    output = ''
  if "error" in output.lower() or "bug" in output.lower() or status != 0:
    paramFile.close()
    permWrite.close()
    if command != "error":
      tmpWrite.write("Error in %s\n%s: %s"%(command,status,output))
      writeOut("Error in %s\n%s: %s"%(command,status,output))
      writeOut('</pre></body></html>\n')
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()

    sys.exit("Error in %s\n%s: %s"%(command,status,output))
  print(output)
  tmpWrite.write(output+'<br>')
  tmpWrite.close()
  if debug:
    makeCopy(htmlTmp,htmlOut)
    outWrite = open(htmlOut,'a')
    outWrite.write('</pre></body></html>\n')
    outWrite.close()
  tmpWrite = open(htmlTmp,'a')
  return status,output

#function used to copy htmlTmp to htmlOut so client can see progress of geofold
def makeCopy(inFile,outFile):
  try:
    out = open(outFile,'w+')
  except IOError:
    return [1,'Could not open file %s'%(out)]
  try:
    readIn = open(inFile,'r')
  except IOError:
    out.close()
    return [1,'Could not open file %s'%(inFile)]
  for line in readIn:
    out.write(line)
  out.close()
  readIn.close()
  return [0,'']


#Takes the imagemap generated by dot and fixed by fixMap and creates the final DAG html file
def makeMap(htmlDir,LName):
  mapFile = "%s/%s.dag.map"%(htmlDir,LName) #client-side imagemap created by dot
  dagFile = "%s.dag"%(LName)  #dag generated by geofold
  try:
    mapIn = open(mapFile,'r')
  except IOError:
    return [1, "Couldn't open file %s"%(mapFile)]
  try:
    mapOut = open("%s/%s.html"%(htmlDir,dagFile),'w+')
  except IOError:
    mapIn.close()
    return [1,"Couldn't open file %s/%s.html"%(htmlDir,dagFile)]
  #replace the first line of the file
  mapOut.write('<html><body><img src="%s.png" border=0 usemap="#dagmap">\n'%(dagFile))
  mapOut.write('<map name="dagmap">\n')
  for line in mapIn:
    if line.rstrip() != '<map id="FILENAME" name="FILENAME">':
      mapOut.write(line)
  mapOut.write('</body></html>')
  mapOut.close()
  mapIn.close()
  return [0,'']

#find parameters in the parameters file
def findParam (paramFile, target):
  global parameters
  paramFile.seek(0,0) #Go back to the beginning of the file
  if target == "ORANGE":  #Special case for ORANGE
    for line in paramFile:
      if target in line:
        s = line.split()
        s.pop(0)
        i = 0
        while i <len(s):
          s[i]=s[i].strip(' ,\n\r') #get rid of non-numerical characters that may be present
          if s[i]=='':  #remove empty components
            s.pop(i)
            i-=1
          i+=1
        if len(s) == 0:
          s.append(findParam(paramFile,'OMEGA')[1])
        parameters[target]=s
        return [0,s]
    s = []
    s.append(findParam(paramFile,'OMEGA')[1]) #if ORANGE not in parameters, set it to OMEGA
    parameters[target]=s
    return [0,s]
  else:
    for line in paramFile:
      if target in line:
        s= line.split()
        if len(s) <= 1:
          return [1,'']
        else:
          print(s)
          parameters[target]=s[1]
          return [0,s[1]]

  return [1, '']

# Setting global variables
global tmpWrite
global outWrite
global htmlTmp
global htmlOut
global parameters
global thermal
global debug
global gnuplot
debug = False
parameters = {}

conf = "default.conf"

########################## DIRECTORIES #############################
### SET THESE DIRECTORIES AS FOLLOWS:
# gDir is the geofold directory created by tar -zxvf geofold.tgz
# tmpDir is a temporary file directory. Everything in that
# directory can be deleted later. The PDB file should be copied
# to tmpDir to start.
# logDir is a place to hold the log file.
# htmlDir is a public_html directory where HTML
# file may be viewed.
# jobDir holds job files that are read by geod.csh
# (not used in this script except to clean up)
# thisDir is the directory where you are running this script.
if len(sys.argv) == 3:
  conf = sys.argv[2]

configuration = readConf(conf)
thisDir = os.getcwd()
baseDir = configuration['baseDir']
gDir = configuration['gDir']
bDir = configuration['bDir']
maskerDir = configuration['maskerDir']
tmpDir = configuration['tmpDir']
pdbDir = configuration['pdbDir']
logDir = configuration['logDir']
htmlDir = configuration['htmlDir']
jobDir = configuration['jobDir']
paramTemplate = configuration['paramTemplate']
baseURL = configuration['baseURL']
outputURL = configuration['outputURL']
dot = configuration['dot']
convert = configuration['convert']
gnuplot = configuration['gnuplot']


if len(sys.argv) < 2 or len(sys.argv) > 4:
  #### JOB level variable from parameter file ####
  # pdbCode is the basename of a PDB file to be found in TMPDIR
  # chain is a set of chain IDs within the PDB file, or _ for space.
  # LName is a unique name for this job or the process ID
  # OName is a unique name for a previous job, to be used to skip GEOFOLD.
  # UName is a unique name for the current job (or username?)
  print("USAGE: rungeofold.py parametersFile [configuration.conf]")
  sys.exit()
else:
  #Read parameters from the parameters file
  paramFilename = sys.argv[1]
  if len(sys.argv) >= 3:
    Username = sys.argv[2]
  try:
    paramFile = open(paramFilename, 'r')
  except IOError:
    sys.exit(paramFilename+" could not be opened.")
  status,pdbCode= findParam(paramFile, "PDBCODE")
  if status != 0:
    paramFile.close()
    sys.exit("PDBCODE line not found.")
  print("PDBCODE %s" %(pdbCode))
  try:
    status, thermal = findParam(paramFile,"THERMAL")
    if thermal == "1":
      thermal = True
    else:
      thermal = False
  except LookupError:
      thermal = False
  try:
    keyword = findParam(paramFile, "KEYWORD")[1]
  except LookupError:
    keyword = ""
  status,chain=  findParam(paramFile, "CHAIN ")
  if status!=0:
    paramFile.close()
    sys.exit(paramFilename+" chainid not found.")
  print("Chain "+chain)
  status,LName=  findParam(paramFile, "LNAME")
  if status != 0:
    print(paramFilename+" LNAME value not found.")
    LName= str(os.getpid())
    if keyword != '':
      LName = '%s.%s'%(keyword,pdbCode)
  #Strip LName of any suspicious characters
  LName = LName.strip("/*{}[]!@#$%^&();:<>,?\\~\"\'")
  print("LName " + LName)
  findParam(paramFile, "EMAIL")
  status,OName=  findParam(paramFile, "ONAME")
  if status != 0:
    OName=LName
  rerun = OName == LName
  UName = pdbCode + chain + LName
  print("ONAME "+OName)
  status,email = findParam(paramFile,"EMAIL")
  #parameters needed to be added to parameters dicitonary
  findParam(paramFile,"REDUCING")
  status,debug = findParam(paramFile,"DEBUG")
  if status != 0: debug = False
  else: debug = debug == "1"
  print("debug: %s"%debug)
  status,barrels = findParam(paramFile,"BARRELMOVES")
  barrels = barrels == "1"
  if status != 0:
    barrels = True

  #job directory setup
  tmpDir = "%s/%s"%(tmpDir,LName)
  htmlDir = "%s/%s"%(htmlDir,LName)
  try:
    os.makedirs(tmpDir,0755)
    os.makedirs(htmlDir,0755)
  except OSError as e:
    pass
  #dagDir setup
  #htmlDirs = htmlDir.split('/')
  #gDirs = gDir.split('/')
  #dagDir = '/'.join(htmlDirs[len(gDirs)+1:len(htmlDirs)])
  os.environ["dagDir"] = htmlDir
  cp = "cp -p %s/isegment.cgi %s/isegment.cgi"%(gDir,htmlDir)
  commands.getstatusoutput(cp)

  #Added for do-over script interface
  status,redo = findParam(paramFile,"REDO")
  if status == 0 and len(redo) != 0:
    print tmpDir
    cmd = "cp -v %s%s/* %s/"%(tmpDir.split(LName)[0],redo,tmpDir)
    status,output = commands.getstatusoutput(cmd)
    print cmd
    print status
    print output
    cmd = "find %s -name %s* -print"%(tmpDir,redo)
    status,output = commands.getstatusoutput(cmd)
    print cmd
    print status
    print output
    print "line parsing"
    output = output.split('\n')
    print output
    for line in output:
      line = line.split(redo)
      if len(line) > 1:
        print line
        cmd = "mv -v %s/%s%s %s/%s%s"%(tmpDir,redo,line[1],tmpDir,LName,line[1])
        print cmd
        status,output = commands.getstatusoutput(cmd)
        print status
        print output


  #### FILES ####
  htmlRefresh= "%s/header_refresh.html" %(bDir)
  htmlHead= "%s/header.html" %(bDir)
  htmlTmp= "%s/%s.html" %(tmpDir,LName)
  htmlLog= "%s/%s.log" %(logDir,LName)
  htmlOut= "%s/%s.html" %(htmlDir,LName)
  htmlPerm = "%s/%s.perm" %(tmpDir,LName)

  ###PROGRAMS####
  maxTraffic= gDir+"/maxTraffic"
  mtCut= 0.1
  #convert = "/usr/bin/convert"
  #dot = "/usr/bin/dot"

  try:
    tmpWrite = open(htmlTmp, 'w+')
  except IOError:
    paramFile.close()
    sys.exit("Couldn't open file %s"%(htmlTmp))
  try:
    permWrite = open(htmlPerm,'w+')
  except IOError:
    paramFile.close()
    sys.exit("Couldn't Open file %s"%(htmlPerm))

  #Add refresh info to tmpWrite + htmlOut
  tmpWrite.write('<html><head><meta http-equiv="refresh" content="2"></head><body>')
  outWrite = open(htmlOut,'w+')
  outWrite.write('<html><head><meta http-equiv="refresh" content="2"></head><body><pre>')
  outWrite.close()


  #### VARIABLES ####
  # omega is the desolvation free energy
  # maxSplit is the maximum number of splits for each intermediate
  #   of unfolding.
  # folding (folding) = 1 is folding, 0 is unfolding
  # halflife = 1 for stopping at the half-life, 0 to go to
  #equilibrium in Unfoldsim

  #variable parameters
  status,folding= findParam (paramFile, "FING")
  if status != 0:
    folding = 0
  status, halflife= findParam (paramFile, "HLFE")
  if status !=0:
    halflife = 0
  status,bCut= findParam (paramFile, "BREAKCUT")
  if status != 0:
    paramFile.close()
    sys.exit("BREAKCUT not found")
  status, pCut= findParam (paramFile, "PIVOTCUT")
  if status != 0:
    paramFile.close()
    sys.exit("PIVOTCUT not found")
  status,hCut= findParam (paramFile, "HINGECUT")
  if status != 0:
    paramFile.close()
    sys.exit("HINGECUT not found")
  status,wat= findParam (paramFile, "WAT")
  if status != 0:
    paramFile.close()
    sys.exit("WAT not found")
  status,voids=findParam (paramFile, "VOIDENTROPY")
  if status != 0:
    paramFile.close()
    sys.exit("VOIDENTROPY not found")
  voids = float(voids)
  status,omegaRange= findParam (paramFile, "ORANGE")
  print("ORANGE ",omegaRange)
  status,maxSplit = findParam (paramFile, "MAXSPLIT")
  if status != 0:
    maxSplit = 4
  status, doIt = findParam (paramFile, "RUNGEOFOLD")
  if status != 0:
    doIt = 1
  doIt = int(doIt)
  print("doIt: %s"%(doIt))
  os.environ['LNAME'] = LName
  if doIt == 1 :
    print "Running  GEOFOLD"
  else:
    print "Skipping GEOFOLD"
    refresh = "cat %s > %s" %(htmlRefresh,htmlTmp)
    commands.getstatusoutput(refresh)
    tmpWrite.write("<h3>Skipping GeoFOLD. UnfoldSim results for " + UName + "</h3><br>")
    tmpWrite.write("<h4>Using previously calculated DAG file " + OName + " </h4><br>")
    tmpWrite.write("<pre><br>")
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    tmpWrite.write("</body></html>")
    if OName == UName:
      print("ERROR Need old job name as arg 4 to skip geofold.")
      paramFile.close()
      permWrite.close()
      tmpWrite.close()
      sys.exit()
    cp = "cp %s/%s.dag %s/%s.dag" %(tmpDir,OName,tmpDir,LName)
    commands.getstatusoutput(cp)
  if doIt == 1 :
    print("============= PARAMETERS =============")
    tmpWrite.write("============= PARAMETERS =============<br>")
    writeOut("============= PARAMETERS =============\n")
    if os.path.getsize(paramFilename)==0:
      tmpWrite.write("Parameters file is empty!!\n")
      sys.exit("Parameters file is empty!!!")
    tmpWrite.write("Using these parameters<br>")
    for line in paramFile:
      tmpWrite.write(line+'<br>')
    paramFile.close()
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()

    print("============= GETCHAIN =============")
    tmpWrite.write("============= GETCHAIN =============<br>")
    writeOut("============= GETCHAIN =============\n")
    pdbFile = "%s/%s.pdb" %(pdbDir,pdbCode)
    if not os.path.isfile(pdbFile):
      sys.exit("File not found: %s\n" %(pdbFile))
    tmpWrite.write("Extracting protein atoms from %s<br>" %(pdbCode))
    getProtein = "%s/xgetchain + < %s/%s.pdb > %s/%s.tmp" %(gDir,pdbDir,pdbCode,tmpDir,LName)
    runProgram(getProtein)
    tmpWrite.write("Extracting chains %s<br>" %(chain))
    getChains = "%s/xgetchain %s < %s/%s.tmp > %s/%s.pdb" %(gDir,chain,tmpDir,LName,tmpDir,LName)
    runProgram(getChains)
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    print("============= RENUMBER =============" )
    tmpWrite.write("============= RENUMBER =============<br>" )
    writeOut("============= RENUMBER =============\n")
    LNamePDB = "%s/%s.pdb" %(tmpDir,LName)
    if not os.path.isfile(LNamePDB):
      sys.exit("File not found: %s\n" %(LNamePDB))
    renumberOne = "%s/xrenumber_one %s %s/%s.tmp" %(gDir,LNamePDB,tmpDir,LName)
    runProgram(renumberOne)
    mv = "mv %s/%s.tmp %s" %(tmpDir,LName,LNamePDB)
    print(mv)
    commands.getstatusoutput(mv)
    cp = "cp %s %s/%s.pdb" %(LNamePDB,htmlDir,LName)
    print(cp)
    commands.getstatusoutput(cp)
    if os.path.getsize(LNamePDB)==0:
      tmpWrite.write("Error. empty file.\n")
      if debug:
        tmpWrite.close()
        makeCopy(htmlTmp,htmlOut)
        tmpWrite = open(htmlTmp,'a')
        outWrite = open(htmlOut,'a')
        outWrite.write('</pre></body></html>\n')
        outWrite.close()
    print("============= 3to1 (extract sequence) =============")
    tmpWrite.write("============= 3to1 (extract sequence) =============<br>")
    writeOut("============= 3to1 (extract sequence) =============\n")
    x3to1 = "%s/x3to1 + < %s > %s/%s.seq" %(gDir,LNamePDB,tmpDir,LName)
    runProgram(x3to1)
    tmpWrite.write("Complete sequence (all chains ) is<br>")
    seqFile = open("%s/%s.seq"%(tmpDir,LName))
    for line in seqFile:
      tmpWrite.write(line+'<br>')
    seqFile.close()
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    print("============= PDB2CIJ (extract contacts) =============")
    tmpWrite.write("============= PDB2CIJ (extract contacts) =============<br>\n")
    writeOut("============= PDB2CIJ (extract contacts) =============\n")
    pdb2cij = "%s/xpdb2cij %s 8. > %s/%s.cij" %(gDir,LNamePDB,tmpDir,LName)
    runProgram(pdb2cij)
    tmpWrite.write("Number of contacts found: <br>")
    writeNumContacts = "wc -l %s/%s.cij | awk '{print $1}'>> %s" %(tmpDir,LName,htmlTmp)
    status,output = commands.getstatusoutput(writeNumContacts)
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    print("============= PDB2HB (extract H-bonds SS-bonds) =============")
    tmpWrite.write("============= PDB2HB (extract H-bonds SS-bonds) =============<br>")
    writeOut("============= PDB2HB (extract H-bonds SS-bonds) =============\n")
    pdb2hb = "%s/xpdb2hb %s %s > %s/%s.hb" %(gDir,paramFilename,LNamePDB,tmpDir,LName)
    runProgram(pdb2hb)
    writeParam = "echo '\nHBONDS %s/%s.hb' >> %s" %(tmpDir,LName,paramFilename)
    status,output = commands.getstatusoutput(writeParam)
    print("============= CONTACTMASK =============")
    tmpWrite.write("============= CONTACTMASK =============<br>\n")
    writeOut("============= CONTACTMASK =============\n")
    writeTime = "Time before running CONTACTMASK "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)

    ##masker setup
    os.environ['MASKERDIR']=maskerDir
    os.environ['DESIGN_HOME']=os.environ['MASKERDIR']
    os.environ['MSIZE']="512"
    os.environ['DTHETA']="9.0"
    os.environ['MASKLIB']="%s/%s.mask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['VMASK']="%s/%s.vmask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['VTRIANGLE']="%s/%s.tri" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['MASKTEMPLATE']="%s/%s.sort" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['MASKDAT']="%s/%s.dat" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['ATOMLIB']="%s/atoms.lib" %(os.environ['MASKERDIR'])
    os.environ['COUNTBITS']="%s/cbits.bin" %(os.environ['MASKERDIR'])
    os.environ['COLLARS']="%s/collars%s.bin" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['SLOPES']="%s/slopes%s.dat" %(os.environ['MASKERDIR'],os.environ['DTHETA'])
    os.environ['FORCEFIELD']="%s/gromacs_ff.prm" %(os.environ['MASKERDIR'])

    contactmask = "%s/xcontactmask %s %s/%s.sas 1.4" %(maskerDir,LNamePDB,tmpDir,LName)
    runProgram(contactmask)
    writeTime = "Time after running CONTACTMASK "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)
    cp = "cp %s/%s.sas %s/%s.sas" %(tmpDir,LName,htmlDir,LName)
    commands.getstatusoutput(cp)
    line = 'Pairwise contact surfaces <a href = "%s.sas">file.</a><br>' %(LName)
    tmpWrite.write(line)
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    writeParam = "echo CONTACTS %s/%s.sas >> %s" %(tmpDir,LName,paramFilename)
    commands.getstatusoutput(writeParam)
    #================================ HB2CIJ ================================#
    print("============= HB2CIJ =============")
    tmpWrite.write("============= HB2CIJ =============<br>\n")
    writeOut("============= HB2CIJ =============\n")
    hb2cij = "%s/seams/xhb2cij %s/%s.hb %s/%s.hbcij"%(gDir,tmpDir,LName,tmpDir,LName)
    runProgram(hb2cij)
    #=======================================================================#
    print("============= PDB2SEAMS (extract seams) =============")
    tmpWrite.write("============= PDB2SEAMS (extract seams) =============<br>\n")
    writeOut("============= PDB2SEAMS (extract seams) =============\n")
    if barrels:
      pdb2seams = "%s/seams/xpdb2seams %s %s/%s.sas %s/%s.hbcij > %s/%s.seams" %(gDir,LNamePDB,tmpDir,LName,tmpDir,LName,tmpDir,LName)
    #   pdb2seams = "%s/seams/xpdb2seams %s %s/%s.sas > %s/%s.seams" %(gDir,LNamePDB,tmpDir,LName,tmpDir,LName)
      os.environ['GDIR']=gDir
      os.environ['TMPDIR']=tmpDir
      runProgram(pdb2seams)
    else:
      hbFile = open("%s/%s.hb"%(tmpDir,LName),'r')
      nres = 0
      found = False
      while not found:
        for line in hbFile:
          if "residues" in line:
            line = line.split()
            nres = line[1]
            found = True
            break
      hbFile.close()
      seamsFile = open("%s/%s.seams"%(tmpDir,LName),'w+')
      seamsFile.write(" Size of PDB file = %s\n"%(nres))
      seamsFile.write("# ========== Results from PDB2SEAMS =========\n")
      seamsFile.write("# NBARRELS <Number of Barrels>\n")
      seamsFile.write("# BARREL : <id> <Number of Seams>\n")
      seamsFile.write("# SEAM : <id> <Number of Buttons> <Total Energy> <x1-beta1> <x2-beta1> <x1-beta2> <x2-beta2>\n")
      seamsFile.write("# BUTTON : <id> <Energy-Beta1> <Energy-beta2>\n")
      seamsFile.write("NBARRELS      0\n")
      seamsFile.close()
    writeParam = "echo SEAMS %s/%s.seams >> %s" %(tmpDir,LName,paramFilename)
    commands.getstatusoutput(writeParam)
    print("============= SPLITSEAMS =============")
    tmpWrite.write("<br>============= SPLITSEAMS =============<br>\n")
    writeOut("============= SPLITSEAMS =============\n")
    if barrels:
      splitseams = "%s/xsplitseams %s/%s.pdb %s/%s.seams %s/%s.sas %s/%s.hb %s %s/%s.split > %s/%s.splitlog" %(gDir,tmpDir,LName,tmpDir,LName,tmpDir,LName,tmpDir,LName,paramFilename,tmpDir,LName,tmpDir,LName)
      runProgram(splitseams)
      writeParam = "echo SEAMS %s/%s.seams >> %s" %(tmpDir,LName,paramFilename)
    print("============= VOIDMASK =============")
    tmpWrite.write("============= VOIDMASK =============<br>\n")
    writeOut("============= VOIDMASK =============\n")


    writeTime = "Time before running VOIDMASK "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)

    ## This should be read from a file......
    ##masker setup
    os.environ['MASKERDIR']=maskerDir
    os.environ['DESIGN_HOME']=os.environ['MASKERDIR']
    os.environ['MSIZE']="512"
    os.environ['DTHETA']="9.0"
    os.environ['MASKLIB']="%s/%s.mask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['VMASK']="%s/%s.vmask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['VTRIANGLE']="%s/%s.tri" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['MASKTEMPLATE']="%s/%s.sort" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['MASKDAT']="%s/%s.dat" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['ATOMLIB']="%s/atoms.lib" %(os.environ['MASKERDIR'])
    os.environ['COUNTBITS']="%s/cbits.bin" %(os.environ['MASKERDIR'])
    os.environ['COLLARS']="%s/collars%s.bin" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
    os.environ['SLOPES']="%s/slopes%s.dat" %(os.environ['MASKERDIR'],os.environ['DTHETA'])
    os.environ['FORCEFIELD']="%s/gromacs_ff.prm" %(os.environ['MASKERDIR'])

    if voids != 0  :
      voidmask = "%s/xvoidmask %s %s/%s.void.pdb 1.4 1.2 1.4" %(maskerDir,LNamePDB,tmpDir,LName)
      runProgram(voidmask)
      cp = "cp %s/%s.void.pdb %s/%s.void.pdb" %(tmpDir,LName,htmlDir,LName)
      commands.getstatusoutput(cp)
    else:
      cp = "cp %s/%s.pdb %s/%s.void.pdb" %(tmpDir,LName,tmpDir,LName)
      commands.getstatusoutput(cp)
      cp = "cp %s/%s.pdb %s/%s.void.pdb" %(tmpDir,LName,htmlDir,LName)
      commands.getstatusoutput(cp)
    writeTime = "Time after running VOIDMASK "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)
    if os.path.getsize("%s/%s.void.pdb"%(tmpDir,LName))==0:
      tmpWrite.write("ERROR in VOIDMASK. Empty file.")
      if debug:
        tmpWrite.close()
        makeCopy(htmlTmp,htmlOut)
        tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write("ERROR in VOIDMASK. Empty file.")
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
      sys.exit("ERROR in VOIDMASK. Empty file.")
    line = 'Coordinates with void positions: <a href="%s.void.pdb"> %s.void.pdb </a><br>' %(LName,LName)
    tmpWrite.write(line)
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
    print("============= GEOFOLD =============")
    tmpWrite.write("============= GEOFOLD =============<br>")
    writeOut("============= GEOFOLD =============\n")
    writeTime = "Time before running GEOFOLD "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)
    geofold = "%s/xgeofold %s/%s.void.pdb %s/%s.dag %s > %s/%s.dag.log" %(gDir,tmpDir,LName,tmpDir,LName,paramFilename,tmpDir,LName)
    runProgram(geofold)
    writeTime = "Time after running GEOFOLD "+time.strftime("%c")+'<br>'
    tmpWrite.write(writeTime)
    writeOut(writeTime)
    dagFile = "%s/%s.dag" %(tmpDir,LName)
    if os.path.getsize(dagFile)==0:
      tmpWrite.write("ERROR in GEOFOLD. Empty file.\n")
      if debug:
        tmpWrite.close()
        makeCopy(htmlTmp,htmlOut)
      outWrite = open(htmlOut,'a')
      outWrite.write("ERROR in GEOFOLD. Empty file.")
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
      sys.exit("ERROR in GEOFOLD. Empty file.")
    cp = "cp %s/%s.dag %s/%s.dag" %(tmpDir,LName,htmlDir,LName)
    commands.getstatusoutput(cp)
    line = 'Directed acyclic graph (DAG) <a href="%s.dag">file.</a><br>' %(LName)
    tmpWrite.write(line)
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
  outWrite.close()
  ##SKIPGEOFOLD
  LNamePDB = "%s/%s.pdb" %(tmpDir,LName)
  if doIt != 3:
      print(doIt == 3)
      print(doIt)
      print("============= UNFOLDSIM =============")
      tmpWrite.write("Unfolding %s%s<br>" %(pdbCode,chain))
      tmpWrite.write("============= UNFOLDSIM =============<br>")
      writeOut("============= UNFOLDSIM =============\n")
      nn = 0
      for value in omegaRange:
        nn += 1
        cp = "cp %s/%s.dag %s/%s_%s.dag" %(tmpDir,LName,tmpDir,LName,nn)
        commands.getstatusoutput(cp)
        if not thermal:
          sed = "sed -e \"s/^OMEGA .*/OMEGA %s/\" %s > %s.1" %(value,paramFilename,paramFilename)
        else:
          sed = "sed -e \"s/^TEMPERATURE .*/TEMPERATURE %s/\" %s > %s.1"%(value, paramFilename,paramFilename)
        status,output=commands.getstatusoutput(sed)
        logFile = "%s/%s_%s.log" %(tmpDir,LName,nn)
        if not thermal:
          print("============= run %s omega = %s =============" %(nn,value))
          tmpWrite.write("============= run %s omega = %s =============<br>" %(nn,value))
          writeOut("============= run %s omega = %s =============\n" %(nn,value))
        else:
          print("============= run %s temp = %s K =============" %(nn,value))
          tmpWrite.write("============= run %s temp = %s K =============<br>" %(nn,value))
          writeOut("============= run %s temp = %s K =============\n" %(nn,value))
        writeTime = "Time before running UNFOLDSIM "+time.strftime("%c") +'<br>'
        tmpWrite.write(writeTime)
        writeOut(writeTime)
        unfoldsim = "%s/xunfoldsim %s/%s_%s.dag %s.1 > %s" %(gDir,tmpDir,LName,nn,paramFilename,logFile)
        tmpWrite.write(unfoldsim+'<br>')
        runProgram(unfoldsim)
        writeTime = "Time after running UNFOLDSIM "+time.strftime("%c")+'<br>'
        tmpWrite.write(writeTime)
        writeOut(writeTime)
        tmpWrite.write("<p><pre><br>")
        # grep = "grep ^\"TIMECOURSE\" %s | tail -1 >> %s" %(logFile,htmlTmp)
        log = open(logFile,'r')
        lines = []
        for line in log:
          linesplit = line.split()
          if len(linesplit) != 0 and linesplit[0]=='TIMECOURSE':
            lines.append(line)
        if len(lines)==0:
          sys.stderr.write("No timecourse data\n")
          tmpWrite.write("No timecourse data\n")
          runProgram("error")
        tmpWrite.write(lines[len(lines)-1])
        log.close()
        if debug:
          tmpWrite.close()
          makeCopy(htmlTmp,htmlOut)
          tmpWrite = open(htmlTmp,'a')
          outWrite = open(htmlOut,'a')
          outWrite.write('</pre></body></html>\n')
          outWrite.close()
  print("============= PATHWAY2PS =============")
  tmpWrite.write("============= PATHWAY2PS =============<br>")
  writeOut("============= PATHWAY2PS =============\n")
  nn = 0
  for value in omegaRange:
    nn+=1
    pathway2ps = "%s/xpathway2ps %s/%s.seq %s/%s_%s.dag.age %s/%s.cij %s/%s_%s.ps 4" %(gDir,tmpDir,LName,tmpDir,LName,nn,tmpDir,LName,tmpDir,LName,nn)
    tmpWrite.write(pathway2ps+'<br>')
    status,output = commands.getstatusoutput(pathway2ps)
    tmpWrite.write(output+'<br>')
    if debug:
      tmpWrite.close()
      makeCopy(htmlTmp,htmlOut)
      tmpWrite = open(htmlTmp,'a')
      outWrite = open(htmlOut,'a')
      outWrite.write('</pre></body></html>\n')
      outWrite.close()
  print("============= CONVERT =============")
  tmpWrite.write("============= CONVERT =============<br>")
  writeOut("============= CONVERT =============\n")
  nn=0
  for value in omegaRange:
    nn+=1
    runConvert = "%s -trim -geometry 100 -background white %s/%s_%s.ps %s/%s_%s_thumb.png" %(convert,tmpDir,LName,nn,tmpDir,LName,nn)
    tmpWrite.write(runConvert+'<br>')
    commands.getstatusoutput(runConvert)
    runConvert = "%s -background white %s/%s_%s_thumb.png %s/%s_%s_thumb.png" %(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    commands.getstatusoutput(runConvert)
    runConvert = "%s -trim -background white %s/%s_%s.ps %s/%s_%s.png" %(convert,tmpDir,LName,nn,tmpDir,LName,nn)
    commands.getstatusoutput(runConvert)
    runConvert = "%s -background white %s/%s_%s.png %s/%s_%s.png" %(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    commands.getstatusoutput(runConvert)

  print("============= MAXTRAFFIC =============")
  tmpWrite.write("============= MAXTRAFFIC =============<br>")
  writeOut("============= MAXTRAFFIC =============\n")
  nn = 0
  for value in omegaRange:
    nn+=1
    runMaxTraffic = "%s %s/%s_%s.dag.out %s > %s/%s_%s.dot" %(maxTraffic,tmpDir,LName,nn,mtCut,tmpDir,LName,nn)
    tmpWrite.write(runMaxTraffic+'<br>')
    runProgram(runMaxTraffic)
  print("============= DOT =============")
  tmpWrite.write("============= DOT =============<br>")
  writeOut("============= DOT =============\n")
  nn = 0
  for value in omegaRange:
    nn+=1
    runDot = "%s -Gdpi=100 -Gsize=6,100 -Tcmapx  -o %s/%s_%s.dag.map %s/%s_%s.dot" %(dot,htmlDir,LName,nn,tmpDir,LName,nn)
    tmpWrite.write(runDot+'<br>')
    runProgram(runDot)
    #fixMap('%s/%s_%s.dag.map'%(htmlDir,LName,nn))
    runDot = "%s -Gdpi=100 -Gsize=6,100 -Tpng -o %s/%s_%s.dag.png %s/%s_%s.dot" %(dot,htmlDir,LName,nn,tmpDir,LName,nn)
    tmpWrite.write(runDot+'<br>')
    runProgram(runDot)
  print("============= CONVERT =============")
  tmpWrite.write("============= CONVERT =============<br>")
  writeOut("============= CONVERT =============\n")
  nn = 0
  for value in omegaRange:
    nn+=1
    runConvert = "%s -trim -geometry 100 %s/%s_%s.dag.png  %s/%s_%s.dag_thumb.png" %(convert,htmlDir,LName,nn,tmpDir,LName,nn)
    tmpWrite.write(runConvert+'<br>')
    runProgram(runConvert)
    runConvert = "%s %s/%s_%s.dag_thumb.png %s/%s_%s.dag_thumb.png" %(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    tmpWrite.write(runConvert+'<br>')
    runProgram(runConvert)
  #################### CREATE FINAL OUTPUT HTML FILE ##########################
  writeOut("Generating final output HTML....")
  nn = 0
  for value in omegaRange:
    nn+=1
    LNameValue = "%s_%s"%(LName,nn)
    makeMap(htmlDir,LNameValue)
  #Get h/s-bond info
  h,s = getHS("%s/%s"%(tmpDir,LName))
  ## Main output page
  #header of document
  permWrite.write('<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">\n')
  permWrite.write('<html>\n<head><title>GeoFold Server</title>\n')
  permWrite.write('<link rev=made href="mailto:bystrc@rpi.edu">\n</head>\n')
  permWrite.write('<body>\n<img src="../../banner.gif">\n')
  #table
  permWrite.write('<table width = "80%">\n<tr><td colspan="2">\n')
  #molscript gif
  permWrite.write('<img src="%s.gif" alt="image not found for %s.pdb. Try re-loading.">\n'%(LName,LName))
  permWrite.write('</td>\n')
  #basic information
  permWrite.write('<td colspan="4">\n<h1>Unfolding...</h1>\n')
  permWrite.write('<h4>%s</h4>\n'%(LName))
  permWrite.write('<h4><a href="%s.pdb">Coordinate file (%s)</a></h4>\n'%(LName,LName))
  permWrite.write('<h4><a href="./%s_1.dag.out">Unfolding graph for %s</a></h4>\n'%(LName,LName))
  #h/s-bond info
  permWrite.write('<h4><a href="./%s.zip" download>Download as zip file</a></h4>\n'%(LName))
  permWrite.write('<h5>Number of H-bonds fond: %s</h5>\n'%(h))
  permWrite.write('<h5>Number of SS-bonds found: %s</h5>\n'%(s))
  #Sequence info
  permWrite.write('<pre>\n')
  try:
    seqFile = open("%s/%s.seq"%(tmpDir,LName))
  except IOError:
    print("Couldn't get sequence information")
  with seqFile:
    for line in seqFile:
      permWrite.write(line)
    seqFile.close()
  permWrite.write('</pre>\n</td>\n</tr>\n<tr>\n')

  #Timecourse info
  nn = 0
  for value in omegaRange:
    nn+=1
    if not thermal:
      permWrite.write('<td><a href="%s_%s.log">timecourse</a> for &omega; = %s<br><a href="%s_%s.tc.png"><img src="%s_%s.tc_thumb.png" width=100 alt="timecourse data missing"></a></td>\n'%(LName,nn,value,LName,nn,LName,nn))
    else:
      permWrite.write('<td><a href="%s_%s.log">timecourse</a> for temp = %s K<br><a href="%s_%s.tc.png"><img src="%s_%s.tc_thumb.png" width=100 alt="timecourse data missing"></a></td>\n'%(LName,nn,value,LName,nn,LName,nn))
  permWrite.write('</tr>\n<tr>\n')

  #Ageplots
  nn = 0
  for value in omegaRange:
    nn+=1
    permWrite.write('<td>Age plot for %s<br><a href="%s_%s.png"><img src="%s_%s_thumb.png" width=100 alt="age plot data missing"></a>\n</td>\n'%(value,LName,nn,LName,nn))
  permWrite.write('</tr>\n<tr>\n')

  #DAGs
  nn = 0
  for value in omegaRange:
    nn+=1
    permWrite.write('<td>Pathways at %s<br><a href="%s_%s.dag.html"><img src="%s_%s.dag_thumb.png" width=100 alt="pathway data missing"></a></td>\n'%(value,LName,nn,LName,nn))
  permWrite.write('</tr>\n<tr>\n')

  #Energy Profiles
  nn = 0
  for value in omegaRange:
    nn+=1
    permWrite.write('<td>Energy profile for %s<br><a href="%s_%s.nrg.html"><img src="%s_%s.nrg_thumb.png" width=100 alt="energy data missing"></a></td>\n'%(value,LName,nn,LName,nn))
  permWrite.write('</tr>\n</table>\n')

  permWrite.write('<h3><A HREF="http://www.bioinfo.rpi.edu/geofold/howtoreadit.htm">How to Read GeoFOLD Output</A></h3>\n')
  #Equilibrium info
  permWrite.write('<hr width="80%">\n')
  if halflife == 1:
    if folding == 1:
      permWrite.write("<h4>Times and concentrations at 1/2 folded")
    else:
      permWrite.write("<h4>Times and concentrations at 1/2 unfolded")
  else:
    if folding == 1:
      permWrite.write("<h4>Equilibrium concentrations and rates under folding conditions")
    else:
      permWrite.write("<h4>Equilibrium concentrations and rates under unfolding conditions")
  #output this info as a table
  permWrite.write('</h4>\n<table width = "80%">\n')
  if not thermal:
    permWrite.write('<tr><td>&Delta;G<sub>solvation</sub></td><td>time(sec)</td><td>final [F]</td><td>final [U]</td><td>final[I]</td><td>half-life (sec)</td></tr>\n')
  else:
    permWrite.write('<tr><td>Temperature, K</td><td>time(sec)</td><td>final [F]</td><td>final [U]</td><td>final[I]</td><td>half-life (sec)</td></tr>\n')
  plotFile= open("%s/%s.plot"%(tmpDir,LName),'w+')
  nn = 0
  for value in omegaRange:
    nn+=1
    logFile = open("%s/%s_%s.log" %(tmpDir,LName,nn),'r')
    permWrite.write('<tr><td>%s </td>'%(value))
    #find the last TIMECOURSE line in logFile and writeout its contents
    timecourses = []
    for line in logFile:
      line = line.split()
      if len(line)!=0:
        if line[0]=='TIMECOURSE':
          timecourses.append(line)
    line = timecourses[len(timecourses)-1]
    logFile.close()
    # permWrite.write('<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>'%(line[1],line[2],line[3],line[4],line[5]))
    permWrite.write('<td>%s</td><td>%s</td><td>%s</td><td>%s</td>'%(line[1],line[2],line[3],line[4]))
    #calculate ln(ku) stuff for plot
    i = 0
    #find the last timecourse where the unfolded state is < 50 if unfolding, folded state if otherwise
    if folding != 1:
      while i<len(timecourses) and float(timecourses[i][3])<50:
        i+=1
    else:
      while i < len(timecourses) and float(timecourses[i][2])<50:
        i+=1
    i-=1
    # hl = math.log(math.log(2)/float(timecourses[i][1]))
    try:
      lnku = math.log(math.log(2)/float(timecourses[i][5]))
    except ZeroDivisionError:
      lnku = 0
      print("\n\n\nDivided by zero\n\n\n")
    permWrite.write('<td>%s</td>'%(timecourses[i][5]))
    if lnku != 0:
      if thermal:
        print(1./float(value))
        plotFile.write("%s %s\n"%(1./float(value),lnku))
      else:
        plotFile.write("%s %s\n"%(value,lnku))
  plotFile.close()
  permWrite.write('</table>\n')

  print("============= FIT_POLY =============")
  writeOut("\n============= FIT_POLY =============\n")
  tmpWrite.write('<pre>\n')
  ##INSERT ADJUSTMENT STUFF HERE
  try:
    if not thermal:
      fit("%s/%s.plot"%(tmpDir,LName),"%s/%s_fit.plot"%(tmpDir,LName))
    else:
      commands.getstatusoutput("cp %s/%s.plot %s/%s_fit.plot"%(tmpDir,LName,tmpDir,LName))
  except Exception:
    commands.getstatusoutput("cp %s/%s.plot %s/%s_fit.plot"%(tmpDir,LName,tmpDir,LName))
  fit_poly = "%s/xfit_poly %s/%s_fit.plot 1 %s/%s.fit >> %s" %(gDir,tmpDir,LName,tmpDir,LName,htmlTmp)
  runProgram(fit_poly)
  #Predicting half-life in water (omega = 30.) 95?
  fitFile = open("%s/%s.fit"%(tmpDir,LName),'r')
  line = fitFile.readline()
  lines = []
  lsq = False
  sl = 1
  ic = 0
  while not lsq:
    line = fitFile.readline()
    if "Least-squares" in line:
      lsq = True
      lines.append(fitFile.readline())
      lines.append(fitFile.readline())
      lines.append(fitFile.readline())
      print lines
      ic = lines[1]
      sl = lines[2]
      ic = ic.split()
      sl = sl.split()
      ic = float(ic[1])
      sl = float(sl[1])
    lnku = (sl*(1./float(wat)))+ic
  hl = math.log(2)/math.exp(lnku)
  fitFile.close()
  if not thermal:
    permWrite.write('<h5>Projected half-life of unfolding in pure water (solvation weight = %s): %s seconds</h5><br>'%(wat,hl))
  else:
    permWrite.write('<h5>Projected half-life of unfolding in pure water (temperature = %s K): %s seconds</h5><br>'%(wat,hl))

  #final graph outputs
  permWrite.write('<table>\n<tr>\n<td>\n')
  if not thermal:
    permWrite.write('<h4>ln(k<sub>u</sub>) vs &omega; for %s </h4>'%(LName))
  else:
    permWrite.write('<h4>ln(k<sub>u</sub>) vs temp for %s </h4>'%(LName))
  permWrite.write('<a href="%s.gp.png"><img src="%s.gp_thumb.png" alt="data missing"></a><br>\n</td>\n'%(LName,LName))
  permWrite.write('<td>\n<h4>Superimposed energy profiles for %s </h4>'%(LName))
  permWrite.write('<a href="%s_all.nrg.png"><img src="%s_all.nrg_thumb.png" alt="data missing"></a><br>\n</td></tr></table>'%(LName,LName))
  #Completion time
  permWrite.write('<br>Completed %s<br>'%(time.strftime("%c")))
  permWrite.write('<p>Modify and do over:\n')

  #### Create gnuplot image
  print("============= GNUPLOT =============")
  writeOut("============= GNUPLOT =============\n")
  nn = 0
  for value in omegaRange:
    nn+=1
    plotTimeCourse("%s/%s"%(tmpDir,LName),nn)
    runConvert = '%s -trim -geometry 100 %s/%s_%s.tc.ps %s/%s_%s.tc_thumb.png'%(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    runProgram(runConvert)
    runConvert = '%s -trim %s/%s_%s.tc.ps %s/%s_%s.tc.png'%(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    runProgram(runConvert)
    status,output = writeEnergyProfile(tmpDir,htmlDir,LName, nn)
    if status != 0:
      sys.exit(output)
    picLink(htmlDir,LName,nn)
    runConvert = '%s -trim -geometry 100 %s/%s_%s.nrg.ps %s/%s_%s.nrg_thumb.png'%(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    runProgram(runConvert)
    runConvert = '%s -trim %s/%s_%s.nrg.ps %s/%s_%s.nrg.png'%(convert,tmpDir,LName,nn,htmlDir,LName,nn)
    runProgram(runConvert)
  status,output = energyProfileAll("%s/%s"%(tmpDir,LName),omegaRange)
  if status != 0:
    print(status)
    print(output)
  else:
    runConvert = '%s -trim -geometry 200 %s/%s_all.nrg.ps %s/%s_all.nrg_thumb.png'%(convert,tmpDir,LName,htmlDir,LName)
    runProgram(runConvert)
    runConvert = '%s -trim %s/%s_all.nrg.ps %s/%s_all.nrg.png'%(convert,tmpDir,LName,htmlDir,LName)
    runProgram(runConvert)
  status,output = createGnuplot("%s/%s"%(tmpDir,LName),wat)
  #status = 1
  if status != 0:
    print(status)
    print(output)
  else:
    runConvert = "%s -trim %s/%s.gp.ps %s/%s.gp.png" %(convert,tmpDir,LName,tmpDir,LName)
    runProgram(runConvert)
    runConvert = "%s %s/%s.gp.png %s/%s.gp.png" %(convert,tmpDir,LName,htmlDir,LName)
    runProgram(runConvert)
    runConvert = "%s -trim -geometry 200 %s/%s.gp.ps %s/%s.gp_thumb.png" %(convert,tmpDir,LName,tmpDir,LName)
    runProgram(runConvert)
    runConvert = "%s %s/%s.gp_thumb.png %s/%s.gp_thumb.png" %(convert,tmpDir,LName,htmlDir,LName)
    runProgram(runConvert)
  ######Create "do over" button
  createForm(permWrite,paramFilename)
  #####BACK button
  permWrite.write('<h3><a href="../../geofold.php">Back to GeoFold server</a></h3><br>\n')
  permWrite.write('</body></html>\n')
  permWrite.close()
  #generate molscript
  paramFile = open(paramFilename,'r')
  status, script = findParam(paramFile,'MOLSCRIPT')
  paramFile.close()
  if script == '1' or status != 0:
    writeOut("========== MOLSCRIPT ==========\n")
    status,output = makeMolScript('%s/%s'%(htmlDir,LName))
    if status != 0:
      print(status)
      print(output)
  cp = "cp -v %s %s" %(htmlPerm,htmlOut)
  status,output=commands.getstatusoutput(cp)
  tmpWrite.close()
  nn=0
  for value in omegaRange:
    nn+=1
    commands.getstatusoutput("cp %s/%s_%s.dag.out %s/%s_%s.dag.out"%(tmpDir,LName,nn,htmlDir,LName,nn))
    commands.getstatusoutput("cp %s/%s_%s.log %s/%s_%s.log" %(tmpDir,LName,nn,htmlDir,LName,nn))
  makeZip(htmlDir,LName)
