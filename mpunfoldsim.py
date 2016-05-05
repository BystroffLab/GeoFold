# Wenyin San
# Feb 2016
import os
import os.path
import sys
import commands
import time
import traceback
import math
import numpy as np
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
from georansac import fit


# set comm, rank & size

try :
	comm  = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	name = MPI.Get_processor_name()

except Exception as e :
    traceback.print_exc ()
    print "error : %s" % e
    sys.exit (1)


conf = "bach_sanw.conf"

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

"""
It processes a given file and returns a HMM object ready to roll.
"""
# initialize status and status iterator

try:
    argWrite = open(argFile, 'w+')
  except IOError:
    paramFile.close()
    sys.exit("Couldn't open file %s"%(argFile))
  argFile.write('omegaRange:\n')
  for value in omegaRange:
  	argFile.write("%s\n" %(value))
  argFile.write('LName:\n')
  argFile.write("%s\n" %(LName))
  argFile.write('paramFilename:\n')
  argFile.write("%s\n" %(paramFilename))
  argFile.write('thermal:\n')
  argFile.write("%s\n" %(thermal))
  argFile.write('doIT:\n')
  argFile.write("%s\n" %(doIT))
  argFile.close() 

def readArg(argFile)
    try:
    	argRead = open(argFile, 'r')
    except IOError:
    	argFile.close()
    	sys.exit("Couldn't open file %s"%(argFile))

    output = {}
	STATUSES = ('omegaRange:', 'LName:', 'paramFilename:', 'thermal:', 'doIT:')
	statuses = iter(STATUSES)
	states = {}
    omegaRange = []
    for line in argRead:
        # match all words in each line
        columns = re.findall("[^\s]+", line.strip())

        if columns[0] in ['omegaRange:', 'LName:', 'paramFilename:', 'thermal:', 'doIT:']:
            status = next(statuses)
            continue

        elif status == 'omegaRange:':
    		omegaRange.append(float(columns[0]))
            continue

        elif status == 'LName:':
    		Lname = columns[0]
    		continue

        elif status == 'paramFilename:':
    		paramFilename = columns[0]
    		continue

    	elif status == 'thermal:':
    		thermal = bool(columns[0])
    		continue

    	elif status == 'doIT:':
    		doIT = int(columns[0])
    		continue

	return omegaRange, LName, paramFilename, thermal, doIT

##++++++++++++++++++++++++++++++++++++++++++
##++++++++++++++++++++++++++++++++++++++++++
argFile = "%s/argFile.txt" %(tmpDir)

omegaRange, LName, paramFilename, thermal, doIT = readArg(argFile)

oR_size = math.ceil(float(omegaRange/n))

def Range():
	if rank != n-1:
		# 0 < rank < n-1
		start = (rank-1)*size
		end = rank*size
	else:
		# rank == n-1
		start = (rank-1)*size
		end = len(omegaRange)
	return start, end

if rank == 0:
	print(doIt == 3)
	print(doIt)
	print("============= UNFOLDSIM =============")
	tmpWrite.write("Unfolding %s%s<br>" %(pdbCode,chain))
	tmpWrite.write("============= UNFOLDSIM =============<br>")
	writeOut("============= UNFOLDSIM =============\n")

else： 
	start, end = Range()
	nn = start
    for value in range(start, end):
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

	  