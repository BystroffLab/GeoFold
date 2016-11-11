# Wenyin San
# Feb 2016
import re
import os
import os.path
import sys
import commands
import time
import traceback
import math
import time
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
"""
readConf(confFile) method from original version of GeoFold
bach_sanw.conf saved at the home directory

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
#added by SAN
#B. Walcott publicly accessible results
thisDir   = os.getcwd()
baseDir   = thisDir
gDir      = baseDir
bDir      = baseDir
maskerDir = gDir+'/masker'
pdbDir    = gDir+'/pdbs'
tmpDir    = baseDir+"/tmp"
logDir    = baseDir+"/log"
htmlDir   = baseDir+"/html"
jobDir    = baseDir+"/job"
paramTemplate=gDir+'/parameters'
#baseURL='http://www.bioinfo.rpi.edu/bystrc/geofold'
outputURL='output'

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
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def writeOut(x):
  global htmlOut
  try:
    wout = open(htmlOut,'a')
  except IOError:
    sys.exit("Couldn't open file: %s"%(htmlOut))
  wout.write(x)
  wout.close()

def runProgram(command):
  global tmpWrite
  global writeOut
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

# initialize status and status iterator
def readArg(argFile):
    global htmlOut
    try:
    	argRead = open(argFile, 'r')
    except IOError:
    	sys.exit("Couldn't open file %s"%(argFile))
    
    """
    STATUSES = ('omegaRange:', 'LName:', 'paramFilename:', 'thermal:', 'doIt:', 'debug:', 'gDir:', 'htmlTmp:', 'htmlOut:')
    statuses = iter(STATUSES)
    states = {}
    """
    parList = {}
    omegaRange = [] 
    for line in argRead:
        if line[0] != "#":
            line = line.split()
            if len(line) > 2:
                if line[0] == 'omegaRange':
                    curr_line = line.split(' ')
                    for i in range(1, len(curr_line)):
                        omegaRange.append(int(curr_line[i]))
                else:
                    line[1] = " ".join(line[1:])
            if len(line) == 2:
                parList[line[0]] = line[1]
    argRead.close()
    """
    for line in argRead:
        # match all words in each line
        columns = re.findall("[^\s]+", line.strip())

	if columns[0] in ['omegaRange:', 'LName:', 'paramFilename:', 'thermal:', 'doIt:', 'debug:', 'gDir:', 'htmlTmp:', 'htmlOut:']:
            status = next(statuses)
            continue

        elif status == 'omegaRange:':
    	    omegaRange.append(float(columns[0]))
            continue

        elif status == 'LName:':
	    LName = columns[0]
    	    continue

        elif status == 'paramFilename:':
	    paramFilename = columns[0]
	    continue

    	elif status == 'thermal:':
            thermal = str2bool(columns[0])
	    #thermal = bool(columns[0])
	    continue

    	elif status == 'doIt:':
	    doIt = int(columns[0])
	    continue
 
        elif status == 'debug:':
	    doIt = str2bool(columns[0])
	    continue
        
        elif status == 'gDir:':
	    gDir = columns[0]
	    continue

        elif status == 'htmlTmp:':
	    htmlTmp = columns[0]
	    continue

        elif status == 'htmlOut:':
	    htmlOut = columns[0]
	    continue
    """
    LName = parList['LName']
    paramFilename = parList['paramFilename']
    thermal = parList['thermal']
    doIt = parList['doIt']
    debug = parList['debug']
    gDir = parList['gDir']
    htmlTmp = parList['htmlTmp']
    htmlOut = parList['htmlOut']
    return omegaRange, LName, paramFilename, thermal, doIt, debug, gDir, htmlTmp, htmlOut

##++++++++++++++++++++++++++++++++++++++++++
##++++++++++++++++++++++++++++++++++++++++++

# chunk size calculated of omegaRange
#
#chunk_size = math.ceil(float(omegaRange/n))

"""
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
"""
def main():
    global tmpWrite
    global writeOut
    global htmlOut
    global htmlTmp
    global debug
    tmpDir = sys.argv[1]
    tmp_par_file = "%s/tmp_par_file.txt" %(tmpDir)
    omegaRange = []
    LName = " "
    paramFilename = " "
    thermal = False
    doIt = 0
    debug = False
    gDir = " "
    htmlTmp = " "
    htmlOut = " "
    omegaRange, LName, paramFilename, thermal, doIt, debug, gDir, htmlTmp, htmlOut = readArg(tmp_par_file)

    if rank == 0:
	"""
	Wait for gather
	"""
        print "Root rank started"


    else:
        
        if rank <= len(omegaRange):
            tmpWrite = open(htmlTmp,'a')
            value = omegaRange[rank-1]
            if not thermal:
                sed = "sed -e \"s/^OMEGA .*/OMEGA %s/\" %s > %s.1" %(value,paramFilename,paramFilename)
            else:
                sed = "sed -e \"s/^TEMPERATURE .*/TEMPERATURE %s/\" %s > %s.1"%(value, paramFilename,paramFilename)
            status,output=commands.getstatusoutput(sed)
            logFile = "%s/%s_%s.log" %(tmpDir,LName,rank)
            #print("Thermal is %s:" %(thermal))
            if not thermal:
                print("============= run %s omega = %s =============" %(rank,value))
                tmpWrite.write("============= run %s omega = %s =============<br>" %(rank,value))
                writeOut("============= run %s omega = %s =============\n" %(rank,value))
            else:
                print("============= run %s temp = %s K =============" %(rank,value))
                tmpWrite.write("============= run %s temp = %s K =============<br>" %(rank,value))
                writeOut("============= run %s temp = %s K =============\n" %(rank,value))
            writeTime = "Time before running UNFOLDSIM "+time.strftime("%c") +'<br>'
            tmpWrite.write(writeTime)
            writeOut(writeTime)
            unfoldsim = "%s/xunfoldsim %s/%s_%s.dag %s.1 > %s" %(gDir,tmpDir,LName,rank,paramFilename,logFile)
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

if __name__ == "__main__":
    main()
	  