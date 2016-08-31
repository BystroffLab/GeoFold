#!/usr/bin/python

"""
Used to confirm the detection of beta sheets
"""

import os, sys

#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------
USAGE = "\
Find residues in beta strands using Stride \n\n\
USAGE: seams_stride.py <input PDB filename> <output filename>"

def main (args):
    if len (args) < 2:
        print USAGE
        sys.exit (0)

    pdb = args [1]
    try:
      stridePath = args [2]
    except IndexError:
      print("Error reading args : %s" %args)
      sys.exit(1)
    betas =  getBetasStride (pdb, stridePath)
    for b in betas:
        print b
    
#-------------------------------------------------------------------------------
# getBetasStride: 
#-------------------------------------------------------------------------------
def getBetasStride (pdbFilename, stridePath):
    sys.stderr.write("getBetasStride\n")
    strideFilename = runStride (pdbFilename, stridePath)

    betas = processStride (strideFilename)

    return betas

#-------------------------------------------------------------------------------
# writeBetas:
#-------------------------------------------------------------------------------
def writeBetas (betas, outFilename):
    sys.stderr.write("writeBetas\n")
    outFile = open (outFilename, "w")

    outFile.write (str (len (betas)))
    for i in betas:
        outFile.write (" " + str (i))

    outFile.close()
    
#-------------------------------------------------------------------------------
# runStride
#-------------------------------------------------------------------------------
def runStride (pdbFilename, stridePath):
    import os
    sys.stderr.write("runStride\n")
    SEAMS_PATH = os.environ ["GDIR"] + "/seams"
    gDir = os.environ["GDIR"]+'/'
    stride = SEAMS_PATH + "/stride "
    stemName = name (pdbFilename)
    strideFilename = stridePath + stemName + "-tmp"  + ".stride"

    #log (["!!!!!!!!!!!!!", "runStride", pdbFilename, strideFilename, currentDir])
    #strValue = runProgram (["stride", pdbFilename, strideFilename], currentDir)


    cmm = stride + pdbFilename + " > " + strideFilename
    sys.stderr.write(cmm+'\n')
    
#"/home/garrel/svn/geofold/seams/stride /home/garrel/svn/geofold/tmp/rpl14.pdb"
    #os.system (cmm)
    import commands
    status, output = commands.getstatusoutput(cmm)
    sys.stderr.write("%s: %s\n"%(status,output))

    #import subprocess
    #output = subprocess.Popen (["stride", pdbFilename], cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]
    #output = subprocess.Popen ([SEAMS_PATH, "stride", pdbFilename], cwd="/tmp", stdout=subprocess.PIPE, shell=True).communicate()[0]

    #outFile = open (strideFilename, "w")
    #for line in output:
    #    outFile.write (line)
    # outFile.close()

    return strideFilename

#-------------------------------------------------------------------------------
# processStride
#-------------------------------------------------------------------------------
def processStride (strideFilename):
    sys.stderr.write("processStride\n")
    strideFile = open (strideFilename)

    # Skips header lines until reach the values of aminos and structures and others

    seqAABeta = []

    for line in strideFile:
        
        lst = line.split ()
        if lst[0] != "ASG": 
          continue

        try:
          residue = int (lst[3])
        except LookupError:
          try:
            residue = int(lst[3][0:len(lst[3])])
          except LookupError:
            print("Error parsing lst[3]: %s" %lst)
        except ValueError:
          continue

        try:
          structure = lst [5]
        except LookupError:
          print("Error parsing lst[5]: %s" %lst)
          sys.exit(1)
        except ValueError:
          print("Error invalid value lst[5]: %s" %lst[5])
          sys.exit(1)

        if structure in ["E", "B"]:
            sys.stderr.write("Strand found!\n")
            seqAABeta.append (residue)

    return seqAABeta

#-------------------------------------------------------------------------------
# runProgram: call a external program that returns a value running on "workingDir"
#-------------------------------------------------------------------------------
def runProgram (listOfParams, workingDir):
    import subprocess
    sys.stderr.write("runProgram\n")
    output = subprocess.Popen (listOfParams[:2], cwd=workingDir, stdout=subprocess.PIPE).communicate()[0]

    return output
##################################################################
# Delete temporal files created during the process
##################################################################
def clean (listOfFiles):
    for file in listOfFiles:
        os.system ("rm " + file)

###############################################################################
# Returns the stem BASENAME  of a relative or full name with some extension
###############################################################################
def name (namefile):
    baseName = os.path.basename (namefile)
    newNamefile = baseName.split(".") [0]

    return newNamefile

###############################################################################
# Print a message to the error output stream
###############################################################################
def log (message):
    string=""
    for i in message:
        string += str (i) + " "
    
    #sys.stderr.write (">>>> SecStr: " + string+"\n")

###############################################################################
# Call the main function 
###############################################################################
if __name__ == "__main__": 

    main (sys.argv)

