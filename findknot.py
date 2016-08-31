# A program designed to detect knots in protein
# Based on W Taylor Nature 406(6798):916-919

import argparse
import sys
import os
import os.path

#Reads a PDB file and outputs CA coordinates as a list
def readPDB(pdbFile):
  #output matrix
  calphas = [[]]
  #open file
  try:
    infile = open(pdbFile,'r')
  except IOError:
    sys.exit("Failed to open file %s"%(pdbFile))
  #parse the file
  for line in infile:
    line = line.split()
    if(line[0]=='ATOM'):
      if(line[2]=='CA'):
        tmp = [line[6],line[7],line[8]]
        calphas.append(tmp)
  return calphas
  
#runs the smoothing algorithm
def smooth(chain):
  None
 
#main program

#argument parser
parser = argparse.ArgumentParser(description='Check PDB files for knots.')
parser.add_argument('pdbFile',metavar='p',nargs='+',help='PDB file(s) to be processed')
parser.add_argument('-v',dest='verbose',action='store_true',default=False,help='Run this program verbosely.')
args = parser.parse_args()

#read PDB files into CA matrices
pdbcas=[[[]]]
for file in args.pdbFile:
  pdbcas.append(readPDB(file))
  if args.verbose:
    for ca in pdbcas[len(pdbcas)-1]:
      print(ca)
    
#run smoothing algorithm