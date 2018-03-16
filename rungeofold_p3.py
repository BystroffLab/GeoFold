
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
B. Walcott 11 June 2014 
B. Walcott 16 March 2018"""

# imports
import os
import os.path
import sys
import subprocess
import time
import math
import argparse
import asyncio

# async stuff
# https://kevinmccarthy.org/2016/07/25/streaming-subprocess-stdin-and-stdout-with-asyncio-in-python/
async def _read_stream(stream, cb):  
    while True:
        line = await stream.readline()
        if line:
            cb(line)
        else:
            break

async def _stream_subprocess(cmd, stdout_cb, stderr_cb):  
    process = await asyncio.create_subprocess_exec(*cmd,
            stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)

    await asyncio.wait([
        _read_stream(process.stdout, stdout_cb),
        _read_stream(process.stderr, stderr_cb)
    ])
    return await process.wait()


def execute(cmd, stdout_cb, stderr_cb):  
    loop = asyncio.get_event_loop()
    rc = loop.run_until_complete(
        _stream_subprocess(
            cmd,
            stdout_cb,
            stderr_cb,
    ))
    loop.close()
    return rc

# https://stackoverflow.com/questions/41283595/how-to-redirect-python-subprocess-stderr-and-stdout-to-multiple-files
class MultiOut(object):
    def __init__(self, *args):
        self.handles = args

    def write(self, s):
        for f in self.handles:
            f.write(s)

def readParameters(parFile):
    output = {}
    try:
        parIn = open(parFile)
    except IOError:
        sys.exit("parameter file %s could not be opened"%(parFile))
    lines = [line for line in parIn]
    parIn.close()
    # pdbcode
    line = [x for x in lines if x.split()[0] == 'PDBCODE']
    if len(line) == 0: 
        sys.exit("ERROR: Parameters file missing PDBCODE")
    else:
        output['pdbcode'] = line[-1].split()[1]
    # thermal
    line = [x for x in lines if x.split()[0] == 'THERMAL']
    if len(line) == 0: 
        output['thermal'] = False
    else:
        output['thermal'] = line[-1].split()[1] == '1'
    # keyword
    line = [x for x in lines if x.split()[0] == 'KEYWORD']
    if len(line) == 0:
        output['keyword'] = ''
    else:
        output['keyword'] = line[-1].split()[1]
    # chain
    line = [x for x in lines if x.split()[0] == 'CHAIN']
    if len(line) == 0:
        sys.exit("ERROR: Parameters file missing CHAIN")
    else:
        output['chain'] = line[-1].split()[1]
    # lname
    line = [x for x in lines if x.split()[0] == 'LNAME']
    if len(line) == 0:
        output['lname'] = "%s.%s"%(output['keyword'],output['pdbcode'])
    else:
        output['lname'] = line[-1].split()[1]
    # strip lname of suspicious characters
    output['lname'] = output['lname'].strip("/*{}[]!@#$%^&();:<>,?\\~\"\'")
    # email
    line = [x for x in lines if x.split()[0] == 'EMAIL']
    if len(line) == 0:
        output['email'] = None
    else:
        output['email'] = line[-1].split()[1]
    # oname
    line = [x for x in lines if x.split()[0] == 'ONAME']
    if len(line) == 0:
        output['rerun'] = False
    else:
        output['rerun'] = True
        output['oname'] = output['lname']
    output['uname'] = '%s%s%s'%(output['pdbcode'],output['chain'],output['lname'])
    # reducing
    line = [x for x in lines if x.split()[0] == 'REDUCING']
    if len(line) == 0: 
        output['reducing'] = False
    else:
        output['reducing'] = line[-1].split()[1] == '1'
    # debug
    line = [x for x in lines if x.split()[0] == 'DEBUG']
    if len(line) == 0: 
        output['debug'] = False
    else:
        output['debug'] = line[-1].split()[1] == '1'
    # barrelmoves
    line = [x for x in lines if x.split()[0] == 'BARRELMOVES']
    if len(line) == 0: 
        output['barrelmoves'] = False
    else:
        output['barrelmoves'] = line[-1].split()[1] == '1'
    # OMEGA
    line = [x for x in lines if x.split()[0] == 'OMEGA']
    if len(line) == 0:
        output['omega'] = None
    else:
        output['omega'] = line[-1].split()[1]
    # ORANGE
    line = [x for x in lines if x.split()[0] == 'ORANGE']
    if len(line) == 0:
        if output['omega'] is None:
            sys.exit("ERROR: ORANGE and OMEGA missing in parameters file")
        else:
            output['orange'] = [output['omega']]
    else:
        output['orange'] = line[-1].split()[1:]
    # REDO
    line = [x for x in lines if x.split()[0] == 'REDO']
    if len(line) == 0:
        output['redo'] = None
    else:
        output['redo'] = line[-1].split()[1]
    
    return output

def main():
    # parser
    parser = argparse.ArgumentParser()
    # config file
    parser.add_argument("-conf",default="default.conf",help="config file to use for GeoFold")
    # parameters file
    parser.add_argument("-par",help="parameters file for this run of GeoFold")
    args = parser.parse_args()
    # commandline parameters
    
    # parameters
    parameters = readParameters(args.par)
    # getchain
    # renumber
    # 3to1
    # pdb2cij
    # pdb2hb
    # contactmask
    # pdb2seams2
    # splitseams
    # voidmask
    # geofold
    # unfoldsim
    # pathway2ps
    # convert
    # maxtraffic
    # dot
    # convert
    # fit_poly
    # gnuplot
    # molscript
    

if __name__ == '__main__': main()