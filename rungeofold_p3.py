
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

def main():
    parser = argparse.ArgumentParser()
    # commandline parameters
    
    # parameters
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