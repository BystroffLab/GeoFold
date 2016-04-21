#!/usr/bin/env python
# Runs the example and produce output
# Execute using mpirun
# Wenyin San
# March 2016
from shell import shell
from launch import launch
from runmpi import get_mpirun
from get_cpu_num import get_cpu_num
import os
import sys

def test_mpunfoldsim(argFile):
	MPIRUN = get_mpirun()
	nproc = get_cpu_num()# shell("make > make.log")
	this_dir = os.getcwd()
	"""
		omegaRange = omegaRange
	    tmpDir = arg['tmpDir']
	    LName = arg['LName']
	    nn = arg['nn'] 
	    paramFilename = arg['paramFilename'] 
	    thermal = arg['thermal']
	    htmlTmp = arg['htmlTmp']
	    doIT = arg['doIT']
    """
	command = this_dir + "/hosts" + this_dir + "/py_hello_world.py" #+ " " + "argFile"
	status, output = launch(command, runcmd=MPIRUN, nproc=nproc, pipe=False)
	# Save the output in a log file
	f = open("run.log", "w")
	f.seek(0)
	f.write(str(output))
	f.close()




	