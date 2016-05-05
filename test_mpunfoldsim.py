#!/usr/bin/env python
# Runs the example and produce output
# Execute using mpirun
# Wenyin San
# March 2016
# Case for saving files
from shell import shell
from launch import launch
from runmpi import get_mpirun
from get_cpu_num import get_cpu_num
import os
import sys
this_dir = "/bach1"+os.getcwd()
def test_mpunfoldsim():
        MPIRUN = get_mpirun()
        #nproc = get_cpu_num()
        # shell("make > make.log")
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
        command = this_dir + "/py_mpfiles.py" #+ " " + "argFile"
	status, output = launch(command, runcmd=MPIRUN, hostfile = None, nproc=5, pipe=False)
	# Save the output in a log file
	f = open("run.log", "w")
	f.seek(0)
	f.write(str(output))
	f.close()

if __name__ == '__main__':
        test_mpunfoldsim()




	
