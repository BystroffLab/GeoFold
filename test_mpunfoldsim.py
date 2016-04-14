#!/usr/bin/env python
# Runs the example and produce output
# Execute using mpirun
# Wenyin San
# March 2016
if __name__ == '__main__':
	from shell import shell
	from launch import launch
	from runmpi import get_mpirun
	from get_cpu_num import get_cpu_num
	import os
	import sys

	MPIRUN = get_mpirun()
	nproc = get_cpu_num()
	# shell("make > make.log")
	status, output = launch("./mpunfoldsim.py", runcmd=MPIRUN, nproc=nproc, pipe=False)
	# Save the output in a log file
	f = open("run.log", "w")
	f.seek(0)
	f.write(str(output))
	f.close()