# Launch parallellization
# If run on a known machine, will generate a run script and submit a job for execution. 
# Uses 'mpirun' by default
# Wenyin San
# March 2016
import os
import sys
from shell import shell
import get_cpu_num

this_dir = "/bach1"+os.getcwd()

def launch(command, runcmd = "mpirun -np", hostfile = None, nproc=None, output=None, pipe=False, verbose=True):
    """
    Launch MPI jobs
    status, output = launch(command, nproc, output=None)
    runcmd     Command for running parallel job; defaults to "mpirun -np"    
    command    The command to run (string)
    nproc      Number of processors (integer)
    output     Optional name of file for output
    """
    # Determine number of CPUs on this machine
    if nproc == None:
        nproc = get_cpu_num()
    if hostfile == None:
        hostfile = this_dir+"/hosts" 

    #cmd = runcmd + " " + str(nproc) +' --hostfile ' + hostfile + " python " + command
    cmd = runcmd + " " + str(nproc) + " python " + command
    if output != None:
        cmd = cmd + " > "+output

    if verbose == True:
         print cmd    

    return shell(cmd, pipe=pipe)
