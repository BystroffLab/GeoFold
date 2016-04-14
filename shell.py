# Run shell command
# Wenyin San
# March 2016

import sys
import os
import time
import subprocess

try:
    from subprocess import STDOUT, PIPE, call, Popen
    lib = "call"
except:
    # Use os.system (depreciated)
    from os import popen4, system
    lib = "system"
    
def shell(command, pipe=False):
    start = time.time()
    output = None
    status = 0
    if lib == "system":
        if pipe:
            handle = popen4(command)
            output = handle[1].read()
        else:
            status = system(command)
    else:
        # lib = "call"
        if pipe:
            child = Popen(command, stderr=STDOUT, stdout=PIPE, shell=True)
            output = child.stdout.read()
            status = child.returncode
        else:
            status = call(command, shell=True)
    end = time.time()
    return status, output