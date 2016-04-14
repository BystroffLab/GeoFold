#!/usr/bin/python
import os
import sys

import subprocess

programs = {
    'hello_world': ('hello_world', 4)
    #,
    #'py_hello_world': ('py_hello_world', 4)
}

if __name__ == '__main__':
    if len(sys.argv)>1:
        program_name = sys.argv[1] 
    else:
        program_name = None

    if not program_name in programs:
        print 'Enter program name to run. Possible programs are: {0}'.format(programs.keys())

    else:
        # Compile before running
        with open(os.devnull, 'wb') as devnull:
            subprocess.call(
                ['cd ./{0} && make'.format(programs[program_name][0])],
                stdout=devnull, stderr=subprocess.STDOUT, shell=True)

        mpirun = os.environ.get('MPIRUN', 'mpirun')
        if not os.environ.get('MPI_HOSTS'):
            hosts = '' 
        else:
            hosts = '-f {0}'.format(os.environ.get('MPI_HOSTS'))

        sys_call = '{0} -n {1} {2} ./{3}/{4}'.format(mpirun, programs[program_name][1], hosts, programs[program_name][0], program_name)

        if len(programs[program_name]) > 2:
            sys_call = '{0} {1}'.format(sys_call, ' '.join(programs[program_name][2]))

        print sys_call
        subprocess.call([sys_call], shell=True)