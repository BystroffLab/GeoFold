# Get MPIRUN
# Wenyin San
# March 2016
from distutils import spawn

def get_mpirun( default="mpirun -np" ):
  # if it exists environment variable named MPIRUN, return MPIRUN; 
  # else return a default mpirun command

  MPIRUN = None
  if spawn.find_executable("mpirun") is not None:
  	MPIRUN  = "mpirun -np"

  elif spawn.find_executable("mpiexec") is not None:
  	MPIRUN  = "mpiexec"
    
  # return a default mpirun command	
  if MPIRUN is None:
    MPIRUN = default
    print "getmpirun: using the default ", default

  return MPIRUN