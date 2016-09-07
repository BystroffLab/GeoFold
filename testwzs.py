"""runs xcontactmask for the testing of the wzs solvation model"""
import time
import os
import commands
def main():
  print("============= CONTACTMASK =============")
  writeTime = "Time before running CONTACTMASK "+time.strftime("%c")
  print(writeTime)
  maskerDir = "masker"
  LNamePDB = "tmp/2b3p.4602_3.pdb"
  tmpDir = "tmp"
  LName = "wzs"
  ##masker setup
  os.environ['MASKERDIR']=maskerDir
  os.environ['DESIGN_HOME']=os.environ['MASKERDIR']
  os.environ['MSIZE']="512"
  os.environ['DTHETA']="9.0"
  os.environ['MASKLIB']="%s/%s.mask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['VMASK']="%s/%s.vmask" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['VTRIANGLE']="%s/%s.tri" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['MASKTEMPLATE']="%s/%s.sort" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['MASKDAT']="%s/%s.dat" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['ATOMLIB']="%s/atoms.lib" %(os.environ['MASKERDIR'])
  os.environ['COUNTBITS']="%s/cbits.bin" %(os.environ['MASKERDIR'])
  os.environ['COLLARS']="%s/collars%s.bin" %(os.environ['MASKERDIR'],os.environ['MSIZE'])
  os.environ['SLOPES']="%s/slopes%s.dat" %(os.environ['MASKERDIR'],os.environ['DTHETA'])
  os.environ['FORCEFIELD']="%s/gromacs_ff.prm" %(os.environ['MASKERDIR'])
  
  contactmask = "%s/xcontactmask %s %s.sas 1.4" %(maskerDir,LNamePDB,LName)
  status,output = commands.getstatusoutput(contactmask)
  print("%s: %s"%(status,output))
  writeTime = "Time after running CONTACTMASK "+time.strftime("%c")
  print(writeTime)
  print("%s: %s"%(status,output))
    

if __name__=="__main__":
  main()
