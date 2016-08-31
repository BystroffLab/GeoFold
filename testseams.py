"""testing seams"""
import commands
import os

gDir = "./"
LNamePDB = "tmp/2b3p.4602_3.pdb"
tmpDir = "./tmp"
LName = "wzs"

print("============= PDB2SEAMS (extract seams) =============")
pdb2seams = "%s/seams/xpdb2seams %s %s/%s.sas" %(gDir,LNamePDB,tmpDir,LName)
os.environ['GDIR']=gDir
os.environ['TMPDIR']=tmpDir
status,output = commands.getstatusoutput(pdb2seams)
print("%s: %s"%(status,output))
