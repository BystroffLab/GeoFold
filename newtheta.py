import commands
pdbs = ['1e65','1aps','1avz','1fkf','1imq','1jyg','1k0s','1m9s','1n88','1nti','1o6x','1rfa','1ris','1rlq','1spr','1ubq','1urn','2ci2','2ptl','3gb1']
for pdb in pdbs:
  cmd = "python rungeofold.py spsa/%sp.tmp.par walcob"%(pdb)
  status,output = commands.getstatusoutput(cmd)
  print("%s: %s"%(status,output))
