import commands
pdbs = ['1C9O','1CSP','1G6P','1MJC','1M9ST','1PNJ','1SHF','1SHG','1SRL','2VKN','1JO8']

for pdb in pdbs:
  cmd = "python rungeofold.py seamtest/%s.par walcob & python rungeofold.py seamtest/%ss.par walcob"%(pdb,pdb)
  status,output = commands.getstatusoutput(cmd)
  print("%s: %s"%(status,output))
