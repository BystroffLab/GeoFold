import commands
for n in range(300,410,10):
  cmd = "python rungeofold.py 2lp5.%s.par walcob"%(n)
  status,output = commands.getstatusoutput(cmd)
  print("%s: %s"%(status,output))
