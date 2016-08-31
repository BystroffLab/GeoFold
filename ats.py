import commands

names = ['noseam','seam2','seam','seamx2','seamx']

for value in names:
  cmd = 'python rungeofold.py aTS_%s.par walcob'%(value)
  status,output = commands.getstatusoutput(cmd)
  print("%s: %s"%(status,output))
