import commands

cmd = []
for n in range(300,410,10):
  file1 = "2lp5.%s.par"%(n)
  file2 = "2lp5.%s.nos.par"%(n)
  read = open(file1,'r')
  write = open(file2,'w+')
  for line in read:
    if "BARRELMOVES" in line:
      write.write("BARRELMOVES 0\n")
    if "LNAME" in line:
      write.write("LNAME 2lp5.%s.nos\n"%(n))
    else:
      write.write(line)
  read.close()
  write.close()
  cmd.append("python rungeofold.py %s walcob"%(file2))

for command in cmd:
  status,output = commands.getstatusoutput(command)
  print("%s: %s"%(status,output))
