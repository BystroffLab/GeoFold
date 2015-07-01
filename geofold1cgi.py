#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import commands
import os
import urllib
import os.path
import cgi
import cgitb
cgitb.enable

def getchains(pdbfile):
  """Extracts chain information from PDB file
  returns a list of the chain characters"""
  try:
    read = open(pdbfile,'r')
  except IOError:
    print("getchains: Couldn't open file %s"%(pdbfile))
    raise
  else:
    result = []
    for line in read:
      line = line.split()
      if line[0]=='ATOM':
        if line[4] not in result and line[4].isalpha():
          result.append(line[4])
        elif "_" not in result and not line[4].isalpha():
          result.append("_")
    read.close()
    return result

def printStyle():
  """print the style info necessary for collapsible portion"""
  print('<style type="text/css">')
  print('.row { vertical-align: top; height:auto !important; }')
  print('.list {display:none; }')
  print('.show { display: none; }')
  print('.hide:target + .show {display: inline; }')
  print('.hide:target {display: none; }')
  print('.hide:target ~ .list {display:inline; }')
  print('@media print { .hide, .show { display: none; } }')
  print('</style>')
  
def collapseHead():
  """print the first part of the collapsible portion"""
  print('<div class="row">')
  print(' <a href="#hide1" class="hide" id="hide1">Show advanced options</a>')
  print(' <a href="#show1" class="show" id="show1">Hide advanced options</a>')
  print(' <div class="list">')
  

#HTML header
print "Content-Type: text/html;charset=utf-8"
print
#Test statement
print("<html><head>")
printStyle()
#set directories
basedir = "/bach1/home/flex"
urldir = "/bach1/home/flex/public_html/geofold"
tmpdir = "%s/server/geofold/tmp"%(basedir)
dbgfile = "%s/output/debug.out"%(urldir)
#default parameters file
paramfile = "%s/server/geofold/bin/parameters"%(basedir)
#pdb repository
pdbdir = "%s/server/data/pdb"%(basedir)
pdbunit = "%s/server/data/pdb1"%(basedir)
settings= "settings.html"
pid = os.getpid()
outdir = "output/"
jobdir = "jobs/"
#read form data
form = cgi.FieldStorage()
#set variables from form data
email = form["email_address"].value
try:
  pdbid = form["pdbid"].value
  pdbfile = form["pdbfile"].value
  code = False
except KeyError:
  pdbcode = form["pdbcode"].value
  code = True
keyword = form["keyword"].value
keyword = keyword.split()
keyword = "_".join(keyword)
#remove the words bug and error from keyword
keyword = keyword.lower().split("bug")
keyword = "bg".join(keyword)
keyword = keyword.lower().split("error")
keyword = "err".join(keyword)
keyword = keyword.strip("/*{}[]!@#$%^&();:<>,?\\~\"\'")
lname="%s.%s"%(keyword,pid)
#check if pdb is from code or from file-upload
if not code:
  if pdbid == "": file = True
  else: file = False
  #create uploaded pdbfile
  if file:
    if pdbfile != "" :
      try:
        pdbwrite = open("%s/%s.pdb"%(tmpdir,pid),'w+')
      except IOError:
        print("<br>Couldn't open pdbfile")
        raise IOError
      else:
        for line in pdbfile:
          pdbwrite.write(line)
        pdbwrite.close()
        commands.getstatusoutput("chmod 0777 %s/%s.pdb"%(tmpdir,pid))
    else:
      print("<br>No pdbfile found.")
      raise IOError
  else:
    url = "http://pdb.org/pdb/files/%s.pdb"%(pdbid.upper()[0:4])
    urllib.urlretrieve(url,"%s/%s.pdb"%(tmpdir,pid))
else:
  print("moving file")
  status,output = commands.getstatusoutput("cp %s/%s.pdb %s/%s.pdb"%(tmpdir,pdbcode,tmpdir,pid))
  print("file moved")
  if status !=0: print("%s: %s"%(status,output))
#extract chains from pdb file
status,output = commands.getstatusoutput("cp %s/%s.pdb %s/%s.pdb"%(tmpdir,pid,pdbdir,pid))
if status != 0:
  print("%s: %s"%(status,output))
chains = getchains("%s/%s.pdb"%(tmpdir,pid))
#write HTML form to go to next cgi script
print("<title>GeoFOLD: %s</title><h3>GeoFOLD: %s</h3></head>"%(lname,lname))
print("<body>")
#chains
if len(chains)==0:
  print("Error, no chains found.")
  raise IOError
print("<br>Select chains for <strong>%s</strong>"%(lname))
print('<form method="post" action="geofold2cgi.py">')
for chain in chains:
  print('<br><input type="checkbox" name="chain_%s" value="%s">%s'%(chain,chain,chain))
#hidden inputs
#LNAME##
print('<input type="hidden" name="lname" value="%s">'%(lname))
#EMAIL##
print('<input type="hidden" name="email" value="%s">'%(email))
#PDBCODE##
print('<input type="hidden" name="pdbcode" value="%s">'%(pid))
#OMEGA##
print('<input type="hidden" name="omega" value="5.">')
#INTERMEDIATE##
print('<input type="hidden" name="intermediate" value="0">')
#BARRELMOVES##
print('<input type="hidden" name="barrelmoves" value="1">')
print('<br>mouseover for more info')
#ORANGE
print('<br><span title="Enter a range of surface tension values (space-separated with decimal points ie. 5. 10. 20.)">&omega; range</span><input type="text" name="orange" placeholder="ex: 5. 10. 20." value="" size=40>')
#checkboxes
#FOLDING checkbox
print('<br><input type="checkbox" name="fing" value="1">Folding simulation')
#HALFLIFE checkbox
print('<br><input type="checkbox" name="hlfe" value="1">Run until half folded/unfolded')
#RUNGEOFOLD checkbox hidden
#print('<br><input type="checkbox" name="rungeofold" value="1" checked> Run Geofold (only uncheck if you are re-running a job)')
print('<input type="hidden" name="rungeofold" value="1">')
#REDUCING checkbox
print('<br><input type="checkbox" name="reducing" value="1">Reduce cysteine disulfides')
#MOLSCRIPT checkbox
print('<br><input type="checkbox" name="molscript" value="1" checked>Create gif of protein')
#Show debug info?
print('<br><input type="checkbox" name="debug" value="1"> Show detailed output while running')
#Advanced options
print('<br><strong>Advanced Options</strong>')
collapseHead()
print('<table><tr><td>')
#BREAKCUT# 0.05
print('<br><span title="Cut-off percentage (expressed as decimal) for break move">Break cut</span>:<br> <input type="text" name="breakcut" value="0.05" size=15>')
#PIVOTCUT# 0.01
print('<br><span title="Cut-off percentage (expressed as decimal) for pivot move">Pivot cut</span>:<br> <input type="text" name="pivotcut" value="0.01" size=15>')
#HINGECUT# 0.5
print('<br><span title="Hinge moves must be able rotate 60&deg;*this number">Hinge cut</span>:<br> <input type="text" name="hingecut" value="0.5" size=15>')
#SEAMCUT# 10.
print('<br><span title="Cut off for seam move">Seam cut</span>:<br><input type="text" name="seamcut" value="10" size=15>')
#BREAKPOINTENTROPY# 90.
print('<br><span title="Entropy gained upon a break, in kJ/mol/K.  A break is the separation of two unconnected segments or chains.">Break entropy</span>:<br><input type="text" name="breakpointentropy" value="90." size=15>')
print('</td><td>')
#HINGEPOINTENTROPY# 30.
print('<br><span title="Entropy gained upon a hinge move, in kJ/mol/K.  A hinge is the rotation of a segment around two points located at either end of the rotating segment.">Hinge entropy</span>:<br><input type="text" name="hingepointentropy" value="30." size=15>')
#TEMPERATURE# 300.
print('<br><span title="in Kelvin">Temperature</span>:<br><input type="text" name="temperature" value="300." size=15>')
#CONCENTRATION# 1.
print('<br><span title="Total protein concentration in mol/L to be assigned ot the system.  If the protein is an oligomer, then this concentration may affect the unfolding pathway.  Otherwise, it has no effect.">Concentration</span>:<br><input type="text" name="concentration" value="1." size=15>')
#VOIDENTROPY# 0.
print('<br><span title="Entropy gained when a buried void space is exposed by unfolding, in kJ/mol/K.">Void entropy</span>:<br><input type="text" name="voidentropy" value="0." size=15>')
#SOLIDITY# 1000.
print('<br><span title="When the total buried surface of an intermediate is less than SOLIDITY, the intermediate is treated as a liquid, otherwise it is a solid. For a solid, the entropy of unfolding is expressed after the transition state. For a liquid, it is expressed before the transition state.">Solidity</span>:<br><input type="text" name="solidity" value="1000." size=15>')
print('</td><td>')
#HBONDENERGY# 1000. -> 1.
print('<br><span title="J/mol per hydrogen bond.">H-bond energy</span>:<br><input type="text" name="hbondenergy" value="100." size=15>')
#HAMMONDSCALE# 1000.
print('<br><span title="Hammond postulated that the tranition state lies closer to the higher energy ground state on the reaction coordinate. The reaction coordinate in GeoFold is the solvent accessible surface area (SAS). The transition state position (&Theta;-m) is defined for each unfolding step as the fraction of &Delta;SAS expressed at the t-state. HAMMONDSCALE (H) is in units of surface area &Aring;&sup2;. The location of the transitions state, &Theta;-m = ( tanh( dNRG/ H) +1 ) /0.5">Hammond scale</span>:<br><input type="text" name="hammondscale" value="1000." size=15>')
#SIDECHAINENTROPY# 1.
print('<br><span title="Scale factor for sidechain entropy. Each amino acid has a different sidechain entropy, which is the difference in entropy between bound and fully solvent-exposed states. In this case, the amount of sidechain entropy expressed in a given unfolding step is equal to its intrisic entropy times the fraction of its total surface area exposed in that step.">Sidechain entropy</span>:<br><input type="text" name="sidechainentropy" value="1." size=15>')
#HINGEBARRIER# 0. or 50.
print('<br><span title="Optional additional barrier for hinge moves">Hinge barrier</span>:<br><input type="text" name="hingebarrier" value="0." size=15>')
#PIVOTBARRIER# 0. or 20.
print('<br><span title="Optional additional barrier for pivot moves.">Pivot barrier</span>:<br><input type="text" name="pivotbarrier" value="0." size=15>')
print('</td><td>')
#WATER# 30.
print('<br><span title="Virtual denaturant level equivalent to water.  kJ/mol/&Aring;&sup2;">Water &omega;</span>:<br><input type="text" name="water" value="30." size=15>')
#MAXSPLIT# 4
print('<br><span title="Maximum number of elemental unfolding steps applied to one intermediate. Optional values are 2,4,8,16, and 32. Speed is effected by higher settings.">Max split</span>:<br><input type="text" name="maxsplit" value="4" size=15>')
#MAXTIME# 10.
print('<br><span title="Maximum time to run simulation (s)">Max time</span>:<br><input type="text" name="maxtime" value="10." size=15>')
#MINSEG# 4
print('<br><span title="Smallest length of chain for pivot move.">Minimum segment</span>:<br><input type="text" name="minseg" value="4" size=15>')
#CAVITATION# 0.000001
print('<br><span title="Cavitation.">Cavitation</span>:<br><input type="text" name="cavitation" value="0.000001" size=15>')
print('</td></tr></table>')
#FLORY# 1
print('<br>Flory effect Modeling:')
print('<br><input type="radio" name="FLORY" value="0">Don\'t model')
print('<br><input type="radio" name="FLORY" value="1" checked>Take lowest path')
print('<br><input type="radio" name="FLORY" value="2">Take highest path')
print('<br><input type="radio" name="FLORY" value="3">Take lowest weighted avg path')
print('<br><span title="weight to be used for non-covalent contacts.">Non-covalent weight</span>:<br><input type="text" name="FLORYW" value="8." size=15>')
print('</div></div>')
#submit button
print('<br><input type="submit" name = "submit" value="submit"></form>')
  
print "</body></html>"
