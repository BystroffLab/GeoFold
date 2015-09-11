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

def toggle(name,checked,msg):
    if checked:
        check = " checked"
    else:
        check = ""
    print('<br>%s<div class="onoffswitch">'%(msg))
    print('<input type="checkbox" name="%s" class="onoffswitch-checkbox" id="%s"%s>'%(name,name,check))
    print('<label class="onoffswitch-label" for="%s">'%(name))
    print('<span class="onoffswitch-inner"></span>')
    print('<span class="onoffswitch-switch"></span>')
    print('</label></div>')

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
        if line[4] not in result and line[4].isalnum():
          result.append(line[4])
        elif "_" not in result and not line[4].isalnum():
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
  
def makeParameters(form,out):
  """creates the parameters file from the data read from the form"""
  chainset = False
  chains = 'CHAIN '
  checkboxes = ['reducing','molscript','hlfe','fing']
  for value in form:
    if "chain_" in value:
      chains = "%s%s"%(chains,form[value].value)
      chainset = True
    else:
      try:
        out.write("%s %s\n"%(value.upper(),form[value].value))
      except KeyError:
        out.write("%s 0\n"%(value.upper()))
  for value in checkboxes:
    try:      
      if form[value].value == 'on':
          out.write('%s 1\n'%(value.upper()))
      else:
          out.write("%s 0\n"%(value.upper())
      #out.write("%s %s\n"%(value.upper(),form[value].value))
    except KeyError:
      out.write("%s 0\n"%(value.upper()))
  if not chainset:
    chains += '.'
  out.write("%s\n"%(chains))
  out.close()  
  
def firstScript(form):
    """
    #Walcob settings
    #set directories
    basedir = "/bach1/home/walcob"
    urldir = "/bach1/home/walcob/public_html/GeoFold"
    tmpdir = "%s/GeoFold/tmp"%(basedir)
    dbgfile = "%s/output/debug.out"%(urldir)
    #default parameters file
    paramfile = "%s/GeoFold/parameters"%(basedir)
    #pdb repository
    pdbdir = "%s/GeoFold/pdbs"%(basedir)
    pdbunit = pdbdir
    settings= "settings.html"
    pid = os.getpid()
    outdir = "output/"
    jobdir = "jobs/"
    """
    #Flex settings    
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
      if pdbid == "": 
        file = True
      else: 
        file = False
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
      status,output = commands.getstatusoutput("cp -v %s/%s.pdb %s/%s.pdb"%(tmpdir,pdbcode,tmpdir,pid))
      if status !=0: print("%s: %s"%(status,output))
      print("%s: %s"%(status,output))
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
    print('<form method="post" action="geocgi.cgi">')
    print('<input type="hidden" name="script" value=2>')
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
    print('<input type="hidden" name="omega" value="1.">')
    #INTERMEDIATE##
    print('<input type="hidden" name="intermediate" value="0">')
    #BARRELMOVES##
    print('<input type="hidden" name="barrelmoves" value="1">')
    print('<br>mouseover for more info')
    #ORANGE
    print('<br><span title="Enter a range of surface tension values or temperatures if denaturing thermally (space-separated with decimal points ie. 5. 10. 20.)">&omega; range</span><input type="text" name="orange" placeholder="ex: 5. 10. 20." value="" size=40>')
    #checkboxes/toggles
    #FOLDING checkbox
    #print('<br><input type="checkbox" name="fing" value="1">Folding simulation')
    toggle('fing',False,'Folding Simulation')
    #HALFLIFE checkbox
    #print('<br><input type="checkbox" name="hlfe" value="1">Run until half folded/unfolded')
    toggle('hlfe',False,'Run until half folded/unfolded')
    #RUNGEOFOLD checkbox hidden
    #print('<br><input type="checkbox" name="rungeofold" value="1" checked> Run Geofold (only uncheck if you are re-running a job)')
    print('<input type="hidden" name="rungeofold" value="1">')
    #REDUCING checkbox
    #print('<br><input type="checkbox" name="reducing" value="1">Reduce cysteine disulfides')
    toggle('reducing',False,'Reduce cysteine disulfides')
    #MOLSCRIPT checkbox
    #print('<br><input type="checkbox" name="molscript" value="1" checked>Create gif of protein')
    toggle('molscript',True,'Create GIF of protein')
    #thermal
   # print('<br><input type="checkbox" name="thermal" value="1">Use thermal denaturation')
    toggle('thermal',False,'Use thermal denaturation')
    #Show debug info?
    #print('<br><input type="checkbox" name="debug" value="1">Show detailed output while running')
    toggle('debug',False,'Show detailed output while running')
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
    print('<br><span title="Virtual denaturant level equivalent to water.  kJ/mol/&Aring;&sup2; If using thermal denaturation, this will be the temperature to which unfolding rates will be extrapolated.">Water &omega;</span>:<br><input type="text" name="water" value="1." size=15>')
    #MAXSPLIT# 4
    print('<br><span title="Maximum number of elemental unfolding steps applied to one intermediate. Optional values lp5.7310_nos.htmlare 2,4,8,16, and 32. Speed is effected by higher settings.">Max split</span>:<br><input type="text" name="maxsplit" value="4" size=15>')
    #MAXTIME# 10.
    print('<br><span title="Maximum time to run simulation (s)">Max time</span>:<br><input type="text" name="maxtime" value="10." size=15>')
    #MINSEG# 4
    print('<br><span title="Smallest length of chain for pivot move.">Minimum segment</span>:<br><input type="text" name="minseg" value="4" size=15>')
    #CAVITATION# 0.000001
    print('<br><span title="Cavitation.">Cavitation</span>:<br><input type="text" name="cavitation" value="0.000001" size=15>')
    print('</td></tr></table>')
    #FLORY# 1
    print('<br>Flory effect Modeling:')
    print('<br><input type="radio" name="flory" value="0">Don\'t model')
    print('<br><input type="radio" name="flory" value="1" checked>Take lowest path')
    print('<br><input type="radio" name="flory" value="2">Take highest path')
    print('<br><input type="radio" name="flory" value="3">Take lowest weighted avg path')
    print('<br><span title="weight to be used for non-covalent contacts.">Non-covalent weight</span>:<br><input type="text" name="floryw" value="8." size=15>')
    print('</div></div>')
    #submit button
    print('<br><input type="submit" name = "submit" value="submit"></form>')
      
    print "</body></html>"
    
def secondScript(form):
    #directories and basic settings
    waitimage = "http://www.bioinfo.rpi.edu/bystrc/pub/mpeg/peptide.gif"
    #set directories
    basedir = "/bach1/home/flex"
    gdir = "%s/server/geofold"%(basedir)
    urldir = "/bach1/home/flex/public_html/geofold"
    tmpdir = "%s/server/geofold/tmp"%(basedir)
    #default parameters file
    paramfile = "%s/server/geofold/bin/parameters"%(basedir)
    #pdb repository
    pdbdir = "%s/server/data/pdb"%(basedir)
    pdbunit = "%s/server/data/pdb1"%(basedir)
    settings= "settings.html"
    pid = os.getpid()
    outdir = "output/"
    jobdir = "%s/jobs"%(gdir)
    #set url name
    try:
      lname = form["lname"].value
    except KeyError:
      lname = "%s.%s"%(form["keyword"].value,form["pdbcode"].value)
    url = "%s/%s%s/%s.html"%(urldir,outdir,lname,lname)
    urlwrite = "http://www.bioinfo.rpi.edu/geofold/%s%s/%s.html"%(outdir,lname,lname)
    #urlwrite = "http://bach1.bio.rpi.edu/walcob/GeoFold/%s%s/%s.html"%(outdir,lname,lname)
    
    #write HTML output
    print('<html><head>')
    print('<meta http-equiv="refresh" content="30;url=%s">'%(urlwrite))
    print('</head><body><h4>Creating GeoFOLD job. Please wait...</h4><br>')
    #print('</body></html>')
    #new URL settings
    """
    try:
      out = open(url,'w+')
    except IOError:
      print 'IOERROR'
      """
    #reset url name
    url = "%s.html"%(lname)
    print url
    #write-out parameters file
    file = open("%s/%s.par"%(tmpdir,lname),'w+')
    makeParameters(form,file)
    print 'parameters made in %s/%s.par'%(tmpdir,lname)
    #get chain info for html results
    chains = ''
    for value in form:
      if "chain_" in value:
        chains+=form[value].value
      else:
        chains = '.'
    print('<h4>Your GeoFold job is in the queue but has not yet started.</h4>\n')
    print('<p>Your input coordinates were uploaded as filename %s chains %s\n'%(form["pdbcode"].value,chains))
    print('<p>Notifications will be sent to %s\n'%(form["email"].value))
    print('<p><a href="%s">BOOKMARK THIS PAGE</a>:\n'%(url))
    print('<p>Stay on this page, or return to this page to see your results,')
    print('which will include the following:<br><ul>\n')
    print('<li>The timecourse of unfolding as concentrations of Folded, Unfolded, and Intermediate states.</li>\n')
    print('<li>An unfolding pathway in the form of a clickable pathway tree.</li>\n')
    print('<li>An Age Plot expressing the order of contact loss during unfolding.</li>\n')
    print('<li>Unfolding/folding kinetics simulations and associated plots.</li>\n')
    print('</ul>')
    print('<p>You will see results for multiple values of &omega; (virtual denaturant).\n')
    print('<p>You will have the option to Do Over, changing the &omega; values or any other parameters.\n')
    print('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/settings.html">Click on this page')
    print('to see a brief explanation of the parameters.</a>\n')
    print('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/howtoreadit.htm">Click on this page')
    print('to see a guide to GeoFold output.</a>\n')
    print('<p><img src="%s"><br>\n'%(waitimage))
    print('</body></html>')
    #write job file
    print 'attempting to write job to %s/%s.job'%(jobdir,form['lname'].value)
    job = open("%s/%s.job"%(jobdir,form["lname"].value),'w+')
    job.write('%s %s %s'%(form["pdbcode"].value,chains,form["lname"].value))
    try:
      oname = form["oname"].value
    except KeyError:
      oname = ''
    job.write(' %s'%(oname))
    job.close()
    print 'job written to %s/%s.job'%(jobdir,form['lname'].value)
    #deprotect files
    commands.getstatusoutput('chmod 0777 %s/%s'%(outdir,url))
    commands.getstatusoutput('chmod 0777 %s/%s.job'%(jobdir,form["lname"].value))
    commands.getstatusoutput('chmod 0777 %s/%s.par'%(tmpdir,form["lname"].value))
    
def redo(form):
    #based on firstScript    
    checked = {'1':' checked','0':''} #['',' checked']
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
    #read info from oldparfile
    oldParameterFile = form['oldParameters'].value
    readOldParameters = open(oldParameterFile,'r')
    oldParameters = {}
    for line in readOldParameters:
      if len(line) != 1:
        line = line.split()
        oldParameters[line[0]]=' '.join(line[1:])
    readOldParameters.close()
    #keyword is changed
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
    #Note to users about redo script
    print('<title>GeoFOLD: %s rerun</title><h3>GeoFOLD: %s</h3></head>\n'%(lname,lname))
    print('<body><br><p>Please configure your rerun of GeoFOLD.  Below are the options selected from the previous run.  Please adjust them accordingly.')    
    print("<br>Select chains for <strong>%s</strong>"%(lname))
    #recreate form from firstScript using default info from oldparfile as needed
    print('<form method="post" action="geocgi.cgi">')
    #hidden input script = 4
    print('<input type="hidden" name="script" value=4>')
    #chains
    chains = getchains("%s/%s.pdb"%(tmpdir,oldParameters['PDBCODE']))
    if len(chains)==0:
      print("Error, no chains found.")
      raise IOError 
    for chain in chains:
      if chain in oldParameters['CHAIN']:
        print('<br><input type="checkbox" name="chain_%s" value="%s" checked>%s'%(chain,chain,chain))
      else:
        print('<br><input type="checkbox" name="chain_%s" value="%s">%s'%(chain,chain,chain))
    #hidden inputs
    print('<input type="hidden" name="lname" value="%s">'%(lname))
    print('<input type="hidden" name="email" value="%s">'%(oldParameters['EMAIL']))
    print('<input type="hidden" name="redo" value="%s">'%(oldParameters['LNAME']))
    print('<input type="hidden" name="pdbcode" value="%s">'%(oldParameters['PDBCODE']))
    print('<input type="hidden" name="omega" value="1.">')
    print('<input type="hidden" name="intermediate" value="0">')
    print('<input type="hidden" name="barrelmoves" value="1">')
    print('<input type="hidden" name="oldParameters" value="%s">'%(oldParameterFile))
    print('<br>mouseover for more info')
    #ORANGE
    print('<br><span title="Enter a range of solvation weights or temperatures (space-separated with decimal points ie. 5. 10. 20.)">&omega; range</span><input type="text" name="orange" value="%s" size=40>'%(oldParameters['ORANGE']))
    #Folding
    #print('<br><input type="checkbox" name="fing" value="1"%s>Folding simulation'%(checked[oldParameters['FING']]))
    toggle('fing',oldParameters['FING']==1,'Folding simulation')
    #Half-life
    #print('<br><input type="checkbox" name="hlfe" value="1"%s>Run until half folded/unfolded'%(checked[oldParameters['HLFE']]))
    toggle('hlfe',oldParameters['HLFE']==1,'Run until half folded/unfolded')    
    #reducing
    #print('<br><input type="checkbox" name="reducing" value="1"%s>Reduce cysteine disulfides'%(checked[oldParameters['REDUCING']]))
    toggle('reducing',oldParameters['REDUCING']==1,'Reduce cysteine disulfides')    
    #molscript
    #print('<br><input type="checkbox" name="molscript" value="1"%s>Create gif of protein'%(checked[oldParameters['MOLSCRIPT']]))
    toggle('molscript',oldParameters['MOLSCRIPT']==1,'Create GIF of protein')
    #thermal
    try:
    	#print('<br><input type="checkbox" name="thermal" value="1"%s>Thermal Denaturaton'%(checked[oldParameters["THERMAL"]]))
      toggle('thermal',oldParameters['THERMAL']==1,'Thermal Denaturation')
    except:
        #print('<br><input type="checkbox" name="thermal" value="1">Thermal Denaturaton')
        toggle('thermal',False,'Thermal Denaturation')
    #debug
    try:
    	#print('<br><input type="checkbox" name="debug" value="1"%s>Show detailed output while running'%(checked[oldParameters["DEBUG"]]))
      toggle('debug',oldParameters['DEBUG']==1,'Show detailed output while running')
    except:
      #print('<br><input type="checkbox" name="debug" value="1">Show detailed output while running')
      toggle('debug',False,'Show detailed output while running')
    #Advanced options
    print('<br><strong>Advanced Options</strong>')
    collapseHead()
    print('<table><tr><td>')
    #breakcut
    print('<br><span title="cut-off percentage (expressed as decimal) for break move">Break cut</span>:<br><input type="text" name="breakcut" value="%s" size=15>'%(oldParameters['BREAKCUT']))
    #pivotcut
    print('<br><span title="Cut-off percentage (expressed as decimal) for pivot move">Pivot cut</span>:<br><input type="text" name="pivotcut" value="%s" size=15>'%(oldParameters["PIVOTCUT"]))
    #hingecut
    print('<br><span title="Cut-off percentage (expressed as decimal) for hinge move">Hinge cut</span>:<br><input type="text" name="hingecut" value="%s" size=15>'%(oldParameters["HINGECUT"]))
    #seamcut
    print('<br><span title="Cut off for seam move">Seam cut</span>:<br><input type="text" name="seamcut" value="%s" size=15>'%(oldParameters["SEAMCUT"]))
    #breakpointentropy
    print('<br><span title="Entropy gained upon a break, in kJ/mol/K.  A break is the separation of two unconnected segments or chains.">Break entropy</span>:<br><input type="text" name="breakpointentropy" value="%s" size=15>'%(oldParameters["BREAKPOINTENTROPY"]))
    print('</td><td>')    
    #hingepointentropy
    print('<br><span title="Entropy gained upon a hinge move, in kJ/mol/K.  A hinge is the rotation of a segment around two points located at either end of the rotating segment.">Hinge entropy</span>:<br><input type="text" name="hingepointentropy" value="%s" size=15>'%(oldParameters["HINGEPOINTENTROPY"]))
    #temperature
    print('<br><span title="in Kelvin">Temperature</span>:<br><input type="text" name="temperature" value="%s" size=15>'%(oldParameters["TEMPERATURE"]))
    #concentration
    print('<br><span title="Total protein concentration in mol/L to be assigned ot the system.  If the protein is an oligomer, then this concentration may affect the unfolding pathway.  Otherwise, it has no effect.">Concentration</span>:<br><input type="text" name="concentration" value="%s" size=15>'%(oldParameters["CONCENTRATION"]))    
    #voidentropy
    print('<br><span title="Entropy gained when a buried void space is exposed by unfolding, in kJ/mol/K.">Void entropy</span>:<br><input type="text" name="voidentropy" value="%s" size=15>'%(oldParameters["VOIDENTROPY"]))
    #solidity
    print('<br><span title="When the total buried surface of an intermediate is less than SOLIDITY, the intermediate is treated as a liquid, otherwise it is a solid. For a solid, the entropy of unfolding is expressed after the transition state. For a liquid, it is expressed before the transition state.">Solidity</span>:<br><input type="text" name="solidity" value="%s" size=15>'%(oldParameters["SOLIDITY"])) 
    print('</td><td>')
    #h-bond energy
    print('<br><span title="J/mol per hydrogen bond.">H-bond energy</span>:<br><input type="text" name="hbondenergy" value="%s" size=15>'%(oldParameters["HBONDENERGY"]))    
    #Hammond
    print('<br><span title="Hammond postulated that the tranition state lies closer to the higher energy ground state on the reaction coordinate. The reaction coordinate in GeoFold is the solvent accessible surface area (SAS). The transition state position (&Theta;-m) is defined for each unfolding step as the fraction of &Delta;SAS expressed at the t-state. HAMMONDSCALE (H) is in units of surface area &Aring;&sup2;. The location of the transitions state, &Theta;-m = ( tanh( dNRG/ H) +1 ) /0.5">Hammond scale</span>:<br><input type="text" name="hammondscale" value="%s" size=15>'%(oldParameters["HAMMONDSCALE"]))
    #sidechainentropy
    print('<br><span title="Scale factor for sidechain entropy. Each amino acid has a different sidechain entropy, which is the difference in entropy between bound and fully solvent-exposed states. In this case, the amount of sidechain entropy expressed in a given unfolding step is equal to its intrisic entropy times the fraction of its total surface area exposed in that step.">Sidechain entropy</span>:<br><input type="text" name="sidechainentropy" value="%s" size=15>'%(oldParameters["SIDECHAINENTROPY"]))
    #hingebarrier
    print('<br><span title="Optional additional barrier for hinge moves">Hinge barrier</span>:<br><input type="text" name="hingebarrier" value="%s" size=15>'%(oldParameters["HINGEBARRIER"]))
    #pivotbarrier
    print('<br><span title="Optional additional barrier for pivot moves.">Pivot barrier</span>:<br><input type="text" name="pivotbarrier" value="%s" size=15>'%(oldParameters["PIVOTBARRIER"]))
    print('</td><td>')    
    #water
    print('<br><span title="Virtual denaturant level equivalent to water.  kJ/mol/&Aring;&sup2;.  If using thermal denaturation, this will be the temperature to which unfolding rates will be extrapolated.">Water &omega;</span>:<br><input type="text" name="water" value="%s" size=15>'%(oldParameters["WATER"]))
    #maxsplit
    print('<br><span title="Maximum number of elemental unfolding steps applied to one intermediate. Optional values are 2,4,8,16, and 32. Speed is effected by higher settings.">Max split</span>:<br><input type="text" name="maxsplit" value="%s" size=15>'%(oldParameters["MAXSPLIT"]))
    #maxtime
    print('<br><span title="Maximum time to run simulation (s)">Max time</span>:<br><input type="text" name="maxtime" value="%s" size=15>'%(oldParameters["MAXTIME"]))
    #minseg
    print('<br><span title="Smallest length of chain for pivot move.">Minimum segment</span>:<br><input type="text" name="minseg" value="%s" size=15>'%(oldParameters["MINSEG"]))
    #cavitation
    print('<br><span title="Cavitation.">Cavitation</span>:<br><input type="text" name="cavitation" value="%s" size=15>'%(oldParameters["CAVITATION"]))    
    #end table
    print('</td></tr></table>')
    #FLORY
    flory = ['','','','']
    flory[int(oldParameters["FLORY"])] = ' checked'
    print('<br>Flory effect Modeling:')
    num = 0
    print('<br><input type="radio" name="flory" value="0"%s>Don\'t model'%(flory[0]))
    print('<br><input type="radio" name="flory" value="1"%s>Take lowest path'%(flory[1]))
    print('<br><input type="radio" name="flory" value="2"%s>Take highest path'%(flory[2]))
    print('<br><input type="radio" name="flory" value="3"%s>Take lowest weighted avg path'%(flory[3]))
    print('<br><span title="weight to be used for non-covalent contacts.">Non-covalent weight</span>:<br><input type="text" name="floryw" value="%s" size=15>'%(oldParameters["FLORYW"]))
    #submit
    print('</div></div>\n<br><input type="submit" name="submit" value="submit"></form>\n</body></html>')
    
    #script 4 will determine which parts of geofold will need to be rerun if any and direct things accordingly

def fourthScript(form):
    """This script determines what parts of GeoFold need to be rerun and will send the results to secondScript"""
    oldParameters = {}
    changeParameters = ['REDUCING','BREAKCUT','PIVOTCUT','HINGECUT','SEAMCUT','BREAKPOINTENTROPY','HINGEPOINTENTROPY','TEMPERATURE','VOIDENTROPY','SOLIDITY','HBONDENERGY','HAMMONDSCALE','SIDECHAINENTROPY','HINGEBARRIER','PIVOTBARRIER','WATER','MAXSPLIT','MINSEG','CAVITATION','FLORY','FLORYW']
    Parameters = {}
    runGeoFold = '0'
    #Read in old parameters
    oldParFile = open(form['oldParameters'].value,'r')
    for line in oldParFile:
        if len(line) != 1:
          line = line.split()
          oldParameters[line[0]] = " ".join(line[1:])
    oldParFile.close()
    #Check if any parameters have changed
    for parameter in changeParameters:
        if oldParameters[parameter] != form.getvalue(parameter.lower(),'0'):
            runGeoFold = '1'
    #form.add_field('rungeofold', runGeoFold)
    #Everything else is the same as secondScript, so copypasta
    
    #directories and basic settings
    waitimage = "http://www.bioinfo.rpi.edu/bystrc/pub/mpeg/peptide.gif"
    #set directories
    basedir = "/bach1/home/flex"
    gdir = "%s/server/geofold"%(basedir)
    urldir = "/bach1/home/flex/public_html/geofold"
    tmpdir = "%s/tmp"%(gdir)
    #default parameters file
    paramfile = "%s/server/geofold/bin/parameters"%(basedir)
    #pdb repository
    pdbdir = "%s/server/data/pdb/"%(basedir)
    pdbunit = "%s/server/data/pdb1/"%(basedir)
    settings= "settings.html"
    pid = os.getpid()
    outdir = "output/"
    jobdir = "%s/jobs"%(gdir)
    #set url name
    try:
      lname = form["lname"].value
    except KeyError:
      lname = "%s.%s"%(form["keyword"].value,form["pdbcode"].value)
    url = "%s/%s%s/%s.html"%(urldir,outdir,lname,lname)
    
    urlwrite = "http://www.bioinfo.rpi.edu/geofold/%s%s/%s.html"%(outdir,lname,lname)
    
    #write HTML output
    os.environ["JOBURL"] = urlwrite
    print('<meta http-equiv="refresh" content="4;url=%s">'%(urlwrite))
    print('</head><body><h4>Creating GeoFOLD job. Please wait...</h4><br>')
    
    #new URL settings
    #reset url name
    url = "%s/%s.html"%(lname,lname)
    #write-out parameters file
    file = open("%s/%s.par"%(tmpdir,lname),'w+')
    makeParameters(form,file)
    write_param = open('%s/%s.par'%(tmpdir,lname),'a')
    write_param.write("RUNGEOFOLD "+runGeoFold)
    write_param.close()
    #get chain info for html results
    chains = ''
    for value in form:
      if "chain_" in value:
        chains+=form[value].value
      else:
        chains = '.'
    #Write out initial HTML results page
    print('<h4>Your GeoFold job is in the queue but has not yet started.</h4>\n')
    print('<p>Your input coordinates were uploaded as filename %s chains %s\n'%(form["pdbcode"].value,chains))
    print('<p>Notifications will be sent to %s\n'%(form["email"].value))
    print('<p><a href="%s">BOOKMARK THIS PAGE</a>:\n'%(url))
    print('<p>Stay on this page, or return to this page to see your results,')
    print('which will include the following:<br><ul>\n')
    print('<li>The timecourse of unfolding as concentrations of Folded, Unfolded, and Intermediate states.</li>\n')
    print('<li>An unfolding pathway in the form of a clickable pathway tree.</li>\n')
    print('<li>An Age Plot expressing the order of contact loss during unfolding.</li>\n')
    print('<li>Unfolding/folding kinetics simulations and associated plots.</li>\n')
    print('</ul>')
    print('<p>You will see results for multiple values of &omega; (virtual denaturant).\n')
    print('<p>You will have the option to Do Over, changing the &omega; values or any other parameters.\n')
    print('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/settings.html">Click on this page')
    print('to see a brief explanation of the parameters.</a>\n')
    print('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/howtoreadit.htm">Click on this page')
    print('to see a guide to GeoFold output.</a>\n')
    print('<p><img src="%s"><br>\n'%(waitimage))
    print('</body></html>')
    #out.close()
    #write job file
    job = open("%s/%s.job"%(jobdir,lname),'w+')
    job.write('%s %s %s'%(form["pdbcode"].value,chains,lname))
    try:
      oname = form["oname"].value
    except KeyError:
      oname = ''
    job.write(' %s'%(oname))
    job.close()
    #deprotect files
    commands.getstatusoutput('chmod 0777 %s/%s'%(outdir,url))
    commands.getstatusoutput('chmod 0777 %s/%s.job'%(jobdir,lname))
    commands.getstatusoutput('chmod 0777 %s/%s.par'%(tmpdir,lname))
    

#HTML header
print "Content-Type: text/html;charset=utf-8\n\n"
#Test statement
print("<html><head>")
print('<link type="text/css" rel="stylesheet" href="toggle.css"/>')
printStyle()

#print("</head><body>")
form = cgi.FieldStorage()
query = int(form['script'].value)
if query == 1:
    firstScript(form)
elif query == 2:
    secondScript(form)
elif query == 3:
    redo(form)
elif query == 4:
    fourthScript(form)
else:
     print("</head><body><pre>Invalid query: %s"%(query))
     for entry in form:
       print entry
       print form[entry].value
       print
     print "</pre></body></html>"
       
