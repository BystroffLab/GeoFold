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

"""Reads input from geofold1cgi.py and converts it into a parameters file.
This script will also generate all other associated job files and add it to
the Geofold daemon queue.  Will redirect to the results page
FORM DATA:
lname email pdbcode omega intermediate barrelmoves orange folding halflife
rungeofold reducing molscript breakcut pivotcut hingecut seamcut
breakpointentropy hingepointentropy temperature concentration
voidentropy solidity hbondenergy hammondscale sidechainenropy
hingebarrier pivotbarrier water maxsplit maxtime minseg cavitation"""

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
      out.write("%s %s\n"%(value.upper(),form[value].value))
    except KeyError:
      out.write("%s 0\n"%(value.upper()))
  if not chainset:
    chains += '.'
  out.write("%s\n"%(chains))
  out.close()
  
#HTML header
print "Content-Type: text/html;charset=utf-8"
print
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
#read form data
form = cgi.FieldStorage()
#set url name
try:
  lname = form["lname"].value
except KeyError:
  lname = "%s.%s"%(form["keyword"].value,form["pdbcode"].value)
url = "%s/%s%s.html"%(urldir,outdir,lname)

urlwrite = "http://www.bioinfo.rpi.edu/geofold/%s%s.html"%(outdir,form["lname"].value)

#write HTML output
print('<html><head>')
print('<meta http-equiv="refresh" content="2;url=%s">'%(urlwrite))
print('</head><body><h4>Creating GeoFOLD job. Please wait...</h4><br>')
print('</body></html>')

#new URL settings
out = open(url,'w+')
#reset url name
url = "%s.html"%(form["lname"].value)
#write-out parameters file
file = open("%s/%s.par"%(tmpdir,form["lname"].value),'w+')
makeParameters(form,file)
#get chain info for html results
chains = ''
for value in form:
  if "chain_" in value:
    chains+=form[value].value
  else:
    chains = '.'
#Write out initial HTML results page
out.write('<html><head><title>%s</title>\n'%(url))
out.write('<meta http-equiv="refresh" content="4;url=%s">\n'%(url))
out.write('</head><body>\n')
out.write('<h4>Your GeoFold job is in the queue but has not yet started.</h4>\n')
out.write('<p>Your input coordinates were uploaded as filename %s chains %s\n'%(form["pdbcode"].value,chains))
out.write('<p>Notifications will be sent to %s\n'%(form["email"].value))
out.write('<p><a href="%s">BOOKMARK THIS PAGE</a>:\n'%(url))
out.write('<p>Stay on this page, or return to this page to see your results,')
out.write('which will include the following:<br><ul>\n')
out.write('<li>The timecourse of unfolding as concentrations of Folded, Unfolded, and Intermediate states.</li>\n')
out.write('<li>An unfolding pathway in the form of a clickable pathway tree.</li>\n')
out.write('<li>An Age Plot expressing the order of contact loss during unfolding.</li>\n')
out.write('<li>Unfolding/folding kinetics simulations and associated plots.</li>\n')
out.write('</ul>')
out.write('<p>You will see results for multiple values of &omega; (virtual denaturant).\n')
out.write('<p>You will have the option to Do Over, changing the &omega; values or any other parameters.\n')
out.write('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/settings.html">Click on this page')
out.write('to see a brief explanation of the parameters.</a>\n')
out.write('<p><a href="http://www.bioinfo.rpi.edu/bystrc/geofold/howtoreadit.htm">Click on this page')
out.write('to see a guide to GeoFold output.</a>\n')
out.write('<p><img src="%s"><br>\n'%(waitimage))
out.write('</body></html>')
out.close()
#write job file
job = open("%s/%s.job"%(jobdir,form["lname"].value),'w+')
job.write('%s %s %s'%(form["pdbcode"].value,chains,form["lname"].value))
try:
  oname = form["oname"].value
except KeyError:
  oname = ''
job.write(' %s'%(oname))
job.close()
#deprotect files
commands.getstatusoutput('chmod 0777 %s/%s'%(outdir,url))
commands.getstatusoutput('chmod 0777 %s/%s.job'%(jobdir,form["lname"].value))
commands.getstatusoutput('chmod 0777 %s/%s.par'%(tmpdir,form["lname"].value))
