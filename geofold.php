<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<html>
<head><title>GEOFOLD Protein Unfolding Pathways</title>
<link rev=made href="mailto:bystrc@rpi.edu">
</head>
<BODY background="geofold_backdrop.jpg">
<form method="post" action="./geocgi.cgi?1" enctype="multipart/form-data">
<center>
<table width=100%>
  <tr>
    <td>
      <center><h1><font color="#23AA13"><i>GeoFold server</i></font></h1></center>
      <font color="#FF5555">BETA version!!!</font>
<br>Send bug reports to bystrc@rpi.edu
<br>SYSTEM UPDATE IN PROGRESS. Switching to python.... Please be patient.
<br><font size=1>
Tue Jun 17 14:49:59 EDT 2014
</font>
    </td>
    <td>
      <p><img alt="GeoFold Logo" width=180 src="GFlogo.jpg" align="bottom"><br>
      <a href="http://www.rpi.edu">
      <font color="#FF0000">
      RPI
      </font>
      </a>
    </td>
  </tr>
</table>
</center>
<p> Email address:
<input type="text" name="email_address" value="" size=40>
<p>Key word(s):
<input type="text" name="keyword" value="" placeholder="Enter a unique id for this job" size=80>
<p>Enter a 4-letter <a href="http://ww.pdb.org">PDB</a> identifier
<input type="text" name="pdbid" value="" placeholder="ex. 2b3p" size=4>
...or upload a file in PDB format
<input type="file" name="pdbfile" value="" maxlength=10000 size=20>
<INPUT TYPE="hidden" NAME="params" VALUE="/bach1/home/bystrc/server/geofold/bin/parameters">
<!--<INPUT TYPE="hidden" NAME="plus" VALUE="+">
<INPUT TYPE="hidden" NAME="slash" VALUE="/">
<INPUT TYPE="hidden" NAME="backslash" VALUE="\">
<INPUT TYPE="hidden" NAME="colon" VALUE=":">
<INPUT TYPE="hidden" NAME="space" VALUE=" ">
<INPUT TYPE="hidden" NAME="at" VALUE="@">-->
<p><INPUT TYPE="submit" NAME=".submit" VALUE="Run GeoFold">
</FORM>
<table>
  <tr>
    <td>
          <p><h5>Jobs completed since Tue Jul 23, 2008:&nbsp;<font color='#AA1111'>
	  597
        </font></b></h5>
          </td>
  </tr>
  <tr>
    <td>
      <FORM METHOD="POST" ACTION="./showfinished.cgi" ><INPUT TYPE="submit" NAME=".submit" VALUE="View Recent"></FORM>
    </td>
  </tr>
</table>
<p><A HREF="download.html">Download GeoFOLD programs</A><br>
<p><A HREF="howtoreadit.htm">How to read GeoFOLD output</A><br>
<p><A HREF="settings.html">Expert settings</A><br>
<p>Please cite 
<A HREF="http://onlinelibrary.wiley.com/doi/10.1002/prot.23249/full">
Ramakrishnan V, Srinivasan S, Salem SM, Zaki MJ , Matthews SJ, Colon W, & Bystroff C. (2012) GeoFold: Topology-based protein unfolding pathways capture the effects of engineered disulfides on kinetic stability. Proteins 80(3):920-934. 
</A><br>
Questions?  Contact 
<a href="mailto:bystrc@rpi.edu">bystrc@rpi.edu</a>
<br>
<font size=5>
   <a href="../">Bystroff Lab Home Page</a>
</font>
<br>
<font size=2>Last updated: 
Wed Aug 24 12:08:33 EDT 2011
</font>
<p>Funding from <p><A HREF="grant.html"><img src="nsf.gif" border=0></A>
</BODY>
</HTML>
