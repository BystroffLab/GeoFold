#!/bin/csh
setenv OUTDIR GeoFold/output

## echo '<\!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">'
echo 'Content-type: text/html'
echo ""
echo ""
echo '<HTML>'
echo '<HEAD><TITLE>Recent GeoFold simulations</TITLE>'
echo '<LINK REV=MADE HREF="mailto:bystrc@rpi.edu">'
echo '<BASE HREF="http://www.bioinfo.rpi.edu/geofold/">'
echo '</HEAD>'
echo '<BODY background="geofold_backdrop.jpg">'
echo '<center>'
echo '<table><tr>'
/bin/ls -1t $OUTDIR/*[0-9][0-9][0-9][0-9]/*.html | \
sed -e "s#"$OUTDIR"/##" |  \
awk -F / 'BEGIN{n=0}{n++;if(n==5){printf "</tr>\n<tr>\n";n=0;};printf "<td><A HREF=\"'$OUTDIR'/%s/%s.html\">%s</A></td>\n", $1, $1, $1}' 
echo '</tr></table>'
echo '<h3><A HREF="index.html">Back to GeoFold server</A></h3>'
echo '</center>'
echo '</BODY>'
echo '</HTML>'
