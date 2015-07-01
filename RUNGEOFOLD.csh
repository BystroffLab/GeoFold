#### RUNGEOFOLD server copy ####
# Last modified: Wed Jan 11 10:27:23 EST 2012
# This script is run by the geod.csh server script.
# By modifying the directories below, you can
# run this script and generate HTML output.
# C.Bystroff  29-NOV-2008
## MODIFIED: this version runs dot command remotely on bach1,
## because the version on bach1 is more up to date, allowing -Timap_np
## cb 6-aug-10
########################################################
#####           COMMAND LINE PARAMETERS           ######
########################################################
# If no arguments are provided, then the script expects
# PDBCODE CHAIN LNAME and ONAME to be already set,
# i.e. by prepending setenv lines. see makeqcom
if ($#argv == 3) then
  # PDBCODE is the basename of a PDB file to be found in TMPDIR
  setenv PDBCODE $1
  # CHAIN is a set of chain IDs within the PDB file, or _ for space.
  setenv CHAIN $2
  # LNAME is a unique name for this job. Used as the root of file names.
  setenv LNAME $3
  setenv ONAME $LNAME
else if ($#argv >= 4 ) then
  # ONAME is a unique name for a previous job, to 
  # be used to skip dag calculation step.
  setenv ONAME $4
else if ($#argv > 0 ) then
  echo "Usage: RUNGEOFOLD.CSH <PDB file basename> <chain identifiers> <unique name> <old unique name>"
  exit
endif
# PID is the process ID, unique to this job
setenv PID $$
## Report the command line or prepend data. Die if not set.
echo "pid=$PID"
echo "$PDBCODE $CHAIN $LNAME $ONAME"
if ($status) exit
setenv TARG $PDBCODE$CHAIN
########################################################
#####              DIRECTORIES                    ######
#####   THESE NEED TO BE SET UPON INSTALLATION    ######
########################################################
# URLDIR
setenv URLDIR /bach1/home/bystrc/public_html
# BASEDIR
setenv BASEDIR /bach1/home/bystrc/server
# THISDIR is the directory where this script originated.
setenv THISDIR $BASEDIR/geofold/bin
# SRCDIR is the directory under which are installed packages
setenv SRCDIR $BASEDIR/geofold/src
# GDIR is the geofold directory created by tar -zxvf geofold.tgz
setenv GDIR $SRCDIR/geofold
# LOGDIR is a place to hold the log file.
setenv LOGDIR $BASEDIR/geofold/log
# HTMLDIR is a public_html directory where HTML
# file may be viewed.
setenv HTMLDIR $URLDIR/geofold/output
# JOBDIR holds job files that are read by geod.csh
# (not used in this script except to clean up)
setenv JOBDIR $URLDIR/geofold/jobs
# TMPDIR is a temporary file directory. Everything in that
# directory can be deleted later. The PDB file should be copied
# to TMPDIR to start.
setenv TMPDIR $BASEDIR/geofold/tmp
# OUTPUTURL is the directory under BASEURL where all server output is sent.
setenv OUTPUTURL "output/"
# BASEURL is the full url of the site containing the output/ directory
setenv BASEURL "http://www.bioinfo.rpi.edu/geofold/"
########################################################################
#####      SCRIPTS AND PROGRAMS REQUIRED BY THE SERVER            ######
#####  install/reset these utilities when installing the server.  ######
########################################################################
## ImageMagick convert
## To install convert and other Imagemagick tools, 
## go to http://www.imagemagick.org/script/index.php
setenv CONVERT /usr/bin/convert
## GraphViz dot program. 
## To install dot, go to http://www.graphviz.org/Download_linux_rhel.php
## setenv DOT 'ssh bach1 /usr/bin/dot'
setenv DOT 'ssh bach1 /usr/bin/dot'
## gnuplot is called within these csh scripts
## To install gnuplot, go to: http://www.gnuplot.info/
setenv CREATEGNUPLOT $THISDIR/creategnuplot.csh
setenv CREATEPROFILE $THISDIR/energyprofile.csh
setenv ALLPROFILES $THISDIR/energyprofile_all.csh
setenv CREATETRACE $THISDIR/createtrace.csh
## Molscript
## To install molscript and molauto, go to: http://www.avatar.se/molscript/
setenv MOLSCRIPT /bach1/home/bystrc/bin/molscript
setenv MOLAUTO /bach1/home/bystrc/bin/molauto
setenv MAKEMOL $THISDIR/makemolscriptmovie.csh
## UNIX (overrides login script settings)
## Please make sure these paths are correct on
## your system. This may be unnecessary in most
## cases, but it will override any re-settings of
## these Unix commands that may be present in the
## startup scripts. I hope "setenv" is never
## re-set! 
setenv MV '/bin/mv -f'
setenv CP '/bin/cp -f'
setenv LS /bin/ls
setenv CAT /bin/cat
setenv GREP /bin/grep
setenv RM '/bin/rm -f'
## note: add awk
#########################################
#####       GEOFOLD PROGRAMS        #####
#########################################
setenv XGEOFOLD $GDIR/xgeofold
setenv XUNFOLDSIM $GDIR/xunfoldsim
setenv GETCHAIN $GDIR/xgetchain
setenv GETSEQ $GDIR/x3to1
setenv MAXTRAFFIC $GDIR/maxTraffic
setenv RENUMBER  $GDIR/xrenumber_one 
setenv PDB2CIJ $GDIR/xpdb2cij
setenv PDB2HB $GDIR/xpdb2hb
setenv VOIDMASK $GDIR/masker/xvoidmask
setenv CONTACTMASK $GDIR/masker/xcontactmask
setenv PATHWAY2PS $GDIR/xpathway2ps
setenv FITPOLY $GDIR/xfit_poly 
setenv PDB2SEAMS  $GDIR/seams/xpdb2seams 
####################################################
##### FILES INCLUDED WITH THE GEOFOLD SERVER. ######
########################################################################
# PARAMTEMPLATE is the template for the parameters file. Don't touch this.
# other files are programs, output files etc.
setenv PARAMTEMPLATE ${THISDIR}/parameters
## HTML headers, w/ and w/o refresh. The refresh version
## keeps refreshing, until the page is replaced with the
## non-refresh version.
setenv HTMLREFRESH ${THISDIR}/header_refresh.html 
setenv HTMLHEAD ${THISDIR}/header.html 
###########################################
##### FILES GENERATED BY THE SERVER #######
###########################################
setenv PARAMFILE ${TMPDIR}/${LNAME}.par
setenv HTMLTMP ${TMPDIR}/${LNAME}.html	
setenv HTMLLOG ${LOGDIR}/${LNAME}.log
setenv HTMLOUT ${HTMLDIR}/${LNAME}.html	
setenv HTMLPERM ${TMPDIR}/${LNAME}.perm
setenv MAILOUT ${TMPDIR}/${LNAME}.eml
#####################
##### VARIABLES #####
#####################
# masker_setup contains filenames and variables
# for voidmask and contactmask
source $GDIR/masker/masker_setup.csh
# MTCUT is the minimum probability of a pathway, for plotting
# by maxTraffic
setenv MTCUT 0.1
# OMEGA is the desolvation free energy
setenv OMEGA 30
# MSPL is the maximum number of splits for each intermediate
#   of unfolding.
setenv MSPL 4
# FING (folding) = 1 is folding, 0 is unfolding
setenv FING 0
# HLFE (halflife) = 0 to go to equilibrium in UnfoldSim.
# HLFE (halflife) = 1 for stopping at the half-life, measured using [U]
# HLFE (halflife) = 2 for stopping at the half-life, measured using [F]
setenv HLFE 0
# MAXRES is the value of maxres in geofold_global.f90
setenv MAXRES 1600

if !(-e $PARAMFILE) then
  echo "$PARAMFILE not found. Copying $PARAMTEMPLATE instead."
  $CP $PARAMTEMPLATE $PARAMFILE
  chmod +w $PARAMFILE
endif
setenv EMAIL `$GREP ^"EMAIL" $PARAMFILE | awk '{print $2}'`
echo "EMAIL $EMAIL"
setenv KEYWORD `$GREP ^"KEYWORD" $PARAMFILE | awk '{print $2}'`
echo "KEYWORD $KEYWORD"
setenv VOIDS `$GREP ^"VOIDENTROPY" $PARAMFILE | awk '{printf "%d", $2}'`
echo "VOIDENTROPY $VOIDS"
setenv FING `$GREP ^"FOLDING" $PARAMFILE | awk '{print $2}'`
echo "FOLDING $FING"
setenv HLFE `$GREP ^"HALFLIFE" $PARAMFILE | awk '{print $2}'`
echo "HALFLIFE $HLFE"
setenv BCUT `$GREP ^"BREAKCUT" $PARAMFILE | awk '{print $2}'`
echo "BREAKCUT $BCUT"
setenv PCUT `$GREP ^"PIVOTCUT" $PARAMFILE | awk '{print $2}'`
echo "PIVOTCUT $PCUT"
setenv HCUT `$GREP ^"HINGECUT" $PARAMFILE | awk '{print $2}'`
echo "HINGECUT $HCUT"
setenv INTM `$GREP ^"INTERMEDIATE" $PARAMFILE | awk '{print $2}'`
echo "INTERMEDIATE $INTM"
setenv WAT `$GREP ^"WATER" $PARAMFILE | awk '{print $2}'`
echo "WATER $WAT"
setenv SEAMS `$GREP ^"BARRELMOVE" $PARAMFILE | awk '{print $2}'`
echo "BARRELMOVE $SEAMS"
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if ($#ORANGE == 0) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
if (`$GREP -c ^"MOLSCRIPT" $PARAMFILE`) then
  setenv MOVIE `$GREP ^"MOLSCRIPT" $PARAMFILE | awk '{print $2}'`
else
  setenv MOVIE " "
endif
echo "MOLSCRIPT $MOVIE"
if (`$GREP -c ^"MAXSPLIT" $PARAMFILE`) then
  setenv MSPL `$GREP ^"MAXSPLIT" $PARAMFILE | awk '{print $2}'`
  if (-e $GDIR/xgeofold_split$MSPL) then
    setenv XGEOFOLD $GDIR/xgeofold_split$MSPL
  endif
endif
###### Do again based on previous run
if !( $ONAME == $LNAME ) then
  $CP ${TMPDIR}/$ONAME.dag ${TMPDIR}/$LNAME.dag
  $CP ${TMPDIR}/$ONAME.pdb ${TMPDIR}/$LNAME.pdb
  $CP ${TMPDIR}/$ONAME.cij ${TMPDIR}/$LNAME.cij
  $CP ${TMPDIR}/$ONAME.seq ${TMPDIR}/$LNAME.seq
  $CP ${TMPDIR}/$ONAME.hb  ${TMPDIR}/$LNAME.hb
  $CP ${HTMLDIR}/$ONAME.gif ${HTMLDIR}/$LNAME.gif
endif
if (`$GREP -c ^"RUNGEOFOLD" $PARAMFILE`) then
  setenv DOIT `$GREP ^"RUNGEOFOLD" $PARAMFILE | awk '{print $2}'`
  if ( $DOIT ) then
    echo "Running GeoFOLD"
  else
    echo "Skipping GeoFOLD"
    $CAT $HTMLREFRESH > $HTMLTMP
    ## echo "<h5>RUNGEOFOLD=0</h5>" >> $HTMLTMP
    echo "<h3>Skipping GeoFOLD DAG calculation.</h3>" >> $HTMLTMP
    echo "<h4>Using previously calculated DAG file $ONAME </h4>" >> $HTMLTMP
    echo "<p>Changes in GeoFOLD parameters (POINTENTROPY, REDUCING) ignored. " >> $HTMLTMP
    echo "<p>Running on "`hostname`"  "`date` >> $HTMLTMP
    echo "<pre>" >> $HTMLTMP
    setenv REDUCING 0
    if (`$GREP -c ^"REDUCING" $PARAMFILE`) then
      setenv REDUCING `$GREP ^"REDUCING" $PARAMFILE | awk '{print $2}'`
    endif
    $CAT $HTMLTMP > $HTMLOUT
    if ( $ONAME == $LNAME ) then
      echo "geod.csh DEAMON ERROR Need old job name as arg 4 to skip geofold."
      echo "geod.csh DEAMON ERROR Need old job name as arg 4 to skip geofold."  >> $HTMLOUT
      echo '</body></html>' >> $HTMLOUT
      exit
    endif
    echo '</body></html>' >> $HTMLOUT
    goto SKIPGEOFOLD
  endif
else
  echo 'RUNGEOFOLD line is missing!'
  echo 'RUNGEOFOLD line is missing!' >> $HTMLOUT
  echo '</body></html>' >> $HTMLOUT
  exit
endif
$CAT $HTMLREFRESH > $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT

echo "<h3>GeoFOLD results for $EMAIL </h3>" >> $HTMLTMP
echo "<pre>" >> $HTMLTMP
echo "Job $LNAME" >> $HTMLTMP
if !( $ONAME == $LNAME ) then
  echo 'Based on <A HREF="'${BASEURL}${OUTPUTURL}${ONAME}.html'">'${ONAME}'</A>' >> $HTMLTMP
endif
echo "Unfolding $TARG" >> $HTMLTMP
echo "Running on "`hostname` >> $HTMLTMP
echo "============= PARAMETERS =============" 
echo "============= PARAMETERS =============" >> $HTMLTMP
if !(-e $PARAMFILE) then
  $CP $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file not found"
  echo "Parameters file not found $PARAMFILE" >> $HTMLTMP
else if (-z $PARAMFILE) then
  $CP $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file is empty"
  echo "Parameters file is empty $PARAMFILE" >> $HTMLTMP
else
  echo "OK"
endif
echo "Using these parameters" >> $HTMLTMP
$CAT $PARAMFILE >> $HTMLTMP
echo "============= checking for PDB file ============="
if (-e ${TMPDIR}/$PDBCODE.error ) then
    echo "<p>" >>  $HTMLTMP
    $CAT ${TMPDIR}/$PDBCODE.error >>  $HTMLTMP
    $CAT $HTMLTMP > $HTMLOUT
    echo '</PRE></body></html>' >> $HTMLOUT
    $RM ${TMPDIR}/$PDBCODE.error
    echo "Problem."
    exit
else if !(-e ${TMPDIR}/$PDBCODE.pdb) then
    echo "<p>ERROR. PDB file is missing. " >>  $HTMLTMP
    $CAT $HTMLTMP > $HTMLOUT
    echo '</PRE></body></html>' >> $HTMLOUT
    echo "Problem."
    exit
else
    echo "OK"
endif
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= GETCHAIN =============" 
echo "============= GETCHAIN =============" >> $HTMLTMP
echo "Extracting protein atoms from $PDBCODE " >> $HTMLTMP
echo "REMARK Extracting protein atoms from $PDBCODE " > ${TMPDIR}/$LNAME.pdb
(	$GETCHAIN + < ${TMPDIR}/$PDBCODE.pdb \
        > ${TMPDIR}/$LNAME.tmp	|| \
        ( echo "Error in GETCHAIN" ; exit )  ) >>& $HTMLTMP
echo "Extracting chains $CHAIN " >> $HTMLTMP
(	$GETCHAIN $CHAIN < ${TMPDIR}/$LNAME.tmp \
        >> ${TMPDIR}/$LNAME.pdb	|| \
        ( echo "Error in GETCHAIN" ; exit )  ) >>& $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
if !($status) echo "OK"
echo "============= RENUMBER =============" 
echo "============= RENUMBER =============" >> $HTMLTMP
(  $RENUMBER ${TMPDIR}/$LNAME.pdb ${TMPDIR}/$LNAME.tmp || \
   ( echo "error in renumber "; exit ) ) >>& $HTMLTMP
  $MV ${TMPDIR}/$LNAME.tmp ${TMPDIR}/$LNAME.pdb   
  $CP ${TMPDIR}/$LNAME.pdb  $HTMLDIR/$LNAME.pdb
  if (-z ${TMPDIR}/$LNAME.pdb ) then
    echo "ERROR. empty file. " >>  $HTMLTMP
    $CAT $HTMLTMP > $HTMLOUT
    echo '</PRE></body></html>' >> $HTMLOUT
    exit
  else
    echo "OK"
  endif
  echo 'Using these <A HREF="'$LNAME.pdb'">coordinates.</A>' >> $HTMLTMP
  $CAT $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
#################### CREATE ROTATING GIF IMAGE ##########################
if !( $MOVIE ) then
  echo "============= not running molscript =============" 
  echo "============= not running molscript =============" >> $HTMLTMP
  $CP ${HTMLDIR}/$MOVIE.gif $HTMLDIR/$LNAME.gif
else
  echo "============= run molscript using $LNAME =============" 
  echo "============= run molscript using $LNAME =============" >> $HTMLTMP
  if ( $ONAME == $LNAME ) then
    $MAKEMOL ${TMPDIR}/$LNAME $HTMLDIR &
  ## $MV ${TMPDIR}/$ONAME.gif $HTMLDIR/
  ## touch $HTMLDIR/$ONAME.gif
  endif
endif
echo "============= 3to1 (extract sequence) =============" 
echo "============= 3to1 (extract sequence) =============" >> $HTMLTMP
(	$GETSEQ "+" < ${TMPDIR}/$LNAME.pdb > ${TMPDIR}/$LNAME.seq || \
        ( echo "3to1 ended with errors" ; exit )    ) >>& $HTMLTMP
        echo "Complete sequence (all chains ) is"  >>& $HTMLTMP
        $CAT ${TMPDIR}/$LNAME.seq  >>& $HTMLTMP
  $CAT $HTMLTMP > $HTMLOUT
  setenv RRR `grep ^"#" ${TMPDIR}/$LNAME.seq | head -1 | awk '{print $(NF-4)}'`
  if ($RRR > $MAXRES) then
    echo '</PRE><h2>Maximum number of residues exceeded: '$RRR' out of '$MAXRES'"</h2></body></html>' >> $HTMLOUT
    exit 1
  endif
  echo '</PRE></body></html>' >> $HTMLOUT

echo "============= PDB2CIJ (extract contacts) =============" 
echo "============= PDB2CIJ (extract contacts) =============" >> $HTMLTMP
(	$PDB2CIJ ${TMPDIR}/$LNAME.pdb 8. \
        > ${TMPDIR}/$LNAME.cij || \
        ( echo "Error in PDB2CIJ " ; exit )    ) >>& $HTMLTMP
 echo "OK"
 echo -n "Number of contacts found: " >>& $HTMLTMP
 wc -l ${TMPDIR}/$LNAME.cij | awk '{print $1}' >>& $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= PDB2HB (extract H-bonds SS-bonds) =============" 
echo "============= PDB2HB (extract H-bonds SS-bonds) =============" >> $HTMLTMP
  setenv REDUCING 0
  if (`$GREP -c ^"REDUCING" $PARAMFILE`) then
    setenv REDUCING `$GREP ^"REDUCING" $PARAMFILE | awk '{print $2}'`
  endif
(	$PDB2HB ${TMPDIR}/$LNAME.pdb \
        > ${TMPDIR}/$LNAME.hb || \
        ( echo "Error in PDB2HB " ; exit )    ) >>& $HTMLTMP
 echo "OK"
 echo -n "Number of H-bonds found: " >>& $HTMLTMP
 $GREP -c " H " ${TMPDIR}/$LNAME.hb >>& $HTMLTMP
 echo -n "Number of SS-bonds found: " >>& $HTMLTMP
 $GREP -c " S " ${TMPDIR}/$LNAME.hb >>& $HTMLTMP
 $GREP " S " ${TMPDIR}/$LNAME.hb | awk '{printf "(%d,%d)",$1,$2}' >>& $HTMLTMP
 echo "HBONDS ${TMPDIR}/$LNAME.hb " >> $PARAMFILE
 if ($REDUCING) then
   awk '$(NF)!="S"' ${TMPDIR}/$LNAME.hb > ${TMPDIR}/$LNAME.tmp
   $MV ${TMPDIR}/$LNAME.tmp ${TMPDIR}/$LNAME.hb
   echo "<br>REDUCING conditions. SS-bonds removed. " >>& $HTMLTMP
 endif
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= CONTACTMASK =============" 
$CAT $HTMLTMP > $HTMLOUT
echo "============= currently running CONTACTMASK =============" >> $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
echo "============= CONTACTMASK =============" >> $HTMLTMP
echo "	$CONTACTMASK ${TMPDIR}/$LNAME.pdb ${TMPDIR}/$LNAME.sas 1.4 " >>& $HTMLTMP
(	echo "Time before running CONTACTMASK `date`"    ) >>& $HTMLTMP
(	$CONTACTMASK ${TMPDIR}/$LNAME.pdb \
        ${TMPDIR}/$LNAME.sas 1.4 || \
        ( echo "Error in CONTACTMASK " ; exit )   ) >>& $HTMLTMP
(	echo "Time after running CONTACTMASK `date`" ) >>& $HTMLTMP
  echo "OK"
  $CP ${TMPDIR}/$LNAME.sas $HTMLDIR/$LNAME.sas
  echo 'Pairwise contact surfaces <A HREF="'$LNAME.sas'">file.</A>' >> $HTMLTMP
  echo "CONTACTS ${TMPDIR}/$LNAME.sas " >> $PARAMFILE
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT


if ($SEAMS) then
  ### Seams preprocessing. Added by Luis Garreta (June/04/2013)
  ### Seams with energy. Added by Luis Garreta (June/04/2013)
  echo "============= BARRELMOVES =============" 
  # Three parameters: pdb, sas, path geofold (stride)
  echo "============= BARRELMOVES =============" >> $HTMLTMP
  (	echo "Time before running BARRELMOVES `date`") >>& $HTMLTMP
  echo "$PDB2SEAMS ${TMPDIR}/$LNAME.pdb ${TMPDIR}/$LNAME.sas"' > '"${TMPDIR}/$LNAME.seams" >>& $HTMLTMP
  (	$PDB2SEAMS ${TMPDIR}/$LNAME.pdb ${TMPDIR}/$LNAME.sas > ${TMPDIR}/$LNAME.seams || \
       ( echo "Error in PDB2SEAMS " ; exit )    ) >>& $HTMLTMP
  echo "OK"
  echo "Checking for a beta barrel " >>& $HTMLTMP
  cat $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "SEAMS ${TMPDIR}/$LNAME.seams" >> $PARAMFILE
else
  echo "============= NO BARRELMOVES calculated =============" 
  grep -v ^"SEAMS" $PARAMFILE > /tmp/junk
  mv /tmp/junk $PARAMFILE
endif

if ($VOIDS) then
  echo "============= VOIDMASK =============" 
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running VOIDMASK =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= VOIDMASK =============" >> $HTMLTMP
  (	echo "Time before running VOIDMASK `date`") >>& $HTMLTMP
  (	$VOIDMASK ${TMPDIR}/$LNAME.pdb \
        ${TMPDIR}/$LNAME.void 1.4 1.2 1.4 || \
        ( echo "Error in VOIDMASK " ; exit )   ) >>& $HTMLTMP
  (	echo "Time after running VOIDMASK `date`" ) >>& $HTMLTMP
  echo "OK"
  if (-z ${TMPDIR}/$LNAME.void ) then
      (	echo "ERROR in VOIDMASK . Empty file. " ) >>& $HTMLTMP
    $CAT $HTMLTMP > $HTMLOUT
    echo '</PRE></body></html>' >> $HTMLOUT
    exit
  endif
  $CP ${TMPDIR}/$LNAME.void  ${HTMLDIR}/$LNAME.void.pdb
  echo 'Coordinates with void positions: <A HREF="'$LNAME.void.pdb'">'$LNAME.void.pdb'</A>' \
      >>& $HTMLTMP
else
  $CP ${TMPDIR}/$LNAME.pdb ${TMPDIR}/$LNAME.void
  echo 'Void positions not calculated' >>& $HTMLTMP
endif
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= GEOFOLD =============" 
$CAT $HTMLTMP > $HTMLOUT
echo "============= currently running GEOFOLD, This may take time, please be patient =============" >> $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
echo "============= GEOFOLD =============" >> $HTMLTMP
echo "Time before running GEOFOLD `date`"  >>& $HTMLTMP
echo "$XGEOFOLD  ${TMPDIR}/$LNAME.void ${TMPDIR}/$LNAME.dag ${PARAMFILE} "  >>& $HTMLTMP
(	$XGEOFOLD ${TMPDIR}/$LNAME.void ${TMPDIR}/$LNAME.dag \
        ${PARAMFILE} || \
        ( echo "Error in GEOFOLD " ; exit )   ) >>& $HTMLTMP
echo "OK"
echo "Time after running GEOFOLD `date`" >>& $HTMLTMP

if (-z ${TMPDIR}/$LNAME.dag ) then
      (	echo "ERROR in GEOFOLD . Empty file. " ) >>& $HTMLTMP
   $CAT $HTMLTMP > $HTMLOUT
   echo '</PRE></body></html>' >> $HTMLOUT
   exit
endif
if !(-e ${TMPDIR}/$LNAME.dag ) then
   echo "ERROR in GEOFOLD . No output file. " >>& $HTMLTMP
   $CAT $HTMLTMP > $HTMLOUT
   echo '</PRE></body></html>' >> $HTMLOUT
   exit
endif
  $CP ${TMPDIR}/$LNAME.dag  ${HTMLDIR}/$LNAME.dag
  ## $CP ${TMPDIR}/$LNAME.dag  ${TMPDIR}/$ONAME.dag 
  echo 'Directed acyclic graph (DAG) <A HREF="'$LNAME.dag'">file.</A>' >>& $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

SKIPGEOFOLD:
echo "============= UNFOLDSIM =============" 
echo "Unfolding $TARG" >> $HTMLTMP
echo "============= UNFOLDSIM =============" >> $HTMLTMP
@ NN = 0
while ( $#ORANGE )
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running UNFOLDSIM, omega=$ORANGE[1], please wait a few minutes =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  sed -e "s/^OMEGA .*/OMEGA $ORANGE[1]/" $PARAMFILE > ${PARAMFILE}.1
  @ NN ++
  setenv LOGFILE $TMPDIR/${LNAME}_$NN.log
  echo "============= run $NN omega=$ORANGE[1] =============" 
  echo "============= run $NN omega=$ORANGE[1] =============" >> $HTMLTMP
  echo "Time before running UNFOLDSIM $NN `date`"    >>& $HTMLTMP
  echo "$XUNFOLDSIM ${TMPDIR}/${LNAME}.dag ${PARAMFILE}.1 > $LOGFILE " >> $HTMLTMP
  $XUNFOLDSIM ${TMPDIR}/${LNAME}.dag ${PARAMFILE}.1 > $LOGFILE || \
    ( echo "ERROR when running UNFOLDSIM " ; exit )
  echo "Time after running UNFOLDSIM OMEGA=$ORANGE[1] : `date`"         >>& $HTMLTMP
  echo "<p><pre>"  >>& $HTMLTMP
  $GREP ^"TIMECOURSE" $LOGFILE | tail -1 >>& $HTMLTMP
  $CAT $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT

  echo "============= PATHWAY2PS =============" 
  $CAT $HTMLTMP > $HTMLOUT
  ## Skip the rest if the unfold simulation failed.
  if !(-e ${TMPDIR}/${LNAME}.dag.path) then
    echo "ERROR when running UNFOLDSIM. No pathway data for OMEGA=$ORANGE[1]"
    echo "ERROR when running UNFOLDSIM. No pathway data for OMEGA=$ORANGE[1]" >>& $HTMLTMP
    shift ORANGE
    continue
  else
    $CP ${TMPDIR}/${LNAME}.dag.path ${HTMLDIR}/${LNAME}_${NN}.path
  endif
  if (-e ${TMPDIR}/${LNAME}.dag.age) then
    $CP ${TMPDIR}/${LNAME}.dag.age ${HTMLDIR}/${LNAME}_${NN}.age
  endif
  echo "============= currently running PATHWAY2PS, should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  ## Age plot: colored contact map
  echo "============= PATHWAY2PS =============" >> $HTMLTMP
  echo "$PATHWAY2PS ${TMPDIR}/${LNAME}.seq ${TMPDIR}/${LNAME}.dag.path ${TMPDIR}/${LNAME}.cij ${TMPDIR}/${LNAME}.ps 4"	 >> $HTMLTMP
  $PATHWAY2PS ${TMPDIR}/${LNAME}.seq ${TMPDIR}/${LNAME}.dag.path ${TMPDIR}/${LNAME}.cij ${TMPDIR}/${LNAME}.ps 4	 >>& $HTMLTMP
  $CAT $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= CONVERT =============" 
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running CONVERT , should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= CONVERT =============" >> $HTMLTMP
  echo "$CONVERT -trim -geometry 100 ${TMPDIR}/$LNAME.ps  ${TMPDIR}/${LNAME}_age_${NN}_thumb.gif"  >> $HTMLTMP
  $CONVERT -trim -geometry 100 ${TMPDIR}/$LNAME.ps  ${TMPDIR}/${LNAME}_age_${NN}_thumb.gif
  $CONVERT -trim -border 5 -bordercolor black ${TMPDIR}/${LNAME}_age_${NN}_thumb.gif ${HTMLDIR}/${LNAME}_age_${NN}_thumb.jpg
  $CONVERT -trim -geometry 600 ${TMPDIR}/${LNAME}.ps  ${TMPDIR}/${LNAME}_age_${NN}.gif
  $CONVERT -trim ${TMPDIR}/${LNAME}_age_${NN}.gif  ${HTMLDIR}/${LNAME}_age_${NN}.jpg
  ###### Draw DAG image ######
  echo "============= MAXTRAFFIC =============" 
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running MAXTRAFFIC , should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= MAXTRAFFIC ============="  >> $HTMLTMP
  echo "$MAXTRAFFIC ${TMPDIR}/${LNAME}.dag.out $MTCUT ${OUTPUTURL}${LNAME}_${NN}.dag.out > ${TMPDIR}/$LNAME.dot"  >> $HTMLTMP
  $MAXTRAFFIC ${TMPDIR}/${LNAME}.dag.out $MTCUT ${OUTPUTURL}${LNAME}_${NN}.dag.out > ${TMPDIR}/$LNAME.dot
  $CP ${TMPDIR}/${LNAME}.dag.out ${HTMLDIR}/${LNAME}_${NN}.dag.out
  echo "============= DOT =============" 
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running DOT , should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= DOT ============="  >> $HTMLTMP
  echo "$DOT -Tpng  -o ${HTMLDIR}/${LNAME}.dag_${NN}.png -Timap_np -o ${HTMLDIR}/${LNAME}.dag_${NN}.map ${TMPDIR}/${LNAME}.dot"  >> $HTMLTMP
  $DOT -Tpng  -o ${HTMLDIR}/${LNAME}.dag_${NN}.png -Timap_np -o ${HTMLDIR}/${LNAME}.dag_${NN}.map ${TMPDIR}/${LNAME}.dot
  echo "$DOT -Tps  -o ${TMPDIR}/${LNAME}.dot_${NN}.ps  ${TMPDIR}/${LNAME}.dot   "  >> $HTMLTMP
  $DOT -Tps  -o ${TMPDIR}/${LNAME}.dot_${NN}.ps  ${TMPDIR}/${LNAME}.dot   
  echo "============= CONVERT =============" 
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= currently running CONVERT , should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= CONVERT ============="  >> $HTMLTMP
  echo "$CONVERT -trim -geometry 100 ${TMPDIR}/${LNAME}.dot_${NN}.ps ${HTMLDIR}/${LNAME}.dag_${NN}_thumb.png" >> $HTMLTMP
        $CONVERT -trim -geometry 100 ${TMPDIR}/${LNAME}.dot_${NN}.ps ${HTMLDIR}/${LNAME}.dag_${NN}_thumb.png
  ##echo "$CONVERT -trim ${TMPDIR}/${LNAME}.dag_${NN}_thumb.gif ${HTMLDIR}/${LNAME}.dag_${NN}_thumb.jpg" >> $HTMLTMP
  ##      $CONVERT -trim -border 5 -bordercolor green ${TMPDIR}/${LNAME}.dag_${NN}_thumb.gif ${HTMLDIR}/${LNAME}.dag_${NN}_thumb.jpg
  $CAT $HTMLTMP > $HTMLOUT
  echo "============= done running CONVERT =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  ###### NEW subscript ######
  echo "============= energyprofile ============="  
  $CAT $HTMLTMP > $HTMLOUT  
  echo "============= currently running energyprofile, should be quick =============" >> $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= energyprofile =============" >> $HTMLTMP
  echo "$CREATEPROFILE ${TMPDIR}/${LNAME}"   >> $HTMLTMP
  $CREATEPROFILE ${TMPDIR}/${LNAME} ${TMPDIR}/${LNAME}.nrg_${NN}
  $CAT $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  echo "============= CONVERT ============="
  echo "============= CONVERT =============" >> $HTMLTMP
  echo "$CONVERT -trim -geometry 100 ${TMPDIR}/$LNAME.nrg_${NN}.ps  ${TMPDIR}/${LNAME}_nrg_${NN}_thumb.gif"  >> $HTMLTMP
  $CONVERT -trim -geometry 100 ${TMPDIR}/$LNAME.nrg_${NN}.ps  ${TMPDIR}/${LNAME}_nrg_${NN}_thumb.gif
  $CONVERT -trim -border 5 -bordercolor blue ${TMPDIR}/${LNAME}_nrg_${NN}_thumb.gif ${HTMLDIR}/${LNAME}_nrg_${NN}_thumb.jpg
  $CONVERT -trim ${TMPDIR}/${LNAME}.nrg_${NN}.ps  ${TMPDIR}/${LNAME}_nrg_${NN}.gif
  $CONVERT -trim ${TMPDIR}/${LNAME}_nrg_${NN}.gif  ${HTMLDIR}/${LNAME}_nrg_${NN}.jpg
  ##### done with unfolding simulation
  shift ORANGE
end
@ MM = $NN
#################### CREATE FINAL OUTPUT HTML FILE ##########################
## Main output page
$CAT $HTMLHEAD > $HTMLPERM
$CP ${TMPDIR}/$LNAME.pdb  $HTMLDIR/
$CP ${TMPDIR}/$LNAME.seq  $HTMLDIR/
################ CREATE a TABLE ROW of OUTPUT data #####################
echo '<table width="80%">' >> $HTMLPERM
echo '  <tr>' >> $HTMLPERM
echo '    <td colspan="2">' >> $HTMLPERM
echo '      <img src="'${LNAME}.gif'" alt="image not found for '${LNAME}.pdb'. Try re-loading.">' >> $HTMLPERM
echo '    </td>' >> $HTMLPERM
echo '    <h3><A HREF="'${LNAME}.pdb'">Coordinate file ('$TARG') used for this simulation</A></h3>' >> $HTMLTMP
echo '    <td colspan="4">' >> $HTMLPERM
if ( $FING ) then
  echo '    <h1>Folding...<h1> ' >>& $HTMLPERM
  echo '    <h1>Folding...<h1> ' >>& $HTMLTMP
else
  echo '    <h1>Unfolding...<h1> ' >>& $HTMLPERM
  echo '    <h1>Unfolding...<h1> ' >>& $HTMLTMP
endif
echo '      <h4>'$LNAME'</h4>' >> $HTMLPERM
if !( $ONAME == $LNAME ) then
  echo '      <h4>Based on <A HREF="'${BASEURL}${OUTPUTURL}${ONAME}.html'">'${ONAME}'</A></h4>' >> $HTMLPERM
endif
echo '      <h4><A HREF="'${LNAME}.pdb'">Coordinate file ('$TARG')</A></h4>' >> $HTMLPERM
echo '      <h4><A HREF="'${LNAME}_1.dag.out'">Unfolding graph for '$TARG' </A></h4>' >> $HTMLPERM
echo '      <h4><A HREF="'${LNAME}_1.dag.out'">Unfolding graph for '$TARG' </A></h4>' >> $HTMLTMP
echo -n    "<h5>Number of H-bonds found: " >>& $HTMLPERM
$GREP -c    " H " ${TMPDIR}/$LNAME.hb >>& $HTMLPERM
echo       "</h5>" >>& $HTMLPERM
echo -n    "<h5>Number of SS-bonds found: " >>& $HTMLPERM
$GREP -c    " S " ${TMPDIR}/$LNAME.hb >>& $HTMLPERM
if ($REDUCING) then
   echo " REDUCING conditions. SS-bonds removed. " >>& $HTMLPERM
endif
echo       "</h5>" >>& $HTMLPERM
echo '      <h4>Sequence (all chains)</h4>' >> $HTMLTMP
echo '      <pre>' >> $HTMLPERM
echo '      <pre>' >> $HTMLTMP
$CAT ${TMPDIR}/${LNAME}.seq >> $HTMLPERM
$CAT ${TMPDIR}/${LNAME}.seq >> $HTMLTMP
echo '      </pre>' >> $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '      </PRE></body></html>' >> $HTMLOUT
echo '      </pre>' >> $HTMLPERM
echo '    </td>' >> $HTMLPERM
echo '  </tr>' >> $HTMLPERM
################ CREATE a TABLE ROW of TIMECOURSES #####################
echo '################ CREATE a TABLE ROW of gnuplot TIMECOURSES #####################'
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
echo '  <tr>' >> $HTMLPERM
@ NN = 0
while ( $#ORANGE )
  @ NN ++
  setenv LOGFILE $TMPDIR/${LNAME}_$NN.log
  $CP $LOGFILE $HTMLDIR/
  echo "$CREATETRACE $LOGFILE $INTM" >> $HTMLTMP
  $CREATETRACE $LOGFILE $INTM
  ls -l $LOGFILE*
  echo "$CONVERT -trim -geometry 100 $TMPDIR/${LNAME}_$NN.ps $HTMLDIR/${LNAME}_$NN.thumb.gif" >> $HTMLTMP
  $CONVERT -trim -geometry 100 $TMPDIR/${LNAME}_$NN.ps $HTMLDIR/${LNAME}_$NN.thumb.gif
  echo "$CONVERT -trim $TMPDIR/${LNAME}_$NN.ps $HTMLDIR/${LNAME}_$NN.gif" >> $HTMLTMP
  $CONVERT -trim $TMPDIR/${LNAME}_$NN.ps $HTMLDIR/${LNAME}_$NN.gif
  echo '  <td><A href="'${LNAME}_$NN.log'">timecourse</A> for omega='$ORANGE[1]'<br><A HREF="'${LNAME}_$NN.gif'"><img src="'${LNAME}_$NN.thumb.gif'" width=100 alt="timecourse data missing"></A></td>' >> $HTMLPERM
  shift ORANGE
end
echo '  </tr>' >> $HTMLPERM
################ CREATE a TABLE ROW of age plots #####################
echo '################ CREATE a TABLE ROW of age plots ####################'
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
@ NN = 0
echo '  <tr>' >> $HTMLPERM
while ( $#ORANGE )
  @ NN ++
  echo '<td>Age plot for '$ORANGE[1]' <br><A HREF="'${LNAME}_age_${NN}.jpg'"><img src="'${LNAME}_age_${NN}_thumb.jpg'" width=100 alt="age plot data missing" ></A>' >> $HTMLPERM
  echo '<br><A HREF="'${LNAME}_${NN}.path'">Concentrations at 1/2 folded</A>' >> $HTMLPERM
  echo '<br><A HREF="'${LNAME}_${NN}.age'">Contact ages in seconds</A></td>' >> $HTMLPERM
  shift ORANGE
end
echo '  </tr>' >> $HTMLPERM
################ CREATE a TABLE ROW of DAGs #####################
echo '################ CREATE a TABLE ROW of DAGs ####################'
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
@ NN = 0
echo '  <tr>' >> $HTMLPERM
while ( $#ORANGE )
  @ NN ++
  echo '<td>Pathways at '$ORANGE[1]' <br><A HREF="'${LNAME}.dag_${NN}.html'"><img src="'${LNAME}.dag_${NN}_thumb.png'" width=100 alt="pathway data missing"></A></td>' >> $HTMLPERM
  ## DAG image map file
  echo '<html><body><img src="'${LNAME}.dag_${NN}.png'" border=0 usemap="#dagmap">' >  ${HTMLDIR}/${LNAME}.dag_${NN}.html
  echo "============= convert ismap to usemap =============" 
  ## maxTraffic outputs a map file to be used with the 'ismap' attribute. However, ismap is not compatible with
  ## the current HTTP server. So, the following extracts lines from the map file to construct a usemap, which is appended
  ## to the HTML file.
  echo '<map name="dagmap">' >>  ${HTMLDIR}/${LNAME}.dag_${NN}.html
  grep -v "base referer" ${HTMLDIR}/${LNAME}.dag_${NN}.map | \
  awk '{printf "<AREA HREF=\"%s\" shape=\"rect\" coords=\"%s,%s\"/>\n",$2,$3,$4}' >> ${HTMLDIR}/${LNAME}.dag_${NN}.html
  echo '</map></body></html>' >>  ${HTMLDIR}/${LNAME}.dag_${NN}.html

  shift ORANGE
end
echo '  </tr>' >> $HTMLPERM
################ CREATE a TABLE ROW of energy profiles #####################
echo '################ CREATE a TABLE ROW of energy profiles ####################'
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
@ NN = 0
echo '  <tr>' >> $HTMLPERM
while ( $#ORANGE )
  @ NN ++
  echo '<td>Energy profile for '$ORANGE[1]' <br><A HREF="'${LNAME}_nrg_${NN}.jpg'"><img src="'${LNAME}_nrg_${NN}_thumb.jpg'" width=100 alt="energy data missing"></A></td>' >> $HTMLPERM
  shift ORANGE
end
echo '  </tr>' >> $HTMLPERM
echo '</table>' >> $HTMLPERM
echo '<h3><A HREF="../howtoreadit.htm">How to Read GeoFOLD Output</A></h3>' >> $HTMLPERM

####### Output unfolding rates at different omega values #######
echo '<hr width="80%">' >>& $HTMLPERM
echo -n '<h4>' >>& $HTMLPERM
echo -n '<h4>' >>& $HTMLTMP
if ( $HLFE > 0 ) then
  if ( $FING ) then
    echo 'Times and concentrations at 1/2 folded' >>& $HTMLPERM
    echo 'Times and concentrations at 1/2 folded' >>& $HTMLTMP
  else
    if ( $HLFE > 1 ) then
      echo 'Times and concentrations at <1/2 folded' >>& $HTMLPERM
      echo 'Times and concentrations at <1/2 folded' >>& $HTMLTMP
    else
      echo 'Times and concentrations at >1/2 unfolded' >>& $HTMLPERM
      echo 'Times and concentrations at >1/2 unfolded' >>& $HTMLTMP
    endif
  endif
else
  if ( $FING ) then
    echo 'Equilibrium concentrations and rates starting from unfolded' >>& $HTMLPERM
    echo 'Equilibrium concentrations and rates starting from unfolded' >>& $HTMLTMP
  else
    echo 'Equilibrium concentrations and rates starting from folded' >>& $HTMLPERM
    echo 'Equilibrium concentrations and rates starting from folded' >>& $HTMLTMP
  endif
endif
echo "<p><pre>"  >>& $HTMLPERM
echo "<p><pre>"  >>& $HTMLTMP
echo -n "dGsolvation   "  >>& $HTMLPERM
echo -n "dGsolvation   "  >>& $HTMLTMP
echo "   time(sec)   final[F]  final[U]  final[I]  sum[c](ppb) half_life(sec)" >>& $HTMLPERM
echo "   time(sec)   final[F]  final[U]  final[I]  sum[c](ppb) half_life(sec)" >>& $HTMLTMP
@ NN = 0
set ORANGE=(`$GREP ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`$GREP ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
echo -n "" > $TMPDIR/${LNAME}.plot
while ( $#ORANGE )
  @ NN ++
  setenv LOGFILE $TMPDIR/${LNAME}_$NN.log
  echo -n "$ORANGE[1] kJ/mol/A^2   "  >>& $HTMLPERM
  echo -n "$ORANGE[1] kJ/mol/A^2   "  >>& $HTMLTMP
  $GREP ^"TIMECOURSE" $LOGFILE | tail -1 | sed -e "s/TIMECOURSE//" >>& $HTMLPERM
  $GREP ^"TIMECOURSE" $LOGFILE | tail -1 | sed -e "s/TIMECOURSE//" >>& $HTMLTMP
  if ( $FING ) then
    ## don't make a chevron plot
    ## we don't trust folding rates from geoFold.
  else
    ## HLFE <= 1, find half-life of [U] (increasing) from the trajectory file.
    ## HLFE >= 2, find half-life of [F] (decreasing) from the trajectory file.
    if ($HLFE == 2) then
      ## look at [F] 
      $GREP ^"TIMECOURSE" $LOGFILE | awk '$3>50.' | tail -1 \
      | awk -v omega=$ORANGE[1] '{print omega, log(log(2)/$2)}' >> $TMPDIR/${LNAME}.plot
    else if ($HLFE == 1) then
      ## look at [U] 
      $GREP ^"TIMECOURSE" $LOGFILE | awk '$4<50.' | tail -1 \
      | awk -v omega=$ORANGE[1] '{print omega, log(log(2)/$2)}' >> $TMPDIR/${LNAME}.plot
    else 
      ## folding to equilibrium . Follow [F].
      setenv EQF `$GREP ^"TIMECOURSE" $LOGFILE | tail -1 | awk '{print $3}'`                    
      $GREP ^"TIMECOURSE" $LOGFILE | awk -v eqf=$EQF '$3>(50.+(eqf/2.))' | tail -1 \
      | awk -v omega=$ORANGE[1] '{print omega, log(log(2)/$2)}' >> $TMPDIR/${LNAME}.plot
    endif
  endif
  shift ORANGE
end 
echo "</pre>"  >>& $HTMLPERM
echo "</pre>"  >>& $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT
setenv FIT $TMPDIR/${LNAME}.fit
echo "============= FIT_POLY =============" 
echo "<pre>"  >>& $HTMLTMP
$FITPOLY $TMPDIR/${LNAME}.plot 1 $FIT  >>& $HTMLTMP
setenv IC `$GREP -A 2 "Least-squares" $FIT | tail -1 | awk '{print $2}'`
setenv SL `$GREP -A 3 "Least-squares" $FIT | tail -1 | awk '{print $2}'`
setenv LNKU `awk -v OM=$WAT -v sl=$SL -v ic=$IC 'BEGIN{print sl*OM + ic}'`
echo "</pre>"  >>& $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT
echo -n "<h5>Projected half-life of unfolding " >>& $HTMLPERM
echo -n "<h5>Projected half-life of unfolding " >>& $HTMLTMP
echo "in pure water (dGsolvation=$WAT kJ/mol/A^2) : " >>& $HTMLPERM
echo "in pure water (dGsolvation=$WAT kJ/mol/A^2) : " >>& $HTMLTMP
awk -v lnku=$LNKU 'BEGIN{printf "%e", log(2)/exp(lnku)}' >>& $HTMLPERM
awk -v lnku=$LNKU 'BEGIN{printf "%e", log(2)/exp(lnku)}' >>& $HTMLTMP
echo " seconds</h5><br> " >>& $HTMLPERM
echo " seconds</h5><br> " >>& $HTMLTMP
##### Create gnuplot image for 1/2 chevron plot
echo "Running $CREATEGNUPLOT"
$CREATEGNUPLOT $TMPDIR/${LNAME}  $WAT
$CONVERT -trim $TMPDIR/${LNAME}.gp.ps $TMPDIR/${LNAME}.gp.gif
$CONVERT -trim $TMPDIR/${LNAME}.gp.gif $HTMLDIR/${LNAME}.gp.jpg
$CONVERT -trim -geometry 200 $TMPDIR/${LNAME}.gp.ps $TMPDIR/${LNAME}.gp_thumb.gif
$CONVERT -trim -border 5 -bordercolor yellow $TMPDIR/${LNAME}.gp_thumb.gif $HTMLDIR/${LNAME}.gp_thumb.jpg
##### Create superimposed energy profiles
echo "============= energyprofile =============" >> $HTMLTMP
echo "Running $ALLPROFILES"
$ALLPROFILES $TMPDIR/${LNAME}.nrg  $TMPDIR/${LNAME}.all
$CONVERT -trim $TMPDIR/${LNAME}.all.ps $TMPDIR/${LNAME}.all.gif
$CONVERT -trim $TMPDIR/${LNAME}.all.gif $HTMLDIR/${LNAME}.all.jpg
$CONVERT -trim -geometry 200 $TMPDIR/${LNAME}.all.ps $TMPDIR/${LNAME}.all_thumb.gif
$CONVERT -trim -border 5 -bordercolor red $TMPDIR/${LNAME}.all_thumb.gif $HTMLDIR/${LNAME}.all_thumb.jpg
##### Make a table row with chevron and energy profiles
echo '<table>' >> $HTMLPERM
echo '  <tr>' >> $HTMLPERM
echo '    <td>' >> $HTMLPERM
echo '<h4>ln(ku) vs omega for '$TARG' </h4><A HREF="'${LNAME}.gp.jpg'"><img src="'${LNAME}.gp_thumb.jpg'" alt="data missing"></A><br>' >> $HTMLTMP
echo '<h4>ln(ku) vs omega for '$TARG' </h4><A HREF="'${LNAME}.gp.jpg'"><img src="'${LNAME}.gp_thumb.jpg'" alt="data missing"></A><br>' >> $HTMLPERM
echo '    </td>' >> $HTMLPERM
echo '    <td>' >> $HTMLPERM
echo '<h4>Superimposed energy profiles for '$TARG' </h4><A HREF="'${LNAME}.all.jpg'"><img src="'${LNAME}.all_thumb.jpg'" width=200 alt="data missing"></A><br>' >> $HTMLTMP
echo '<h4>Superimposed energy profiles for '$TARG' </h4><A HREF="'${LNAME}.all.jpg'"><img src="'${LNAME}.all_thumb.jpg'" width=200 alt="data missing"></A><br>' >> $HTMLPERM
echo '    </td>' >> $HTMLPERM
echo '  </tr>' >> $HTMLPERM
echo '</table>' >> $HTMLPERM
echo "<br>Completed `date` <br>" >> $HTMLPERM
echo "<br>Completed `date` <br>" >> $HTMLTMP
$CAT $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT
###### Create "do over" button.
## geofold2.cgi is a cgi fortran program that reads POST stream
## and leaves a job file for the geod.csh deamon to run.
## Source for that program is gf2cgi.f90
###### 
echo '<p>Modify and do over: '                                                     >> $HTMLPERM
echo '<FORM METHOD="POST" ACTION="../geofold2.cgi">'                               >> $HTMLPERM
echo '(<A HREF="../settings.html">Click for more info.</A>):<br>'                  >> $HTMLPERM
echo '<textarea NAME="params" rows=20 cols=60 WRAP="off">'                         >> $HTMLPERM
## LNAME is the unique name of the current job 
echo 'LNAME '${LNAME}                                                              >> $HTMLPERM
$GREP -v ^"LNAME"  $PARAMFILE                                                       >> $HTMLPERM
echo '</textarea>'                                                                 >> $HTMLPERM
echo '<INPUT type="hidden" name="email_address" value="'$EMAIL'">'                 >> $HTMLPERM
echo '<INPUT type="hidden" name="keyword" value="'$KEYWORD'">'                     >> $HTMLPERM
echo '<INPUT type="hidden" name="pdbid" value="'$PDBCODE'">'                       >> $HTMLPERM
echo '<INPUT type="hidden" name="chains" value="'$CHAIN'">'                        >> $HTMLPERM
echo '<INPUT type="hidden" name="slash" value="/">'                                >> $HTMLPERM
echo '<INPUT type="hidden" name="space" value=" ">'                                >> $HTMLPERM
echo '<INPUT type="hidden" name="plus" value="+">'                                 >> $HTMLPERM
echo '<INPUT type="hidden" name="at" value="@">'                                   >> $HTMLPERM
echo '<INPUT type="hidden" name="percent" value="%">'                              >> $HTMLPERM
echo '<INPUT type="hidden" name="lname" value="'$LNAME'">'                         >> $HTMLPERM
echo '<INPUT TYPE="hidden" NAME="unit" VALUE="0">'                                 >> $HTMLPERM
echo '<INPUT TYPE="hidden" NAME="pid" VALUE="'$PID'">'                             >> $HTMLPERM
echo '<p><INPUT TYPE="submit" NAME=".submit" VALUE="Do over">'                     >> $HTMLPERM
echo '</FORM>' >> $HTMLPERM
##### BACK button
echo '<h3><A HREF="../server.php">Back to GeoFold server</A></h3><br>' >> $HTMLPERM
echo '</body></html>' >> $HTMLPERM
$CP $HTMLPERM $HTMLOUT
###### Send notification by email ######
echo "Your GeoFOLD results are available at the following URL:" > $MAILOUT
echo "${BASEURL}${OUTPUTURL}${LNAME}.html " >> $MAILOUT
$CAT $MAILOUT | mail -s "Your GeoFOLD results are ready" $EMAIL 
###### CLEAN UP OLD FILES #######
#find $TMPDIR/ -mtime +240 -exec rm -vfr {} \;
find $JOBDIR/ -mtime +240 -exec rm -vfr {} \;
#find $HTMLDIR/ -mtime +240 -exec rm -vfr {} \;
