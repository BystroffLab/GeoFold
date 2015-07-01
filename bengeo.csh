#### RUNGEOFOLD ####
# This script runs the GeoFOLD suite of programs
# The input files are (1) a parameters file and (2) a PDB file
# All you need to do is modify the generic "paramaters"
# file so that it has the basename for the PDB file (i.e. usually the
# 4-letter PDB code), the chain letter(s) or "_" for
# blank, and any modifications you might want.
# Name the new file appropriately and then run this script:
# ./bengeo.csh mypdb.pdb myfile.par 
# output files will appear in subdirectories output, tmp, log.
# If these directories do not already exist, they are created.
# The current directory must contain the geofold package,
# installed by unpacking geofold.tgz
# tar -zxvf geofold.tgz
# ...then compiled...
# cd geofold
# make clean
# make all
# ... that's it.
##
#--------------------------------------------------------
# By modifying the directories below, you can
# run this script and generate HTML output.
# C.Bystroff  29-NOV-2008
#B. Walcott 14-MAY-2014
##### COMMAND LINE PARAMETERS ######
setenv PARAMFILE $2
if !(-e $PARAMFILE) then
  echo "$PARAMFILE not found. "
  exit
endif
setenv PDBFILE $1
if !(-e $PDBFILE) then
  echo "$PDBFILE not found. "
  exit
endif
##### JOB level variable from parameter file ######
# PDBCODE is the basename of a PDB file to be found in TMPDIR
# CHAIN is a set of chain IDs within the PDB file, or _ for space.
# LNAME is a unique name for this job or the process ID
# ONAME is a unique name for a previous job, to be used to skip GEOFOLD.
# UNAME is a unique name for the current job
setenv PDBCODE `grep ^"PDBCODE" $PARAMFILE | awk '{print $2}'`
if ($PDBCODE == "") then
  echo "PDBCODE line not found."
  exit
else
  echo "PDBCODE $PDBCODE"
endif
setenv CHAIN `grep ^"CHAIN" $PARAMFILE | awk '{print $2}'`
if ($CHAIN == "") then
  echo "CHAIN line not found."
  exit
else
  echo "CHAIN $CHAIN"
endif
setenv LNAME `grep ^"LNAME" $PARAMFILE | awk '{print $2}'`
if ($LNAME == "") then
  echo "LNAME line not found."
  setenv LNAME $$
endif
echo "LNAME $LNAME"
setenv TARG $PDBCODE$CHAIN
setenv UNAME $TARG$LNAME
setenv ONAME `grep ^"ONAME" $PARAMFILE | awk '{print $2}'`
if ($ONAME == "") then
  echo "ONAME line not found."
  setenv ONAME $UNAME
endif
echo "ONAME $ONAME"
## setenv EMAIL `echo $LNAME | sed -e "s/000.*//"`
########################## DIRECTORIES #############################
### SET THESE DIRECTORIES AS FOLLOWS:	
# GDIR is the geofold directory created by tar -zxvf geofold.tgz
# TMPDIR is a temporary file directory. Everything in that
# directory can be deleted later. The PDB file should be copied
# to TMPDIR to start.
# LOGDIR is a place to hold the log file.
# HTMLDIR is a public_html directory where HTML
# file may be viewed.
# JOBDIR holds job files that are read by geod.csh
# (not used in this script except to clean up)
# THISDIR is the directory where you are running this script.
setenv THISDIR `pwd`
setenv GDIR $THISDIR
setenv LOGDIR $THISDIR/log
setenv HTMLDIR $THISDIR/output
setenv JOBDIR $THISDIR/jobs
setenv TMPDIR $THISDIR/tmp
##setenv BASEURL "http://www.bioinfo.rpi.edu/bystrc/geofold/"
setenv BASEURL "file:///home/bystroff3/public_html/"
setenv OUTPUTURL "output/"
########################### PROGRAMS ###############################
## setenv CONVERT /ext2/www/html/applications/mfold/bin32/convert
setenv CONVERT convert
setenv DOT dot
setenv MAXTRAFFIC $GDIR/maxTraffic
setenv MTCUT 0.1
setenv CREATEGNUPLOT $GDIR/creategnuplot.csh
##### FILES ######
# PARAMTEMPLATE is the template for the parameters file. Don't touch this.
# PARAMFILE is a copy of PARAMTEMPLATE, and may be modified.
# other files are programs, output files etc.
setenv PARAMTEMPLATE ${THISDIR}/parameters
setenv HTMLREFRESH ${THISDIR}/header_refresh.html 
setenv HTMLHEAD ${THISDIR}/header.html 
## setenv PARAMFILE ${TMPDIR}/${LNAME}.par
#setenv XGEOFOLD $GDIR/xgeofold_split4
setenv XGEOFOLD $GDIR/xgeofold
setenv HTMLTMP ${TMPDIR}/${UNAME}.html	
setenv HTMLLOG ${LOGDIR}/${UNAME}.log
setenv HTMLOUT ${HTMLDIR}/${UNAME}.html	
setenv HTMLPERM ${TMPDIR}/${UNAME}.perm
##### VARIABLES #####
# OMEGA is the desolvation free energy
# MSPL is the maximum number of splits for each intermediate
#   of unfolding.
# FING (folding) = 1 is folding, 0 is unfolding
# HLFE (halflife) = 1 for stopping at the half-life, 0 to go to equilibrium in UnfoldSim.
# 
setenv OMEGA 30
setenv MSPL 4
setenv FING 0
setenv HLFE 0

if(`grep -c ^"SIDECHAINHB" $PARAMFILE`) then
  setenv SIDECHAINHB `grep ^"SIDECHAINHB" $PARAMFILE | awk '{print $2}'`
else
  setenv SIDECHAINHB 1
endif
echo "SIDECHAINHB $SIDECHAINHB"
setenv FING `grep ^"FOLDING" $PARAMFILE | awk '{print $2}'`
echo "FOLDING $FING"
setenv HLFE `grep ^"HALFLIFE" $PARAMFILE | awk '{print $2}'`
echo "HALFLIFE $HLFE"
setenv BCUT `grep ^"BREAKCUT" $PARAMFILE | awk '{print $2}'`
echo "BREAKCUT $BCUT"
setenv PCUT `grep ^"PIVOTCUT" $PARAMFILE | awk '{print $2}'`
echo "PIVOTCUT $PCUT"
setenv HCUT `grep ^"HINGECUT" $PARAMFILE | awk '{print $2}'`
echo "HINGECUT $HCUT"
setenv WAT `grep ^"WATER" $PARAMFILE | awk '{print $2}'`
echo "WATER $WAT"
set ORANGE=(`grep ^"ORANGE" $PARAMFILE | cut -c 7- `)
if ($ORANGE == 0) then
  set ORANGE=(`grep ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
if (`grep -c ^"MAXSPLIT" $PARAMFILE`) then
  setenv MSPL `grep ^"MAXSPLIT" $PARAMFILE | awk '{print $2}'`
  if (-e $GDIR/xgeofold_split$MSPL) then
    setenv XGEOFOLD $GDIR/xgeofold_split$MSPL
  endif
endif
if (`grep -c ^"RUNGEOFOLD" $PARAMFILE`) then
  setenv DOIT `grep ^"RUNGEOFOLD" $PARAMFILE | awk '{print $2}'`
  if ( $DOIT ) then
    echo "Running GeoFOLD"
  else
    echo "Skipping GeoFOLD"
    cat $HTMLREFRESH > $HTMLTMP
    echo "<h3>Skipping GeoFOLD. UnfoldSim results for $UNAME </h3>" >> $HTMLTMP
    echo "<h4>Using previously calculated DAG file $ONAME </h4>" >> $HTMLTMP
    echo "<pre>" >> $HTMLTMP
    cat $HTMLTMP > $HTMLOUT
    echo '</body></html>' >> $HTMLOUT
    if ( $ONAME == $UNAME ) then
      echo "ERROR Need old job name as arg 4 to skip geofold."
      exit
    endif
    cp ${TMPDIR}/$ONAME.dag ${TMPDIR}/$UNAME.dag
    goto SKIPGEOFOLD
  endif
endif
cat $HTMLREFRESH > $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT

echo "<h3>GeoFOLD results for $UNAME </h3>" >> $HTMLTMP
echo "<pre>" >> $HTMLTMP
echo "Job $UNAME" >> $HTMLTMP
echo "Unfolding $TARG" >> $HTMLTMP
echo "============= PARAMETERS =============" 
echo "============= PARAMETERS =============" >> $HTMLTMP
if !(-e $PARAMFILE) then
  cp $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file not found!!" >> $HTMLTMP
else if (-z $PARAMFILE) then
  cp $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file is empty!!" >> $HTMLTMP
endif
echo "Using these parameters" >> $HTMLTMP
cat $PARAMFILE >> $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= GETCHAIN =============" 
echo "============= GETCHAIN =============" >> $HTMLTMP
cp ./$PDBFILE ${TMPDIR}/$PDBCODE.pdb
if !(-e ${TMPDIR}/$PDBCODE.pdb ) then
  echo "File not found: ${TMPDIR}/$PDBCODE.pdb "
  exit
endif
echo "Extracting protein atoms from $PDBCODE " >> $HTMLTMP
(	$GDIR/xgetchain + < ${TMPDIR}/$PDBCODE.pdb \
        > ${TMPDIR}/$UNAME.tmp	|| \
        ( echo "Error in GETCHAIN" ; exit )  ) >>& $HTMLTMP
echo "Extracting chains $CHAIN " >> $HTMLTMP
(	$GDIR/xgetchain $CHAIN < ${TMPDIR}/$UNAME.tmp \
        > ${TMPDIR}/$UNAME.pdb	|| \
        ( echo "Error in GETCHAIN" ; exit )  ) >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= RENUMBER =============" 
echo "============= RENUMBER =============" >> $HTMLTMP
if !(-e ${TMPDIR}/$UNAME.pdb ) then
  echo "File not found: ${TMPDIR}/$UNAME.pdb "
  exit
endif
(  $GDIR/xrenumber_one ${TMPDIR}/$UNAME.pdb ${TMPDIR}/$UNAME.tmp || \
   ( echo "error in renumber "; exit ) ) >>& $HTMLTMP
  mv ${TMPDIR}/$UNAME.tmp ${TMPDIR}/$UNAME.pdb   
  cp ${TMPDIR}/$UNAME.pdb  $HTMLDIR/$UNAME.pdb
  if (-z ${TMPDIR}/$UNAME.pdb ) then
    echo "ERROR. empty file. " >>  $HTMLTMP
    cat $HTMLTMP > $HTMLOUT
    echo '</PRE></body></html>' >> $HTMLOUT
    exit
  endif
  echo 'Using these <A HREF="'$UNAME.pdb'">coordinates.</A>' >> $HTMLTMP
  cat $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT

echo "============= 3to1 (extract sequence) =============" 
echo "============= 3to1 (extract sequence) =============" >> $HTMLTMP
(	$GDIR/x3to1 "+" < ${TMPDIR}/$UNAME.pdb > ${TMPDIR}/$UNAME.seq || \
        ( echo "3to1 ended with errors" ; exit )    ) >>& $HTMLTMP
        echo "Complete sequence (all chains ) is"  >>& $HTMLTMP
        cat ${TMPDIR}/$UNAME.seq  >>& $HTMLTMP
  cat $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT

echo "============= PDB2CIJ (extract contacts) =============" 
echo "============= PDB2CIJ (extract contacts) =============" >> $HTMLTMP
(	$GDIR/xpdb2cij ${TMPDIR}/$UNAME.pdb 8. \
        > ${TMPDIR}/$UNAME.cij || \
        ( echo "Error in PDB2CIJ " ; exit )    ) >>& $HTMLTMP
 echo -n "Number of contacts found: " >>& $HTMLTMP
 wc -l ${TMPDIR}/$UNAME.cij | awk '{print $1}' >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= PDB2HB (extract H-bonds SS-bonds) =============" 
echo "============= PDB2HB (extract H-bonds SS-bonds) =============" >> $HTMLTMP
(	$GDIR/xpdb2hb $PARAMFILE ${TMPDIR}/$UNAME.pdb \
        > ${TMPDIR}/$UNAME.hb || \
        ( echo "Error in PDB2HB " ; exit )    ) >>& $HTMLTMP
 echo -n "Number of H-bonds found: " >>& $HTMLTMP
 grep -c " H " ${TMPDIR}/$UNAME.hb >>& $HTMLTMP
 echo -n "Number of SS-bonds found: " >>& $HTMLTMP
 grep -c " S " ${TMPDIR}/$UNAME.hb >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
echo "HBONDS ${TMPDIR}/$UNAME.hb" >> $PARAMFILE

echo "============= CONTACTMASK =============" 
echo "============= CONTACTMASK =============" >> $HTMLTMP
(	echo "Time before running CONTACTMASK `date`"                                                               ) >>& $HTMLTMP
(	source $GDIR/masker/masker_setup.csh; \
        $GDIR/masker/xcontactmask ${TMPDIR}/$UNAME.pdb \
        ${TMPDIR}/$UNAME.sas 1.4 || \
        ( echo "Error in CONTACTMASK " ; exit )   ) >>& $HTMLTMP
(	echo "Time after running CONTACTMASK `date`" ) >>& $HTMLTMP
  cp ${TMPDIR}/$UNAME.sas $HTMLDIR/$UNAME.sas
  echo 'Pairwise contact surfaces <A HREF="'$UNAME.sas'">file.</A>' >> $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
echo "CONTACTS ${TMPDIR}/$UNAME.sas" >> $PARAMFILE

echo "============= VOIDMASK =============" 
echo "============= VOIDMASK =============" >> $HTMLTMP
(	echo "Time before running VOIDMASK `date`") >>& $HTMLTMP
(	source $GDIR/masker/masker_setup.csh; \
        $GDIR/masker/xvoidmask ${TMPDIR}/$UNAME.pdb \
        ${TMPDIR}/$UNAME.void 1.4 1.2 1.4 || \
        ( echo "Error in VOIDMASK " ; exit )   ) >>& $HTMLTMP
(	echo "Time after running VOIDMASK `date`" ) >>& $HTMLTMP
cp ${TMPDIR}/$UNAME.void  ${HTMLDIR}/$UNAME.void.pdb
if (-z ${TMPDIR}/$UNAME.void ) then
      (	echo "ERROR in VOIDMASK . Empty file. " ) >>& $HTMLTMP
  cat $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  exit
endif
echo 'Coordinates with void positions: <A HREF="'$UNAME.void.pdb'">'$UNAME.void.pdb'</A>' \
      >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

echo "============= GEOFOLD =============" 
echo "============= GEOFOLD =============" >> $HTMLTMP
	echo "Time before running GEOFOLD `date`"  >>& $HTMLTMP
	echo "$XGEOFOLD  ${TMPDIR}/$UNAME.void ${TMPDIR}/$UNAME.dag $PARAMFILE "  >>& $HTMLTMP
	( $XGEOFOLD  ${TMPDIR}/$UNAME.void ${TMPDIR}/$UNAME.dag $PARAMFILE  || \
        ( echo "Error in GEOFOLD " ; exit )   ) >>& $HTMLTMP
	echo "Time after running GEOFOLD `date`" >>& $HTMLTMP
if (-z ${TMPDIR}/$UNAME.dag ) then
      (	echo "ERROR in GEOFOLD . Empty file. " ) >>& $HTMLTMP
   cat $HTMLTMP > $HTMLOUT
   echo '</PRE></body></html>' >> $HTMLOUT
   exit
endif
  cp ${TMPDIR}/$UNAME.dag  ${HTMLDIR}/$UNAME.dag
  cp ${TMPDIR}/$UNAME.dag  ${TMPDIR}/$ONAME.dag 
  echo 'Directed acyclic graph (DAG) <A HREF="'$UNAME.dag'">file.</A>' >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

SKIPGEOFOLD:
echo "============= UNFOLDSIM =============" 
echo "Unfolding $TARG" >> $HTMLTMP
echo "============= UNFOLDSIM =============" >> $HTMLTMP
@ NN = 0
while ( $#ORANGE )
  sed -e "s/^OMEGA .*/OMEGA $ORANGE[1]/" $PARAMFILE > ${PARAMFILE}.1
  @ NN ++
  setenv LOGFILE $TMPDIR/${UNAME}_$NN.log
  echo "============= run $NN omega=$ORANGE[1] =============" 
  echo "============= run $NN omega=$ORANGE[1] =============" >> $HTMLTMP
  echo "Time before running UNFOLDSIM $NN `date`"    >>& $HTMLTMP
  echo "$GDIR/xunfoldsim ${TMPDIR}/$ONAME.dag ${PARAMFILE}.1 > $LOGFILE " >> $HTMLTMP
  $GDIR/xunfoldsim ${TMPDIR}/$ONAME.dag ${PARAMFILE}.1 > $LOGFILE || \
    ( echo "ERROR when running UNFOLDSIM " ; exit )
  echo "Time after running UNFOLDSIM `date`"         >>& $HTMLTMP
  echo "<p><pre>"  >>& $HTMLTMP
  grep ^"TIMECOURSE" $LOGFILE | tail -1 >>& $HTMLTMP
  cat $HTMLTMP > $HTMLOUT
  echo '</PRE></body></html>' >> $HTMLOUT
  shift ORANGE
end
@ MM = $NN

echo "============= PATHWAY2PS =============" 
echo "============= PATHWAY2PS =============" >> $HTMLTMP
echo "$GDIR/xpathway2ps ${TMPDIR}/${ONAME}.seq ${TMPDIR}/${ONAME}.dag.path ${TMPDIR}/${ONAME}.cij ${TMPDIR}/${UNAME}.ps 4"	 >> $HTMLTMP
$GDIR/xpathway2ps ${TMPDIR}/${ONAME}.seq ${TMPDIR}/${ONAME}.dag.path ${TMPDIR}/${ONAME}.cij ${TMPDIR}/${UNAME}.ps 4	 >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT
echo "============= CONVERT =============" 
echo "============= CONVERT =============" >> $HTMLTMP
echo "$CONVERT -trim -geometry 100 ${TMPDIR}/$UNAME.ps  ${TMPDIR}/${UNAME}_thumb.gif"  >> $HTMLTMP
$CONVERT -trim -geometry 100 ${TMPDIR}/$UNAME.ps  ${TMPDIR}/${UNAME}_thumb.gif
$CONVERT ${TMPDIR}/${UNAME}_thumb.gif ${HTMLDIR}/${UNAME}_thumb.jpg
$CONVERT -trim -geometry 600 ${TMPDIR}/${UNAME}.ps  ${TMPDIR}/${UNAME}.gif
$CONVERT ${TMPDIR}/${UNAME}.gif  ${HTMLDIR}/${UNAME}.jpg
echo 'ageplots made'
##### DOT must be installed for the following to work
$DOT >& /dev/null
if ( $status ) then
  echo "Dot must be installed in order to generate the DAG image."
else
  echo "============= MAXTRAFFIC =============" 
  echo "============= MAXTRAFFIC ============="  >> $HTMLTMP
  echo "$MAXTRAFFIC ${TMPDIR}/${ONAME}.dag.out $MTCUT > ${TMPDIR}/$UNAME.dot"  >> $HTMLTMP
  $MAXTRAFFIC ${TMPDIR}/${ONAME}.dag.out $MTCUT > ${TMPDIR}/$UNAME.dot
  echo "============= DOT =============" 
  echo "============= DOT ============="  >> $HTMLTMP
  echo "$DOT -Tjpg  -o ${HTMLDIR}/${UNAME}.dag.jpg -Timap -o ${HTMLDIR}/${UNAME}.dag.map ${TMPDIR}/${UNAME}.dot"  >> $HTMLTMP
  $DOT -Tjpg  -o ${HTMLDIR}/${UNAME}.dag.jpg -Timap -o ${HTMLDIR}/${UNAME}.dag.map ${TMPDIR}/${UNAME}.dot
  echo "$DOT -Tps  -o ${TMPDIR}/${UNAME}.dot.ps  ${TMPDIR}/${UNAME}.dot   "  >> $HTMLTMP
  $DOT -Tps  -o ${TMPDIR}/${UNAME}.dot.ps  ${TMPDIR}/${UNAME}.dot   
  echo "============= CONVERT =============" 
  echo "============= CONVERT ============="  >> $HTMLTMP
  echo "$CONVERT -trim -geometry 100 ${TMPDIR}/${UNAME}.dot.ps ${TMPDIR}/${UNAME}.dag_thumb.gif" >> $HTMLTMP
  $CONVERT -trim -geometry 100 ${TMPDIR}/${UNAME}.dot.ps ${TMPDIR}/${UNAME}.dag_thumb.gif
  echo "$CONVERT ${TMPDIR}/${UNAME}.dag_thumb.gif ${HTMLDIR}/${UNAME}.dag_thumb.jpg" >> $HTMLTMP
  $CONVERT ${TMPDIR}/${UNAME}.dag_thumb.gif ${HTMLDIR}/${UNAME}.dag_thumb.jpg
endif
#################### CREATE FINAL OUTPUT HTML FILE ##########################
## DAG image map file
echo '<html><body><A HREF="'${UNAME}.dag.map'" border="0"><img src="'${UNAME}.dag.jpg'" ismap></A></body></html>' >  ${HTMLDIR}/${UNAME}.dag.html
## Main output page
cat $HTMLHEAD > $HTMLPERM
cp ${TMPDIR}/$ONAME.pdb  $HTMLDIR/
cp ${TMPDIR}/$ONAME.seq  $HTMLDIR/
echo '<h3>Unfolding Age plot for '$TARG' </h3><A HREF="'${UNAME}.jpg'"><img src="'${UNAME}_thumb.jpg'" alt="'${UNAME}_thumb'"></A><br>' >> $HTMLPERM
echo '<h3>Unfolding Age plot for '$TARG' </h3><A HREF="'${UNAME}.jpg'"><img src="'${UNAME}_thumb.jpg'" alt="'${UNAME}_thumb'"></A><br>' >> $HTMLTMP
echo '<h3>Unfolding Pathway plot for '$TARG' </h3><A HREF="'${UNAME}.dag.html'"><img src="'${UNAME}.dag_thumb.jpg'" alt="'${UNAME}.dag'"></A><br>' >> $HTMLPERM
echo '<h3>Unfolding Pathway plot for '$TARG' </h3><A HREF="'${UNAME}.dag.html'"><img src="'${UNAME}.dag_thumb.jpg'" alt="'${UNAME}.dag'"></A><br>' >> $HTMLTMP
echo '<h3><A HREF="'${ONAME}.dag'">Unfolding graph (text file) for '$TARG' </A></h3>' >> $HTMLPERM
echo '<h3><A HREF="'${ONAME}.dag'">Unfolding graph (text file) for '$TARG' </A></h3>' >> $HTMLTMP
echo '<h3><A HREF="'${ONAME}.pdb'">Coordinate file ('$TARG') used for this simulation</A></h3>' >> $HTMLPERM
echo '<h3><A HREF="'${ONAME}.pdb'">Coordinate file ('$TARG') used for this simulation</A></h3>' >> $HTMLTMP
echo '<h3>Sequence (all chains)</h3>' >> $HTMLPERM
echo '<h3>Sequence (all chains)</h3>' >> $HTMLTMP
echo '<pre>' >> $HTMLPERM
echo '<pre>' >> $HTMLTMP
cat ${TMPDIR}/${ONAME}.seq >> $HTMLPERM
cat ${TMPDIR}/${ONAME}.seq >> $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</PRE></body></html>' >> $HTMLOUT

#### Output unfolding rates at different omega values ####
echo '</pre>' >> $HTMLPERM
echo '</pre>' >> $HTMLTMP
echo -n '<h4>' >>& $HTMLPERM
echo -n '<h4>' >>& $HTMLTMP
if ( $HLFE ) then
  if ( $FING ) then
    echo 'Times and concentrations at 1/2 folded' >>& $HTMLPERM
    echo 'Times and concentrations at 1/2 folded' >>& $HTMLTMP
  else
    echo 'Times and concentrations at 1/2 unfolded' >>& $HTMLPERM
    echo 'Times and concentrations at 1/2 unfolded' >>& $HTMLTMP
  endif
else
  if ( $FING ) then
    echo 'Equilibrium concentrations and rates under folding conditions' >>& $HTMLPERM
    echo 'Equilibrium concentrations and rates under folding conditions' >>& $HTMLTMP
  else
    echo 'Equilibrium concentrations and rates under unfolding conditions' >>& $HTMLPERM
    echo 'Equilibrium concentrations and rates under unfolding conditions' >>& $HTMLTMP
  endif
endif
echo "<p><pre>"  >>& $HTMLPERM
echo "<p><pre>"  >>& $HTMLTMP
echo -n "dGsolvation   "  >>& $HTMLPERM
echo -n "dGsolvation   "  >>& $HTMLTMP
echo "   time(sec)   final[F]  final[U]  final[I]  sum[c](ppb) half_life(sec)" >>& $HTMLPERM
echo "   time(sec)   final[F]  final[U]  final[I]  sum[c](ppb) half_life(sec)" >>& $HTMLTMP
@ NN = 0
set ORANGE=(`grep ^"ORANGE" $PARAMFILE | cut -c 7- `)
if !( $#ORANGE ) then
  set ORANGE=(`grep ^"OMEGA " $PARAMFILE | awk '{print $2}'`)
endif
echo -n "" > $TMPDIR/${UNAME}.plot
while ( $#ORANGE )
  @ NN ++
  setenv LOGFILE $TMPDIR/${UNAME}_$NN.log
  echo -n "$ORANGE[1] kJ/mol/A^2   "  >>& $HTMLPERM
  echo -n "$ORANGE[1] kJ/mol/A^2   "  >>& $HTMLTMP
  grep ^"TIMECOURSE" $LOGFILE | tail -1 | sed -e "s/TIMECOURSE//" >>& $HTMLPERM
  grep ^"TIMECOURSE" $LOGFILE | tail -1 | sed -e "s/TIMECOURSE//" >>& $HTMLTMP
  grep ^"TIMECOURSE" $LOGFILE | awk '$4<50.' | tail -1 \
    | awk -v omega=$ORANGE[1] '{print omega, log(log(2)/$2)}' >> $TMPDIR/${UNAME}.plot
  shift ORANGE
end 
echo "</pre>"  >>& $HTMLPERM
echo "</pre>"  >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT
setenv FIT $TMPDIR/${UNAME}.fit
echo "============= FIT_POLY =============" 
echo "<pre>"  >>& $HTMLTMP
$GDIR/xfit_poly $TMPDIR/${UNAME}.plot 1 $FIT  >>& $HTMLTMP
setenv IC `grep -A 2 "Least-squares" $FIT | tail -1 | awk '{print $2}'`
setenv SL `grep -A 3 "Least-squares" $FIT | tail -1 | awk '{print $2}'`
setenv LNKU `awk -v OM=$WAT -v sl=$SL -v ic=$IC 'BEGIN{print sl*OM + ic}'`
echo "</pre>"  >>& $HTMLTMP
cat $HTMLTMP > $HTMLOUT
echo '</body></html>' >> $HTMLOUT

echo -n "<h5>Projected half-life of unfolding " >>& $HTMLPERM
echo -n "<h5>Projected half-life of unfolding " >>& $HTMLTMP
echo "in pure water (dGsolvation=$WAT kJ/mol/A^2) : " >>& $HTMLPERM
echo "in pure water (dGsolvation=$WAT kJ/mol/A^2) : " >>& $HTMLTMP
awk -v lnku=$LNKU 'BEGIN{printf "%e", log(2)/exp(lnku)}' >>& $HTMLPERM
awk -v lnku=$LNKU 'BEGIN{printf "%e", log(2)/exp(lnku)}' >>& $HTMLTMP
echo " seconds</h5><br> " >>& $HTMLPERM
echo " seconds</h5><br> " >>& $HTMLTMP
##### Create gnuplot image
$CREATEGNUPLOT $TMPDIR/${UNAME}  $WAT
$CONVERT -trim $TMPDIR/${UNAME}.gp.ps $TMPDIR/${UNAME}.gp.gif
$CONVERT $TMPDIR/${UNAME}.gp.gif $HTMLDIR/${UNAME}.gp.jpg
$CONVERT -trim -geometry 100 $TMPDIR/${UNAME}.gp.ps $TMPDIR/${UNAME}.gp_thumb.gif
$CONVERT $TMPDIR/${UNAME}.gp_thumb.gif $HTMLDIR/${UNAME}.gp_thumb.jpg
echo '<h4>ln(ku) vs omega for '$TARG' </h4><A HREF="'${UNAME}.gp.jpg'"><img src="'${UNAME}.gp_thumb.jpg'" alt="'${UNAME}.gp.jpg'"></A><br>' >> $HTMLTMP
echo '<h4>ln(ku) vs omega for '$TARG' </h4><A HREF="'${UNAME}.gp.jpg'"><img src="'${UNAME}.gp_thumb.jpg'" alt="'${UNAME}.gp.jpg'"></A><br>' >> $HTMLPERM
echo "<br>Completed `date` <br>" >> $HTMLPERM
echo "<br>Completed `date` <br>" >> $HTMLTMP
###### Create "do over" button.
# echo '<FORM METHOD="POST" ACTION="../geofold.cgi" >' >> $HTMLTMP
# echo '<INPUT type="hidden" name="email_address" value="'$UNAME'@rip.ude">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="pdbid" value="'$PDBCODE'">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="chains" value="'$CHAIN'">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="params" value="'$PARAMFILE'">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="slash" value="/">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="space" value=" ">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="plus" value="+">' >> $HTMLTMP
# echo '<INPUT type="hidden" name="lname" value="'$ONAME'">' >> $HTMLTMP
# echo '<p><INPUT TYPE="submit" NAME=".submit" VALUE="Do over">' >> $HTMLTMP
# echo '</FORM>' >> $HTMLTMP
# cat $HTMLTMP > $HTMLOUT
# echo '</body></html>' >> $HTMLOUT
# echo '<FORM METHOD="POST" ACTION="../geofold.cgi" >' >> $HTMLPERM
# echo '<INPUT type="hidden" name="email_address" value="'$UNAME'@rip.ude">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="pdbid" value="'$PDBCODE'">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="chains" value="'$CHAIN'">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="params" value="'$PARAMFILE'">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="slash" value="/">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="space" value=" ">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="plus" value="+">' >> $HTMLPERM
# echo '<INPUT type="hidden" name="lname" value="'$ONAME'">' >> $HTMLPERM
# echo '<p><INPUT TYPE="submit" NAME=".submit" VALUE="Do over">' >> $HTMLPERM
# echo 'based on <A HREF="'${BASEURL}${OUTPUTURL}${ONAME}.html'">'$ONAME'</A>' >> $HTMLPERM
# echo '</FORM>' >> $HTMLPERM
##### BACK button
# echo '<h3><A HREF="../server.php">Back to GeoFold server</A></h3><br>' >> $HTMLPERM
echo '</body></html>' >> $HTMLPERM
cp $HTMLPERM $HTMLOUT

###### CLEAN UP OLD FILES #######
# find $TMPDIR/ -mtime +60 -exec rm -vfr {} \;
# find $JOBDIR/ -mtime +60 -exec rm -vfr {} \;
# find $HTMLDIR/ -mtime +120 -exec rm -vfr {} \;
