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
##### Getting info from parameter file #####
# PDBCODE is PDB ID of PDBFILE
# CHAIN is a set of chain IDs within the PDB file, or _ for space.
# LNAME is a unique name for this job or the process ID
# ONAME is a unique name for a previous job, to be used to skip GEOFOLD.
# UNAME is a unique name for the current job
#
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
##### Directories and Programs #####
setenv thisdir `pwd`
setenv seams $thisdir/seams
setenv xgeofold $thisdir/xgeofold
mkdir $thisdir/out/$UNAME
setenv outdir $thisdir/out/$UNAME
setenv convert /usr/bin/convert

##### VARIABLES #####
# OMEGA is the desolvation free energy
# MSPL is the maximum number of splits for each intermediate
#   of unfolding.
# FING (folding) = 1 is folding, 0 is unfolding
# HLFE (halflife) = 1 for stopping at the half-life, 0 to go to equilibrium in UnfoldSim.
# SCHB = 1 if sidechain hydrogen bonds are to be considered, otherwise 0
# 
setenv OMEGA 30
setenv MSPL 4
setenv FING 0
setenv HLFE 0
if(`grep -c ^"SIDECHAINHB" $PARAMFILE`) then
  setenv SCHB `grep ^"SIDECHAINHB" $PARAMFILE | awk '{print $2}'`
else
  setenv SCHB 1
endif
echo "SIDECHAINHB $SCHB"
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
  if (-e $thisdir/xgeofold_split$MSPL) then
    setenv xgeofold $thisdir/xgeofold_split$MSPL
  endif
endif
if (`grep -c ^"RUNGEOFOLD" $PARAMFILE`) then
  setenv DOIT `grep ^"RUNGEOFOLD" $PARAMFILE | awk '{print $2}'`
  if ( $DOIT ) then
    echo "Running GeoFOLD"
  else
    echo "Skipping GeoFOLD"
    if ( $ONAME == $UNAME ) then
      echo "ERROR Need old job name as arg 4 to skip geofold."
      exit
    endif
    cp ${thisdir}/$ONAME.dag ${thisdir}/$UNAME.dag
    goto SKIPGEOFOLD
  endif
endif
echo "============= PARAMETERS =============" 
if !(-e $PARAMFILE) then
  cp $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file not found!!"
else if (-z $PARAMFILE) then
  cp $PARAMTEMPLATE $PARAMFILE
  echo "Parameters file is empty!!"
endif
echo "Using these parameters"
echo "============= GETCHAIN ============="
if !(-e ${thisdir}/$PDBFILE ) then
  echo "File not found: ${thisdir}/$PDBFILE"
  exit
endif
cp ${thisdir}/$PDBFILE ${outdir}/$PDBCODE.pdb
(	$thisdir/xgetchain + < ${outdir}/$PDBCODE.pdb \
        > ${outdir}/$UNAME.tmp	|| \
        ( echo "Error in GETCHAIN" ; exit )  )
echo "Extracting protein chain atoms from $PDBFILE "
(	$thisdir/xgetchain $CHAIN < ${outdir}/$UNAME.tmp \
        > ${outdir}/$UNAME.pdb	|| \
        ( echo "Error in GETCHAIN" ; exit )  )
echo "============= RENUMBER ============="
if !(-e ${outdir}/$UNAME.pdb ) then
  echo "File not found: ${outdir}/$UNAME.pdb "
  exit
endif
(  $thisdir/xrenumber_one ${outdir}/$UNAME.pdb ${outdir}/$UNAME.tmp || \
   ( echo "error in renumber "; exit ) )
mv ${outdir}/$UNAME.tmp ${outdir}/$UNAME.pdb
if (-z ${outdir}/$UNAME.pdb ) then
  echo "ERROR. empty file. "
  exit
endif
echo "============= 3to1 (extract sequence) =============" 
(	$thisdir/x3to1 "+" < ${outdir}/$UNAME.pdb > ${outdir}/$UNAME.seq || \
        ( echo "3to1 ended with errors" ; exit )    )
echo "============= PDB2HB (extract H-bonds SS-bonds) =============" 
(	$thisdir/xpdb2hb $PARAMFILE ${outdir}/$UNAME.pdb \
        > ${outdir}/$UNAME.hb || \
        ( echo "Error in PDB2HB " ; exit )    )
echo "HBONDS ${outdir}/$UNAME.hb" >> $PARAMFILE
echo "============= PDB2SEAMS2 (extract seams) ============="
(	$seams/xpdb2seams2 ${outdir}/$UNAME.hb \
        > ${outdir}/$UNAME.seams || \
        ( echo "Error in PDB2SEAMS2 " ; exit )    )
echo "SEAMS ${outdir}/$UNAME.seams" >> $PARAMFILE
echo "============= PDB2CIJ (extract contacts) =============" 
(	$thisdir/xpdb2cij ${outdir}/$UNAME.pdb 8. \
        > ${outdir}/$UNAME.cij || \
        ( echo "Error in PDB2CIJ " ; exit )    )
echo "============= CONTACTMASK =============" 
(	source $thisdir/masker/masker_setup.csh; \
        $thisdir/masker/xcontactmask ${outdir}/$UNAME.pdb \
        ${outdir}/$UNAME.sas 1.4 || \
        ( echo "Error in CONTACTMASK " ; exit )   )
echo "CONTACTS ${outdir}/$UNAME.sas" >> $PARAMFILE
echo "============= VOIDMASK =============" 
(	source $thisdir/masker/masker_setup.csh; \
        $thisdir/masker/xvoidmask ${outdir}/$UNAME.pdb \
        ${outdir}/$UNAME.void 1.4 1.2 1.4 || \
        ( echo "Error in VOIDMASK " ; exit )   ) 
if (-z ${outdir}/$UNAME.void ) then
      (	echo "ERROR in VOIDMASK . Empty file. " )
  exit
endif
echo "============= GEOFOLD =============" 
( $xgeofold  ${outdir}/$UNAME.void ${outdir}/$UNAME.dag $PARAMFILE  || \
        ( echo "Error in GEOFOLD " ; exit )   )
if (-z ${outdir}/$UNAME.dag ) then
      (	echo "ERROR in GEOFOLD . Empty file. " )
  exit
endif

SKIPGEOFOLD:
echo "============= UNFOLDSIM =============" 
@ NN = 0
while ( $#ORANGE )
  sed -e "s/^OMEGA .*/OMEGA $ORANGE[1]/" $PARAMFILE > ${PARAMFILE}.1
  @ NN ++
  setenv LOGFILE $outdir/${UNAME}_$NN.log
  echo "============= run $NN omega=$ORANGE[1] =============" 
  $thisdir/xunfoldsim ${outdir}/$ONAME.dag ${PARAMFILE}.1 > $LOGFILE || \
    ( echo "ERROR when running UNFOLDSIM " ; exit )
  shift ORANGE
end
@ MM = $NN
echo "============= PATHWAY2PS =============" 
$thisdir/xpathway2ps ${outdir}/${ONAME}.seq ${outdir}/${ONAME}.dag.path ${outdir}/${ONAME}.cij ${outdir}/${UNAME}.ps 4
echo "============= $convert ============="
#$convert -trim -geometry 100 ${outdir}/$UNAME.ps  ${outdir}/${UNAME}_thumb.gif
#$convert ${outdir}/${UNAME}_thumb.gif ${outdir}/${UNAME}_thumb.jpg
$convert  -geometry 600x600 ${outdir}/${UNAME}.ps  ${outdir}/${UNAME}.jpg
