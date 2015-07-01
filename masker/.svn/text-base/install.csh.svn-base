#!/bin/csh
##### README for MASKER #####
#(0) Unpack masker
#
# tar -zxvf masker.tgz
# cd masker
#
#(1)
# install masker programs.
#
# Current settings:
# source setup.csh
##################################################################################
## set the number on the following line to: 128, 256, 512, 1024, 2048 or 4096   ##
setenv MASKSIZE  512
##################################################################################
if ($MASKSIZE > 2048) then
  setenv DTHETA 4.5
  @ SMALL = $MASKSIZE / 16
  @ NTHETA = 20
else
if ($MASKSIZE > 512) then
  setenv DTHETA 4.5
  @ SMALL = $MASKSIZE / 8
  @ NTHETA = 20
else
  setenv DTHETA 9.0
  @ SMALL = $MASKSIZE / 4
  ## @ SMALL = 7  ## for debugging xtriangle
  @ NTHETA = 10
endif
endif
#
# make clean
# make set$MASKSIZE
# make all
#
#(2)
# make mask library, and associated files.
#
# create new setup file
setenv CURRDIR `pwd`
sed -e "s#setenv MASKERDIR.*#setenv MASKERDIR "$CURRDIR"#" setup.csh | \
sed -e "s#setenv DESIGN_HOME.*#setenv DESIGN_HOME "$CURRDIR"#"           | \
sed -e "s/MSIZE   ..../MSIZE   "$MASKSIZE"/"           | \
sed -e "s/DTHETA  ..../DTHETA  "$DTHETA" /"           > masker_setup.csh
chmod +x masker_setup.csh
source masker_setup.csh
cat masker_setup.csh | awk '{printf "export %s=%s\n",$2,$3}' > masker_setup.sh
cp slopes4.5.dat.old $SLOPES
# create a template PDB file with evenly-spaced points on a sphere (i.e. 1024.pdb)
echo "Creating masks for masker: MASKSIZE=$MSIZE SMALL=$SMALL DTHETA=$DTHETA "
if !(-e $MSIZE.pdb) then
  @ NCYCLE = $MSIZE * 2
  echo "=============== running xmakemask ==============="
  echo "./xmakemask  $MSIZE $MSIZE.pdb 987654321 $NCYCLE"
  ./xmakemask  $MSIZE $MSIZE.pdb 987654321 $NCYCLE
  if ($status) exit
endif
# Sort those points using polar coordinates
  echo "=============== running xsortpdbmask ==============="
  echo "./xsortpdbmask $DTHETA $MSIZE.pdb > $MSIZE.sort"
./xsortpdbmask $DTHETA $MSIZE.pdb > $MSIZE.sort
if ($status) exit
# Cut the template mask all possible ways, creating a library of masks (1024.mask)
  echo "=============== running xbinarymask ==============="
  echo "./xbinarymask $MSIZE.sort $MSIZE.mask $MSIZE.dat $DTHETA"
./xbinarymask $MSIZE.sort $MSIZE.mask $MSIZE.dat $DTHETA
if ($status) exit
# create a template PDB file with evenly-spaces points on a sphere, a smaller set for tiling.
if !(-e $SMALL.pdb) then
  echo "=============== running xmakemask ==============="
  echo "./xmakemask $SMALL $SMALL.pdb 12345678 $MSIZE"
  ./xmakemask $SMALL $SMALL.pdb 12345678 $MSIZE
  if ($status) exit
endif
# Sort those points using polar coordinates
#./xsortpdbmask 10.0 $SMALL.pdb > $SMALL.sort
#if ($status) exit
# Assign points in the larger template ($MSIZE.pdb) to points in the smaller one ($SMALL.pdb)
  echo "=============== running xviewmask ==============="
  echo "./xviewmask $MSIZE.sort $SMALL.pdb $MSIZE.dat ${MSIZE}v.pdb $MSIZE.vmask"
./xviewmask $MSIZE.sort $SMALL.pdb $MSIZE.dat ${MSIZE}v.pdb $MSIZE.vmask
if ($status) exit
# Generate a file of all 2-byte integers and the number of 1-bits in each.
echo "99999999" | ./xcountbits
if ($status) exit
# Generate collars file, premasked edges.
  echo "=============== running xmakecollars ==============="
  echo "./xmakecollars"
./xmakecollars
if ($status) exit
# Generate slopes and intercepts to do linear interpolations
awk -v n=$NTHETA 'BEGIN{for(i=0;i<n;i++) print 0,0;}' > $SLOPES
  echo "=============== running xplottheta ==============="
  echo "./xplottheta 10. 5. | awk 'NF==8 {print $7, $8}' > junk.dat"
./xplottheta 10. 5. | awk 'NF==8 {print $7, $8}' > junk.dat
if ($status) exit
mv -fv junk.dat $SLOPES
# Make triangle file for tiling images  (BUG)
#  echo "=============== running xtriangles ==============="
#  echo "./xtriangles junk.r3d $VTRIANGLE"
#./xtriangles junk.r3d $VTRIANGLE
#if ($status) exit
# test program
source masker_setup.csh
  echo "=============== running xpdbmask ==============="
  echo "./xpdbmask test.pdb "
./xpdbmask test.pdb 
if ($status) then
  echo "ERROR encountered when running xpdbmask"
  exit
endif
#
echo  ">>>>> CHECK YOUR RESULTS: Here are the archived results for the test set:"
cat test.results
setenv BINDIR `echo $PATH | awk -F : '{print $(NF)}'`
echo  -n ">>>>> CHOOSE DIRECTORY FOR EXECUTABLES : [$BINDIR] "
setenv EXEDIR $<
if ($EXEDIR == "") then
  setenv EXEDIR $BINDIR
endif
cp -fv masker_setup.csh pdbmask voidmask contactmask $EXEDIR/
pushd $EXEDIR
ln -sfv voidmask VoidMask
ln -sfv voidmask Voidmask
ln -sfv pdbmask PdbMask
ln -sfv pdbmask PDBMask
ln -sfv contactmask unfold_energy
ln -sfv contactmask ContactMask
popd
echo  ">>>>> All done. MASKER programs are ready to go: pdbmask, voidmask, and contactmask"

