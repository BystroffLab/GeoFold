#!/bin/bash

export stf=/home/walcob/d/Dropbox/Documents/Bystroff/start2fold/$1
export tmp=/home/walcob/d/GeoFold/tmp/$2
export gDir=/home/walcob/d/GeoFold
export pdb=$2

cd $stf
wget https://files.rcsb.org/view/$pdb.pdb
cp -v $stf/$pdb.pdb $gDir/pdbs/
cp -v $gDir/parameters/template.par $gDir/parameters/$pdb.par
emacs $gDir/parameters/$pdb.par
emacs $stf/$pdb.hx

python $gDir/rungeofold.py $gDir/parameters/$pdb.par $gDir/walcob.conf
cp -v $tmp/*.dag.age $stf/
cp -v $tmp/*.cij $stf/
cp -v $tmp/*.seq $stf/
python $gDir/hx2cij.py -hx $stf/$pdb.hx -cij $stf/$pdb.cij -o $stf/$pdb.hx2cij
for i in {1..21}
do
  $gDir/hxpathway2ps $stf/$pdb.seq $stf/${pdb}_${i}.dag.age $stf/$pdb.cij $stf/${pdb}_${i}.hx.ps $stf/$pdb.hx2cij
  convert -trim -background white $stf/${pdb}_${i}.hx.ps $stf/${pdb}_${i}.hx.png
done
