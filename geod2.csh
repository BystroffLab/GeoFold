#!/bin/csh  
#### GeoFold deamon
echo "geod2.csh deamon started `date` on `hostname` in `pwd`"
echo "This deamon drives the python script"
echo "Running on `hostname` `date` in `pwd`" > Running
## /home/bystrc/server/geofold
## NOTE: server/geofold is located on brahms
############ change SHOME to the directory geofold is installed in ##########
setenv SHOME /bach1/home/walcob
setenv LHOME /home/walcob
setenv GDIR $SHOME/GeoFold
setenv BDIR $GDIR
setenv PYTHON /usr/bin/python
echo "Jobs currently running on `hostname` 0" > Load
echo "Running in $GDIR" >> Running
setenv PDBDIR $GDIR/pdbs
if !(-e $PDBDIR ) setenv PDBDIR $LHOME/GeoFold/pdbs
setenv LOGDIR $GDIR/log
setenv BIOLUNIT $PDBDIR
if !(-e $BIOLUNIT ) setenv BIOLUNIT $PDBDIR
setenv GETCHAIN $GDIR/xgetchain
setenv STMP $GDIR/tmp
setenv SJOB $GDIR/jobs
setenv WGETPDB 'ssh bach1 cd server/geofold; /bach1/home/bystrc/bin/wgetpdb'
foreach RQRD ( $GDIR $BDIR $PDBDIR $LOGDIR $BIOLUNIT $GETCHAIN $SHOME $STMP $SJOB )
  if !(-e $RQRD ) then
    echo "ERROR: $RQRD is required for geod2 to run"
    echo "ERROR: $RQRD is required for geod2 to run" >> Running
    echo "Either fix the error or write an exception into geod2.csh"
    echo "Either fix the error or write an exception into geod2.csh" >> Running
    echo "Server exiting................." >> Running
    exit
  endif
end
setenv SKILL kill
@ ALLOWABLE = 2
@ NDONE = 0
if (-e Done) then
  @ NDONE = `cat Done`
endif
echo "JOBS ACCUMULATED SO FAR: $NDONE "

@ LHR = 0
while ( 1 )
  if (-e $GDIR/STOP) exit
  while (-e $GDIR/STOPPED)
    echo "geod.csh queue paused `date`"
    while (-e $GDIR/STOPPED)
     sleep 1
    end
    echo "geod.csh queue unpaused `date`"
  end
  ## COUNT GEOFOLD PROCESSES CURRENTLY RUNNING
  set N = `ps ax -o user -o pid -o command | awk -v w=walcob '$1 == w' | grep "${STMP}/.*csh" | grep -v "grep" | wc -l | awk '{print $1}'`
  if ( $N < $ALLOWABLE ) then
    ## if we can start another job, then do the following:
    ##  1. submit a geofold job using the basename of the file as the basename of the output,
    ##  and argument 1 of the content of the file as the "PDB code".
    foreach CFILE (`/bin/ls -1tra ${SJOB}/*.job | & grep -v "No match"`)
      echo "GEOFOLD JOB ========> $CFILE  `date`"
      if ( -e ${CFILE} ) then
        if ( -z ${CFILE} ) then
          echo "Zero length input file. ${CFILE} ......removing it." 
          rm -vf ${CFILE}
        else
          ## The .job file exists and non-zero .  Set up the job.
          setenv PDBCODE `awk '{print $1}' ${CFILE} | head -1`
          setenv CHAIN `awk '{print $2}' ${CFILE} | head -1`
          setenv BUNIT `awk '{print $3}' ${CFILE} | head -1`
          if (`awk '{print NF}'  ${CFILE} | head -1` > 3 ) then
            setenv ONAME `awk '{print $4}' ${CFILE} | head -1`
            ## NOTE: Do-over. PDB file should already exist in ${STMP}/
          else
            setenv ONAME ""
            ## NOTE: PDB file should already exist in pdb/ or pdb1/.
            ## If it does not, use wget to get it.
            if ( $BUNIT ) then
              if ( -e ${STMP}/$PDBCODE.pdb ) then
              else if ( -e $BIOLUNIT/$PDBCODE.pdb ) then
                cp $BIOLUNIT/$PDBCODE.pdb  ${STMP}/
              else
                $WGETPDB $PDBCODE 1 > wgetlog
                if (`grep -c "Connection refused" wgetlog` > 0) then
                   echo "ERROR: cannot communicate with PDB ftp site" > ${STMP}/$PDBCODE.error
                else
                  mv $PDBCODE.pdb1 $BIOLUNIT/$PDBCODE.pdb
                  cp $BIOLUNIT/$PDBCODE.pdb  ${STMP}/
                endif
              endif
            else
              if ( -e ${STMP}/$PDBCODE.pdb ) then
              else if ( -e $PDBDIR/$PDBCODE.pdb ) then
                cp $PDBDIR/$PDBCODE.pdb  ${STMP}/
              else
                $WGETPDB $PDBCODE > wgetlog
                if (`grep -c "Connection refused" wgetlog` > 0) then
                   echo "ERROR: cannot communicate with PDB ftp site" > ${STMP}/$PDBCODE.error
                else
                  mv $PDBCODE.pdb $PDBDIR/$PDBCODE.pdb
                  cp $PDBDIR/$PDBCODE.pdb  ${STMP}/
                endif
              endif
            endif
            #  chmod 0777 ${STMP}/$PDBCODE.pdb
          endif
          setenv LNAME `basename ${CFILE} .job`
          setenv PARFILE ${LNAME}.par
          mv -fv $CFILE ${CFILE}.x
          # if ($CHAIN == "" ) setenv CHAIN  _
          ## mv ${GDIR}/${SJOB}/${LNAME}.par $STMP/
          $PYTHON $GDIR/rungeofold.py ${STMP}/$PARFILE > $LOGDIR/${LNAME}.log
          # chmod +x $GDIR/tmp/${LNAME}.csh
          ############ NEXT LINE SHOULD BE A JOB SUBMISSION -- MACHINE DEPENDENT
          ## qsub -b n -j y -cwd -o $LOGDIR/${LNAME}.log tmp/${LNAME}.csh $PDBCODE $CHAIN ${LNAME}  ${ONAME}
          # $GDIR/tmp/${LNAME}.csh $PDBCODE $CHAIN ${LNAME}  ${ONAME} >& $LOGDIR/${LNAME}.log &
          ## Old way...
          ## $GDIR/tmp/${LNAME}.csh >& $LOGDIR/${LNAME}.log &
          ## New python script (Erin Gilbert and Ben Walcott)
          # $GDIR/tmp/${LNAME}.csh >& $LOGDIR/${LNAME}.log &
          echo "job submitted: ${CFILE}"
          @ NDONE ++
          echo $NDONE > Done
        endif
      else
        echo "File found, and not found\!\!\!: $CFILE "  
        ## starts with space?
      endif
    end
  endif
  ## As long as the server is running, write a message to the Running file.
  echo "Running on `hostname` `date`" > Running
  sleep 5
  @ HR = `date +%k`
  if !($HR == $LHR) then
    if ($HR == 0) then
      echo "Another day has passed: `date`"
    endif
    @ LHR = $HR
  endif
## check the IP file, put heavy users in the BANNED list. They get a message that says "quota exceeded."
## sort $RECENTUSERS | uniq -c | awk -v nallow=$QUOTA '$1 > nallow {print $2}' > $BANNED
end
