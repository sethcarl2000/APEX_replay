#!/bin/bash

#make sure you have set: 
# setenv RAWFILE_DIR /cache/mss/halla/apex/raw 

echo "Checking from run $1 to run $2"

run=$(($1-1))
fileNames=$()

while [ $run -lt $2 ]
do
    run=$(( $run + 1 ))
    #see if we can find any rawfiles in the cache for this run
    
    if ! test -f $VOLATILE_DIR/production/decode.$run.0.root; 
    then 
	echo "decode-file for $run not found in path $VOLATILE_DIR/production." 
	continue
    fi 
    
    echo "submitting replay job for run $run (right-arm=$arm)"
    
    sbatch --job-name=APEX_replay_$run script-run-replay $run $arm
        
done

echo "done."  
