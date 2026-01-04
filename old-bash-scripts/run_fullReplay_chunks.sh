#!/bin/bash

#make sure you have set: 
# setenv RAWFILE_DIR /cache/mss/halla/apex/raw 

echo "Checking from run $1 to run $2"

run=$(($1-1))

while [ $run -lt $2 ]
do
    run=$(( $run + 1 ))
    #see if we can find any rawfiles in the cache for this run
    
    if ! test -f $RAWFILE_DIR/apex_$run.dat.0; 
    then 
	continue
    fi 
    
    #check to make sure there's at least a bit of data
    n_rawBytes=$(ls -l $RAWFILE_DIR/apex_$run.* | awk '{sum += $5} END{ print sum}')    
    if [[ $n_rawBytes -gt 1000000 ]]; then
	
	continue
    fi 
    
    if [[ $n_rawBytes -lt 1000 ]]; then 
	
	continue 
    fi 
    
    jobStr=$(sbatch --job-name=APEX_fullReplay_smallRun_$run --time=25 --mem-per-cpu=1000 script-replay-chunk $run 2>&1)
    
    run_jobID=${jobStr##* }    

    echo "Run $run -- Submitted. jobID=$run_jobID"
        
done

echo "done."  
