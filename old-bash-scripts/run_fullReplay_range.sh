#!/bin/bash

#make sure you have set: 
# setenv RAWFILE_DIR /cache/mss/halla/apex/raw 

echo "Checking from run $1 to run $2"

arm=$3

run=$(($1-1))
fileNames=$()

while [ $run -lt $2 ]
do
    run=$(( $run + 1 ))
    #see if we can find any rawfiles in the cache for this run
    
    if ! test -f $RAWFILE_DIR/apex_$run.dat.0; 
    then 
	echo "Run $run not found in cache." 
	continue
    fi 
    
    n_parts_raw=$(ls -l $RAWFILE_DIR/apex_$run.* | wc -l)
    
    #do we at least have 1 MB of raw data (if not, let's not waste time/RAM!) 
    n_rawBytes=$(ls -l $RAWFILE_DIR/apex_$run.* | awk '{sum += $5} END{ print sum}')    
    if [ $n_rawBytes -lt 1000000 ] 
    then 
	echo "Skipping run $run, less than 1 MB associated raw data.." 
	continue
    fi 
    
    #first, submit the 'decode' job. we'll need this before we can start the 
    #  track replay using the new code. 
    n_raw_MB=$(( $n_rawBytes / 1000000 )) 
    
    #time, in minutes, that we estimate this job will need (this is a very high
    #  guess, to be safe.)
    DECODE_MB_per_min=400
        
    minutes_decode=$(( $n_raw_MB / $DECODE_MB_per_min + 20 )) 
    
    echo "Run $run, $n_parts_raw subfiles, $n_raw_MB MB raw file size (est. max $minutes_decode mins. decoding)..."
    
    jobStr_decode=$(sbatch --job-name=APEX_decode_$run --time=$minutes_decode script-run-decode $run 2>&1)
    
    jobID_decode=${jobStr_decode##* }
    
    # Note on some of the parameters:
    #  
    #  --dependency=afterok:${jobID_decode}  ~  don't start this job unless the 
    #                                           last one finished ok 
    #  --kill-on-invalid-dep=yes             ~  if the last job fails, then kill 
    #                                           this one, too. 
    # now, if the 'decode' was successful, we can actuall replay the data. 
    sbatch --dependency=afterok:${jobID_decode} --time=100 --kill-on-invalid-dep=yes --job-name=APEX_genPts_$run script-run-replay $run 
    
done

echo "done."  
