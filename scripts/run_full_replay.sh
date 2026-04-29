#!/bin/bash

# set environment vars, and replay data 
source set_apex_replay_env.sh 
cd ${PATH_APEX_REPLAY}

# min raw file size, in MB
min_raw_data_MB=1.0

run0=${1}
run1=${2}

run=$(( ${run0} - 1 ))

job_id_string=""

slurm_array_ids=""

#output format (CSV)
# run# , raw data (MB), [submitted/skipped],
echo "run,raw_data (MB),status"

while [[ $run -lt $run1 ]]
do
    run=$(( ${run} + 1 ))
    echo -n "Run $run --- "
    
    #check if there are any raw files

    n_raw_files=$(ls -1 ${PATH_APEX_CACHE}/apex_$run.* 2>/dev/null | wc -l)
    
    if [[ $n_raw_files -lt 1 ]]
    then
	echo "no raw files. skipped"
	continue
    fi

    echo -n "$n_raw_files raw files found. " 

    #check if the number of raw files is at least 1 MB in total 
    raw_MB=$(./scripts/get_rawdata_MB $run)
    
    if [[ $(echo "$raw_MB < $min_raw_data_MB" | bc -l) == 1 ]]
    then
        echo "less than 1 MB associated raw data ($raw_MB MB).. skipped"
        continue
    fi

    raw_GB=$(echo "scale=3; $raw_MB/1024" | bc) 
    echo -n "$raw_GB GB of associated raw data. " 
      
    # high estimate for time needed for full process
    raw_MB_per_min=180
    
    # estimate how much time we need for this run 
    mins=$( echo "scale=0; $raw_MB/$raw_MB_per_min + 20" | bc) 

    echo -n "$mins minutes alloted to run. "
    
    if [[ "~${slurm_array_ids}" == "~" ]]
    then
	slurm_array_ids="${run}"
    else
	slurm_array_ids="${slurm_array_ids},${run}"
    fi
    
    echo "submitted"
    
    job_id_string="${job_str##* } ${job_id_string}"
    
done

cmd=" --partition=production --array=${slurm_array_ids} --time=240 --job-name=apex_replay_${run0}_${run1} scripts/run-full-replay-array ${PATH_APEX_VOLATILE}/production/replay"

sbatch --test-only $cmd

echo "cmd: 'sbatch $cmd'" 

echo "${job_id_string}" > misc/job_ids_${run0}_${run1}.log 

