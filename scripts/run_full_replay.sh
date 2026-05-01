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

#the output file in which to put slurm array info.
#the structure of the '|' delimited file is as follows:
# array-id | run-number | raw-file-number | raw-file-path
slurm_array_file="${PATH_APEX_VOLATILE}/slurm/array_${run0}_${run1}.list"

#wipe the data file, if it exists
echo -n > "${slurm_array_file}" 

echo "slurm array data file: '${slurm_array_file}'" 

array_id=0

#output format (CSV)
# run# , raw data (MB), [submitted/skipped],
echo "run,n. raw files,raw_data (GB),submitted/skipped"

while [[ $run -lt $run1 ]]
do
    run=$(( ${run} + 1 ))
    
    #check if there are any raw files

    n_raw_files=$(ls -1 ${PATH_APEX_CACHE}/apex_$run.* 2>/dev/null | wc -l)

    #check if there are any raw files
    if [[ $n_raw_files -lt 1 ]]
    then
	continue
    fi

    echo -n "${run},$n_raw_files," 

    #check if the number of raw files is at least 1 MB in total 
    raw_MB=$(./scripts/get_rawdata_MB $run)
    raw_GB=$(echo "scale=6; $raw_MB/1024" | bc) 

    echo -n "${raw_GB},"
    
    if [[ $(echo "$raw_MB < $min_raw_data_MB" | bc -l) == 1 ]]
    then
        echo "skipped"
        continue
    fi

    #echo -n "$raw_GB GB of associated raw data. " 
      
    # high estimate for time needed for full process
    raw_MB_per_min=180
    
    # estimate how much time we need for this run 
    mins=$( echo "scale=0; $raw_MB/$raw_MB_per_min + 20" | bc) 

    #echo -n "$mins minutes alloted to run. "

    rawfile_num=0
    # put all of the info we need in the slurm array data file 
    while read -r line;
    do
	echo "${array_id}|${run}|${rawfile_num}|${line}" >> ${slurm_array_file}

	if [[ "~${slurm_array_ids}" == "~" ]]
	then
	    slurm_array_ids="${array_id}"
	else
	    slurm_array_ids="${slurm_array_ids},${array_id}"
	fi

	array_id=$(( ${array_id} + 1 ))
	rawfile_num=$(( ${rawfile_num} + 1 ))
	
    done < <(ls -1 ${PATH_APEX_CACHE}/apex_$run.*)
    
    
    echo "submitted"
    
    job_id_string="${job_str##* } ${job_id_string}"
    
done

cmd=" --partition=production --array=${slurm_array_ids} --time=240 --job-name=apex_replay_${run0}_${run1} scripts/run-full-replay-array ${slurm_array_file} ${PATH_APEX_VOLATILE}/production/replay"

sbatch $cmd

echo "cmd: 'sbatch $cmd'" 

echo "${job_id_string}" > misc/job_ids_${run0}_${run1}.log 

