#!/bin/bash

source set_apex_replay_env.sh


PATH_LOG="$PATH_APEX_VOLATILE/slurm/replay_log"
cd $PATH_LOG 

output_log="$PATH_APEX_REPLAY/logs/full_replay.log" 

#combines all into one file
echo "slurm-array-job-id|slurm-array-task-id|path-decode-input|path-replay-output|run-number|rawfile-num|segment-num|flags" > $output_log

n_files=0
n_lines=0

while read -r filename; do
    
    #strip the run number off
    slurm_array_job_id=${filename:7:7}

    #strip off everything but the end, and then strip the extension
    slurm_array_task_id=${filename:15}
    slurm_array_task_id=${slurm_array_task_id%.*}

    while read -r line; do

	echo "${slurm_array_job_id}|${slurm_array_task_id}|${line}" >> $output_log
	n_lines=$(( $n_lines + 1 ))
	
    done < <(cat ${filename})

    n_files=$(( $n_files + 1 ))

done < <(ls -1 *)
    
echo "done.\nn. files: ${n_files}, n. lines: ${n_lines}" 
