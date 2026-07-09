#!/bin/bash

# set environment vars, and replay data 
source set_apex_replay_env.sh 
cd ${PATH_APEX_REPLAY}

# min raw file size, in MB
min_raw_data_MB=1.0

run0=${1}
run1=${2}

run=$(( ${run0} - 1 ))

# maximum-ish events to process
raw_event_cap=2500000

path_rawfile_data="logs/rawfile-sizes.list"

array_file_stem="${PATH_APEX_VOLATILE}/slurm/array_lists/array_${run0}_${run1}"

array_events=0
array_index=0
array_file="${array_file_stem}_${array_index}.list" 

if [[ -f "${array_file}" ]]; then rm "${array_file}"; fi 
touch "${array_file}"

min_events_per_file=750

while [[ ${run} -lt ${run1} ]]
do

    run=$(( $run + 1 ))
    
    while read -r line;
    do
	# skip the line if its empty 
	if [[ -z "${line}" ]]; then continue; fi 
	
	# now, split the line into an array
	IFS='|' read -ra line_array <<< "${line}"
	
	line_run="${line_array[0]}"
	
	if [[ "${line_run}" -lt "${run}" ]]; then continue; fi
	
	if [[ "${line_run}" -gt "${run}" ]]; then
	    run=$(( $run + 1 ))
	    if [[ "${run}" -gt "${run1}" ]]; then break; fi 
	fi
	
	rawfile_num="${line_array[1]}"
	start_event="${line_array[2]}"
	end_event="${line_array[3]}"
	events_in_file="${line_array[4]}"
	path="${line_array[5]}"
	
	if [[ $events_in_file -lt $min_events_per_file ]]
	then
	    continue;
	fi
	
	# add this file to the list of files to process
	echo "${line_run}|${rawfile_num}|${start_event}|${end_event}|${path}" >> "${array_file}" 

	array_events=$(( ${array_events} + ${events_in_file} ))
	if [[ ${array_events} -ge ${raw_event_cap} ]]
	then
	    echo "<${0}>: array task ${array_index} assigned ${array_events} raw events."
	    array_events=0
	    array_index=$(( $array_index + 1 ))
	    array_file="${array_file_stem}_${array_index}.list"
	    if [[ -f "${array_file}" ]]; then rm "${array_file}"; fi 
	    touch "${array_file}" 
	fi
        
    done < <(cat "${path_rawfile_data}")
    
done 

#now, let's submit the jobs.

#check if the last array file is empty
n_rows_last=$(cat ${array_file} | wc -l) 
if [[ $n_rows_last -lt 1 ]]
then
    array_last=$(( $array_index - 1 ))
else
    array_last=$array_index
fi


cmd_string="sbatch --job-name=apex_replay_${run0}_${run1} --array=0-${array_last} scripts/run-full-replay-array ${array_file_stem}" 
echo "${cmd_string}" 

eval ${cmd_string}
