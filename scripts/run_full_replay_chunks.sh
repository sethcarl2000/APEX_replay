#!/bin/bash

# set environment vars, and replay data 
source set_apex_replay_env.sh 
cd ${PATH_APEX_REPLAY}

# min raw file size, in MB
min_raw_data_MB=1.0

run0=${1}
run1=${2}

# max number of CODA events to shoot for. If adding one more run overflows this value, then skip it. 
max_events_per_task=5000000

# this gives us the number of CODA events in each raw file 
path_coda_events_CSV="logs/n-CODA-events-per-file.csv" 

# this file will tell each array which jobs it must complete. 
array_file="${PATH_APEX_REPLAY}/array-tasks.csv"

# make an array file for the csv 
root -l -b -q "macros/create_task_csv.C(\"${path_coda_events_CSV}\", \"${array_file}\", ${run0}, ${run1}, ${max_events_per_task})" 

if [[ ! -d "$(pwd)/slurm_payloads" ]]
then
	mkdir "$(pwd)/slurm_payloads"
fi

timestamp="$(date +'%Y-%b-%d_%H.%M.%S')"
slurm_timestamp="$(date +'%Y-%m-%dT%H:%M:%S')"

path_tarball="${PATH_APEX_REPLAY}/slurm_payloads/payload__${timestamp}__run-range_${run0}-${run1}.tar.gz" 

tar -czf "${path_tarball}" build utils decode replay DB outDefs array-tasks.csv execute_array_task.C set_apex_replay_env.sh

# find the last array in the csv 
last_array_id=0
while read -r line 
do
	IFS=' ' read -ra line_array <<< "${line}"
	last_array_id="${line_array[0]}"

done < <(cat ${array_file})

echo "last array id: ${last_array_id}"

cmd_string="sbatch --job-name=apex_replay_${run0}_${run1} --array=0-${last_array_id} scripts/run-full-replay-array ${path_tarball}" 
echo "${cmd_string}" 

logfile="logs/replay_submit_${run0}_${run1}_${timestamp}.log"
errfile="logs/replay_submit_${run0}_${run1}_${timestamp}.log"

# jobstr=$(eval ${cmd_string} 2> ${errfile})
jobstr="Submitted batch job 9960067"
echo "" > ${errfile}

if [[ ! -z "$(cat ${errfile})" ]] 
then
    echo "Caught error trying to submit slurm job:"
    cat ${errfile}
fi



# get the job ID
if [[ "${jobstr}" =~ ([0-9]+) && -n "${BASH_REMATCH[0]}" ]]
then
    slurm_job_id="${BASH_REMATCH[0]}"
else
    slurm_job_id="null"
fi

echo "" > "${logfile}"

echo "# submitted full replay at: ${slurm_timestamp}" >> "${logfile}" 
echo "# run-range: [${run0}, ${run1}]" >> "${logfile}"
echo "# tarball: ${path_tarball}" >> "${logfile}"
echo "# slurm job id: ${slurm_job_id}" >> "${logfile}"
echo "# n. array tasks: $(( ${last_array_id} + 1 ))" >> "${logfile}"
echo "# array task assignments: " >> "${logfile}"
cat ${array_file} >> "${logfile}"


