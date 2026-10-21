# Bash scripts

## run_full_replay_chunks.sh

### Array tasks
First, it creates a list of **array tasks**. The entire replay is submitted as one slurm array; which means that each job in the array is assigned a unique index. this script decides which specific runs each 'task' in the array will execute: 

```bash 
# max number of CODA events to shoot for. If adding one more run overflows this value, then skip it. 
max_events_per_task=4000000

# this gives us the number of CODA events in each raw file 
path_coda_events_CSV="logs/n-CODA-events-per-file.csv" 

# this file will tell each array which jobs it must complete. 
array_file="${PATH_APEX_REPLAY}/array-tasks.csv"

# make an array file for the csv 
root -l -b -q "macros/create_task_csv.C(\"${path_coda_events_CSV}\", \"${array_file}\", ${run0}, ${run1}, ${max_events_per_task})" 
``` 

The `max_events_per_task` is the maximum number of raw CODA events that any task will be assigned; this figure is used as a rough proxy for the amount of time a task will take to execute. The output csv `array-tasks.csv` is a space-delimited csv file, with 3 columns: `task-id|run-number|rawfile-number`. Each task reads a copy of this file, so it knows which rawfiles to process. 

### Tarball code payload
Next, a copy of all the code which will be needed to execute the replay is compressed into a compressed _tarball_: 

```bash
if [[ ! -d "$(pwd)/slurm_payloads" ]]
then
    mkdir "$(pwd)/slurm_payloads"
fi

timestamp="$(date +'%Y-%b-%d_%H.%M.%S')"
slurm_timestamp="$(date +'%Y-%m-%dT%H:%M:%S')"

path_tarball="${PATH_APEX_REPLAY}/slurm_payloads/payload__${timestamp}__run-range_${run0}-${run1}.tar.gz" 

tar -czf "${path_tarball}" build utils decode scripts replay DB outDefs array-tasks.csv execute_specific_run.C set_apex_replay_env.sh
```

Each task will make a copy of this tarball and unpack it to execute the code. This way, a copy of the exact code used to produce a given full-replay pass is available on the ifarm filesystem (and will be saved to tape). 

### Submitting the Slurm task array, and logging the submission: 
Now, the slurm array is submitted (using [run-full-replay-array](#run-full-replay-array), and the command used to submit it is copied to an auto-generated logfile (along with other misc. parameters, such as the current timestamp of submission): 

```bash
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

jobstr=$(eval ${cmd_string} 2> ${errfile})
# jobstr="Submitted batch job 9960067"
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
```

## run-full-replay-array
the script `run-full-replay-array` is what is executed by each slurm task. 
First, it copies the tarball to its allocated `/scratch` directory, and unpacks it: 
```bash
# move to the scratch dir 
PATH_SCRATCH="/scratch/slurm/${SLURM_JOB_ID}" 

FUNC="run-full-replay-array" 

# check to make sure this directory exists
if [[ ! -d "${PATH_SCRATCH}" ]]
then
    echo "scratch directory does not exist! path: ${PATH_SCRATCH}" >&2 
    exit 1 
fi

cd ${PATH_SCRATCH}

path_tarball=${1}
echo "<${FUNC}> [$(date +'%y-%m-%dT%H:%M:%S')]: In script." 
echo "<${FUNC}>: unpacking tarball: ${path_tarball}" 

# unpack the tarball of data / executables
tar -xzf ${path_tarball}

echo "<${FUNC}>: done."
```
Some miscellaneous information is printed, to be saved in the slurm `.out` logfile: 

```bash
# set up APEX env
echo "<${FUNC}>: setting env"

source set_apex_replay_env.sh
echo "<${FUNC}>: pwd: '$(pwd)'" 

OUTPUT_DIR="${PATH_APEX_VOLATILE}/production" 

# load the proper version of the analyzer
#module unload analyzer 
#module load analyzer/1.7.15 

echo "<${FUNC}>: analyzer version:" 
analyzer --version

echo "<${FUNC}>: root version:" 
root --version

echo "<${FUNC}>: uname -a (OS info):"
uname -a
```

We check to make sure the csv of runs/files for this array exists, and then we get a list of all runs we need to process: 

```bash
# check to make sure we have an 'array tasks' file 

PATH_TASK_LIST="array-tasks.csv"

if [[ ! -f "${PATH_TASK_LIST}" ]]
then    
    echo "Task file does not exist! expected rel. path: ${PATH_TASK_LIST}; directory contents: " >&2 
    ls -1 * >&2 
    exit 1 
fi

echo "<${FUNC}>: Starting timestamp: $(date +'%Y.%m.%d %H:%M:%S')" 

# now, execute all tasks: 
echo "<${FUNC}>: Executing task macro:"

# get list of all runs
runlist=$(./scripts/get_runs_in_task.sh "${PATH_TASK_LIST}" ${SLURM_ARRAY_TASK_ID})

if [[ -z "${runlist}" ]]
then
    echo "<${FUNC}>: no runs found for task: ${SLURM_ARRAY_TASK_ID}!!" >&2 
    exit 
fi
```

The macro `execute_specific_run.C` is executed. This macro, well, executes a specific run. Specifically, it looks for all rawfiles associated with a given run in the task array file, and processes them in order. See that file's documentation for details. 
