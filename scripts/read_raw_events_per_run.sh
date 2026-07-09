#!/bin/bash

run_0=${1}
run_1=${2}

# this will be where the output goes
path_outfile=${3}

# wipe this output file
# format will be:
#
#  run-number|rawfile-number|start-event|end-event|total-events|rawfile-path
echo "" > "${path_outfile}"


run=${run_0}

while [[ ${run} -lt ${run_1} ]]
do

    event=0
    # go through each raw-file for this run.
    while read -r path_rawfile;
    do
	
	
    
    
done
