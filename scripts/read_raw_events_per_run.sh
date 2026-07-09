#!/bin/bash

run_0=${1}
run_1=${2}

source set_apex_replay_env.sh
cd ${PATH_APEX_REPLAY}

# this will be where the output goes
path_outfile=${3}

# wipe this output file
# format will be:
#
#  run-number|rawfile-number|start-event|end-event|total-events|rawfile-path
if [[ -f "${path_outfile}" ]]; then rm ${path_outfile}; touch ${path_outfile}; fi


run=$(( ${run_0} - 1 ))

while [[ ${run} -lt ${run_1} ]]
do

    run=$(( ${run} + 1 ))
    #m
    #m
    #m
    #m
    #m
    #m
    # - muon's comment (8 Jul 26)

    rawfile_num=0
    
    path_rawfile="${PATH_APEX_CACHE}/apex_${run}.dat.${rawfile_num}" 
    
    files=""
    
    while [[ -f "${path_rawfile}" ]]
    do
	if [[ $rawfile_num -gt 0 ]]; then
	    files="${files}, \"${path_rawfile}\""
	else
	    files="\"${path_rawfile}\""
	fi
	
	rawfile_num=$(( ${rawfile_num} + 1 ))	
	path_rawfile="${PATH_APEX_CACHE}/apex_${run}.dat.${rawfile_num}"

    done
    
    analyzer -l -b -q "decode/count_coda_events_array.C(${run}, {${files}}, \"${path_outfile}\")"
        
done

