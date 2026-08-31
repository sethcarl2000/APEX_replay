#!/bin/bash

infile=${1}

echo "run-number|rawfile-number|n-CODA-events" 

while read -r line
do
    # the rawfile-sizes.list has the following format:
    # [run-#] | [rawfile-#] | [first-CODA-event] | [last-CODA-event] | [CODA-events-in-file] | [path-to-CODA-file]
    
    # skip blank lines
    if [[ -z "${line}" ]]; then continue; fi


    # parse the line into different bits
    IFS='|' read -ra line_array <<< "${line}"
    
    run_number=${line_array[0]}

    rawfile_number=${line_array[1]}

    n_CODA_events=${line_array[4]}

    # skip invalid files
    if [[ "${n_CODA_events}" == "null" ]]; then continue; fi 
    
    echo "${run_number}|${rawfile_number}|${n_CODA_events}"

done < <(cat ${infile}) 
    
    
    
    
