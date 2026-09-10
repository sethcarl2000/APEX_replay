#!/bin/bash

runslist=""

task_id="${2}" 

while read -r line; 
do

    if [[ "${line}" =~ "task-id"* ]]
    then 
        continue
    fi  

    IFS=' '; read -ra line_array <<< "${line}" 

    task_i="${line_array[0]}"
    if [[ "${task_i}" != "${task_id}" ]]; then continue; fi 

    run="${line_array[1]}"

    if [[ "${runlist}" =~ ${run} ]]; then continue; fi 

    if [[ -z "${runlist}" ]]; then runlist="${run}"; else runlist="${runlist} ${run}"; fi
    
done < <(cat "${1}")

echo ${runlist}
