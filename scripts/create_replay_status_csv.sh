#!/bin/bash

start_time="${1}" 

if [[ -z "${start_time}" || "${start_time}" == "--help" || "${start_time}" == "-h" ]]
then
    echo "usage: ./create_replay_status_csv.sh [start-time]"
    echo ""
    echo " i.e.:"
    echo "" 
    echo "  ./create_replay_status_csv.sh now-12hours" 
    echo ""
    echo " use 'sacct --help' for other valid formats."

    exit
fi

path_output="logs/replay_status_$(date +'%Y-%b-%d_%H.%M.%S').csv"

sacct --user=$(whoami) --starttime=${start_time} --parsable --noheader --format=JobName,JobID,Start,Elapsed,State,MaxRSS,TotalCPU > "${path_output}"

echo "created output csv under path: ${path_output}" 
