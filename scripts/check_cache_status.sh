#!/bin/bash

job_id=${1}

n_total=0
n_pending=0
n_hit=0
n_done=0
n_running=0

while read -r line
do
    if [[ "${line}" == "get request: "* ]];
    then
	continue
    fi
    if [[ "${line}" == "user: "* ]];
    then
	continue
    fi
    if [[ "${line}" == "status: "* ]];
    then
	continue
    fi
    
    if [[ "${line}" == *"done" ]]; then n_done=$(( $n_done + 1 )); fi
    if [[ "${line}" == *"pending" ]]; then n_pending=$(( $n_pending + 1 )); fi
    if [[ "${line}" == *"hit" ]]; then n_hit=$(( $n_hit + 1 )); fi
    if [[ "${line}" == *"running" ]]; then n_running=$(( $n_running + 1 )); fi
    
    
    n_total=$(( $n_total + 1 ))

done < <(jcache status ${job_id}) 



echo    "total    : $n_total "
echo    "_________:__________________________"
echo -n "hit      : $n_hit   ("
echo "$(echo "scale=1; 100*$n_hit/$n_total" | bc) %)" 
echo -n "pending  : $n_pending ("
echo "$(echo "scale=1; 100*$n_pending/$n_total" | bc) %)" 
echo -n "running  : $n_running ("
echo "$(echo "scale=1; 100*$n_running/$n_total" | bc) %)" 
echo -n "done     : $n_done (" 
echo "$(echo "scale=1; 100*$n_done/$n_total" | bc) %)" 
