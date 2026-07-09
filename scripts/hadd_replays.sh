#!/bin/bash

source set_apex_replay_env.sh
cd $PATH_APEX_VOLATILE/production

run_min=$1
run_max=$2

step=25

run=$run_min
while [[ $run -lt $run_max ]];
do
    run0=$run
    run1=$(( $run + $step - 1 ))

    files=""

    run_i=$run0
    while [[ $run_i -lt $run1 ]]
    do
	#echo "run_i $run_i" 
	files_i=$(ls -w 1 replay.$run_i* 2>/dev/null | tr '\n' ' ')
	if [[ -z "$files" ]]
	then 
	    files="${files_i}"
	else
	    files="${files} ${files_i}"
	fi
	#echo "files: $files" 
	run_i=$(( $run_i + 1 ))
    done

    hadd -k -j 8 replay-$run0-$run1.root ${files}

    run=$(( $run + $step ))
done
