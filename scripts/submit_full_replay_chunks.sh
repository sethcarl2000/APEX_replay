#!/bin/bash

step=25

start=${1}
end=${2}

run0=$start

while [[ $run0 -lt $end ]]
do

    run1=$(( $run0 + $step - 1 ))

    ./scripts/run_full_replay.sh ${run0} ${run1}

    run0=$(( $run0 + $step ))

done
