#!/bin/bash

#SBATCH --output=/work/halla/apex/disk1/sethhall/apex_replay/logs/%x.out
#SBATCH --error=/work/halla/apex/disk1/sethhall/apex_replay/logs/%x.err

source set_apex_replay_env.sh
cd ${PATH_APEX_REPLAY}
echo "<${0}>: pwd: '$(pwd)'" 

time root -l -b -q macros/draw_dt.C 
