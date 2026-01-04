#!/bin/bash

###########################################################
#
#  Replays optics runs. optional arguments can be gvien,
#   use '-h' option to learn more.
#   Check below for a list of default parameters. 
# 
###########################################################

# Place where input .root files will be 
stem_input="$VOLATILE_DIR/optics/decode"
#.[run].root

# place where output root files will be
stem_output="$VOLATILE_DIR/optics/replay"
#.[run].root

# name of TTree to be found in raw files
decode_TTree_name="T"

# memory per cpu
mem_per_cpu=2500

# time per event
events_per_min=22500

####################
#  Main program
####################

###########################################################
# Help / Invalid option
###########################################################
Help()
{
    #Displays options
    echo 
    echo "run_decode.sh "
    echo "  Runs decode script on desired cached run. "
    echo "  Will use default options in decode_run.C, and submit them as jobs."
    echo ""
    echo " Usage: run_decode.sh -a [R/L] -s [start-run] {-e [end-run]} [options]"
    echo
    echo " Examples: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
    echo "   Replay runs from 4766-4777, only with LHRS:"
    echo
    echo "     $> ./run_replay_optics.sh -a L -s 4766 -l 4777"
    echo
    echo "   Replay run 4767 (L), requesting 3750 MB of Ram per slrum job:"
    echo
    echo "     $> ./run_replay_optics.sh -a L -s 4767 -m 3750"
    echo
    echo
    echo " Options:  "
    echo "  -a  -  R/L "
    echo
    echo "  -s  -  Run number to replay (required)"
    echo
    echo "  -e  -  Last run number to replay (all runs in-between will be checked)"
    echo
    echo "  -h  -  displays this help message, then quits"
    echo
    echo "  -t  -  Submit jobs for test only, so that the queue can be tested."
    echo
    echo "  -z  -  Submit jobs in the priority queue (default is production)"
    echo
    echo "  -p  -  Print default options, then exit."
    echo
    echo "  -m  -  Request MB of ram per cpu. see example above for usage."
    echo;
}
TestOnly="false"
Run_priority="false"


PrintParams()
{
    echo
    echo "Parameters:"
    echo 
    echo " Input_path     = '$stem_input.[run].root' "
    echo "    | Path to output .root files"
    echo 
    echo " Output_path    = '$stem_output.[run].root' "
    echo "    | Path to output .root files"
    echo
    echo " mem_per_cpu    = $mem_per_cpu MB"
    echo "    | RAM requested per job (each job is single-threadded)"
    echo
    echo " events_per_min = $events_per_min"
    echo "    | time per event. Used to budget time for slurm request"
    echo;
}

InvalidOpt()
{
    echo 
    echo "Error: Invalid option. use '-h' to get a list of valid options."
    echo;
}

FirstRun=0
LastRun=-1

###########################################################
# CheckNumFiles
###########################################################
N_files_raw=-1
GetNRaw()
{
    N_files_raw=$(ls -l $stem_input.$Run* 2>/dev/null | wc -l)
}


###########################################################
# SetArm
###########################################################
Arm="0"
GetArm()
{
    Arm=$1
    if [ $Arm != "R" ] && [ $Arm != "L" ]; then Arm="Null"; fi
}


###########################################################
# Main program
###########################################################

while getopts "htza:s:e:m:" option; do
    case ${option} in
	\?) #Invalid option
	    echo "Invalid option '-$OPTARG'; use '-h' for help"
	    exit;;
	h) # Display help
	    Help
	    exit;;
	p) # Print default options, then exit
	    PrintParams
	    exit;;
	t) # Submit jobs for test only
	    TestOnly="true";;
	z) # Partition-priority
	    Run_priority="true";;
	s) # Get run number
	    FirstRun=${OPTARG};;
	a) # Get arm (required)
	    GetArm ${OPTARG}
	    if [ $Arm = "Null" ];
	    then
		echo "arm '${OPTARG}' is invalid, must be either L or R"
		exit
	    fi
	    ;;
	e) # Get last run number (if none is given, only 'first' run will be processed)
	    LastRun=${OPTARG};;
	m) # Request MB of ram per job
	    mem_per_cpu=${OPTARG};;
    esac
done

if [ $Arm = "0" ]
then
    echo "Error: No arm selected. Use '-a L' or '-a R' to pick arm "
    exit;
fi


#check to make sure there's at least one raw file to decode

echo
if [ $LastRun -le $FirstRun ];
then
    echo "Processing run $FirstRun"
    LastRun=$FirstRun;
else
    echo "Processing runs $FirstRun thru $LastRun"
fi

#print which partition to use
echo
echo -n "Running on partition "
if [ $Run_priority = "true" ];
then echo -n "priority."; else echo -n "production."; fi

echo " (Use flag '-t' to check expected start time)" 


#inform the user if we're in test mode
if [ $TestOnly = "true" ]
then
    echo
    echo " ~~ Test-only mode. ~~"
fi

#Print parameters, and ask to proceed.
PrintParams

read -p "Submit jobs with current settings? [y/N]" ans

case "$ans" in
    [yY] ) echo "ok, proceeding...";;
    [nN] ) echo "exiting..."; exit;;
    * ) echo "Invalid response."; exit;;
esac


Run=$(($FirstRun - 1))


nFiles_total=0

while [ $Run -lt $LastRun ]
do
    let Run++ 
    
    #check to see how many raw files there are
    GetNRaw
    if [ $N_files_raw -lt 1 ]; then continue; fi
    
    #GetRawMB
    echo "Run $Run, $N_files_raw raw file/s..." 

    nFiles_run=0

    nEvents_run=0
    
    #loop through each sub-file which matches this pattern
    for PART_FILE in $(ls -1 $stem_input.$Run.*)
    do
	let nFiles_total++
	
	echo -n "   subfile $nFiles_run, " 
	let nFiles_run++

	#first, check that the file can be opened without issue
	err=$(bash check_TTree.sh $PART_FILE $decode_TTree_name; echo $?)

	#echo -n " err=$err, "
	
	if [ $err -eq 2 ]; then echo "File open error"; continue; fi
	if [ $err -eq 3 ]; then echo "TTree open error"; continue; fi
	
	if [ $err -ne 0 ]; then echo "Unknown error"; continue; fi
	
	#now, get the number of raw events
	nEvents_file=$(bash n_events_in_rootfile.sh $PART_FILE $decode_TTree_name)
	
	echo "$nEvents_file events,"
	
	nEvents_run=$(($nEvents_run + $nEvents_file))
    done

    
    echo -n "   ...$nEvents_run total events. "

    time_budget=$(( $nEvents_run / $events_per_min ))
    time_budget=$(( $time_budget + 20 ))
    
    
    command="--job-name=replay_optics_$Run --mem-per-cpu=$mem_per_cpu --time=$time_budget script-replay-singleArm $Run $Arm $stem_input $stem_output"

    if [ $Run_priority = "true" ];
    then
	command="--partition=priority ${command}"
    fi
    
    if [ $TestOnly = "true" ];
    then
	command="--test-only ${command}"
    fi
    
    

eval "sbatch ${command}"

done
