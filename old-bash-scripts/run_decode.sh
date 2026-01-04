#!/bin/bash

###########################################################
#
#  Decodes runs. optional arguments can be gvien,
#   use '-h' option to learn more.
#   Check below for a list of default parameters. 
# 
###########################################################

# Place where output .root files will be 
Output_path="$VOLATILE_DIR/optics/decode"
#.[run].root

# Variable .odef to use
OutDef_path="outDefs/optics.odef"

# Decode MB per min
DECODE_MB_per_min=300

# Mem-per-cpu, in mb
mem_per_cpu=3500

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
    echo " Usage: run_decode.sh -r {1st-run} [-l {last-run}] [options]"
    echo ""
    echo " Options:  "
    echo "  -r  -  Run number to decode (required)"
    echo
    echo "  -l  -  Last run number to decode (all runs in-between will be checked)"
    echo
    echo "  -h  -  displays this help message, then quits"
    echo
    echo "  -n  -  print the number of raw files this run currently has in the cache, "
    echo "         then exit."
    echo
    echo "  -s  -  Print the number of raw files, and their total size in the cache (in MB)"
    echo
    echo "  -p  -  Print the hard-coded parameters you want to run, then exit."
    echo
    echo "  -t  -  Submit jobs for test only, so that the queue can be tested. "
    echo
    echo "  -z  -  Submit jobs in the priority queue (default is production)"
    echo
    echo "  -m  -  Request memory-per-cpu in MB (use '-p' to see default)"
    echo
}
TestOnly="false"
Run_priority="false"


PrintParams()
{
    echo
    echo "Default parameters: (Change in run_decode.sh)"
    echo 
    echo " Output_path = '$Output_path.[run].root' "
    echo "    | Path to output .root files"
    echo 
    echo " OutDef_path = '$OutDef_path' "
    echo "    | Path to variable definitions"
    echo 
    echo " mem_per_cpu = $mem_per_cpu MB"
    echo "    | Path to variable definitions"
    echo;
}

InvalidOpt()
{
    echo 
    echo "Error: Invalid option. use '-h' to get a list of valid options."
    echo;
}

Run=0
LastRun=-1

###########################################################
# CheckNumFiles
###########################################################
N_files_raw=-1
GetNRaw()
{
    N_files_raw=$(ls -l $RAWFILE_DIR/apex_$Run.* 2>/dev/null | wc -l)
}



###########################################################
# CheckNumFiles
###########################################################
RawMB=0
GetRawMB()
{
    local n_rawBytes=0
    n_rawBytes=$(ls -l "$RAWFILE_DIR/apex_$Run."* 2>/dev/null | awk '{sum += $5} END{ print sum}')
    RawMB=$(( $n_rawBytes / 1000000 ))
}

###########################################################
# Main program
###########################################################

while getopts "hptzs:n:r:l:m:" option; do
    case ${option} in
	\?) #Invalid option
	    echo "Invalid option '-$OPTARG'; use '-h' for help"
	    exit;;
	h) # Display help
	    Help
	    exit;;
	t) # Submit jobs for test only
	    echo "t => slurm test-only mode"
	    TestOnly="true";;
	p) # Print parameters and exit
	    PrintParams
	    exit;;
	n) # Display number of raw files, then exit
	    Run=$OPTARG
	    GetNRaw
	    echo "Run $Run has $N_files_raw raw-file/s on cache."
	    exit;;
	r) # Get run number
	    Run=${OPTARG};;
	l) # Get last run number (if none is given, only 'first' run will be processed)
	    LastRun=${OPTARG};;
	z) # Partition-priority
	    echo "z => running jobs in priority partition"
	    Run_priority="true";;
	s) # Get size of all raw files
	    Run=${OPTARG}
	    GetNRaw
	    if [ $N_files_raw -lt 1 ]; then exit; fi
	    GetRawMB
	    echo "Run $Run Has $N_files_raw raw-file/s with a total size of $RawMB MB"
	    exit;;
	m) # set memory per cpu
	    mem_per_cpu="${OPTARG}"
	    echo "$mem_per_cpu MB";;
    esac
done

if [ $Run_priority = "false" ]; then echo "running on production partition."; fi


#check to make sure there's at least one raw file to decode

if [ $LastRun -le $Run ];
then
    echo "Processing run $Run"
    LastRun=$Run;
else
    echo "Processing runs $Run thru $LastRun"
fi
    
Run=$(($Run - 1))

while [ $Run -lt $LastRun ]
do
    let Run++ 

    echo -n "Run $Run: "
    
    #check to see how many raw files there are
    GetNRaw
    if [ $N_files_raw -lt 1 ];
    then
	echo "No files."
	continue
    fi 

    echo -n "$N_files_raw file/s, " 
    
    GetRawMB
    echo -n "$RawMB MB "

    minutes_decode=$(( $RawMB / $DECODE_MB_per_min + 30 ))

    echo -n "(~ $minutes_decode mins). "
    
    command="--job-name=decode_$Run --time=$minutes_decode script-run-decode $Run $Output_path $OutDef_path --mem-per-cpu=$mem_per_cpu"
    
    # decide whether to run on priority or production (-z for priority)
    if [ $Run_priority = "true" ]
    then
	command="--partition=priority ${command}";
    fi
    
    
    if [ $TestOnly = "true" ]; then command="--test-only ${command}"; fi
    
    eval "sbatch ${command}" 
        
done
