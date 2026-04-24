#!/bin/bash

function print-usage () {
    echo " Usage:"
    echo "     cache_rawFile.sh run#0 run#1 [-p days]"
    echo " will check for raw CODA files in direcotry $PATH_APEX_MSS from run number "
    echo " run#0 to run#1 (inclusive). use -p to pin for a spefified number of days." 
    exit 0
};

if [[ "${1}" == "-h" ]]
then
    print-usage
    exit 0
fi

if [[ "${1}~" == "~" ]]
then
   echo "first arg missing."
   print-usage
   exit 1
fi
   
if [[ "${2}~" == "~" ]]
then
   echo "second arg missing."
   print-usage
   exit 1
fi
   
run0=${1}
run1=${2}


pin_days="null"
if [[ "${3}" == "-p" ]]
then
    #check if arg provided
    if [[ "${4}~" == "~" ]]
    then
	echo "option -p requires arg."
	print-usage
	exit 1
    fi
    
    pin_days=${4}

    echo "pinning all files for $pin_days days."
fi    
    

echo "Checking from run ${run0} to run ${run1}"

run=$(($run0-1))
fileNames=""

total_files=0

while [[ $run -lt $run1 ]]
do
    run=$(( $run + 1 ))

    parts=$(ls -l /mss/halla/apex/raw/apex_$run.dat.* | wc -l)

    total_files=$(( $total_files + $parts ))
    
    if [ $parts -gt 0 ] 
    then 
	
        echo "requesting run ${run} (with ${parts} partFiles)... " 
	#jcache get /mss/halla/apex/raw/apex_$run.dat.* 
	files=$(echo -n $(ls -w 0 /mss/halla/apex/raw/apex_$run.dat.*))
#	files+=" "
	
	fileNames="${fileNames} ${files}"
	#jcache get $files 
    fi 
    
done

cmd="jcache get ${fileNames} -e sethhall@jlab.org "
if [[ "${pin_days}" != "null" ]]
then
    cmd="${cmd} -D ${pin_days}"
fi


#echo "${fileNames}" 
echo "cmd: '${cmd}'" 
#echo "total files to be requested: ${total_files}"
eval "${cmd}"  
echo "done."  
