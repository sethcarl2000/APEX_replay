#!/bin/bash

echo "Checking from run $1 to run $2"

run=$(($1-1))
fileNames=$()

while [ $run -lt $2 ]
do
    run=$(( $run + 1 ))

    parts=$(ls -l /mss/halla/apex/raw/apex_$run.dat.* | wc -l)
    
    if [ $parts -gt 0 ] 
    then 
	
        echo "requesting run $run (with $parts partFiles)... " 
	#jcache get /mss/halla/apex/raw/apex_$run.dat.* 
	files=$(ls /mss/halla/apex/raw/apex_$run.dat.*)
	files+=" "
	
	fileNames+=$files
	#jcache get $files 
    fi 
    
done

jcache get $fileNames -e sethhall@jlab.org
echo "done."  
