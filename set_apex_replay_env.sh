#!/bin/bash

source /home/sethhall/.login

PATH_APEX_REPLAY="/work/halla/apex/disk1/sethhall/apex_replay"

# check to see if we have the apex library loaded
echo -n "<${0}>: checking for apex library in LD_LIBRARY_PATH... "

FOUND_LIB=0

while read -d ':' line; do
    if [[ "${line}" == "${PATH_APEX_REPLAY}/src/build" ]]; then
	FOUND_LIB=1 
    fi  
done < <(echo "${LD_LIBRARY_PATH}")

if [[ ${FOUND_LIB} == 1 ]]; then
    echo "found."
else 
    echo "not found. adding apex lib to path."
    LD_LIBRARY_PATH="${PATH_APEX_REPLAY}/src/build:${LD_LIBRARY_PATH}"
fi

# now, let's define some other environment variables:
PATH_APEX_CACHE="/cache/halla/apex/raw"
PATH_APEX_VOLATILE="/volatile/halla/apex/full_replay" 
PATH_APEX_SCRIPTS="${PATH_APEX_REPLAY}/scripts"
PATH_APEX_MACROS="${PATH_APEX_REPLAY}/macros" 
 
