#!/bin/bash

source /home/sethhall/.login



export PATH_APEX_REPLAY="/work/halla/apex/disk1/sethhall/apex_replay"
#export PATH_APEX_ANALYZER="${PATH_APEX_REPLAY}/analyzer" 

# check to see if we have the apex library loaded
# echo -n "<${0}>: checking for apex library in LD_LIBRARY_PATH... "

#ANALYZER_SRC="/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src"
#ANALYZER_LD_PATH="${PATH_APEX_ANALYZER}/lib64"

#
#ANALYZER_EXE="${PATH_APEX_ANALYZER}/bin/analyzer"

# if someone has already defined a version of the analyzer, let's replace it.
#NEW_LD_LIBRARY_PATH=""

# echo ${LD_LIBRARY_PATH} 
#while read -d ':' line; do
#    # only keep items in this list of paths that are NOT the analyzer 
#    if [[ ${line} != *"/analyzer/"* ]]; then
#
#	#check if this is the first lib loaded
#	if [[ "${NEW_LD_LIBRARY_PATH}~" == "~" ]]; then
#	    # if it is, start it 'fresh' 
#	    NEW_LD_LIBRARY_PATH="${line}"
#	else
#	    # otherwise, append this path to the end 
#	    NEW_LD_LIBRARY_PATH="${line}:${NEW_LD_LIBRARY_PATH}" 
#	fi
#   fi 
#done < <(echo "${LD_LIBRARY_PATH}")

#now, add back in *our* version of the analyzer 
#export LD_LIBRARY_PATH="${ANALYZER_LD_PATH}:${NEW_LD_LIBRARY_PATH}" 

LOADED_APEX=0
while read -d ':' line; do
    if [[ "${line}" == "${PATH_APEX_REPLAY}/src/build" ]]; then
	LOADED_APEX=1 
    fi  
done < <(echo "${LD_LIBRARY_PATH}")

if [[ ${LOADED_APEX} == 0 ]]; then
    export LD_LIBRARY_PATH="${PATH_APEX_REPLAY}/src/build:${LD_LIBRARY_PATH}"
fi


# now, let's define some other environment variables:
export PATH_APEX_CACHE="/cache/halla/apex/raw"
export PATH_APEX_VOLATILE="/volatile/halla/apex/full_replay" 
export PATH_APEX_SCRIPTS="${PATH_APEX_REPLAY}/scripts"
export PATH_APEX_MACROS="${PATH_APEX_REPLAY}/macros" 
 
