#!/bin/bash

source /home/sethhall/.login


export PATH_APEX_REPLAY="/work/halla/apex/disk1/sethhall/apex_replay"
#export PATH_APEX_ANALYZER="${PATH_APEX_REPLAY}/analyzer" 

# set path to our DB files 
DB_DIR="${PATH_APEX_REPLAY}/analyzer/DB"

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

function find-path-in-list () {
    while read -d ':' line; do
	if [[ "${line}" == "${2}" ]]; then
	    echo "found"
	    exit 0
	fi
    done < <(echo "${1}")
    echo "not found"
    exit 0
}



LOADED_APEX=$(find-path-in-list "${LD_LIBRARY_PATH}" "${PATH_APEX_REPLAY}/src/build")
#while read -d ':' line; do
#    if [[ "${line}" == "${PATH_APEX_REPLAY}/src/build" ]]; then
#	LOADED_APEX=1 
#	break; 
#    fi  
#done < <(echo "${LD_LIBRARY_PATH}")

if [[ ${LOADED_APEX} == "not found" ]]; then
    echo "<${0}> apex lib not yet loaded" 
    export LD_LIBRARY_PATH="${PATH_APEX_REPLAY}/decode/build:${LD_LIBRARY_PATH}"
else
    echo "<${0}> apex lib already loaded"
fi

APEX_IN_INCLUDE_DIR=$(find-path-in-list "${ROOT_INCLUDE_PATH}" "${PATH_APEX_REPLAY}/src")

if [[ ${APEX_IN_INCLUDE_DIR} == "not found" ]]; then
    echo "<${0}> apex not yet in include path. adding..."
    export ROOT_INCLUDE_PATH="${ROOT_INCLUDE_PATH}:${PATH_APEX_REPLAY}/src"
else 
    echo "<${0}> apex src already in root include path"
fi
    
# now, let's define some other environment variables:
export PATH_APEX_CACHE="/cache/halla/apex/raw"
export PATH_APEX_VOLATILE="/volatile/halla/apex/full_replay" 
export PATH_APEX_SCRIPTS="${PATH_APEX_REPLAY}/scripts"
export PATH_APEX_MACROS="${PATH_APEX_REPLAY}/macros" 

# add our macros to the list of macros
if [ ! -f "~/.rootrc" ]; then
    echo "~/.rootrc file does not exit. making it:"
    echo "Unix.*.Root.MacroPath: .:${PATH_APEX_MACROS}" > ~/.rootrc
else

    macros=$(cat "~/.rootrc")

    if [[ "${macros}~" == "~" ]]; then
	 echo "~/.rootrc file does not exit. making it:"
	 "Unix.*.Root.MacroPath: .:${PATH_APEX_MACROS}" > ~/.rootrc
    fi

    if [[ "${macros}" != *${PATH_APEX_MACROS}* ]]; then
	echo "${macros}:${PATH_APEX_MACROS}" > ~/.rootrc
	echo "adding apex macro to list of macro paths" 
    fi
fi


