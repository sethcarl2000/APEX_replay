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


#APEX_IN_INCLUDE_DIR=$(find-path-in-list "${ROOT_INCLUDE_PATH}" "${PATH_APEX_REPLAY}/src")
#
#if [[ ${APEX_IN_INCLUDE_DIR} == "not found" ]]; then
#    echo "<${0}> apex not yet in include path. adding..."
#    export ROOT_INCLUDE_PATH="${ROOT_INCLUDE_PATH}:${PATH_APEX_REPLAY}/src"
#else 
#    echo "<${0}> apex src already in root include path"
#fi
    
# now, let's define some other environment variables:
export PATH_APEX_CACHE="/cache/halla/apex/raw"
export PATH_APEX_VOLATILE="/volatile/halla/apex/full_replay" 
export PATH_APEX_SCRIPTS="${PATH_APEX_REPLAY}/scripts"
export PATH_APEX_MACROS="${PATH_APEX_REPLAY}/macros" 

# make sure we have this script to add this path to the list
if [[ ! -x "/home/$(whoami)/.local/bin/add-path-to-list" ]]
then
    echo "<${0}>: Error: script 'add-path-to-list' not found in '~/.local/bin'."
    echo "you can add it by running the following bash command:"
    echo ""
    echo "     ln -s \"/home/sethhall/scripts/add-path-to-list\" \"/home/$(whoami)/.local/bin/.\"" 
    echo "" 
    exit
fi

source add-path-to-list "${PATH_APEX_REPLAY}/src" ROOT_INCLUDE_PATH

source add-path-to-list "${PATH_APEX_REPLAY}/src/build" LD_LIBRARY_PATH

exit

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


