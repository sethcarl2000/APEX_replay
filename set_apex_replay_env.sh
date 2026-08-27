#!/bin/bash

export PATH_APEX_REPLAY="/work/halla/apex/disk1/sethhall/apex_replay"
#export PATH_APEX_ANALYZER="${PATH_APEX_REPLAY}/analyzer" 

# set path to our DB files 
export DB_DIR="${PATH_APEX_REPLAY}/DB"
    
# now, let's define some other environment variables:
export PATH_APEX_CACHE="/cache/halla/apex/raw"
export PATH_APEX_VOLATILE="/volatile/halla/apex/full_replay" 
export PATH_APEX_SCRIPTS="${PATH_APEX_REPLAY}/scripts"
export PATH_APEX_MACROS="${PATH_APEX_REPLAY}/macros" 

add_path_to_list () {

    my_path="${1}"
    path_list="${2}" 
    
    # loop over all paths in the list. if we find a match for the given path, then quit (we don't need to add it). 
    local found="false"
    IFS=':' read -ra path_array <<< "${path_list}" 
    for path_item in "${path_array[@]}" 
    do
        if [[ "${path_item}" == "${my_path}" ]]; then 
            found="true"
            break; 
        fi
    done

    if [[ $found == "false" ]]; then 
        if [[ -z "${path_list}" ]]; then path_list="${my_path}"; else path_list="${my_path}:${path_list}"; fi
    fi 
    echo "${path_list}"
}

export ROOT_INCLUDE_PATH=$(add_path_to_list "${PATH_APEX_REPLAY}/decode" "${ROOT_INCLUDE_PATH}")
export LD_LIBRARY_PATH=$(add_path_to_list "${PATH_APEX_REPLAY}/decode/build" "${LD_LIBRARY_PATH}")


