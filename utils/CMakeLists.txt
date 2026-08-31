
#-------------------------------------------------------------------------------------------
#
#   First, let's build the 'utils' library
#
set(src_utils
    fit_gaus_to_hist.cxx
    get_available_cpus.cxx
    get_env_variable_int.cxx
    get_env_variable_string.cxx
    get_slurm_scratch_directory.cxx
    ErrorHandler.cxx
    intpow.cxx
    is_path_directory.cxx
    is_path_regular_file.cxx
    remove_file_from_disk.cxx
    get_file_size_KB.cxx
    get_task_rawfile_list.cxx
  )

set(headers_utils
    include/APEX/utils.h
    include/APEX/utils/ErrorHandler.h
)

add_library(APEX_utils SHARED ${src_utils})

target_include_directories(APEX_utils PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(APEX_utils PUBLIC ROOT::Core ROOT::RIO)

