cmake_minimum_required(VERSION 3.16)

project(ApexReplay 
    VERSION 1.0 
    LANGUAGES CXX
)

find_package(Podd 1.7 REQUIRED)
 
# add a few separate shared libraries

#-------------------------------------------------------------------------------------------
#
#   First, let's build the 'utils' library
#
set(src_utils
    utils/fit_gaus_to_hist.cxx
    utils/get_available_cpus.cxx
    utils/get_env_variable_int.cxx
    utils/get_env_variable_string.cxx
    utils/get_slurm_scratch_directory.cxx
    utils/ErrorHandler.cxx
    utils/intpow.cxx
    utils/is_path_directory.cxx
    utils/is_path_regular_file.cxx
    utils/remove_file_from_disk.cxx
)

set(headers_utils
    utils/include/APEX/utils.h
    utils/include/APEX/utils/ErrorHandler.h
)

add_library(APEX_utils SHARED ${src_utils})

target_include_directories(APEX_utils PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/utils/include")

target_link_libraries(APEX_utils PUBLIC ROOT::Core ROOT::RIO)


#-------------------------------------------------------------------------------------------
#
#   First, we're going to add the 'decode' libaray, which handles operations with CODA files
#
set(src_decode
    decode/count_physics_events.cxx
    decode/init_decoders.cxx
    decode/process_coda_file.cxx 
    decode/Shower.cxx
)

set(headers_decode
    decode/include/APEX/decode.h
    decode/include/APEX/decode/Shower.h
)

add_library(APEX_decode SHARED ${src_decode})

target_include_directories(APEX_decode PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/decode/include")

# link it with Podd
target_link_libraries(APEX_decode PUBLIC Podd::Podd Podd::HallA APEX_utils)

#include(GNUInstallDirs)
#[[
set(PACKAGE APEX_decode)
install(TARGETS ${PACKAGE}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES ${headers_decode} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

build_root_dictionary(G__APEX_decode ${headers_decode} ${src_decode}
    MODULE APEX_decode 
    LINKDEF decode/LinkDef.h
)
]]

#-------------------------------------------------------------------------------------------
#
#   Next, we're going to add the 'decode' libaray, which handles operations with CODA files
#
set(src_replay
# ---   Replay manager class
    replay/Manager.cxx
# ---   User analysis classes (will be passed to a fro by the dataframe)
    replay/user_analysis_classes/EventHandler.cxx
    replay/user_analysis_classes/PidDetector.cxx
    replay/user_analysis_classes/PidManager.cxx
    replay/user_analysis_classes/ReactVertex.cxx
    replay/user_analysis_classes/RMatrixD.cxx
    replay/user_analysis_classes/S2Hit.cxx
    replay/user_analysis_classes/VDC/ChamberPair.cxx
    replay/user_analysis_classes/VDC/VDC.cxx
    replay/user_analysis_classes/VDC/HitCluster.cxx
    replay/user_analysis_classes/VDC/Track.cxx
    replay/user_analysis_classes/VDC/Hit.cxx
    replay/user_analysis_classes/VDC/HitGroup.cxx
# ---   Replay helper classes
    replay/helpers/grid_search.cxx
    replay/helpers/generate_S2_hits.cxx
    replay/helpers/gen_coinc_events.cxx
    replay/helpers/compute_trackdata.cxx
    replay/helpers/gen_coinc_window_events.cxx
    replay/helpers/generate_vdc_tracks.cxx
    replay/helpers/compute_track_error.cxx
    replay/helpers/gen_rawtracks.cxx
    replay/helpers/refine_track.cxx
    replay/helpers/theta_phi_model.cxx
    replay/helpers/group_vdc_hits.cxx
    replay/helpers/gen_pid_data.cxx
    replay/helpers/gen_pairs.cxx
    replay/helpers/gen_react_vertex.cxx
)

set(headers_replay
# ---   replay manager 
    replay/include/APEX/replay/Manager.h
# ---   replay header
    replay/include/APEX/replay.h
# ---   run parameters 
    replay/include/APEX/run_parameters.h
# ---   some unit definitions
    replay/include/APEX/units.h
# ---   replay helpers
    replay/include/APEX/replay/helpers.h
# --    misc. replay helpers
    replay/include/APEX/replay/EventCounter.h
    replay/include/APEX/replay/UniqueIDGenerator.h
# ---   user analysis classes 
    replay/include/APEX/PmtData.h
    replay/include/APEX/EventCounter.h
    replay/include/APEX/S2Hit.h
    replay/include/APEX/ReactVertex.h
    replay/include/APEX/PidDetector.h
    replay/include/APEX/EventHandler.h
    replay/include/APEX/VDC.h
    replay/include/APEX/VDC/Track.h
    replay/include/APEX/VDC/HitGroup.h
    replay/include/APEX/VDC/HitCluster.h
    replay/include/APEX/VDC/ChamberPair.h
    replay/include/APEX/VDC/Hit.h
    replay/include/APEX/PidManager.h
# ---   misc. helpers
    replay/include/RDFNodeAccumulator.h
    replay/include/RMatrixD.h
)

# add the shared library 
add_library(APEX_replay SHARED ${src_replay})

# add the include dir
target_include_directories(APEX_replay PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/replay/include")

# link against ROOT
target_link_libraries(APEX_replay PUBLIC ROOT::Core ROOT::RIO ROOT::ROOTDataFrame ROOT::MathCore ROOT::Physics APEX_utils)
