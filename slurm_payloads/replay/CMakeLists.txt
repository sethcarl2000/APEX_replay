
#-------------------------------------------------------------------------------------------
#
#   Next, we're going to add the 'decode' libaray, which handles operations with CODA files
#
set(src_replay
# ---   Replay manager class
    Manager.cxx
# ---   User analysis classes (will be passed to a fro by the dataframe)
    user_analysis_classes/EventHandler.cxx
    user_analysis_classes/PidDetector.cxx
    user_analysis_classes/PidManager.cxx
    user_analysis_classes/ReactVertex.cxx
    user_analysis_classes/RMatrixD.cxx
    user_analysis_classes/S2Hit.cxx
    user_analysis_classes/VDC/ChamberPair.cxx
    user_analysis_classes/VDC/VDC.cxx
    user_analysis_classes/VDC/HitCluster.cxx
    user_analysis_classes/VDC/Track.cxx
    user_analysis_classes/VDC/Hit.cxx
    user_analysis_classes/VDC/HitGroup.cxx
# ---   Replay helper classes
    helpers/grid_search.cxx
    helpers/generate_S2_hits.cxx
    helpers/gen_coinc_events.cxx
    helpers/compute_trackdata.cxx
    helpers/gen_coinc_window_events.cxx
    helpers/generate_vdc_tracks.cxx
    helpers/compute_track_error.cxx
    helpers/gen_rawtracks.cxx
    helpers/refine_track.cxx
    helpers/theta_phi_model.cxx
    helpers/group_vdc_hits.cxx
    helpers/gen_pid_data.cxx
    helpers/gen_pairs.cxx
    helpers/gen_react_vertex.cxx
)

set(headers_replay
# ---   replay manager 
    include/APEX/replay/Manager.h
# ---   replay header
    include/APEX/replay.h
# ---   run parameters 
    include/APEX/run_parameters.h
# ---   some unit definitions
    include/APEX/units.h
# ---   replay helpers
    include/APEX/replay/helpers.h
# --    misc. replay helpers
    include/APEX/replay/EventCounter.h
    include/APEX/replay/UniqueIDGenerator.h
# ---   user analysis classes 
    include/APEX/PmtData.h
    include/APEX/EventCounter.h
    include/APEX/S2Hit.h
    include/APEX/ReactVertex.h
    include/APEX/PidDetector.h
    include/APEX/EventHandler.h
    include/APEX/VDC.h
    include/APEX/VDC/Track.h
    include/APEX/VDC/HitGroup.h
    include/APEX/VDC/HitCluster.h
    include/APEX/VDC/ChamberPair.h
    include/APEX/VDC/Hit.h
    include/APEX/PidManager.h
# ---   misc. helpers
    include/RDFNodeAccumulator.h
    include/RMatrixD.h
)

# add the shared library 
add_library(APEX_replay SHARED ${src_replay})

# add the include dir
target_include_directories(APEX_replay PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# link against ROOT
target_link_libraries(APEX_replay PUBLIC ROOT::Core ROOT::RIO ROOT::ROOTDataFrame ROOT::MathCore ROOT::Physics APEX_utils)
