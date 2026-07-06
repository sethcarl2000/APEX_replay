
set(headers
    DataLog.h
    interface.h
)

set(sources
    compute_track_error.cxx
    compute_trackdata.cxx
    fit_gaus_to_hist.cxx
    gen_coinc_events.cxx
    gen_pairs.cxx
    gen_pid_data.cxx
    gen_rawtracks.cxx
    gen_react_vertex.cxx
    generate_S2_hits.cxx
    generate_vdc_tracks.cxx
    get_n_cpus.cxx
    grid_search.cxx
    group_vdc_hits.cxx
    refine_track.cxx
    theta_phi_model.cxx
)

# add this static library as a target
add_library(replay_utils_lib STATIC ${sources} ${headers})

# add this as an include directory
target_include_directories(replay_utils_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link the library to other libraries that it will need
target_link_libraries(replay_utils_lib PUBLIC src_lib ROOT::Core ROOT::RIO ROOT::RDataFrame)