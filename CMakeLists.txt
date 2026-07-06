cmake_minimum_required(VERSION 3.16)

project(ApexReplay LANGUAGES CXX)

find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO ROOTDataFrame MathCore Physics)


add_subdirectory(src)
add_subdirectory(replay_utils)

# add the main executable
add_executable(vdc_track_replay vdc_track_replay.cxx)

target_include_directories(
    vdc_track_replay
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# link it against our other libraries in our subdirectories (and ROOT)
target_link_libraries(vdc_track_replay PRIVATE 
    src_lib
    replay_utils_lib
    ROOT::Core ROOT::RIO ROOT::ROOTDataFrame ROOT::MathCore ROOT::Physics
)

