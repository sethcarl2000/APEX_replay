cmake_minimum_required(VERSION 3.16)

project(CountCodaEvents LANGUAGES CXX)

find_package(Podd 1.7.7)

add_subdirectory(src)

# ------------------------------------------------------------------------------
# Add the executable which counts CODA events in a given input file.
add_executable(count_coda_events count_coda_events.cxx)

target_include_directories(
    count_coda_events
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# link it against our other libraries in our subdirectories (and ROOT)
target_link_libraries(count_coda_events PRIVATE src_lib Podd::Podd)


# ------------------------------------------------------------------------------
# Add the executable which actually decodes CODA events, and saves the result to a ROOT file.  
add_executable(decode_run decode_run.cxx)

target_include_directories(
    decode_run
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# link it against our other libraries in our subdirectories (and ROOT)
target_link_libraries(decode_run PRIVATE src_lib Podd::Podd)


