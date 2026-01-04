#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(ApexReplayLibrary VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Define this preprocessor macro to enable debug statements 
option(DEBUG "Enable debug printing" OFF)
if (DEBUG)
    add_compile_definitions(DEBUG)
endif()

#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    17
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics)

#---------------------------------------------------------
#   Define sources and headers
set(sources 
    TapexS2Hit.cxx
)

set(headers 
    TapexS2Hit.h
    TapexS2PMTOffsets.h
)

#---------------------------------------------------------
#   Generate the ROOT dictionary
root_generate_dictionary(G__ApexReplay ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(ApexReplay SHARED ${headers} ${sources} G__ApexReplay.cxx)
target_include_directories(ApexReplay PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(ApexReplay PUBLIC ROOT::Core)

#---------------------------------------------------------
#   Add application executable
#add_executable(pick_sieveholes pick_sieveholes.cxx G__ApexReplay.cxx ${headers} ${sources})

#target_include_directories(pick_sieveholes PRIVATE
#    "${CMAKE_CURRENT_SOURCE_DIR}"
#    "${ROOT_INCLUDE_DIRS}"
#)

#target_link_libraries(pick_sieveholes PRIVATE
#    ${ROOT_LIBRARIES}
#    ApexReplay
#)
