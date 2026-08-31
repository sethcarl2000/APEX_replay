cmake_minimum_required(VERSION 3.16)

project(ApexReplay 
    VERSION 1.0 
    LANGUAGES CXX
)

find_package(Podd 1.7 REQUIRED)
 
# add a few separate shared libraries.
# see the repective CMakeLists.txt in each subdirectory for details

# tell all shared libraries to put their shared-libray files in one place 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 'utils' shared library
add_subdirectory(utils)

# 'decode' shared library
add_subdirectory(decode)

# 'replay' shared library 
add_subdirectory(replay)
