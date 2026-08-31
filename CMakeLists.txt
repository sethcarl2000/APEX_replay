cmake_minimum_required(VERSION 3.16)

project(ApexReplay 
    VERSION 1.0 
    LANGUAGES CXX
)

find_package(Podd 1.7 REQUIRED)
 
# add a few separate shared libraries.
# see the repective CMakeLists.txt in each subdirectory for details

# 'utils' shared library
add_subdirectory(utils)

# 'decode' shared library
add_subdirectory(decode)

# 'replay' shared library 
add_subdirectory(replay)
