cmake_minimum_required(VERSION 3.16)

add_compile_options("-g")

# Sources.
# List all your source files here. They will be put into a shared library
# that can be loaded from a script.
# List only the implementation files (*.cxx). For every implementation file
# there must be a corresponding header file (*.h).
set(sources
  THaApexShower.cxx
)

set(headers
  THaApexShower.h
)

# add this static library as a target
add_library(src_lib STATIC ${sources} ${headers})

# add this as an include directory
target_include_directories(src_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link the library to other libraries that it will need
target_link_libraries(src_lib PUBLIC Podd::Podd)

