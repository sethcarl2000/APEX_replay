<<<<<<< HEAD
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
#   Define this preprocessor macro to range-checks for user-defined collections of objects
#   as well as constexpr namespace arrays.
#   this effects performance, makes errors easier to interpret during development 
option(RANGE_CHECKS "Enable Range-check exceptions on user-defined classes" ON)
if (RANGE_CHECKS)
    add_compile_definitions(RANGE_CHECKS)
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
    TapexEventHandler.cxx
    TapexVDCNamespace.cxx
    TapexVDCHit.cxx
    TapexVDCHitGroup.cxx
#    TapexVDCChamberPair.cxx
)

set(headers 
    EventCounter.h
    TapexS2Hit.h
    TapexS2PMTOffsets.h
    TapexEventHandler.h
    TapexVDCNamespace.h
    TapexVDCHit.h
    TapexVDCHitGroup.h
#    TapexVDCHitCluster.h 
#    TapexVDCChamberPair.h 
#    TapexVDCTrack.h
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
=======
cmake_minimum_required(VERSION 3.5)

project(ApexReplay LANGUAGES CXX)

#----------------------------------------------------------------------------
# Names of source files and target libraries.
# You do want to modify this section.

# The name of your package.
# The generated shared library will be named lib$(PACKAGE).so/dylib
set(PACKAGE ApexHRS)

# Sources.
# List all your source files here. They will be put into a shared library
# that can be loaded from a script.
# List only the implementation files (*.cxx). For every implementation file
# there must be a corresponding header file (*.h).
set(src
  ApexHRS.cxx
#  TNPoly.cxx
#  TXMap.cxx
#  TOpticsModel.cxx
#  OpticsVectors.cxx
#  OpticsPolynomials.cxx
#  TOpticsRays.cxx
  ApexUtils.cxx
  RMatrixD.cxx
  TReactVertex.cxx
  THaApexShower.cxx
)

# Headers.
# Leave the following line untouched.
string(REPLACE .cxx .h headers "${src}")

# If you have any extra headers needed for building the ROOT dictionary,
# add them here. Otherwise leave this line commented out.
#list(APPEND headers UserClassHeader.h)

# If you have other headers NOT needed for the ROOT dictionary but needed for
# the library, add them here as shown in the line commented out below.
set(allheaders ${headers})
#set(allheaders ${headers} StandaloneHeader1.h StandaloneHeader2.h)

#------------------------------------------------------------------------------
# Do not change anything below here unless you know what you are doing

#----------------------------------------------------------------------------
# Find Podd, if necessary, and load its CMake module
if(NOT TARGET Podd::HallA)
  if(DEFINED ENV{ANALYZER})
    set(CMAKE_PREFIX_PATH "$ENV{ANALYZER}" "${CMAKE_PREFIX_PATH}")
  endif()

  find_package(Podd 1.7 REQUIRED)
endif()
include(PoddCMakeEnv)

set_diagnostic_flags(WALL WEXTRA)
#report_build_info()

#----------------------------------------------------------------------------
# Library
add_library(${PACKAGE} SHARED ${src} ${allheaders} ${PACKAGE}Dict.cxx)

target_compile_options(${PACKAGE} PRIVATE ${${PROJECT_NAME_UC}_DIAG_FLAGS_LIST})

target_include_directories(${PACKAGE}
  PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
target_link_libraries(${PACKAGE} PUBLIC Podd::HallA)

include(GNUInstallDirs)

install(TARGETS ${PACKAGE}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(FILES ${allheaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#----------------------------------------------------------------------------
# ROOT dictionary
build_root_dictionary(${PACKAGE} ${headers}
  TARGETS ${PACKAGE}
  LINKDEF ${PACKAGE}_LinkDef.h
  )
>>>>>>> origin/old-replay-software
