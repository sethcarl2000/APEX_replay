cmake_minimum_required(VERSION 3.5)

add_compile_options("-g")

# Sources.
# List all your source files here. They will be put into a shared library
# that can be loaded from a script.
# List only the implementation files (*.cxx). For every implementation file
# there must be a corresponding header file (*.h).
set(src
  ApexPidDetector.cxx
  ApexPidManager.cxx
  ApexUtils.cxx
  ApexVDC.cxx
  ApexVDCChamberPair.cxx
  ApexVDCHit.cxx
  ApexVDCHitCluster.cxx
  ApexVDCHitGroup.cxx
  ApexVDCTrack.cxx
  RMatrixD.cxx
  THaApexShower.cxx
  THaShower.cxx
  TapexEventHandler.cxx
  TapexReactVertex.cxx
  TapexS2Hit.cxx
)

set(headers
  ApexPidDetector.h
  ApexPidManager.h
  ApexUtils.h
  ApexVDC.h
  ApexVDCChamberPair.h
  ApexVDCHit.h
  ApexVDCHitCluster.h
  ApexVDCHitGroup.h
  ApexVDCTrack.h
  ArmMode.h
  EventCounter.h
  PmtData.h
  RMatrixD.h
  THaApexShower.h
  THaShower.h
  TTargetRays.h
  TapexEventHandler.h
  TapexReactVertex.h
  TapexS2Hit.h
  TrdfHelper.h
  UniqueIDGenerator.h
)

# generate the ROOT dictionary
root_generate_dictionary(G__src_lib ${headers} LINKDEF ApexHRS_LinkDef.h)

# add this static library as a target
add_library(src_lib STATIC ${sources} ${headers} G__src_lib.cxx)

# add this as an include directory
target_include_directories(src_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link the library to other libraries that it will need
target_link_libraries(src_lib PUBLIC ROOT::Core ROOT::RIO ROOT::RDataFrame)

