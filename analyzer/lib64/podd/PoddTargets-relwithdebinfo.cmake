#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Podd::Decode" for configuration "RelWithDebInfo"
set_property(TARGET Podd::Decode APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Podd::Decode PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELWITHDEBINFO "EVIO::EVIO"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib64/libdc.so.1.6.6"
  IMPORTED_SONAME_RELWITHDEBINFO "libdc.so.1.6"
  )

list(APPEND _IMPORT_CHECK_TARGETS Podd::Decode )
list(APPEND _IMPORT_CHECK_FILES_FOR_Podd::Decode "${_IMPORT_PREFIX}/lib64/libdc.so.1.6.6" )

# Import target "Podd::HallA" for configuration "RelWithDebInfo"
set_property(TARGET Podd::HallA APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Podd::HallA PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib64/libHallA.so.1.6.6"
  IMPORTED_SONAME_RELWITHDEBINFO "libHallA.so.1.6"
  )

list(APPEND _IMPORT_CHECK_TARGETS Podd::HallA )
list(APPEND _IMPORT_CHECK_FILES_FOR_Podd::HallA "${_IMPORT_PREFIX}/lib64/libHallA.so.1.6.6" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
