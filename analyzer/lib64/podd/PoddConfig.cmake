
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was PoddConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

message(STATUS "Found Podd: ${PACKAGE_PREFIX_DIR} (found version 1.6.6-apex)")

set_and_check(PODD_DIR "${PACKAGE_PREFIX_DIR}")
set_and_check(PODD_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")

if(IS_DIRECTORY "${PACKAGE_PREFIX_DIR}/lib64/podd/Modules")
  list(APPEND CMAKE_MODULE_PATH "${PACKAGE_PREFIX_DIR}/lib64/podd/Modules")
endif()

include(CMakeFindDependencyMacro)
find_dependency(ROOT 5.10)

include("${PACKAGE_PREFIX_DIR}/lib64/podd/PoddTargets.cmake")

check_required_components(PODD)

