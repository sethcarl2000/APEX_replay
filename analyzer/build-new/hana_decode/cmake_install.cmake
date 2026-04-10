# Install script for directory: /work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so.1.6.6"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so.1.6"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode/libdc.so.1.6.6"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode/libdc.so.1.6"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so.1.6.6"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so.1.6"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/u/group/halla/apps/ROOT/6.30-04/el9/RelWithDebInfo/lib:/u/group/halla/apps/evio/5.3/el9/Linux-x86_64/lib:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode/libdc.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so"
         OLD_RPATH "/u/group/halla/apps/ROOT/6.30-04/el9/RelWithDebInfo/lib:/u/group/halla/apps/evio/5.3/el9/Linux-x86_64/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libdc.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Caen1190Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Caen775Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Caen792Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/CodaDecoder.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/F1TDCModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Fadc250Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/FastbusModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/GenScaler.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Lecroy1875Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Lecroy1877Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Lecroy1881Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Module.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/PipeliningModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Scaler1151.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Scaler3800.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Scaler3801.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Scaler560.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/SimDecoder.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaCodaData.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaCodaDecoder.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaCodaFile.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaCrateMap.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaEpics.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaEvData.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaFastBusWord.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaSlotData.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaUsrstrutils.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/VmeModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/THaBenchmark.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/hana_decode/Decoder.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE FILE FILES "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode/libdc_rdict.pcm")
endif()

