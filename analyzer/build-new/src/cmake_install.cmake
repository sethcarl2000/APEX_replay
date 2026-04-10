# Install script for directory: /work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src

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
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so.1.6.6"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so.1.6"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/libHallA.so.1.6.6"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/libHallA.so.1.6"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so.1.6.6"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so.1.6"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode:/u/group/halla/apps/evio/5.3/el9/Linux-x86_64/lib:/u/group/halla/apps/ROOT/6.30-04/el9/RelWithDebInfo/lib:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/libHallA.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so"
         OLD_RPATH "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode:/u/group/halla/apps/evio/5.3/el9/Linux-x86_64/lib:/u/group/halla/apps/ROOT/6.30-04/el9/RelWithDebInfo/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libHallA.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaADCHelicity.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCoincTime.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaG0Helicity.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPhotoReaction.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRunParameters.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackID.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCTrackID.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaAnalysisObject.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCut.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaG0HelicityReader.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPhysicsModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaS2CoincTime.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackInfo.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaAnalyzer.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCutList.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaGoldenTrack.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPidDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaSAProtonEP.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackOut.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCPointPair.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaApparatus.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaDebugModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaHRS.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCPoint.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPostProcess.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackProj.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaEpicsEvtHandler.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCChamber.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaArrayString.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaDecData.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/BdataLoc.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaHelicity.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPrimaryKine.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaScintillator.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackingDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCWire.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaAvgVertex.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaDetMap.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaHelicityDet.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPrintOption.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaSecondaryKine.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackingModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBPM.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaIdealBeam.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaQWEAKHelicity.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaShower.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTriggerTime.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVarList.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBeam.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaDetectorBase.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaInterface.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaQWEAKHelicityReader.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaSpectrometer.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTwoarmVertex.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVertexModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBeamDet.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaElectronKine.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaNamedList.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRTTI.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaSpectrometerDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaUnRasteredBeam.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVform.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBeamEloss.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaElossCorrection.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaNonTrackingDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRaster.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaString.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDC.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVhist.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBeamInfo.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaEpicsEbeam.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRasteredBeam.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaSubDetector.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCAnalyticTTDConv.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VDCeff.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaBeamModule.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaEvent.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/FileInclude.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaReacPointFoil.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTextvars.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCCluster.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCherenkov.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaExtTarCor.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaOutput.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaReactionPoint.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTotalShower.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCHit.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCluster.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaFilter.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaPIDinfo.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRun.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrack.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCPlane.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaCodaRun.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaFormula.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaParticleInfo.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaRunBase.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaTrackEloss.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaVDCTimeToDistConv.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaEvtTypeHandler.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaScalerEvtHandler.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaEvt125Handler.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/Variable.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VariableArrayVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/FixedArrayVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VectorVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/MethodVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/SeqCollectionVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/SeqCollectionMethodVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VectorObjVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VectorObjMethodVar.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/ExtraData.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/THaGlobals.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VarDef.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/src/VarType.h"
    "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/ha_compiledata.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE FILE FILES "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/libHallA_rdict.pcm")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src/analyzer")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer"
         OLD_RPATH "/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/src:/work/halla/apex/disk1/sethhall/analyzer/analyzer-APEX/build-new/hana_decode:/u/group/halla/apps/evio/5.3/el9/Linux-x86_64/lib:/u/group/halla/apps/ROOT/6.30-04/el9/RelWithDebInfo/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/analyzer")
    endif()
  endif()
endif()

