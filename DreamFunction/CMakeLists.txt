cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(GentleFemto)

# IMPORTANT - SET HERE THE PATH TO YOUR BASE DIR
SET(BaseDir "~/Software/GentleFemto/DreamFunction")

# DONT CHANGE THE REST
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(SRCS
  ${BaseDir}/DreamCF.cxx
  ${BaseDir}/DreamData.cxx
  ${BaseDir}/DreamDist.cxx
  ${BaseDir}/DreamKayTee.cxx
  ${BaseDir}/DreamPair.cxx
  ${BaseDir}/DreamPlot.cxx
  ${BaseDir}/ReadDreamFile.cxx
  ${BaseDir}/global.cxx
  )

set(HEADERS
  ${BaseDir}/DreamCF.h
  ${BaseDir}/DreamData.h
  ${BaseDir}/DreamDist.h
  ${BaseDir}/DreamKayTee.h
  ${BaseDir}/DreamPair.h
  ${BaseDir}/DreamPlot.h
  ${BaseDir}/ReadDreamFile.h
  ${BaseDir}/global.h
  )

ROOT_GENERATE_DICTIONARY(G__GentleFemto ${HEADERS} LINKDEF ${BaseDir}/LinkDef.h)
add_library(GentleFemto SHARED ${SRCS} G__GentleFemto)
target_link_libraries(GentleFemto ${ROOT_LIBRARIES})

add_executable(executeDream ${BaseDir}/Scripts/ExecuteCFDream.C)
target_link_libraries(executeDream GentleFemto)

add_executable(executeKayTee ${BaseDir}/Scripts/ExecuteCFmT.C)
target_link_libraries(executeKayTee GentleFemto)

add_executable(executePlot ${BaseDir}/Scripts/ExecutePlotCF.C)
target_link_libraries(executePlot GentleFemto)
