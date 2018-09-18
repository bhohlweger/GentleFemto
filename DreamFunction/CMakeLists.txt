cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(GentleFemto)

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
  ${CMAKE_SOURCE_DIR}/DreamCF.cxx
  ${CMAKE_SOURCE_DIR}/DreamData.cxx
  ${CMAKE_SOURCE_DIR}/DreamDist.cxx
  ${CMAKE_SOURCE_DIR}/DreamKayTee.cxx
  ${CMAKE_SOURCE_DIR}/DreamPair.cxx
  ${CMAKE_SOURCE_DIR}/DreamPlot.cxx
  ${CMAKE_SOURCE_DIR}/ReadDreamFile.cxx
  ${CMAKE_SOURCE_DIR}/global.cxx
  )

set(HEADERS
  ${CMAKE_SOURCE_DIR}/DreamCF.h
  ${CMAKE_SOURCE_DIR}/DreamData.h
  ${CMAKE_SOURCE_DIR}/DreamDist.h
  ${CMAKE_SOURCE_DIR}/DreamKayTee.h
  ${CMAKE_SOURCE_DIR}/DreamPair.h
  ${CMAKE_SOURCE_DIR}/DreamPlot.h
  ${CMAKE_SOURCE_DIR}/ReadDreamFile.h
  ${CMAKE_SOURCE_DIR}/global.h
  )

ROOT_GENERATE_DICTIONARY(G__GentleFemto ${HEADERS} LINKDEF ${CMAKE_SOURCE_DIR}/LinkDef.h)
add_library(GentleFemto SHARED ${SRCS} G__GentleFemto)
target_link_libraries(GentleFemto ${ROOT_LIBRARIES})

add_executable(executeDream ${CMAKE_SOURCE_DIR}/Scripts/ExecuteCFDream.C)
target_link_libraries(executeDream GentleFemto)

add_executable(executeEmmTee ${CMAKE_SOURCE_DIR}/Scripts/ExecuteCFmT.C)
target_link_libraries(executeEmmTee GentleFemto)

add_executable(executeKayTee ${CMAKE_SOURCE_DIR}/Scripts/ExecuteCFkT.C)
target_link_libraries(executeKayTee GentleFemto)

add_executable(executePlot ${CMAKE_SOURCE_DIR}/Scripts/ExecutePlotCF.C)
target_link_libraries(executePlot GentleFemto)
