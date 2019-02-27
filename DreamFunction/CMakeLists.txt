#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(DreamFunctions)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${binFolder}/DreamFunctions)
include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(SRCS
  ${PROJECT_SOURCE_DIR}/DreamCF.cxx
  ${PROJECT_SOURCE_DIR}/DreamData.cxx
  ${PROJECT_SOURCE_DIR}/DreamDist.cxx
  ${PROJECT_SOURCE_DIR}/DreamKayTee.cxx
  ${PROJECT_SOURCE_DIR}/DreamdEtadPhi.cxx
  ${PROJECT_SOURCE_DIR}/DreamPair.cxx
  ${PROJECT_SOURCE_DIR}/DreamPlot.cxx
  ${PROJECT_SOURCE_DIR}/DreamSystematics.cxx
  ${PROJECT_SOURCE_DIR}/ReadDreamFile.cxx
  ${PROJECT_SOURCE_DIR}/global.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/DreamCF.h
  ${PROJECT_SOURCE_DIR}/DreamData.h
  ${PROJECT_SOURCE_DIR}/DreamDist.h
  ${PROJECT_SOURCE_DIR}/DreamKayTee.h
  ${PROJECT_SOURCE_DIR}/DreamdEtadPhi.h
  ${PROJECT_SOURCE_DIR}/DreamPair.h
  ${PROJECT_SOURCE_DIR}/DreamPlot.h
  ${PROJECT_SOURCE_DIR}/DreamSystematics.h
  ${PROJECT_SOURCE_DIR}/ReadDreamFile.h
  ${PROJECT_SOURCE_DIR}/global.h
  )
add_library(DreamFunctions STATIC ${SRCS})
set_target_properties(DreamFunctions PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(DreamFunctions ${ROOT_LIBRARIES})

add_executable(executeDream ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFDream.C)
target_link_libraries(executeDream DreamFunctions)

add_executable(executeEmmTee ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFmT.C)
target_link_libraries(executeEmmTee DreamFunctions)

add_executable(executeKayTee ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFkT.C)
target_link_libraries(executeKayTee DreamFunctions)

add_executable(executePlot ${PROJECT_SOURCE_DIR}/Scripts/ExecutePlotCF.C)
target_link_libraries(executePlot DreamFunctions)

add_executable(executeSigma ${PROJECT_SOURCE_DIR}/Scripts/ExecuteSigma.C)
target_link_libraries(executeSigma DreamFunctions)

add_executable(executeDPhiDEta ${PROJECT_SOURCE_DIR}/Scripts/ExecuteDPhiDEta.C)
target_link_libraries(executeDPhiDEta DreamFunctions)

add_executable(executeDPhiDEtaAtRad ${PROJECT_SOURCE_DIR}/Scripts/ExecuteDPhiDEtaAtRad.C)
target_link_libraries(executeDPhiDEtaAtRad DreamFunctions)
