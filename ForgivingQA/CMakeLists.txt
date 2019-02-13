#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(ForgivingQA)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${binFolder}/ForgivingQA)
include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(SRCS
  ${PROJECT_SOURCE_DIR}/ForgivingReader.cxx
  ${PROJECT_SOURCE_DIR}/EventQA.cxx
  ${PROJECT_SOURCE_DIR}/TrackQA.cxx
  ${PROJECT_SOURCE_DIR}/DecayQA.cxx
  ${PROJECT_SOURCE_DIR}/ForgivingFitter.cxx
  ${PROJECT_SOURCE_DIR}/MakeHistosGreat.cxx
 )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/ForgivingReader.h
  ${PROJECT_SOURCE_DIR}/EventQA.h
  ${PROJECT_SOURCE_DIR}/TrackQA.h
  ${PROJECT_SOURCE_DIR}/DecayQA.h
  ${PROJECT_SOURCE_DIR}/ForgivingFitter.h
  ${PROJECT_SOURCE_DIR}/MakeHistosGreat.h
  )
add_library(ForgivingFnct STATIC ${SRCS})
set_target_properties(ForgivingFnct PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(ForgivingFnct ${ROOT_LIBRARIES})

add_executable(timeForQA ${PROJECT_SOURCE_DIR}/Scripts/MakeQA.C)
target_link_libraries(timeForQA ForgivingFnct)
