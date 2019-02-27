#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Systematics)

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
SET(KITTY_PATH "${CMAKE_SOURCE_DIR}/GentleKitty")#where are all CATS related .h files
SET(DREAM_PATH "${CMAKE_SOURCE_DIR}/DreamFunction")#where are all CATS related .h files
SET(FORGIVEN_PATH "${CMAKE_SOURCE_DIR}/ForgivingQA")#where are all CATS related .h files
include_directories(${KITTY_PATH})
include_directories(${DREAM_PATH})
include_directories(${FORGIVEN_PATH})

add_executable(systematicsSigma ${PROJECT_SOURCE_DIR}/EvalSigmaSystematics.C)
target_link_libraries(systematicsSigma tidyCATS DreamFunctions ForgivingFnct)

add_executable(AddFiles ${PROJECT_SOURCE_DIR}/AddOutputs.C)
target_link_libraries(AddFiles tidyCATS DreamFunctions ForgivingFnct)

add_executable(executeSysteamticsPP ${PROJECT_SOURCE_DIR}/EvalDreamSyspp.C)
target_link_libraries(executeSysteamticsPP tidyCATS DreamFunctions ForgivingFnct)

add_executable(executeSysteamticsPXi ${PROJECT_SOURCE_DIR}/EvalDreamSyspXi.C)
target_link_libraries(executeSysteamticsPXi tidyCATS DreamFunctions ForgivingFnct)