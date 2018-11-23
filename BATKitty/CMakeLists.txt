#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(BATKitty)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E env bat-config --libs OUTPUT_VARIABLE BAT_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CMAKE_COMMAND} -E env bat-config --cflags OUTPUT_VARIABLE BAT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CMAKE_COMMAND} -E env bat-config --libdir OUTPUT_VARIABLE BAT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CMAKE_COMMAND} -E env bat-config --incdir OUTPUT_VARIABLE BAT_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BAT_CFLAGS}")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS} ${BAT_INCDIR})
add_definitions(${ROOT_CXX_FLAGS})


set(SRCS
  ${PROJECT_SOURCE_DIR}/MyMod.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/MyMod.h
  )

add_library(myBAT STATIC ${SRCS})
target_link_libraries(myBAT -L${BAT_LIBDIR} ${BAT_LIB})


add_executable(runTut ${PROJECT_SOURCE_DIR}/Scripts/runMyTut.cxx)
target_link_libraries(runTut myBAT -L${BAT_LIBDIR} ${BAT_LIB})
