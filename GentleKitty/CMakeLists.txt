CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(GentleKitty)
# SET PATHS #
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(ROOTSYS "/home/hohlweger/alice/sw/ubuntu1604_x86-64/ROOT/v6-10-08-2")
SET(ROOT_INCLUDE "/home/hohlweger/alice/sw/ubuntu1604_x86-64/ROOT/v6-10-08-2/include")#where are all ROOT related .h files
SET(ROOT_LIB "/home/hohlweger/alice/sw/ubuntu1604_x86-64/ROOT/v6-10-08-2/lib")#where are the ROOT related .a and .so file
SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
SET(GSL_LIB "/usr/lib/x86_64-linux-gnu/")#where are the GSL .a and .so files
SET(CATS_INCLUDE "/home/hohlweger/CATS/install/include")#where are all CATS related .h files
SET(CATS_LIB "/home/hohlweger/CATS/install/lib")#where are the CATS related .a files
#add_executable(EXECUTABLE ForBernie.cpp)
# DO NOT CHANGE THE REST #

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${ROOT_INCLUDE})

set(SRCS
  ${PROJECT_SOURCE_DIR}/TidyCats.cxx
  ${PROJECT_SOURCE_DIR}/CATSInput.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/TidyCats.h
  ${PROJECT_SOURCE_DIR}/CATSInput.h
  )

add_library(tidyCATS STATIC ${PROJECT_SOURCE_DIR}/TidyCats.cxx ${PROJECT_SOURCE_DIR}/CATSInput.cxx)
target_link_libraries(tidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(EXECUTABLE ${PROJECT_SOURCE_DIR}/Scripts/ForBernie.C)
target_link_libraries(EXECUTABLE tidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(POnly ${PROJECT_SOURCE_DIR}/Scripts/POnly.C)
target_link_libraries(POnly tidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)


#ROOT_GENERATE_DICTIONARY(G__TidyCATS ${HEADERS} LINKDEF ${PROJECT_SOURCE_DIR}/LinkDef.h)
#add_library(TidyCATS SHARED ${SRCS} G__TidyCATS)
#target_link_libraries(TidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${GSL_LIB} ${ROOT_CXX_LIB_FLAGS} libgsl.a libgsl.so)
#
#add_executable(EXECUTABLE ${PROJECT_SOURCE_DIR}/Scripts/ForBernie.C)
#target_link_libraries(EXECUTABLE -L -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${GSL_LIB} ${ROOT_CXX_LIB_FLAGS} libgsl.a libgsl.so)
#

