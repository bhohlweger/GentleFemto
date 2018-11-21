CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

project(GentleKitty)
# SET PATHS #
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
SET(GSL_LIB "/usr/lib/x86_64-linux-gnu/")#where are the GSL .a and .so files
SET(CATS_INCLUDE "$ENV{CATS}/include")#where are all CATS related .h files
SET(CATS_LIB "$ENV{CATS}/lib")#where are the CATS related .a files
SET(DREAM_PATH "${CMAKE_SOURCE_DIR}/DreamFunction")#where are all CATS related .h files
SET(DREAM_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/../DreamFunction")
#add_executable(EXECUTABLE ForBernie.cpp)
# DO NOT CHANGE THE REST #

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG Minuit)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

#find_library(DREAM_LIB libDreamFunctions.so HINTS ${DREAM_LIB_PATH})

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${DREAM_PATH})

set(SRCS
  ${PROJECT_SOURCE_DIR}/TidyCats.cxx
  ${PROJECT_SOURCE_DIR}/CATSInput.cxx
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.cxx
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.cxx
  ${PROJECT_SOURCE_DIR}/SideBandFit.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/TidyCats.h
  ${PROJECT_SOURCE_DIR}/CATSInput.h
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.h
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.h
  ${PROJECT_SOURCE_DIR}/SideBandFit.h
  )

add_library(tidyCATS STATIC ${SRCS})
target_link_libraries(tidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(POnly ${PROJECT_SOURCE_DIR}/Scripts/POnly.C)
target_link_libraries(POnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(XiOnly ${PROJECT_SOURCE_DIR}/Scripts/XiOnly.C)
target_link_libraries(XiOnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(SideBand ${PROJECT_SOURCE_DIR}/Scripts/SideBand.C)
target_link_libraries(SideBand tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(ppSystematics ${PROJECT_SOURCE_DIR}/Scripts/systematics_pp.C)
target_link_libraries(ppSystematics tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DrawCurves ${PROJECT_SOURCE_DIR}/Scripts/DrawCurves.C)
target_link_libraries(DrawCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(EvalpXiCurves ${PROJECT_SOURCE_DIR}/Scripts/EvalpXiCurves.C)
target_link_libraries(EvalpXiCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pSigma0 ${PROJECT_SOURCE_DIR}/Scripts/pSigma0.C)
target_link_libraries(pSigma0 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)


#add_executable(EXECUTABLE ${PROJECT_SOURCE_DIR}/Scripts/ForBernie.C)
#target_link_libraries(EXECUTABLE tidyCATS DREAM_LIB -DreamFunctions -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)
#ROOT_GENERATE_DICTIONARY(G__TidyCATS ${HEADERS} LINKDEF ${PROJECT_SOURCE_DIR}/LinkDef.h)
#add_library(TidyCATS SHARED ${SRCS} G__TidyCATS)
#target_link_libraries(TidyCATS -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${GSL_LIB} ${ROOT_CXX_LIB_FLAGS} libgsl.a libgsl.so)
#
#add_executable(EXECUTABLE ${PROJECT_SOURCE_DIR}/Scripts/ForBernie.C)
#target_link_libraries(EXECUTABLE -L -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${GSL_LIB} ${ROOT_CXX_LIB_FLAGS} libgsl.a libgsl.so)
#

add_subdirectory(Test)
