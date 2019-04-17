CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

project(GentleKitty)
# SET PATHS #
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(GSL_INCLUDE "/usr/local/Cellar/gsl/2.5/include")#where are all GSL related .h files
  SET(GSL_LIB "/usr/local/Cellar/gsl/2.5/lib")#where are the GSL .a and .so files
else()
  SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
  SET(GSL_LIB "/usr/lib/x86_64-linux-gnu/")#where are the GSL .a and .so files
endif()
SET(CATS_INCLUDE "$ENV{CATS}/include")#where are all CATS related .h files
SET(CATS_LIB "$ENV{CATS}/lib")#where are the CATS related .a files
SET(DREAM_PATH "${CMAKE_SOURCE_DIR}/DreamFunction")#where are all CATS related .h files
SET(DREAM_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/../DreamFunction")
SET(FORGIVEN_PATH "${CMAKE_SOURCE_DIR}/ForgivingQA")#where are all CATS related .h files

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${binFolder}/GentleKitty
    )

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
include_directories(${GSL_INCLUDE}/gsl)
include_directories(${DREAM_PATH})
include_directories(${FORGIVEN_PATH})

set(SRCS
  ${PROJECT_SOURCE_DIR}/TidyCats.cxx
  ${PROJECT_SOURCE_DIR}/CATSInput.cxx
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.cxx
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.cxx
  ${PROJECT_SOURCE_DIR}/PlayWithCats.cxx
  ${PROJECT_SOURCE_DIR}/SideBandFit.cxx
  ${PROJECT_SOURCE_DIR}/SidebandSigma.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/TidyCats.h
  ${PROJECT_SOURCE_DIR}/CATSInput.h
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.h
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.h
  ${PROJECT_SOURCE_DIR}/PlayWithCats.h
  ${PROJECT_SOURCE_DIR}/SideBandFit.h
  ${PROJECT_SOURCE_DIR}/SidebandSigma.h
  )

add_library(tidyCATS STATIC ${SRCS})
set_target_properties(tidyCATS PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(tidyCATS ForgivingFnct -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(POnly ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnly.C)
target_link_libraries(POnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlyOnTheGrid ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnlyOnTheGrid.C)
target_link_libraries(POnlyOnTheGrid tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(ppExtractSource ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/pp_mTSourceExtractor.C)
target_link_libraries(ppExtractSource tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(XiOnly ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/XiOnly.C)
target_link_libraries(XiOnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(XiOnlyOnTheGrid ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/XiOnlyOnTheGrid.C)
target_link_libraries(XiOnlyOnTheGrid tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SideBand ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/SideBand.C)
target_link_libraries(SideBand tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(ppSystematics ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/systematics_pp.C)
target_link_libraries(ppSystematics tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DrawCurves ${PROJECT_SOURCE_DIR}/Scripts/DrawCurves.C)
target_link_libraries(DrawCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(EvalpXiCurves ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/EvalpXiCurves.C)
target_link_libraries(EvalpXiCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pSigma0 ${PROJECT_SOURCE_DIR}/Scripts/pSigma0_single.C)
target_link_libraries(pSigma0 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pSigma0exclusion ${PROJECT_SOURCE_DIR}/Scripts/pSigma0exclusion.C)
target_link_libraries(pSigma0exclusion tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PlotSigmaExclusion ${PROJECT_SOURCE_DIR}/Scripts/plotSigmaExclusion.C)
target_link_libraries(PlotSigmaExclusion tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

#add_executable(PlotSigmaModels ${PROJECT_SOURCE_DIR}/Scripts/plotSigmaModels.C)
#target_link_libraries(PlotSigmaModels tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlymT ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnlymT.C)
target_link_libraries(POnlymT tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(CatsPlay ${PROJECT_SOURCE_DIR}/Scripts/CatsPlay.C)
target_link_libraries(CatsPlay tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PotentialPlay ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/PotentialPlay.C)
target_link_libraries(PotentialPlay tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(TestCats ${PROJECT_SOURCE_DIR}/Scripts/TestCats.C)
target_link_libraries(TestCats tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(LambdaPars ${PROJECT_SOURCE_DIR}/Scripts/LambdaPars.C)
target_link_libraries(LambdaPars tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_subdirectory(Test)
