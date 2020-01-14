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

SET(KITTY_PATH "${CMAKE_SOURCE_DIR}/GentleKitty")#where are all CATS related .h files
SET(FORGIVEN_PATH "${CMAKE_SOURCE_DIR}/ForgivingQA")#where are all CATS related .h files
SET(GAMI_PATH "${CMAKE_SOURCE_DIR}/FemtoGami")#where are all CATS related .h files
SET(SYSTEMATIC_PATH "${CMAKE_SOURCE_DIR}/Systematics")#where are all CATS related .h files
SET(ROOFOLD_INCLUDE "$ENV{ROOUNFOLD_ROOT}/include")#where are all CATS related .h files
SET(ROOFOLD_LIB "$ENV{ROOUNFOLD_ROOT}/lib")#where are the CATS related .a files

include_directories(${KITTY_PATH})
include_directories(${FORGIVEN_PATH})
include_directories(${GAMI_PATH})
include_directories(${SYSTEMATIC_PATH})
include_directories(${CATS_INCLUDE})
include_directories("${ROOFOLD_INCLUDE}")
include_directories(${PROJECT_SOURCE_DIR})

set(SRCS
  ${PROJECT_SOURCE_DIR}/DreamCF.cxx
  ${PROJECT_SOURCE_DIR}/DreamData.cxx
  ${PROJECT_SOURCE_DIR}/DreamDist.cxx
  ${PROJECT_SOURCE_DIR}/DreamKayTee.cxx
  ${PROJECT_SOURCE_DIR}/DreamdEtadPhi.cxx
  ${PROJECT_SOURCE_DIR}/DreamPair.cxx
  ${PROJECT_SOURCE_DIR}/DreamPlot.cxx
  ${PROJECT_SOURCE_DIR}/DreamHEP.cxx
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
  ${PROJECT_SOURCE_DIR}/DreamHEP.h
  ${PROJECT_SOURCE_DIR}/DreamSystematics.h
  ${PROJECT_SOURCE_DIR}/ReadDreamFile.h
  ${PROJECT_SOURCE_DIR}/global.h
  )
add_library(DreamFunctions STATIC ${SRCS})
set_target_properties(DreamFunctions PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(DreamFunctions ${ROOT_LIBRARIES} femtoGami)

add_executable(executeDream ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFDream.C)
target_link_libraries(executeDream DreamFunctions)

add_executable(executeDreamBBar ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFDreamBBar.C)
target_link_libraries(executeDreamBBar DreamFunctions)

add_executable(executeEmmTee ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFmT.C)
target_link_libraries(executeEmmTee DreamFunctions)

add_executable(executeKayTee ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCFkT.C)
target_link_libraries(executeKayTee DreamFunctions)

add_executable(executePlot ${PROJECT_SOURCE_DIR}/Scripts/ExecutePlotCF.C)
target_link_libraries(executePlot DreamFunctions)

add_executable(executePlotSigma ${PROJECT_SOURCE_DIR}/Scripts/executePlotSigma.C)
target_link_libraries(executePlotSigma DreamFunctions)

add_executable(executePlotProton ${PROJECT_SOURCE_DIR}/Scripts/executePlotProton.C)
target_link_libraries(executePlotProton DreamFunctions)

add_executable(executePlotBBar ${PROJECT_SOURCE_DIR}/Scripts/ExecutePlotCFBBar.C)
target_link_libraries(executePlotBBar DreamFunctions)

add_executable(executeSigma ${PROJECT_SOURCE_DIR}/Scripts/ExecuteSigma.C)
target_link_libraries(executeSigma DreamFunctions)

add_executable(executeDPhiDEta ${PROJECT_SOURCE_DIR}/Scripts/ExecuteDPhiDEta.C)
target_link_libraries(executeDPhiDEta DreamFunctions)

add_executable(executeDPhiDEtaAtRad ${PROJECT_SOURCE_DIR}/Scripts/ExecuteDPhiDEtaAtRad.C)
target_link_libraries(executeDPhiDEtaAtRad DreamFunctions)

add_executable(executeHEP ${PROJECT_SOURCE_DIR}/Scripts/ConvertToHEP.C)
target_link_libraries(executeHEP DreamFunctions)

add_executable(executeSigmaHEP ${PROJECT_SOURCE_DIR}/Scripts/ConvertSigmaToHEP.C)
target_link_libraries(executeSigmaHEP DreamFunctions)

add_executable(executeOneFile ${PROJECT_SOURCE_DIR}/Scripts/MakeOneFile.C)
target_link_libraries(executeOneFile DreamFunctions)

add_executable(executeMiniJet ${PROJECT_SOURCE_DIR}/Scripts/ExecuteMiniJet.C)
target_link_libraries(executeMiniJet DreamFunctions)

add_executable(executeMiniJetEPOS ${PROJECT_SOURCE_DIR}/Scripts/ExecuteMiniJetEPOS.C)
target_link_libraries(executeMiniJetEPOS DreamFunctions)

add_executable(executeComparisonMC ${PROJECT_SOURCE_DIR}/Scripts/ExecuteComparisonMC.C)
target_link_libraries(executeComparisonMC DreamFunctions)

add_executable(executeCorrectBinCF ${PROJECT_SOURCE_DIR}/Scripts/ExecuteCorrectBinCF.C)
target_link_libraries(executeCorrectBinCF DreamFunctions)

add_executable(executeQADream ${PROJECT_SOURCE_DIR}/Scripts/ExecuteQADream.C)
target_link_libraries(executeQADream DreamFunctions)

add_executable(executePurityDream ${PROJECT_SOURCE_DIR}/Scripts/ExecutePurityDream.C)
target_link_libraries(executePurityDream DreamFunctions)

add_executable(executeDPhiDEtaBBar ${PROJECT_SOURCE_DIR}/Scripts/ExecuteDPhiDEtaBBar.C)
target_link_libraries(executeDPhiDEtaBBar DreamFunctions)

add_executable(executePhi ${PROJECT_SOURCE_DIR}/Scripts/ExecutePhi.C)
target_link_libraries(executePhi DreamFunctions)

add_executable(executePXi ${PROJECT_SOURCE_DIR}/Scripts/GetpXiVariations.C)
target_link_libraries(executePXi DreamFunctions)
