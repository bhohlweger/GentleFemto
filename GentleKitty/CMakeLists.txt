CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

project(GentleKitty)
# SET PATHS #
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(GSL_INCLUDE "/usr/local/Cellar/gsl/2.6/include")#where are all GSL related .h files
  SET(GSL_LIB "/usr/local/Cellar/gsl/2.6/lib")#where are the GSL .a and .so files
else()
  SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
  SET(GSL_LIB "/usr/lib/x86_64-linux-gnu/")#where are the GSL .a and .so files
endif()
SET(CATS_INCLUDE "$ENV{CATS}/include")#where are all CATS related .h files
SET(CATS_LIB "$ENV{CATS}/lib")#where are the CATS related .a files
SET(DREAM_PATH "${CMAKE_SOURCE_DIR}/DreamFunction")#where are all CATS related .h files
SET(DREAM_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/../DreamFunction")
SET(FORGIVEN_PATH "${CMAKE_SOURCE_DIR}/ForgivingQA")#where are all CATS related .h files
SET(GAMI_PATH "${CMAKE_SOURCE_DIR}/FemtoGami")#where are all CATS related .h files
SET(ROOFOLD "$ENV{ROOUNFOLD_ROOT}")#where are all CATS related .h files

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${binFolder}/GentleKitty
    )

# DO NOT CHANGE THE REST #

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG Minuit)
include(${ROOT_USE_FILE})
include_directories("${ROOFOLD}")
include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
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
include_directories(${GAMI_PATH})

set(SRCS
  ${PROJECT_SOURCE_DIR}/TidyCats.cxx
  ${PROJECT_SOURCE_DIR}/CATSInput.cxx
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.cxx
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.cxx
  ${PROJECT_SOURCE_DIR}/PlayWithCats.cxx
  ${PROJECT_SOURCE_DIR}/SideBandFit.cxx
  ${PROJECT_SOURCE_DIR}/SidebandSigma.cxx
  ${PROJECT_SOURCE_DIR}/VariationAnalysis.cxx
  ${PROJECT_SOURCE_DIR}/VariationmTAnalysis.cxx
  ${PROJECT_SOURCE_DIR}/VariationAnalysispAp.cxx
  ${PROJECT_SOURCE_DIR}/VariationAnalysispAL.cxx
  ${PROJECT_SOURCE_DIR}/VariationAnalysisLAL.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/TidyCats.h
  ${PROJECT_SOURCE_DIR}/CATSInput.h
  ${PROJECT_SOURCE_DIR}/CATSInputSigma0.h
  ${PROJECT_SOURCE_DIR}/CATSLambdaParam.h
  ${PROJECT_SOURCE_DIR}/PlayWithCats.h
  ${PROJECT_SOURCE_DIR}/SideBandFit.h
  ${PROJECT_SOURCE_DIR}/SidebandSigma.h
  ${PROJECT_SOURCE_DIR}/VariationAnalysis.h
  ${PROJECT_SOURCE_DIR}/VariationmTAnalysis.h
  ${PROJECT_SOURCE_DIR}/VariationAnalysispAp.h
  ${PROJECT_SOURCE_DIR}/VariationAnalysispAL.h
  ${PROJECT_SOURCE_DIR}/VariationAnalysisLAL.h
  )

add_library(tidyCATS STATIC ${SRCS})
set_target_properties(tidyCATS PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(tidyCATS ForgivingFnct femtoGami -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev ${ROOT_LIBRARIES} -L${GSL_LIB} -lgsl -lgslcblas)

add_executable(POnly ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnly.C)
target_link_libraries(POnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(ppSource ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/ppSourceCFmaker.C)
target_link_libraries(ppSource tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PLOnly ${PROJECT_SOURCE_DIR}/Scripts/pLAnalysis/PLOnly.C)
target_link_libraries(PLOnly tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlyOnTheGrid ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnlyOnTheGrid.C)
target_link_libraries(POnlyOnTheGrid tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlyOnTheGridforCoalescence ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysisForCoalescence/POnlyOnTheGridCoalescence.C)
target_link_libraries(POnlyOnTheGridforCoalescence tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

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

add_executable(ppSystematicsCoalescence ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysisForCoalescence/systematics_ppBhawani.C)
target_link_libraries(ppSystematicsCoalescence tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DrawCurves ${PROJECT_SOURCE_DIR}/Scripts/DrawCurves.C)
target_link_libraries(DrawCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(EvalpXiCurves ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/EvalpXiCurves.C)
target_link_libraries(EvalpXiCurves tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(dumpSigma0 ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/DumpCF.C)
target_link_libraries(dumpSigma0 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PlayCats ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/PlayCats.C)
target_link_libraries(PlayCats tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pSigma0 ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/pSigma0_single.C)
target_link_libraries(pSigma0 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(coupledChannel ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/coupledChannel.C)
target_link_libraries(coupledChannel tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(someLednicky ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/drawLednicky.C)
target_link_libraries(someLednicky tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PlotSigmaModels ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/plotSigmaModels.C)
target_link_libraries(PlotSigmaModels tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SigmaSource ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/FemtoSource.C)
target_link_libraries(SigmaSource tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SigmaRadiusCF ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/radiusCF.C)
target_link_libraries(SigmaRadiusCF tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PlotSigmaSmear ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/pSigma0smear.C)
target_link_libraries(PlotSigmaSmear tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PlotSigmaCoupled ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/pSigma0coupled.C)
target_link_libraries(PlotSigmaCoupled tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(sigmaModels ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/drawSigmaModelPlot.C)
target_link_libraries(sigmaModels tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SigmaSmearingMatrix ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/smearingMatrixSidebands.C)
target_link_libraries(SigmaSmearingMatrix tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(XiSmearingMatrix ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/smearingMatrixSignal.C)
target_link_libraries(XiSmearingMatrix tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SigmaSmearingSherlock ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/smearingLambdaSidebands.C)
target_link_libraries(SigmaSmearingSherlock tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(XiSmearingSherlock ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/smearingLambdaSidebands.C)
target_link_libraries(XiSmearingSherlock tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(Sasaki ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/SasakiPotential.C)
target_link_libraries(Sasaki tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlymT ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnlymT.C)
target_link_libraries(POnlymT tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(CatsPlay ${PROJECT_SOURCE_DIR}/Scripts/CatsPlay.C)
target_link_libraries(CatsPlay tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(PotentialPlay ${PROJECT_SOURCE_DIR}/Scripts/pXiAnalysis/PotentialPlay.C)
target_link_libraries(PotentialPlay tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(TestCats ${PROJECT_SOURCE_DIR}/Scripts/TestCats.C)
target_link_libraries(TestCats tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(TestCats2 ${PROJECT_SOURCE_DIR}/Scripts/TestCats2.C)
target_link_libraries(TestCats2 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(TestCats3 ${PROJECT_SOURCE_DIR}/Scripts/TestCats3.C)
target_link_libraries(TestCats3 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pLPlot ${PROJECT_SOURCE_DIR}/Scripts/pLAnalysis/MakeAPlot.C)
target_link_libraries(pLPlot tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(LambdaPars ${PROJECT_SOURCE_DIR}/Scripts/LambdaPars.C)
target_link_libraries(LambdaPars tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SourceMtMaker ${PROJECT_SOURCE_DIR}/Scripts/pSigma0Analysis/SourceMtMaker.C)
target_link_libraries(SourceMtMaker tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(POnlyExtractor ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/POnlyExtractor.C)
target_link_libraries(POnlyExtractor tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(mTAnalysis ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/mTExtractor.C)
target_link_libraries(mTAnalysis tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(mTAnalysisCoalescence ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysisForCoalescence/mTExtractorCoalescence.C)
target_link_libraries(mTAnalysisCoalescence tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(mTPLAnalysis ${PROJECT_SOURCE_DIR}/Scripts/pLAnalysis/mTExtractorAddPot.C)
target_link_libraries(mTPLAnalysis tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(runPhi ${PROJECT_SOURCE_DIR}/Scripts/pPhi.C)
target_link_libraries(runPhi tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(radPlotter ${PROJECT_SOURCE_DIR}/Scripts/RadiusDrawer.C)
target_link_libraries(radPlotter tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)
add_executable(radPlotterV2 ${PROJECT_SOURCE_DIR}/Scripts/RadiusDrawerPotSep.C)
target_link_libraries(radPlotterV2 tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(b2Plots ${PROJECT_SOURCE_DIR}/Scripts/PlotsB2Coalescence.C)
target_link_libraries(b2Plots tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(radPlotterForCoalescence ${PROJECT_SOURCE_DIR}/Scripts/RadiusDrawerForCoalescence.C)
target_link_libraries(radPlotterForCoalescence tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)
add_executable(plotCFforCoalescence ${PROJECT_SOURCE_DIR}/Scripts/PlotCFsForCoalescence.C)
target_link_libraries(plotCFforCoalescence tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(lambdaPlotter ${PROJECT_SOURCE_DIR}/Scripts/ppAnalysis/LambdaExtractor.C)
target_link_libraries(lambdaPlotter tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SysVarExtractorpAp ${PROJECT_SOURCE_DIR}/Scripts/BBarAnalysis/SysVarExtractorpAp.C)
target_link_libraries(SysVarExtractorpAp tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SysVarExtractorpAL ${PROJECT_SOURCE_DIR}/Scripts/BBarAnalysis/SysVarExtractorpAL.C)
target_link_libraries(SysVarExtractorpAL tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(SysVarExtractorLAL ${PROJECT_SOURCE_DIR}/Scripts/BBarAnalysis/SysVarExtractorLAL.C)
target_link_libraries(SysVarExtractorLAL tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(LednickyPlay ${PROJECT_SOURCE_DIR}/Scripts/LednickyPlay.C)
target_link_libraries(LednickyPlay tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DmesonPlayground ${PROJECT_SOURCE_DIR}/Scripts/pDmeson/DmesonPlayground.C)
target_link_libraries(DmesonPlayground tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DmesonHaidenbauer ${PROJECT_SOURCE_DIR}/Scripts/pDmeson/DmesonHaidenbauer.C)
target_link_libraries(DmesonHaidenbauer tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(pPhiPotentials ${PROJECT_SOURCE_DIR}/Scripts/pPhiPotentials.C)
target_link_libraries(pPhiPotentials tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(fitDminus ${PROJECT_SOURCE_DIR}/Scripts/pDmeson/fitDminus.C)
target_link_libraries(fitDminus tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(correctDminus ${PROJECT_SOURCE_DIR}/Scripts/pDmeson/correctDminus.C)
target_link_libraries(correctDminus tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)

add_executable(DmesonModels ${PROJECT_SOURCE_DIR}/Scripts/pDmeson/DmesonModels.C)
target_link_libraries(DmesonModels tidyCATS -L${DREAM_LIB_PATH} DreamFunctions)
