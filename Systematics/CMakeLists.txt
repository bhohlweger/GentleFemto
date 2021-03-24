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
SET(GAMI_PATH "${CMAKE_SOURCE_DIR}/FemtoGami")#where are all CATS related .h files
SET(ROOFOLD "$ENV{ROOUNFOLD_ROOT}")#where are all CATS related .h files


include_directories(${KITTY_PATH})
include_directories(${DREAM_PATH})
include_directories(${FORGIVEN_PATH})
include_directories(${GAMI_PATH})
include_directories("${ROOFOLD}")

add_executable(systematicsSigma ${PROJECT_SOURCE_DIR}/EvalSigmaSystematics.C)
target_link_libraries(systematicsSigma tidyCATS DreamFunctions ForgivingFnct)

add_executable(sigmaSidebands ${PROJECT_SOURCE_DIR}/EvalSidebands.C)
target_link_libraries(sigmaSidebands tidyCATS DreamFunctions ForgivingFnct)

add_executable(AddFiles ${PROJECT_SOURCE_DIR}/AddOutputs.C)
target_link_libraries(AddFiles tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicsPP ${PROJECT_SOURCE_DIR}/EvalDreamSyspp.C)
target_link_libraries(systematicsPP tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicsPXi ${PROJECT_SOURCE_DIR}/EvalDreamSyspXi.C)
target_link_libraries(systematicsPXi tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicsPL ${PROJECT_SOURCE_DIR}/EvalDreamSyspL.C)
target_link_libraries(systematicsPL tidyCATS DreamFunctions ForgivingFnct)

add_executable(ppmTVariations ${PROJECT_SOURCE_DIR}/GetmTVariations.C)
target_link_libraries(ppmTVariations tidyCATS DreamFunctions ForgivingFnct)

add_executable(pLmTVariations ${PROJECT_SOURCE_DIR}/GetmTVariationspL.C)
target_link_libraries(pLmTVariations tidyCATS DreamFunctions ForgivingFnct)

add_executable(ppmTVariationsForCoalescence ${PROJECT_SOURCE_DIR}/GetmTVariationsCoalescence.C)
target_link_libraries(ppmTVariationsForCoalescence tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicspAp ${PROJECT_SOURCE_DIR}/EvalDreamSyspAp.C)
target_link_libraries(systematicspAp tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicsLAL ${PROJECT_SOURCE_DIR}/EvalDreamSysLAL.C)
target_link_libraries(systematicsLAL tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicspAL ${PROJECT_SOURCE_DIR}/EvalDreamSyspAL.C)
target_link_libraries(systematicspAL tidyCATS DreamFunctions ForgivingFnct)

add_executable(pALmTVariations ${PROJECT_SOURCE_DIR}/GetmTVariationspAL.C)
target_link_libraries(pALmTVariations tidyCATS DreamFunctions ForgivingFnct)

add_executable(pALmTSystematics ${PROJECT_SOURCE_DIR}/GetmTSystematicspAL.C)
target_link_libraries(pALmTSystematics tidyCATS DreamFunctions ForgivingFnct)

add_executable(pALmTSystematics_Update ${PROJECT_SOURCE_DIR}/GetmTSystematicspAL_Update.C)
target_link_libraries(pALmTSystematics_Update tidyCATS DreamFunctions ForgivingFnct)

add_executable(LALmTSystematics_Update ${PROJECT_SOURCE_DIR}/GetmTSystematicsLAL_Update.C)
target_link_libraries(LALmTSystematics_Update tidyCATS DreamFunctions ForgivingFnct)

add_executable(LALmTVariations ${PROJECT_SOURCE_DIR}/GetmTVariationsLAL.C)
target_link_libraries(LALmTVariations tidyCATS DreamFunctions ForgivingFnct)

add_executable(LALmTSystematics ${PROJECT_SOURCE_DIR}/GetmTSystematicsLAL.C)
target_link_libraries(LALmTSystematics tidyCATS DreamFunctions ForgivingFnct)

add_executable(pALmTAncestors ${PROJECT_SOURCE_DIR}/GetmTVariationspALAncestors.C)
target_link_libraries(pALmTAncestors tidyCATS DreamFunctions ForgivingFnct)

add_executable(LALmTAncestors ${PROJECT_SOURCE_DIR}/GetmTVariationsLALAncestors.C)
target_link_libraries(LALmTAncestors tidyCATS DreamFunctions ForgivingFnct)

add_executable(pApmTAncestors ${PROJECT_SOURCE_DIR}/GetmTVariationspApAncestors.C)
target_link_libraries(pApmTAncestors tidyCATS DreamFunctions ForgivingFnct)

add_executable(systematicsD ${PROJECT_SOURCE_DIR}/EvalDmesonSystematics.C)
target_link_libraries(systematicsD tidyCATS DreamFunctions ForgivingFnct)
