CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(FemtoGami)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# SET PATHS #
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

SET(CATS_INCLUDE "$ENV{CATS}/include")#where are all CATS related .h files
SET(CATS_LIB "$ENV{CATS}/lib")#where are the CATS related .a files

SET(ROOFOLD "$ENV{ROOUNFOLD_ROOT}")#where are all CATS related .h files

SET(ROOFOLD_LIB "$ENV{ROOUNFOLD_ROOT}")#where are the CATS related .a files

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${binFolder}/FemtoGami)
include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

SET(KITTY_PATH "${CMAKE_SOURCE_DIR}/GentleKitty")#where are all CATS related .h files
SET(DREAM_PATH "${CMAKE_SOURCE_DIR}/DreamFunction")#where are all CATS related .h files
SET(FORGIVEN_PATH "${CMAKE_SOURCE_DIR}/ForgivingQA")#where are all CATS related .h files
SET(SYSTEMATIC_PATH "${CMAKE_SOURCE_DIR}/Systematics")#where are all CATS related .h files

include_directories(${KITTY_PATH})
include_directories(${DREAM_PATH})
include_directories(${FORGIVEN_PATH})
include_directories(${SYSTEMATIC_PATH})
include_directories(${CATS_INCLUDE})
include_directories("${ROOFOLD}")
include_directories(${PROJECT_SOURCE_DIR})

set(SRCS
  ${PROJECT_SOURCE_DIR}/LambdaGami.cxx
  ${PROJECT_SOURCE_DIR}/MomentumGami.cxx
  ${PROJECT_SOURCE_DIR}/AnalyseProXi.cxx
  )

set(HEADERS
  ${PROJECT_SOURCE_DIR}/LambdaGami.h
  ${PROJECT_SOURCE_DIR}/MomentumGami.h
  ${PROJECT_SOURCE_DIR}/AnalyseProXi.h
  )
  
add_library(femtoGami STATIC ${SRCS})
set_target_properties(femtoGami PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(femtoGami tidyCATS DreamFunctions ForgivingFnct -L${ROOFOLD} RooUnfold)

add_executable(XiGami ${PROJECT_SOURCE_DIR}/Scripts/XiGami.C)
target_link_libraries(XiGami femtoGami)

add_executable(CloseMe ${PROJECT_SOURCE_DIR}/Scripts/XiClosure.C)
target_link_libraries(CloseMe femtoGami)

add_executable(XiDraw ${PROJECT_SOURCE_DIR}/Scripts/XiDrawWithError.C)
target_link_libraries(XiDraw femtoGami)

add_executable(XiComp ${PROJECT_SOURCE_DIR}/Scripts/XiComparison.C)
target_link_libraries(XiComp femtoGami)


add_executable(OmDraw ${PROJECT_SOURCE_DIR}/Scripts/OmegaDrawWithError.C)
target_link_libraries(OmDraw femtoGami)
