cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(binFolder "${CMAKE_CURRENT_BINARY_DIR}/../install")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${binFolder}/lib)
project(GentleFemto CXX C)
enable_testing()

add_subdirectory(DreamFunction)
add_subdirectory(GentleKitty)
add_subdirectory(Systematics)
add_subdirectory(ForgivingQA )
add_subdirectory(FemtoGami )

