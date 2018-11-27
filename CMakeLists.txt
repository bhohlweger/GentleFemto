cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(binFolder "${CMAKE_CURRENT_BINARY_DIR}/../install/")

project(GentleFemto CXX C)
enable_testing()

add_subdirectory(DreamFunction)
add_subdirectory(GentleKitty)

