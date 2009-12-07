# Project name affects name of solution file.
project (arageli CXX)
cmake_minimum_required(VERSION 2.6)

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory (src/arageli)
add_subdirectory (tool/ts/src/ts)

