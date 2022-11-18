cmake_minimum_required (VERSION 3.0)
project (libproject)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

set( headers
${CMAKE_CURRENT_SOURCE_DIR}/util.h
${CMAKE_CURRENT_SOURCE_DIR}/project.h
)

set( sources
${CMAKE_CURRENT_SOURCE_DIR}/project.cpp
${CMAKE_CURRENT_SOURCE_DIR}/util.cpp
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


