cmake_minimum_required (VERSION 3.0)
project (libvpsc)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

set( headers
	${CMAKE_CURRENT_SOURCE_DIR}/solve_VPSC.h
	${CMAKE_CURRENT_SOURCE_DIR}/block.h
	${CMAKE_CURRENT_SOURCE_DIR}/constraint.h
	${CMAKE_CURRENT_SOURCE_DIR}/exceptions.h
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.h
	${CMAKE_CURRENT_SOURCE_DIR}/variable.h
	${CMAKE_CURRENT_SOURCE_DIR}/assertions.h
)

set( sources
	${CMAKE_CURRENT_SOURCE_DIR}/block.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/blocks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constraint.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/solve_VPSC.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/variable.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cbuffer.cpp
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


