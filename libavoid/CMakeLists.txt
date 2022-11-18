cmake_minimum_required (VERSION 3.0)
project (libavoid)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
add_definitions(-DLIBAVOID_EXPORTS)
set( headers
	${CMAKE_CURRENT_SOURCE_DIR}/assertions.h
		${CMAKE_CURRENT_SOURCE_DIR}/connector.h
		${CMAKE_CURRENT_SOURCE_DIR}/connectionpin.h
		${CMAKE_CURRENT_SOURCE_DIR}/connend.h
		${CMAKE_CURRENT_SOURCE_DIR}/debug.h
                        ${CMAKE_CURRENT_SOURCE_DIR}/dllexport.h
		${CMAKE_CURRENT_SOURCE_DIR}/geometry.h 
		${CMAKE_CURRENT_SOURCE_DIR}/geomtypes.h 
		${CMAKE_CURRENT_SOURCE_DIR}/graph.h 
		${CMAKE_CURRENT_SOURCE_DIR}/junction.h 
		${CMAKE_CURRENT_SOURCE_DIR}/libavoid.h 
		${CMAKE_CURRENT_SOURCE_DIR}/makepath.h 
		${CMAKE_CURRENT_SOURCE_DIR}/obstacle.h 
		${CMAKE_CURRENT_SOURCE_DIR}/orthogonal.h 
		${CMAKE_CURRENT_SOURCE_DIR}/router.h 
		${CMAKE_CURRENT_SOURCE_DIR}/shape.h 
		${CMAKE_CURRENT_SOURCE_DIR}/timer.h 
		${CMAKE_CURRENT_SOURCE_DIR}/vertices.h 
		${CMAKE_CURRENT_SOURCE_DIR}/viscluster.h 
		${CMAKE_CURRENT_SOURCE_DIR}/visibility.h 
		${CMAKE_CURRENT_SOURCE_DIR}/hyperedge.h 
		${CMAKE_CURRENT_SOURCE_DIR}/hyperedgeimprover.h 
		${CMAKE_CURRENT_SOURCE_DIR}/mtst.h 
		${CMAKE_CURRENT_SOURCE_DIR}/hyperedgetree.h 
		${CMAKE_CURRENT_SOURCE_DIR}/scanline.h 
		${CMAKE_CURRENT_SOURCE_DIR}/actioninfo.h 
		${CMAKE_CURRENT_SOURCE_DIR}/vpsc.h 
		${CMAKE_CURRENT_SOURCE_DIR}/debughandler.h
)

set( sources
			${CMAKE_CURRENT_SOURCE_DIR}/connectionpin.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/connector.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/connend.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/geometry.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/geomtypes.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/graph.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/junction.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/makepath.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/obstacle.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/orthogonal.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/router.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/shape.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/timer.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/vertices.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/viscluster.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/visibility.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/vpsc.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/hyperedge.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/hyperedgeimprover.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/mtst.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/hyperedgetree.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/scanline.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/actioninfo.cpp
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


