cmake_minimum_required (VERSION 3.0)
project (libcola)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

set( headers
	${CMAKE_CURRENT_SOURCE_DIR}/cola.h
	${CMAKE_CURRENT_SOURCE_DIR}/cluster.h
	${CMAKE_CURRENT_SOURCE_DIR}/commondefs.h
	${CMAKE_CURRENT_SOURCE_DIR}/compound_constraints.h
	${CMAKE_CURRENT_SOURCE_DIR}/pseudorandom.h
	${CMAKE_CURRENT_SOURCE_DIR}/exceptions.h
	${CMAKE_CURRENT_SOURCE_DIR}/gradient_projection.h
	${CMAKE_CURRENT_SOURCE_DIR}/sparse_matrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/straightener.h 
	${CMAKE_CURRENT_SOURCE_DIR}/output_svg.h 
	${CMAKE_CURRENT_SOURCE_DIR}/unused.h 
	${CMAKE_CURRENT_SOURCE_DIR}/cc_clustercontainmentconstraints.h
	${CMAKE_CURRENT_SOURCE_DIR}/cc_nonoverlapconstraints.h
	${CMAKE_CURRENT_SOURCE_DIR}/box.h
	${CMAKE_CURRENT_SOURCE_DIR}/shapepair.h
)

set( sources
	${CMAKE_CURRENT_SOURCE_DIR}/cola.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/colafd.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/conjugate_gradient.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/gradient_projection.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/straightener.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/connected_components.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/convex_hull.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cluster.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/compound_constraints.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pseudorandom.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/output_svg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cc_clustercontainmentconstraints.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cc_nonoverlapconstraints.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/box.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/shapepair.cpp
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


