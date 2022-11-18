cmake_minimum_required (VERSION 3.0)
project (libtopology)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

set( headers
${CMAKE_CURRENT_SOURCE_DIR}/topology_constraints.h
${CMAKE_CURRENT_SOURCE_DIR}/topology_graph.h
${CMAKE_CURRENT_SOURCE_DIR}/cola_topology_addon.h
${CMAKE_CURRENT_SOURCE_DIR}/orthogonal_topology.h
${CMAKE_CURRENT_SOURCE_DIR}/topology_log.h
${CMAKE_CURRENT_SOURCE_DIR}/util.h
)

set( sources
${CMAKE_CURRENT_SOURCE_DIR}/topology_graph.cpp
${CMAKE_CURRENT_SOURCE_DIR}/topology_constraints.cpp
${CMAKE_CURRENT_SOURCE_DIR}/compute_forces.cpp
${CMAKE_CURRENT_SOURCE_DIR}/topology_constraints_constructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/resize.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cola_topology_addon.cpp
${CMAKE_CURRENT_SOURCE_DIR}/orthogonal_topology.cpp
${CMAKE_CURRENT_SOURCE_DIR}/topology_constraints.h
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


