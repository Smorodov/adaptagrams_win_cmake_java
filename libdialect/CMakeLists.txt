cmake_minimum_required (VERSION 3.0)
project (libdialect)


# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

set( headers
${CMAKE_CURRENT_SOURCE_DIR}/aca.h
${CMAKE_CURRENT_SOURCE_DIR}/chains.h
${CMAKE_CURRENT_SOURCE_DIR}/commontypes.h
${CMAKE_CURRENT_SOURCE_DIR}/constraints.h
${CMAKE_CURRENT_SOURCE_DIR}/expansion.h
${CMAKE_CURRENT_SOURCE_DIR}/faces.h
${CMAKE_CURRENT_SOURCE_DIR}/graphs.h
${CMAKE_CURRENT_SOURCE_DIR}/hola.h
${CMAKE_CURRENT_SOURCE_DIR}/io.h
${CMAKE_CURRENT_SOURCE_DIR}/logging.h
${CMAKE_CURRENT_SOURCE_DIR}/nearalign.h
${CMAKE_CURRENT_SOURCE_DIR}/nodeconfig.h
${CMAKE_CURRENT_SOURCE_DIR}/opts.h
${CMAKE_CURRENT_SOURCE_DIR}/ortho.h
${CMAKE_CURRENT_SOURCE_DIR}/peeling.h
${CMAKE_CURRENT_SOURCE_DIR}/planarise.h
${CMAKE_CURRENT_SOURCE_DIR}/quadaction.h
${CMAKE_CURRENT_SOURCE_DIR}/routing.h
${CMAKE_CURRENT_SOURCE_DIR}/treeplacement.h
${CMAKE_CURRENT_SOURCE_DIR}/trees.h
${CMAKE_CURRENT_SOURCE_DIR}/util.h
)

set( sources
    ${CMAKE_CURRENT_SOURCE_DIR}/aca.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bendseqlookup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chains.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constraints.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/edges.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/expansion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/faces.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/graphs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hola.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/logging.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nearalign.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nexes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nodeconfig.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nodes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ortho.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/peeling.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/planarise.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/qalookup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quadaction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/routing.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sides.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/treeplacement.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/trees.cpp
)
include_directories("./")
 
add_library(${PROJECT_NAME} STATIC ${sources} ${headers})


