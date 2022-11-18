cmake_minimum_required (VERSION 3.0)
project (adaptagrams LANGUAGES C CXX)
add_definitions(-DUNICODE -D_UNICODE -DNOMINMAX -D_USE_MATH_DEFINES)
#add_definitions(-DUSE_ASSERT_EXCEPTIONS -D__ASSERT_FUNCTION -DUNICODE -D_UNICODE -DNOMINMAX -D_USE_MATH_DEFINES)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})


set(CMAKE_DEBUG_POSTFIX _d)

#set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};Modules/")
find_package(Java REQUIRED COMPONENTS Runtime Development)
find_package(JNI REQUIRED)
find_package(SWIG REQUIRED)
include(UseJava)
include(UseSWIG)
include_directories(${JAVA_INCLUDE_PATH})

if (WIN32)
include_directories(${JAVA_INCLUDE_PATH}/win32)
endif(WIN32)

set(CMAKE_SWIG_FLAGS -package adaptagrams)

include_directories("./")
add_subdirectory("libvpsc")
add_subdirectory("libcola")
add_subdirectory("libavoid")
add_subdirectory("libdialect")
add_subdirectory("libproject")
add_subdirectory("libtopology")

set_property(SOURCE adaptagrams.i PROPERTY CPLUSPLUS ON)
set_property(SOURCE adaptagrams.i PROPERTY OUTPUT_DIR "${CMAKE_SOURCE_DIR}/java/")
set_property(SOURCE adaptagrams.i PROPERTY OUTPUT_FILE "${CMAKE_SOURCE_DIR}/java/build")

swig_add_library(adaptagrams LANGUAGE java SOURCES adaptagrams.i)
swig_link_libraries(adaptagrams libvpsc  libcola libavoid  libdialect  libproject libtopology)
#SET_TARGET_PROPERTIES(adaptagrams PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(adaptagrams libvpsc  libcola libavoid  libdialect  libproject libtopology)
