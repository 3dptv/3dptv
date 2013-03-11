# Main CMakeLists file
cmake_minimum_required (VERSION 2.8)

project (3dptv_solution)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# Add definitions and Find TIFF

if (WIN32)
  add_definitions( -D_CONSOLE -D_CRT_SECURE_NO_WARNINGS )
  set(TIFF_LIBRARY "${CMAKE_SOURCE_DIR}/src_c/libtiff.lib")
else()
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/opt/ActiveTcl-8.6")
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/usr/lib")
  find_package( TIFF REQUIRED )
endif()

# Find TCL/TK
find_package( TCL REQUIRED 8.6)

if (TCLTK_FOUND )
	message(STATUS "TCL INCLUDE: " ${TCL_INCLUDE_PATH})
	message(STATUS "TCL LIB: " ${TCL_LIBRARY})
	message(STATUS "TK INCLUDE: " ${TCL_INCLUDE_PATH})
	message(STATUS "TK LIB: " ${TK_LIBRARY})
	link_directories( ${TCL_LIBRARY} ${TK_LIBRARY} )
	include_directories( ${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH} )
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")

message(STATUS "gathering ptv tcl files:")
file(GLOB ptv_tcl "${CMAKE_SOURCE_DIR}/src_tcl/*.tcl")

add_subdirectory("${CMAKE_SOURCE_DIR}/src_c")
