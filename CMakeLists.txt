cmake_minimum_required (VERSION 3.2)

message("This CMakeList.txt is written for 64bit system only")

project (yzLibExamples LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(${MSVC}) 
	message(STATUS "Build with MSVC")
elseif(${UNIX})
	message(STATUS "Build with UNIX or LINUX")
else()
	message(WARNING "Build with unknown compiler")
endif()


###############################################
# Options
###############################################
option(ENABLE_OPENGL "Whether enable OpenGL when build the project" ON)
option(ENABLE_GLEW "Whether enable GLEW when build the project" ON)
option(ENABLE_MKL "Whether enable MKL when build the project" OFF)
option(ENABLE_CUDA "Whether enable CUDA when build the project" OFF)

###############################################
# Dependency
###############################################
message(STATUS "==============================")


# ---------- yzLibraryDeps ----------
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/deps/inc)
	message(FATAL_ERROR "deps not found, update submodule to get yzLibraryDeps")
else()
	message(STATUS "yzLibraryDeps found as submodule")
	include_directories(
		${CMAKE_SOURCE_DIR}
		${PROJECT_BINARY_DIR}
		${CMAKE_SOURCE_DIR}/deps/inc
		${PROJECT_BINARY_DIR} )
	if(${MSVC}) 
		set(BIN_PATH "${CMAKE_SOURCE_DIR}/deps/bin/win64;${CMAKE_SOURCE_DIR}/deps/bin/win32")
	endif()
endif()

# ---------- Eigen ----------
find_package(Eigen3 QUIET)

if(Eigen3_FOUND)
	message(STATUS "Eigen3 found: ${Eigen3_VERSION}")
	set(EIGEN_FOUND TRUE)
elseif(EXISTS ${CMAKE_SOURCE_DIR}/3rdParty/eigen/Eigen)
	message(STATUS "Eigen found as submodule")
	set(EIGEN_FOUND TRUE)
	include_directories(3rdParty/eigen)
else()
	message(WARNING "Eigen not found, update submodule to get Eigen")
endif()

message(STATUS "==============================")


###############################################
# Link
###############################################
if(${MSVC}) 	
	link_directories(${CMAKE_SOURCE_DIR}/deps/lib/win64)	
endif()


###############################################
# Output
###############################################
if(NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE Release)
endif() 
set(CMAKE_DEBUG_POSTFIX "d")

if(${MSVC}) 
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_BINARY_DIR}/bin/Release)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_BINARY_DIR}/bin/Debug)
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_BINARY_DIR}/bin/Release)
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_BINARY_DIR}/bin/Debug)
endif()


###############################################
# Config data path
###############################################
set(DATA_PATH "${CMAKE_SOURCE_DIR}/deps/data")
configure_file(
	"${CMAKE_SOURCE_DIR}/examples/data_path.h.in" 
	"${PROJECT_BINARY_DIR}/data_path.h"
	@ONLY)
message(STATUS "${PROJECT_BINARY_DIR}/data_path.h created,\n DATA_PATH=${DATA_PATH}")


###############################################
# Projects
###############################################
macro(SUBDIRLIST result curdir)
	file(GLOB children RELATIVE ${curdir} ${curdir}/*)
	set(dirlist "")
	foreach(child ${children})
	if(IS_DIRECTORY ${curdir}/${child})
		list(APPEND dirlist ${child})
	endif()
	endforeach()
	set(${result} ${dirlist})
endmacro()

SUBDIRLIST(example_groups ${CMAKE_SOURCE_DIR}/examples)

foreach(group ${example_groups})
	if(EXISTS ${CMAKE_SOURCE_DIR}/examples/${group}/CMakeLists.txt)
		add_subdirectory(examples/${group})				
	else()
		SUBDIRLIST(examples ${CMAKE_SOURCE_DIR}/examples/${group})
		foreach(example ${examples})
			if(EXISTS ${CMAKE_SOURCE_DIR}/examples/${group}/${example}/CMakeLists.txt)
				add_subdirectory(examples/${group}/${example})		
			endif()
		endforeach()
	endif()
endforeach()


