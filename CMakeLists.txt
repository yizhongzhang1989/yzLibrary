cmake_minimum_required (VERSION 3.2)

message("This CMakeList.txt is written for 64bit system only")

project (yzLibExamples LANGUAGES CXX)
set(yzLibrary_ROOT_DIR ${CMAKE_SOURCE_DIR})


if(${MSVC}) 
	message(STATUS "Build with MSVC")
elseif(${UNIX})
	message(STATUS "Build with UNIX or LINUX")
else()
	message(WARNING "Build with unknown compiler")
endif()


###############################################
# Settings
###############################################
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
add_definitions(-DUNICODE -D_UNICODE)

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
# Options and config
###############################################
include(yzLib/yzLib.cmake)
include_directories(${yzLib_INCLUDE_DIRS})


###############################################
# Subdirectory and dependency
###############################################
message(STATUS "==============================")

if(${MSVC}) 
	if(NOT EXISTS ${CMAKE_SOURCE_DIR}/deps/inc)
		# yzLibraryDeps not found, only basic examples
		message(WARNING "deps not found, only basic examples. Update submodule yzLibraryDeps to release examples with GUI.")
	else()
		# yzLibraryDeps found, enable features inside yzLibraryDeps
		message(STATUS "yzLibraryDeps found as submodule, enable GLUT, GLEW, FreeImage, and Windows API")

		yzLib_set_ENABLE(GLUT ON)
		yzLib_set_ENABLE(GLEW ON)
		yzLib_set_ENABLE(FreeImage ON)
		yzLib_set_ENABLE(WINDOWS ON)

		include_directories(${CMAKE_SOURCE_DIR}/deps/inc)
		set(BIN_PATH "${CMAKE_SOURCE_DIR}/deps/bin/win64;${CMAKE_SOURCE_DIR}/deps/bin/win32")
	endif()
elseif(${UNIX})
	find_package(OpenGL QUIET)
	if(OpenGL_FOUND)
		message(STATUS "OpenGL found")
		yzLib_set_ENABLE(GLUT ON)
	else()
		message(WARNING "OpenGL not found, try \n\tsudo apt install libglu1-mesa-dev freeglut3-dev mesa-common-dev")
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
	yzLib_set_ENABLE(Eigen ON)
	include_directories(3rdParty/eigen)
else()
	message(WARNING "Eigen not found, update submodule to get Eigen")
endif()

# ---------- CUDA ----------
if(yzLib_ENABLE_CUDA)
	find_package(CUDA  REQUIRED)

	if(CUDA_FOUND)
		message(STATUS "CUDA version: ${CUDA_VERSION}")
		include_directories(${CUDA_INCLUDE_DIRS})  
		link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
	else()
		message(ERROR "CUDA not found, you need to install CUDA, or disable yzLib_ENABLE_CUDA")
	endif()
endif()

message(STATUS "==============================")

yzLib_Print_Info()


###############################################
# Link
###############################################
if(${MSVC}) 	
	link_directories(${CMAKE_SOURCE_DIR}/deps/lib/win64)	
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
# Example Projects
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


###############################################
# Doxygen
###############################################  
message(STATUS "==================================================")

set(Graphviz_ROOT "C:/Program Files/Graphviz")	# default install path
find_path(Graphviz_BIN_PATH NAMES dot.exe PATHS ${Graphviz_ROOT}/bin)
set(HAVE_DOT "YES")
if(NOT EXISTS ${Graphviz_BIN_PATH}/dot.exe)
	set(Graphviz_BIN_PATH "")
	set(HAVE_DOT "NO")
endif()

if(Graphviz_BIN_PATH)
	message(STATUS "dot.exe found, diagrams will be generated by doxygen")
else()
	message(WARNING "cannot find dot.exe. install Graphviz and set the directory containing dot.exe into Graphviz_BIN_PATH, or the document generated by doxygen will not contail diagrams.")
endif()

configure_file(
  "${PROJECT_SOURCE_DIR}/yzLib_Doxy.in"
  "${PROJECT_SOURCE_DIR}/yzLib_Doxy"
  @ONLY
  )
