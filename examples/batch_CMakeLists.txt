# Set project_name as the directory name
file(RELATIVE_PATH project_name ${CMAKE_CURRENT_LIST_DIR}/.. ${CMAKE_CURRENT_LIST_DIR})
file(RELATIVE_PATH folder_name ${CMAKE_CURRENT_LIST_DIR}/../.. ${CMAKE_CURRENT_LIST_DIR}/..)

###############################################
# Check config
###############################################
include(config.txt OPTIONAL)
if(REQUIRE_OPENGL AND NOT ENABLE_OPENGL)
	message("skipped project: ${folder_name}/${project_name}, OpenGL required")
	return()
endif()
if(REQUIRE_GLEW AND NOT ENABLE_GLEW)
	message("skipped project: ${folder_name}/${project_name}, GLEW required")
	return()
endif()
if(REQUIRE_MKL AND NOT ENABLE_MKL)
	message("skipped project: ${folder_name}/${project_name}, MKL required")
	return()
endif()
if(REQUIRE_CUDA AND NOT ENABLE_CUDA)
	message("skipped project: ${folder_name}/${project_name}, CUDA required")
	return()
endif()
if(REQUIRE_EIGEN AND NOT (EIGEN_FOUND OR Eigen3_FOUND))
	message("skipped project: ${folder_name}/${project_name}, EIGEN required")
	return()
endif()


if(REQUIRE_OPENGL)
	if(${MSVC})
		# we already included OpenGL, so don't need to do anything more
	else()
		find_package(OpenGL REQUIRED)
		find_package(GLUT REQUIRED)
	endif()
endif()
if(REQUIRE_GLEW)
	if(${MSVC})
		# we already included OpenGL, so don't need to do anything more
	else()
		find_package(GLEW REQUIRED)
	endif()
endif()


###############################################
# Add the project 
###############################################
project(${project_name})

file(GLOB  SOURCES
	"./*.h"
	"./*.cpp"
)

add_executable(${project_name} ${SOURCES})

if(REQUIRE_OPENGL)
	if(${MSVC})
			# we already included OpenGL, so don't need to do anything more
	else()
		target_link_libraries(${project_name} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
	endif()
endif()
if(REQUIRE_GLEW)
	if(${MSVC})
		# we already included OpenGL, so don't need to do anything more
	else()
		target_link_libraries(${project_name} ${GLEW_LIBRARIES})
	endif()
endif()


###############################################
# Add .user file for Visual Studio to add bin to path
###############################################
if(${MSVC}) 
	configure_file(
		"${yzLibrary_ROOT_DIR}/examples/vs_tmp.vcxproj.user.in" 
		"${PROJECT_BINARY_DIR}/${example}.vcxproj.user"
		@ONLY)
	# message(STATUS "created: ${PROJECT_BINARY_DIR}/${example}.vcxproj.user")
endif()


###############################################
# Add this project to corresponding folder
###############################################
set_target_properties(${project_name} PROPERTIES FOLDER ${folder_name})


message(STATUS "added project: ${folder_name}/${project_name}")
