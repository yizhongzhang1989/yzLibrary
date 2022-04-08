

#################### yzLib functions ####################

function(yzLib_Print_Info)
	message(STATUS "++++++++++++++++++++++++++++++++++++++++++++++++++")
	message(STATUS "yzLib included, add yzLib_INCLUDE_DIR into your include directory to use")
	message(STATUS "yzLib_INCLUDE_DIR=${yzLib_INCLUDE_DIR}")

	message(STATUS "config:")
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
		if(_variableName MATCHES yzLib_ENABLE_*)
			message(STATUS "${_variableName}=${${_variableName}}")
		endif()
	endforeach()

	message(STATUS "--------------------------------------------------")
endfunction()


function(yzLib_set_ENABLE MODULE_NAME VALUE)	# copied from ceres solver
	set(VAR_NAME yzLib_ENABLE_${MODULE_NAME})
	get_property(IS_DEFINED_IN_CACHE CACHE ${VAR_NAME} PROPERTY VALUE SET)
	if (NOT IS_DEFINED_IN_CACHE)
		message(FATAL_ERROR "yzLib_set_ENABLE(${MODULE_NAME} ${VALUE}), ${VAR_NAME} not defined.")
	endif()
	get_property(HELP_STRING CACHE ${VAR_NAME} PROPERTY HELPSTRING)
	get_property(VAR_TYPE CACHE ${VAR_NAME} PROPERTY TYPE)
	set(${VAR_NAME} ${VALUE} CACHE ${VAR_TYPE} "${HELP_STRING}" FORCE)

	configure_file(
	"${yzLib_DIR}/yzLib_config.h.in"
	"${yzLib_DIR}/yzLib_config.h")
endfunction()



#################### yzLib config ####################

set(yzLib_DIR ${CMAKE_CURRENT_LIST_DIR})
get_filename_component(yzLib_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../ ABSOLUTE)

option(yzLib_ENABLE_WINDOWS		"Whether Enable Windows in yzLib"		OFF)
option(yzLib_ENABLE_Eigen		"Whether Enable EIGEN in yzLib"			OFF)
option(yzLib_ENABLE_MKL			"Whether Enable MKL in yzLib"			OFF)
option(yzLib_ENABLE_GLUT		"Whether Enable OpenGL in yzLib"		OFF)
option(yzLib_ENABLE_GLEW		"Whether Enable GLEW in yzLib"			OFF)
option(yzLib_ENABLE_FreeImage	"Whether Enable FreeImage in yzLib"		OFF)
option(yzLib_ENABLE_CUDA		"Whether Enable CUDA in yzLib"			OFF)

configure_file(
	"${yzLib_DIR}/yzLib_config.h.in"
	"${yzLib_DIR}/yzLib_config.h")

