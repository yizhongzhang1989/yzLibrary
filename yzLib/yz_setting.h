/***********************************************************/
/**	\file
	\brief		Setting of yzLib
	\details	All control options of yzLib should be contained
				in this file. All files in the root directory 
				should include this file.
	\author		Yizhong Zhang
	\date		6/19/2012
*/
/***********************************************************/
#ifndef __YZ_SETTING_H__
#define __YZ_SETTING_H__

#include "yzLib/yz_3rd_header_macro.h"

/**
	enable print info flag

	if you don't want to print anything, disable this macro
*/
#define	PRINT_INFO

#ifndef	PRINT_INFO	//	switch the macro because BE_QUIET is already used
#	define BE_QUIET
#endif

/**
	enable openmp

	A lot of functions can be accelerated by openmp, but for a lot 
	of compilers, <omp.h> must be included and omp enabled. But if
	we don't want openmp anyway, disable the macro
*/
#define	ENABLE_OPENMP

/**
	include yz_animation

	disable this macro when you don't want animation related functions
	to be included
*/
#define INCLUDE_YZ_ANIMATION

/**
	include yz_physics

	physics is used for physically based simulation, but if we
	don't want that to be included in yzLib, disable this macro
*/
#define	INCLUDE_YZ_PHYSICS

/**
	include yz_system

	this is used when you don't want system related functions
	to be included
*/
#define INCLUDE_YZ_SYSTEM

/**
	include yz_kinect

	if we don't want to include any kinect related resources, disable
	this macro
*/
#define INCLUDE_YZ_KINECT

/**
	include yz_opengl ignore yz_vector

	This is used when you just want to use yz_opengl, but don't 
	want to include yz_vector.
*/
//#define INCLUDE_YZ_OPENGL_IGNORE_YZ_VECTOR

//	yz_animation will use yz_vector
#ifndef IGNORE_YZ_ANIMATION
#	ifdef INCLUDE_YZ_OPENGL_IGNORE_YZ_VECTOR
#		error yz_animation will use yz_vector, cannot define INCLUDE_YZ_OPENGL_IGNORE_YZ_VECTOR
#	endif
#endif


#endif	//	__YZ_SETTING_H__