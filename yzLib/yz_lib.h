/***********************************************************/
/**	\file
	\brief		the root of yz_lib, include everything in yzLib
	\details	If certain module is not for c++ standard library, then
				the header of that module should be included before this 
				header, then this module will be enabled automatically.

				Include sequence of different modules should not change
				because some files contain conditional include or conditional
				compile information. If these modules don't follow a 
				pre-designed sequence, it is possible that not all modules
				are loaded.
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_LIB_H__
#define __YZ_LIB_H__

//	setting
#include "yzLib/yz_setting.h"

//	math
#include "yzLib/yz_math.h"

//	geometry
#include "yzLib/yz_geometry.h"

//	image
#include "yzLib/yz_image.h"

//	animation
#ifdef	INCLUDE_YZ_ANIMATION
#	include "yzLib/yz_animation.h"
#endif

//	physics
#ifdef	INCLUDE_YZ_PHYSICS
#	include "yzLib/yz_physics.h"
#endif

//	utils
#include "yzLib/yz_utils.h"

//	system
#ifdef	INCLUDE_YZ_SYSTEM
#	include "yzLib/yz_system.h"
#endif

//	opengl, if gl.h is included ahead
#ifdef YZ_gl_h
#	include "yzLib/yz_opengl.h"
#endif

//	cuda, if cuda.h is included ahead
#ifdef YZ_cuda_h
#	include "yzLib/yz_cuda.h"
#endif

//	windows, if windows.h is include ahead
#ifdef YZ_windows_h
#	include "yzLib/yz_windows.h"
#endif

#endif	//	__YZ_LIB_H__