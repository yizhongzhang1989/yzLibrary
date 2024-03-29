/***********************************************************/
/**	\file
	\brief		OpenGL related functions in yzLib
	\details	This file must be included after <gl.h>, or error 
				at compile time. 

				OpenGL module is not an essential part of yzLib,
				so if <gl.h> is not included ahead, this file is 
				simply ignored to insure compile time safety.

				yz_vector_oengl_utils.h will include yz_vector,
				be cautious if your implementation have vector 
				implementation too.
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_OPENGL_H__
#define __YZ_OPENGL_H__

#pragma warning(push)
#pragma warning(disable: 4267)	//	disable warning of size_t int conversion

//	setting
#include "yzLib_config.h"
#include "yzLib/yz_setting.h"

//	for x64 programs
#ifdef _WIN64
#	pragma comment(lib, "glut64.lib")
#endif

//	if glew.h included ahead
#ifdef yzLib_ENABLE_GLEW
#	include "yzLib/yz_opengl/yz_fbo.h"
#	include "yzLib/yz_opengl/yz_vbo.h"
#	include "yzLib/yz_opengl/yz_shader.h"
#endif


//	utils , if glut.h / freeglut.h included ahead
#ifdef yzLib_ENABLE_GLUT

#	include "yzLib/yz_opengl/yz_opengl_utils.h"
#	include "yzLib/yz_opengl/yz_glut_window.h"
#	include "yzLib/yz_opengl/yz_texture.h"
#	include "yzLib/yz_opengl/yz_ascii_displayer.h"

#	ifndef IGNORE_YZ_ANIMATION
#		include "yzLib/yz_opengl/yz_animation_opengl_utils.h"
#	endif

#	ifndef INCLUDE_YZ_OPENGL_IGNORE_YZ_VECTOR
#		include "yzLib/yz_opengl/yz_vector_opengl_utils.h"
#	endif

#endif

#pragma warning(pop)

#endif	//	__YZ_OPENGL_H__