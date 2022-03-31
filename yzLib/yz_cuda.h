/***********************************************************/
/**	\file
	\brief		This file include compoments to make using of CUDA easier
	\details	include the following header ahead if you want to use modules in yzLib

				#include <cuda.h>				//	yz_cuda

				#include <cuda_runtime_api.h>	//	needed all the time

				#include <cuda_gl_interop.h>	//	OpenGL in cuda

				#include <cublas_v2.h>			//	CUBLAS

				#include <imdebug.h>			//	cuda debug

				Since cuda source file are .cu files that need to compile
				using nvcc, they cannot be included in header only, so 
				manually add .cu files to your project if you need them
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_CUDA_H__
#define __YZ_CUDA_H__

//	setting
#include "yzLib/yz_setting.h"


#if !( defined(YZ_cuda_h) && defined(YZ_cuda_runtime_api_h) )
#	error yz_cuda.h must be included after cuda.h and cuda_runtime_api.h
#endif

//	------------------------
//	link cuda lib
//	------------------------
#ifdef YZ_LINK_CUDART_LIB
#	pragma comment( lib, "cudart.lib" )
#endif

#ifdef YZ_LINK_CUBLAS_LIB
#	pragma comment( lib, "cublas.lib" )
#endif

//	------------------------
//	include yzLib
//	------------------------

#include "yzLib/yz_cuda/yz_cuda_basic.h"

#ifdef YZ_imdebug_h	//	cuda debug use 3rd party library Image Debugger, make sure you have that library installed
#	include "yzLib/yz_cuda/yz_cuda_debug.h"
#endif

#ifdef YZ_cuda_gl_interop_h	//	use opengl in cuda
#	include "yzLib/yz_cuda/yz_cuda_opengl.h"
#	ifdef YZ_glew_h
#		include "yzLib/yz_cuda/yz_cuda_fbo.h"
#	endif
#endif	//	#ifdef __CUDA_GL_INTEROP_H__


#endif	//	__YZ_CUDA_H__