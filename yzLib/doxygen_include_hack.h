/***********************************************************/
/**	\file
	\brief		Hack Header for Doxygen Documentation
	\details	The purpose of this file is to hack doxygen.
				We define the macro of these files, so that
				conditional include will be enabled in yzLib.
				In this way, the include graph will be complete.

				Caution! This file is just used for hack doxygen,
				it should never be used in yzLib
	\author		Yizhong Zhang
	\date		6/11/2012
*/
/***********************************************************/
#ifndef __DOXYGEN_INCLUDE_HACK_H__
#define __DOXYGEN_INCLUDE_HACK_H__

//	to enable yz_opengl
#define __GLEW_H__				///<	GL/glew.h
#define __GL_H__				///<	gl.h
#define __glut_h__				///<	GL/glut.h
#define __freeglut_h__			///<	freeglut.h

//	to enable yz_cuda
#define __cuda_cuda_h__			///<	cuda.h
#define __CUDA_RUNTIME_API_H__	///<	cuda_runtime_api.h
#define __CUDA_GL_INTEROP_H__	///<	cuda_gl_interp.h
#define CUBLAS_V2_H_			///<	cublas_v2.h

//	to enable yz_windows
#define _WINDOWS_				///<	windows.h

//	to enable kinect
#define NUIAPI					///<	NuiApi.h

//	to enable mkl
#define _MKL_H_					///<	mkl.h
#define	_MKL_TYPES_H_			///<	mkl_types.h
#define _MKL_LAPACKE_H_			///<	mkl_lapacke.h
#define _MKL_SPBLAS_H_			///<	mkl_spblas.h
#define __MKL_DSS_H				///<	mkl_dss.h
#define _MKL_RCISOLVER_H_		///<	mkl_rci.h

//	to enable visualization
#define IMDEBUG_H				///<	imdebug.h
#define FREEIMAGE_H				///<	FreeImage.h

#include "yzLib/yz_lib.h"



#endif	//	__DOXYGEN_INCLUDE_HACK_H__