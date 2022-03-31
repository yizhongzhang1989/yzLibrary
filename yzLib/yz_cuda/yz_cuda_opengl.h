/***********************************************************/
/**	\file
	\brief		functions to handle opengl in cuda
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_CUDA_OPENGL_H__
#define __YZ_CUDA_OPENGL_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_cuda_h
#	error yz_cuda_opengl.h must be included after cuda.h
#endif
#ifndef YZ_cuda_runtime_api_h
#	error yz_cuda_opengl.h must be included after cuda_runtime_api.h
#endif
#ifndef YZ_cuda_gl_interop_h
#	error yz_cuda_opengl.h must be included after cuda_gl_interop.h
#endif
#ifndef YZ_gl_h
#	error yz_cuda_opengl.h must be included after gl.h
#endif


#include "yzLib/yz_cuda/yz_cuda_basic.h"

namespace yz{
namespace cuda{
/**
	Init OpenGL in CUDA
*/
extern "C"
inline void InitCudaGL(int device_id = 0){
	cudaSafeCall( cudaGLSetGLDevice( device_id ), "set cuda device" );
}

}	//	namespace cuda
}	//	namespace yz

#endif	//	__YZ_CUDA_OPENGL_H__