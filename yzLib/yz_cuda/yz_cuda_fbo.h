/***********************************************************/
/**	\file
	\brief		FBO with CUDA extension
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_CUDA_FBO_H__
#define __YZ_CUDA_FBO_H__

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
#ifndef	YZ_gl_h
#	error yz_cuda_opengl.h must be included after gl.h
#endif
#ifndef YZ_glew_h
#	error yz_cuda_opengl.h must be included after glew.h
#endif

#include <iostream>
#include "yzLib/yz_opengl/yz_fbo.h"
#include "yzLib/yz_cuda/yz_cuda_basic.h"

namespace yz{
namespace cuda{
/**
	CUDA enabled FBO

	If CUDA is enabled, we can get the texture and write to a GPU memory
	space. In this way, we can get an image that displayed using OpenGL.

	After render to FBO, just call GetFBOData() to get the image
*/
class CudaFBO : public opengl::FBO{
public:
	float* tex_ptr_d;

public:
	//	constructor and destructor
	CudaFBO():opengl::FBO(), tex_ptr_d(0), cuda_resource_tex(0){}
	~CudaFBO(){
		Reset();
	}

	/**
		Initialize FBO according to given dimension, replace original InitFBO function

		\return		1: succeed;		0: failed
	*/
	inline int InitFBO(int width, int height){
		if( !opengl::FBO::InitFBO(width, height) ){
			std::cout << "Init FBO in CudaFBO failed" << std::endl;
			return 0;
		}

		//	CUDA alloc memory
		cudaSafeCall( cudaMalloc( (void**)&tex_ptr_d, sizeof(float)*4*tex_width*tex_height ), "CudaFBO alloc texture buffer" );

		//	CUDA
		cudaSafeCall( cudaGraphicsGLRegisterImage( &cuda_resource_tex, tex_id, GL_TEXTURE_2D, cudaGraphicsMapFlagsNone ), "CudaFBO register cuda resources" );

		return 1;
	}

	/**
		Get Data From Texture

		\return		1: succeed;		0: failed
	*/
	inline int GetFBOData(){
		if( ! (status & FBO_STATUS_BIT_INITIALIZED) ){
			std::cout << "fbo : " << fbo_id << " not initialized" << std::endl;
			return 0;
		}
		if( status & FBO_STATUS_BIT_BIND_FBO_ON ){
			std::cout << "fbo : " << fbo_id << " has started render, cannot get fbo data" << std::endl;
			return 0;
		}

		cudaSafeCall( cudaGraphicsMapResources( 1, &cuda_resource_tex, 0 ), "FBO map resources" );
		cudaArray* cuArray;
		cudaSafeCall( cudaGraphicsSubResourceGetMappedArray( &cuArray, cuda_resource_tex, 0, 0 ), "FBO get texture cuda array" );
		cudaMemcpyFromArray(tex_ptr_d, cuArray, 0, 0, sizeof(float)*4*tex_width*tex_height, cudaMemcpyDeviceToDevice);
		cudaSafeCall( cudaGraphicsUnmapResources( 1, &cuda_resource_tex, 0 ), "FBO unmap resources" );

		return 1;
	}

	/**
		Print Help Information
	*/
	inline void Help(){
		FBO::Help();

		std::cout
			<< "besides basic FBO, CudaFBO can get texture to device memory by calling GetFBOData\n"
			<< std::endl;
	}

protected:
	struct cudaGraphicsResource* cuda_resource_tex;

protected:
	/**
		Reset the class
	*/
	inline void Reset(){
		if( cuda_resource_tex ){
			cudaGraphicsUnregisterResource(cuda_resource_tex);
			cuda_resource_tex = NULL;
		}
		if( tex_ptr_d ){
			cudaSafeCall( cudaFree(tex_ptr_d) );
			tex_ptr_d = NULL;
		}

		FBO::Reset();
	}
};

}	//	namespace cuda
}	//	namespace yz

#endif	//	__YZ_CUDA_FBO_H__