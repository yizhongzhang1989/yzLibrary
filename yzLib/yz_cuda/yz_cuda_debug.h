/***********************************************************/
/**	\file
	\brief		Debug GPU Memory
	\details	3rd party library imdebug is used
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_CUDA_DEBUG_H__
#define __YZ_CUDA_DEBUG_H__

#include "yzLib/yz_setting.h"

#if !(defined(YZ_cuda_h) && defined(YZ_cuda_runtime_api_h) && defined(YZ_imdebug_h))
#	error yz_cuda_debug.h must be included after cuda.h, cuda_runtime_api.h and imdebug.h
#endif

#include "yzLib/yz_cuda/yz_cuda_basic.h"

namespace yz{
namespace cuda{
/**
	Watch GPU Memory using luminance

	\param	ptr_d		device memory ptr
	\param	width		width of the memory space
	\param	height		height of the memory space
*/
inline void watch(const float* ptr_d, int width, int height){
	int size = width * height;
	float* ptr = new float[size];
	cpyd2h( ptr, ptr_d, sizeof(float)*size );
	imdebug("lum b=32f w=%d h=%d %p", width, height, ptr);
	delete ptr;
}

/**
	Watch GPU Memory using RGB

	\param	ptr_d		device memory ptr
	\param	width		width of the memory space
	\param	height		height of the memory space
*/
inline void watch(const float3* ptr_d, int width, int height){
	int size = width * height;
	float* ptr = new float[size*3];
	cpyd2h( ptr, ptr_d, sizeof(float3)*size );
	imdebug("rgb b=32f w=%d h=%d %p", width, height, ptr);
	delete ptr;
}

/**
	Watch GPU Memory as int* using luminance

	\param	ptr_d		device memory ptr
	\param	width		width of the memory space
	\param	height		height of the memory space
*/
inline void watch(const int* ptr_d, int width, int height){
	int size = width * height;
	float* float_ptr = new float[size];
	int* int_ptr = (int*)float_ptr;
	cpyd2h( int_ptr, ptr_d, sizeof(int)*size );	//	transform from GPU as int
	for( int i=0; i<width*height; i++ ){		//	transform each element to float
		float_ptr[i] = int_ptr[i];
	}
	imdebug("lum b=32f w=%d h=%d %p", width, height, float_ptr);
	delete float_ptr;
}

/**
	Watch Packed GPU Memory of unsigned char*

	\param	ptr_d		device memory ptr
	\param	pack_size	1:	draw as luminance image; 
						3:	draw as rgb image; 
						4:	draw as rgba image
	\param	width		width of the memory space
	\param	height		height of the memory space
*/
inline void watch(const unsigned char* ptr_d, int pack_size, int width, int height){
	int size = width * height * pack_size;
	unsigned char* ptr = new unsigned char[size];
	cpyd2h( ptr, ptr_d, sizeof(unsigned char)*size );
	if( pack_size == 1 )
		imdebug("lum w=%d h=%d %p", width, height, ptr);
	else if( pack_size == 3 )
		imdebug("rgb w=%d h=%d %p", width, height, ptr);
	else if( pack_size == 4 )
		imdebug("rgba w=%d h=%d %p", width, height, ptr);
	
	delete ptr;
}

/**
	Watch Packed GPU Memory

	\param	ptr_d		device memory ptr
	\param	pack_size	1:	draw as luminance image; 
						2:	draw as luminance alpha image;
						3:	draw as rgb image; 
						4:	draw as rgba image
	\param	width		width of the memory space
	\param	height		height of the memory space
*/
inline void watch(const float* ptr_d, int pack_size, int width, int height){
	int size = width * height * pack_size;
	float* ptr = new float[size];
	cpyd2h( ptr, ptr_d, sizeof(float)*size );
	if( pack_size == 1 )
		imdebug("lum b=32f w=%d h=%d %p", width, height, ptr);
	else if( pack_size == 2 )
		imdebug("luma b=32f w=%d h=%d %p", width, height, ptr);
	else if( pack_size == 3 )
		imdebug("rgb b=32f w=%d h=%d %p", width, height, ptr);
	else if( pack_size == 4 )
		imdebug("rgba b=32f w=%d h=%d %p", width, height, ptr);
	delete ptr;
}

}	//	namespace cuda
}	//	namespace yz

#endif	//	__YZ_CUDA_DEBUG_H__