/***********************************************************/
/**	\file
	\brief		Basic Functions in CUDA
	\details	recommended to include yz_cuda.h as the main entry
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_CUDA_BASIC_H__
#define __YZ_CUDA_BASIC_H__

#include "yzLib/yz_setting.h"

#if !(defined(YZ_cuda_h) && defined(YZ_cuda_runtime_api_h))
#	error yz_cuda_basic.h must be included after cuda.h and cuda_runtime_api.h
#endif

#include <iostream>

namespace yz{
namespace cuda{
/**
	Call CUDA Fucntion with Safe Guard
*/
extern "C"
inline void cudaSafeCall(cudaError err, char* msg=NULL){
	if( cudaSuccess != err ){
		printf( "CUDA error(%s): %s\n", msg, cudaGetErrorString(err) );
		exit(-1);	//	cuda error cannot recover
	}
}
#define	cudaSafeCall yz::cuda::cudaSafeCall

#ifdef YZ_CUBLAS_ENABLE
/**
	Call CUBLAS Fucntion with Safe Guard
*/
extern "C"
inline void cublasSafeCall(cublasStatus_t err, char* msg=NULL){
	if( CUBLAS_STATUS_SUCCESS != err ){
		std::cout << "cublas error(" << msg << "): " << err << std::endl;
		exit(-1);
	}
}
#define cublasSafeCall yz::cuda::cublasSafeCall
#endif

//	alloc memory
#define cudaNew(ptr_d, size)	cudaSafeCall( cudaMalloc((void**)&(ptr_d), (size)), "malloc gpu memory")
#define cudaDelete(ptr_d)		cudaSafeCall( cudaFree(ptr_d), "free gpu memory" )

//	memory copy
#define cpyh2d( mem_d, mem_h, size )	cudaSafeCall( cudaMemcpy( mem_d, mem_h, size, cudaMemcpyHostToDevice ), "cpyh2d" )
#define cpyd2h( mem_h, mem_d, size )	cudaSafeCall( cudaMemcpy( mem_h, mem_d, size, cudaMemcpyDeviceToHost ), "cpyd2h" )
#define cpyd2d( des, src, size )		cudaSafeCall( cudaMemcpy( des, src, size, cudaMemcpyDeviceToDevice ), "cpyd2d" )

/**
	Set Device Memory through Main Memory

	\param	ptr_d		device memory ptr
	\param	val			value to be set
	\param	size		number of values to be set, total memory size: sizeof(val)*size
*/
template<typename T>
inline void setMemDeviceByCPU(T* ptr_d, T val, int size){
	T* ptr = new T[size];
	for(int i=0; i<size; i++)
		ptr[i] = val;

	cpyh2d(ptr_d, ptr, sizeof(T)*size);

	delete[] ptr;
}

}	//	namespace cuda
}	//	namespace yz

#endif	//	__YZ_CUDA_BASIC_H__