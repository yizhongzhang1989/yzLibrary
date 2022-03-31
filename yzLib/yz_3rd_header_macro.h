/***********************************************************/
/**	\file
	\brief		Transform Macro of 3rd header files
	\details	it is possible that different version of 3rd 
				party libraries use different macro to control
				header file. If this happens, yz library include
				will cause error. To solve this, macro is 
				transformed in this file
	\author		Yizhong Zhang
	\date		10/5/2012
*/
/***********************************************************/
#ifndef __YZ_3RD_HEADER_MACRO_H__
#define __YZ_3RD_HEADER_MACRO_H__


//	opengl
#if defined(__gl_h_) || defined(__GL_H__)				//	gl.h
#	define	YZ_gl_h
#endif

#if defined(__glut_h__) || defined(__GLUT_H__)			//	glut.h
#	define	YZ_glut_h
#endif

#if defined(__glew_h__) || defined(__GLEW_H__)			//	glew.h
#	define	YZ_glew_h
#endif

#if defined(__freeglut_h__) || defined(__FREEGLUT_H__)	//	freeglut.h
#	define	YZ_freeglut_h
#endif

//	cuda
#ifdef __cuda_cuda_h__			//	cuda.h
#	define	YZ_cuda_h
#endif

#ifdef __CUDA_RUNTIME_API_H__	//	cuda_runtime_api.h
#	define	YZ_cuda_runtime_api_h
#	define	YZ_LINK_CUDART_LIB
#endif

#ifdef __CUDA_GL_INTEROP_H__	//	cuda_gl_interp.h
#	define	YZ_cuda_gl_interop_h
#endif

#ifdef CUBLAS_V2_H_				//	cublas_v2.h
#	define	YZ_cublas_v2_h
#	define	YZ_LINK_CUBLAS_LIB
#endif


//	windows
#ifdef _WINDOWS_				//	windows.h
#	define	YZ_windows_h
#endif


//	kinect
#ifdef NUIAPI
#	define	YZ_NuiApi_h
#endif


//	eigen
#if (defined(EIGEN_CORE_H) || defined(EIGEN_CORE_MODULE_H)) && defined(EIGEN_LU_MODULE_H) && defined(EIGEN_CHOLESKY_MODULE_H) && defined(EIGEN_QR_MODULE_H) && defined(EIGEN_SVD_MODULE_H) && defined(EIGEN_GEOMETRY_MODULE_H) && defined(EIGEN_EIGENVALUES_MODULE_H)
#	define	YZ_eigen_dense_h
#endif

#ifdef EIGEN_SPARSE_MODULE_H
#	define	YZ_eigen_sparse_h
#endif


//	mkl
#ifdef _MKL_H_					//	mkl.h
#	define	YZ_mkl_h
#	define	YZ_LINK_MKL_RT_LIB
#endif

#ifdef _MKL_TYPES_H_			//	mkl_types.h
#	define	YZ_mkl_types_h
#endif

#ifdef _MKL_LAPACKE_H_			//	mkl_lapacke.h
#	define	YZ_mkl_lapacke_h
#	define	YZ_LINK_MKL_RT_LIB
#endif

#ifdef _MKL_SPBLAS_H_			//	mkl_spblas.h
#	define	YZ_mkl_spblas_h
#	define	YZ_LINK_MKL_RT_LIB
#endif

#ifdef __MKL_DSS_H				//	mkl_dss.h
#	define	YZ_mkl_dss_h
#	define	YZ_LINK_MKL_RT_LIB
#endif

#ifdef _MKL_RCISOLVER_H_		//	mkl_rci.h
#	define	YZ_mkl_rci_h
#	define	YZ_LINK_MKL_RT_LIB
#endif


//	imdebug
#ifdef IMDEBUG_H				//	imdebug.h
#	define	YZ_imdebug_h
#endif


//	free image
#ifdef FREEIMAGE_H				//	FreeImage.h
#	define	YZ_FreeImage_h
#endif


#endif	//	__YZ_3RD_HEADER_MACRO_H__