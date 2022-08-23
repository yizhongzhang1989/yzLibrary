/***********************************************************/
/**	\file
	\brief		a collection of math algorithms in my library
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_MATH_ALGORITHM_H__
#define __YZ_MATH_ALGORITHM_H__

//	config
#include "yzLib_config.h"

//	setting
#include "yzLib/yz_setting.h"

//	------------------------
//	link mkl lib
//	------------------------
#ifdef YZ_LINK_MKL_RT_LIB
#	pragma comment(lib, "mkl_rt.lib")
#endif

//	lookup table
#include "yzLib/yz_math/yz_lookup_table.h"

//	interpolation
#include "yzLib/yz_math/yz_interpolation.h"

//	resampling
#include "yzLib/yz_math/yz_resampling.h"

//	filter
#include "yzLib/yz_math/yz_filter.h"

//	statistics
#include "yzLib/yz_math/yz_statistics/yz_statistics_basic.h"
#include "yzLib/yz_math/yz_statistics/yz_regression_2d.h"
#include "yzLib/yz_math/yz_statistics/yz_regression_3d.h"

//	geometry
#include "yzLib/yz_math/yz_graph_algorithm.h"
#include "yzLib/yz_math/yz_curve.h"

//	blas and sparse
#include "yzLib/yz_math/yz_blas_sparse_operators.h"
#include "yzLib/yz_math/yz_solver.h"

#ifdef yzLib_ENABLE_Eigen
#	include "yzLib/yz_math/yz_eigen_direct_solver.h"
#endif

#ifdef yzLib_ENABLE_MKL
#	include "yzLib/yz_math/yz_mkl_spblas_c.h"
#	include "yzLib/yz_math/yz_pca.h"
#endif


#endif	//	__YZ_MATH_ALGORITHM_H__