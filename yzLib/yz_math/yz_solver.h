/***********************************************************/
/**	\file
	\brief		solver of linear system
	\details	
	\author		Yizhong Zhang
	\date		10/6/2012
*/
/***********************************************************/
#ifndef __YZ_SOLVER_H__
#define __YZ_SOLVER_H__

#ifdef YZ_mkl_dss_h										//	enable mkl direct solver
#	include "yzLib/yz_math/yz_mkl_direct_solver.h"
#endif

#if defined(YZ_mkl_rci_h) && defined(YZ_mkl_spblas_h)	//	enable mkl iterative solver
#	include "yzLib/yz_math/yz_mkl_iter_solver.h"
#endif

namespace yz{





}	//	end namespace yz

#endif	//	__YZ_SOLVER_H__