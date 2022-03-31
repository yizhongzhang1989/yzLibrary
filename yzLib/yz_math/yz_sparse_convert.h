/***********************************************************/
/**	\file
	\brief		Convert Between Different Sparse Matrix Formats
	\author		Yizhong Zhang
	\date		9/9/2012
*/
/***********************************************************/
#ifndef __YZ_SPARSE_CONVERT_H__
#define __YZ_SPARSE_CONVERT_H__

#include <assert.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"

namespace yz{

//	========================================
///@{
/**	@name Coo CSR convert
*/
//	========================================

/**
	convert from Coo to CSR

	csr and coo must both be symmetric or not symmetric
*/
template<typename T>
int convertCooToCSR(SparseMatrixCSR<T>& csr, const SparseMatrixCoo<T>& coo){
	if( (csr.IsSymmetric() && !coo.IsSymmetric() ) ||
		(!csr.IsSymmetric() && coo.IsSymmetric() ) ){
		#ifndef BE_QUIET
			std::cout << "error: convertCooToCSR, cannot convert between symmetric and non-symmetric matrix" << std::endl;
		#endif
		return 0;
	}

	//	we copy coo and perform reorder
	SparseMatrixCoo<T> coo2 = coo;
	coo2.Reorder();

	//	setup csr
	csr.Reset();
	csr.SetDimension(coo2.row_num, coo2.col_num);
	csr.SetIndexing(coo2.indexing);
	csr.value = coo2.value;
	csr.col_id = coo2.col_id;

	//	calculate start index of each row
	csr.row_start[0] = bool(coo2.indexing);
	int end_pos = bool(coo2.indexing);
	for(int i=0; i<coo2.row_num; i++){
		while( end_pos != coo2.row_id.size() && coo2.row_id[end_pos] <= i )
			end_pos ++;
		csr.row_start[i+1] = end_pos;
	}

	return 1;
}

/**
	convert from CSR to Coo

	csr and coo must both be symmetric or not symmetric
*/
template<typename T>
int convertCSRToCoo(SparseMatrixCoo<T>& coo, const SparseMatrixCSR<T>& csr){
	if( (csr.IsSymmetric() && !coo.IsSymmetric() ) ||
		(!csr.IsSymmetric() && coo.IsSymmetric() ) ){
		#ifndef BE_QUIET
			std::cout << "error: convertCSRToCoo, cannot convert between symmetric and non-symmetric matrix" << std::endl;
		#endif
		return 0;
	}

	//	setup csr
	coo.Reset();
	coo.SetDimension(csr.row_num, csr.col_num);
	coo.SetIndexing(csr.indexing);
	coo.value = csr.value;
	coo.col_id = csr.col_id;

	//	calculate row index
	coo.row_id.resize(csr.NNZ());
	for(int i=0; i<csr.row_num; i++){
		for(int j=csr.row_start[i]-bool(csr.indexing); j<csr.row_start[i+1]-bool(csr.indexing); j++){
			coo.row_id[j] = i+bool(csr.indexing);
		}
	}

	return 1;
}


///@}

//	========================================
///@{
/**	@name Symmetric Non-Symmetric convert
*/
//	========================================

/**
	convert symmetric csr matrix to non-symmetric csr matrix

	\todo	implement this function

	\param	csr			the target matrix
	\param	csr_sym		the source matrix, which is symmetric
	\return				whether convert succeed
*/
template<typename T>
int convertCSRSymToCSR(SparseMatrixCSR<T> csr, const SparseMatrixCSR<T>& csr_sym){
	if( csr.IsSymmetric() ){
		#ifndef BE_QUIET
			std::cout << "error: convertCSRSymToCSR, the target is still symmetric matrix" << std::endl;
		#endif
		return 0;
	}

	if( !csr_sym.IsSymmetric() ){	//	the source is not symmetric, no need to convert
		csr = csr_sym;
		return 1;
	}

	std::cout << "error: convertCSRSymToCSR, this function not implemented" << std::endl;

	return 0;
}

///@}

}	//	end namespace yz

#endif	//	__YZ_SPARSE_CONVERT_H__