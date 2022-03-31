/***********************************************************/
/**	\file
	\brief		Visualize Sparse Matrix
	\author		Yizhong Zhang
	\date		9/7/2012
*/
/***********************************************************/
#ifndef __YZ_SPARSE_VISUALIZE_H__
#define __YZ_SPARSE_VISUALIZE_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_imdebug_h
#	error yz_sparse_visualize.h must be included after imdebug.h
#endif

#include <iostream>
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_math/yz_sparse_convert.h"

namespace yz{ namespace utils{

//	========================================
///@{
/**	@name Visualize Sparse Matrix
*/
//	========================================

/**
	plot coordinate format sparse matrix

	currently, use imdebug to visualize non-zero elements
*/
template<typename T>
inline int plotSparseMatrix(const SparseMatrixCoo<T>& spa_mat){
	const SparseMatrixCoo<T>* spa_mat_ptr = &spa_mat;

	//	check the matrix
	if( (*spa_mat_ptr).col_num <=0 || (*spa_mat_ptr).row_num <= 0 ){
		std::cout << "sparse matrix dimensin error" << std::endl;
		return 0;
	}
	if( (*spa_mat_ptr).value.size() != (*spa_mat_ptr).row_id.size() || (*spa_mat_ptr).value.size() != (*spa_mat_ptr).col_id.size() ){
		std::cout << "sparse matrix vector size doesn't match" << std::endl;
		return 0;
	}

	SparseMatrixCoo<T> ghost_spa_mat;
	if( spa_mat.indexing ){
		ghost_spa_mat = spa_mat;
		ghost_spa_mat.SetIndexing(0);
		spa_mat_ptr = &ghost_spa_mat;
	}

	//	visualize the matrix
	if( (*spa_mat_ptr).row_num * (*spa_mat_ptr).col_num < 2048*2048 ){	//	the matrix size dons't exceed imdebug limit
		char* tmp = new char[(*spa_mat_ptr).row_num * (*spa_mat_ptr).col_num];
		memset(tmp, 0, sizeof(char)*(*spa_mat_ptr).row_num*(*spa_mat_ptr).col_num);
		for(int i=0; i<(*spa_mat_ptr).value.size(); i++)
			tmp[(*spa_mat_ptr).row_id[i]*(*spa_mat_ptr).col_num+(*spa_mat_ptr).col_id[i]] = 255;

		imdebug("lum w=%d h=%d %p", (*spa_mat_ptr).col_num, (*spa_mat_ptr).row_num, tmp);

		delete tmp;
	}
	else{	//	the matrix is too big
		int width = (*spa_mat_ptr).col_num, height = (*spa_mat_ptr).row_num;
		while(width*height > 2048*2048){
			width /= 2;
			height /= 2;
		}

		char* tmp = new char[width*height];
		memset(tmp, 0, sizeof(char)*width*height);
		for(int i=0; i<(*spa_mat_ptr).value.size(); i++){
			int m = float((*spa_mat_ptr).row_id[i]) / (*spa_mat_ptr).row_num * width;
			int n = float((*spa_mat_ptr).col_id[i]) / (*spa_mat_ptr).col_num * height;
			tmp[m*width + n] = 255;
		}

		imdebug("lum w=%d h=%d %p", width, height, tmp);

		delete tmp;
	}

	return 1;
}

/**
	plot CSR format sparse matrix

	currently, use imdebug to visualize non-zero elements
*/
template<typename T>
inline int plotSparseMatrix(const SparseMatrixCSR<T>& spa_mat){
	if( spa_mat.IsSymmetric() ){
		SparseMatrixCooSym<T> mat_coo;
		convertCSRToCoo(mat_coo, spa_mat);
		mat_coo.SetIndexing(0);
		return	plotSparseMatrix(mat_coo);
	}
	else{
		SparseMatrixCoo<T> mat_coo;
		convertCSRToCoo(mat_coo, spa_mat);
		mat_coo.SetIndexing(0);
		return	plotSparseMatrix(mat_coo);
	}
}

///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_SPARSE_VISUALIZE_H__