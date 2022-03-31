/***********************************************************/
/**	\file
	\brief		PCA
	\details	functions used for PCA, MKL C interface is used
	\author		Yizhong Zhang
	\date		8/27/2012
*/
/***********************************************************/
#ifndef __YZ_PCA_H__
#define __YZ_PCA_H__

#include "yzLib/yz_setting.h"

#if !defined(YZ_mkl_h) || !defined(YZ_mkl_lapacke_h)
#	error yz_pca.h must be included after mkl.h and mkl_lapacke.h
#endif

#include <stdlib.h>
#include "yzLib/yz_math/yz_vector_utils.h"
#include "yzLib/yz_math/yz_matrix_utils.h"

namespace yz{

//	========================================
///@{
/**	@name PCA related functions
*/
//	========================================

/**
	eigen decomposition, calculate eigen vectors, eigen values and mean

	MKL C interface is used to do SVD decomposition.

	\param	eigen_vector	return eigen vectors of data, space must be big enough
	\param	eigen_value		return eigen values corresponding to eigen vectors, 
							in decreasing order, space must be big enough
	\param	mean			return mean vector, space must be big enough
	\param	data			data to perform PCA, arranged vector by vector. 
							after calling this function, the value of data will be changed
	\param	dim_num			dimensions of the data
	\param	vec_num			how many vectors does the data have
	\param	ret_eigens		how many eigen vectors do we expact to return, if <=0, then return all
	\return					actual eigen number returned. return 0 if decomposition failed
*/
inline int calculateEigens(float* eigen_vector,
						   float* eigen_value, 
						   float* mean, 
						   float* data, 
						   int dim_num, 
						   int vec_num, 
						   int ret_eigens = 0){

	int size = (dim_num<vec_num? dim_num : vec_num);	//	size = min(dim_num, vec_num)

	//	step 1, subtract mean from each dimension
	setVectorValue(mean, 0, dim_num);
	for( int i=0; i<vec_num; i++ )
		setVectorAddVector(mean, data+dim_num*i, dim_num);
	setVectorScale(mean, 1.0/vec_num, dim_num);
	for( int i=0; i<vec_num; i++ )
		setVectorSubVector(data+dim_num*i, mean, dim_num);

	//	step 2, SVD using Lapack C interface
	MKL_INT info;
	float* superb = new float[size];
	if( ret_eigens<=0 || size<ret_eigens ){	//	input space is big enough
		ret_eigens = size;
		info = LAPACKE_sgesvd(
			LAPACK_COL_MAJOR, 'S', 'N', dim_num, vec_num, data, dim_num,
			eigen_value, eigen_vector, dim_num, NULL, vec_num, superb);
		assert( info==0 );
	}
	else{						//	we don't need all the eigen values
		float* u = new float[dim_num*size];
		float* s = new float[size];
		info = LAPACKE_sgesvd(
			LAPACK_COL_MAJOR, 'S', 'N', dim_num, vec_num, data, dim_num,
			s, u, dim_num, NULL, vec_num, superb);
		assert( info==0 );
		memcpy(eigen_vector, u, sizeof(float)*dim_num*ret_eigens);
		memcpy(eigen_value, s, sizeof(float)*ret_eigens);
		delete u;
		delete s;
	}

	delete superb;

	if( info != 0 )
		return 0;
	else
		return ret_eigens;
}

/**
	calculate the projection coefficients of a vector on orthogonal eigen vectors

	v = v_mean + a1 * e1 + a2 * e2 + ... + an * en

	a1, a2, ..., an are the coefficients we want to calculate

	\param	coef			return the coefficient of each eigen vector
	\param	vec				the vector to compute
	\param	mean			mean of original data
	\param	eigen_vector	eigen vectors of original data
	\param	dim_num			the dimension of each vector
	\param	eigen_num		eigen vector number
*/
inline void calculateProjectionOnEigens(float* coef, 
										float* vec, 
										float* mean, 
										float* eigen_vector, 
										int dim_num, 
										int eigen_num){
	float* tmp_vec = new float[dim_num];
	memcpy(tmp_vec, vec, sizeof(float)*dim_num);
	setVectorSubVector(tmp_vec, mean, dim_num);

	for(int i=0; i<eigen_num; i++)
		coef[i] = calculateVectorDot(tmp_vec, eigen_vector+i*dim_num, dim_num);

	delete tmp_vec;
}

/**
	reduce the dimension of a vector using mean and given eigen vectors

	\param	reduced_vec		the result vector after dimension reduction
	\param	original_vec	the original vector
	\param	mean			the mean of original data
	\param	eigen_vector	eigen vectors of original data
	\param	dim_num			dimension of vector
	\param	eigen_num		the number of eigen vectors, should be smaller than dim_num
*/
inline void calculateDimensionReductionOnEigens(float* reduced_vec,
												float* original_vec,
												float* mean,
												float* eigen_vector,
												int dim_num,
												int eigen_num){
	float* coef = new float[dim_num];
	float* tmp_vec = new float[dim_num];

	//	first, calculate coefficient of each eigen vector
	calculateProjectionOnEigens(coef, original_vec, mean, eigen_vector, dim_num, eigen_num);

	//	second, reconstruct the vector
	memcpy(reduced_vec, mean, sizeof(float)*dim_num);
	for(int i=0; i<eigen_num; i++){
		memcpy(tmp_vec, eigen_vector+i*dim_num, sizeof(float)*dim_num);
		setVectorScale(tmp_vec, coef[i], dim_num);
		setVectorAddVector(reduced_vec, tmp_vec, dim_num);
	}
	delete coef;
	delete tmp_vec;
}

/**
	reduce the dimension of a vector using mean and given eigen vectors

	\param	reduced_vec		the result vector after dimension reduction
	\param	coef			return the coefficient of each eigen vector
	\param	original_vec	the original vector
	\param	mean			the mean of original data
	\param	eigen_vector	eigen vectors of original data
	\param	dim_num			dimension of vector
	\param	eigen_num		the number of eigen vectors, should be smaller than dim_num
*/
inline void calculateDimensionReductionOnEigens(float* reduced_vec,
												float* coef,
												float* original_vec,
												float* mean,
												float* eigen_vector,
												int dim_num,
												int eigen_num){
	float* tmp_vec = new float[dim_num];

	//	first, calculate coefficient of each eigen vector
	calculateProjectionOnEigens(coef, original_vec, mean, eigen_vector, dim_num, eigen_num);

	//	second, reconstruct the vector
	memcpy(reduced_vec, mean, sizeof(float)*dim_num);
	for(int i=0; i<eigen_num; i++){
		memcpy(tmp_vec, eigen_vector+i*dim_num, sizeof(float)*dim_num);
		setVectorScale(tmp_vec, coef[i], dim_num);
		setVectorAddVector(reduced_vec, tmp_vec, dim_num);
	}
	delete tmp_vec;
}

///@}

}	//	namespace yz

#endif	//	__YZ_PCA_H__