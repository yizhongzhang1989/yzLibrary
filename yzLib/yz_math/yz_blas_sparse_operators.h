/***********************************************************/
/**	\file
	\brief		Blas and Sparse by Operator Overloading
	\details	We do not use const type, because mkl interface
				doesn't accept const
	\author		Yizhong Zhang
	\date		9/11/2012
*/
/***********************************************************/
#ifndef __YZ_BLAS_SPARSE_H__
#define __YZ_BLAS_SPARSE_H__

#include <assert.h>
#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_dense_vector.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_math/yz_solver.h"
#include "yzLib/yz_math/yz_sparse_convert.h"
#ifdef YZ_mkl_spblas_h	//	if mkl_spblas.h included, we can use mkl functions
#	include "yzLib/yz_math/yz_mkl_spblas_c.h"
#endif

namespace yz{

//	========================================
///@{
/**	@name Sparse Matrix * Dense Vector
*/
//	========================================

/**
	Coordinate Sparse Matrix * Dense Vector

	If dimension of matrix and vector don't match, then print
	debug message and exit

	If mkl_spblas.h included, then we use mkl to calculate
	if both matrix and vector are float or double
*/
template<typename T1, typename T2>
inline DenseVector<PROMOTE_T1_T2> mul(SparseMatrixCoo<T1>&	spa_mat, 
									  DenseVector<T2>&		den_vec){
	if( spa_mat.col_num != den_vec.Dim() ){
		std::cout << "error:mul(coo), sparse matrix * dense vector failed, dimension doesn't match" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_spblas_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
		DenseVector<float> result(den_vec.Dim());

		if( spa_mat.IsSymmetric() ){
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_scoosymv('U', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (float*)&den_vec[0], (float*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_scoosymv('U', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (float*)&den_vec[0], (float*)&result[0]);
			}
		}
		else{
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_scoogemv('N', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (float*)&den_vec[0], (float*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_scoogemv('N', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (float*)&den_vec[0], (float*)&result[0]);
			}
		}

		return result;
	}
	else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
		DenseVector<double> result(den_vec.Dim());

		if( spa_mat.IsSymmetric() ){
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_dcoosymv('U', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (double*)&den_vec[0], (double*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_dcoosymv('U', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (double*)&den_vec[0], (double*)&result[0]);
			}
		}
		else{
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_dcoogemv('N', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (double*)&den_vec[0], (double*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_dcoogemv('N', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
					spa_mat.NNZ(), (double*)&den_vec[0], (double*)&result[0]);
			}
		}

		return result;
	}
#endif	//	if we don't have mkl or the type is not legal, we use the calculation ourselves

	DenseVector<TYPE_PROMOTE(T1, T2)> result(den_vec.Dim());
	if( spa_mat.IsSymmetric() ){	//	symmetric matrix, the lower triangular should also be added
		for(int i=0; i<spa_mat.NNZ(); i++){
			int r_id = spa_mat.row_id[i];
			int c_id = spa_mat.col_id[i];
			result[r_id] += spa_mat.value[i] * den_vec[c_id];
			if( r_id != c_id )
				result[c_id] += spa_mat.value[i] * den_vec[r_id];
		}
	}
	else{							//	non-symmetric matrix, just ordinary use
		for(int i=0; i<spa_mat.NNZ(); i++){
			int r_id = spa_mat.row_id[i];
			int c_id = spa_mat.col_id[i];
			result[r_id] += spa_mat.value[i] * den_vec[c_id];
		}
	}
	return result;
}

/**
	CSR Sparse Matrix * Dense Vector

	If dimension of matrix and vector don't match, then print
	debug message and exit

	If mkl_spblas.h included, then we use mkl to calculate
	if both matrix and vector are float or double
*/
template<typename T1, typename T2>
inline DenseVector<PROMOTE_T1_T2> mul(SparseMatrixCSR<T1>&	spa_mat, 
									  DenseVector<T2>&		den_vec){
	if( spa_mat.col_num != den_vec.Dim() ){
		std::cout << "error:mul(csr), sparse matrix * dense vector failed, dimension doesn't match" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_spblas_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
		DenseVector<float> result(den_vec.Dim());

		if( spa_mat.IsSymmetric() ){
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_scsrsymv('U', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(float*)&den_vec[0], (float*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_scsrsymv('U', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(float*)&den_vec[0], (float*)&result[0]);
			}
		}
		else{
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_scsrgemv('N', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(float*)&den_vec[0], (float*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_scsrgemv('N', spa_mat.row_num,
					(float*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(float*)&den_vec[0], (float*)&result[0]);
			}
		}

		return result;
	}
	else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
		DenseVector<double> result(den_vec.Dim());

		if( spa_mat.IsSymmetric() ){
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_dcsrsymv('U', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(double*)&den_vec[0], (double*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_dcsrsymv('U', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(double*)&den_vec[0], (double*)&result[0]);
			}
		}
		else{
			if( spa_mat.indexing ){	//	one-based indexing
				mkl::c_mkl_dcsrgemv('N', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(double*)&den_vec[0], (double*)&result[0]);
			}
			else{					//	zero-based indexing
				mkl::c_mkl_cspblas_dcsrgemv('N', spa_mat.row_num,
					(double*)&spa_mat.value[0], &spa_mat.row_start[0], &spa_mat.col_id[0], 
					(double*)&den_vec[0], (double*)&result[0]);
			}
		}

		return result;
	}
#endif	//	if we don't have mkl or the type is not legal, we use the calculation ourselves

	DenseVector<TYPE_PROMOTE(T1, T2)> result(den_vec.Dim());
	if( spa_mat.IsSymmetric() ){
		for(int r_id=0; r_id<spa_mat.row_num; r_id++)
			for(int j=spa_mat.row_start[r_id]; j<spa_mat.row_start[r_id+1]; j++){
				int c_id = spa_mat.col_id[j];
				result[r_id] += spa_mat.value[j] * den_vec[c_id];
				if( r_id != c_id )
					result[c_id] += spa_mat.value[j] * den_vec[r_id];
			}
	}
	else{
		for(int r_id=0; r_id<spa_mat.row_num; r_id++)
			for(int j=spa_mat.row_start[r_id]; j<spa_mat.row_start[r_id+1]; j++){
				int c_id = spa_mat.col_id[j];
				result[r_id] += spa_mat.value[j] * den_vec[c_id];
			}
	}
	return result;
}


///@}
//	========================================
///@{
/**	@name Sparse Matrix * value
*/
//	========================================
/**
	CSR Sparse Matrix * value
*/
inline SparseMatrixCSR<double> mul(const SparseMatrixCSR<double>&	spa_mat,
								   double							rhs){
	SparseMatrixCSR<double> result(spa_mat);
	for(int i=0; i<spa_mat.value.size(); i++){
		result.value[i] *= rhs;
	}
	return result;
}

/**
	CSR Sparse Matrix * value
*/
inline SparseMatrixCSR<float> mul(const SparseMatrixCSR<float>&	spa_mat,
								  float							rhs){
	SparseMatrixCSR<float> result(spa_mat);
	for(int i=0; i<spa_mat.value.size(); i++){
		result.value[i] *= rhs;
	}
	return result;
}

/**
	CSR Sparse Matrix * value
*/
template<typename T>
inline SparseMatrixCSR<double> mul(double							lhs,
								   const SparseMatrixCSR<double>&	spa_mat){
	SparseMatrixCSR<double> result(spa_mat);
	for(int i=0; i<spa_mat.value.size(); i++){
		result.value[i] *= lhs;
	}
	return result;
}

/**
	CSR Sparse Matrix * value
*/
template<typename T>
inline SparseMatrixCSR<float> mul(float							lhs,
								  const SparseMatrixCSR<float>&	spa_mat){
	SparseMatrixCSR<float> result(spa_mat);
	for(int i=0; i<spa_mat.value.size(); i++){
		result.value[i] *= lhs;
	}
	return result;
}

///@}

//	========================================
///@{
/**	@name Sparse Matrix * Sparse Matrix
*/
//	========================================
/**
	CSR Sparse Matrix * CSR Sparse Matrix

	we use mkl to calculate the result, so only float and double are allowed
*/
template<typename T1, typename T2>
inline SparseMatrixCSR<PROMOTE_T1_T2> mul(SparseMatrixCSR<T1>&	spa_mat1, 
										  SparseMatrixCSR<T2>&	spa_mat2){
	if( spa_mat1.col_num != spa_mat2.row_num ){
		std::cout << "error:mul(csr), csr*csr failed, dimension doesn't match" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_spblas_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	SparseMatrixCSR<T1> ghost_spa_mat1;
	SparseMatrixCSR<T2> ghost_spa_mat2;
	SparseMatrixCSR<T1>* spa_mat1_ptr = &spa_mat1;
	SparseMatrixCSR<T2>* spa_mat2_ptr = &spa_mat2;
	bool mkl_useable_flag = false;

	if( (Is_float<T1>::check_type && Is_float<T2>::check_type) ||	//	this problem can be solved using mkl only if both 
		(Is_double<T1>::check_type && Is_double<T2>::check_type) ){	//	T1 and T2 are float or double

		mkl_useable_flag = true;
		if( spa_mat1.IsSymmetric() || !spa_mat1.indexing ){	//	we can only handle non-symmetric, one-based indexing matrix
			//	handle symmetric
			if( spa_mat1.IsSymmetric() ){
				if( !convertCSRSymToCSR(ghost_spa_mat1, spa_mat1) ){
					std::cout << "error:mul(csr), csr*csr failed, cannot convert spa_mat1 to non-symmetric" << std::endl;
					exit(0);
				}
			}
			else
				ghost_spa_mat1 = spa_mat1;
			//	handle zero-based indexing
			if( !spa_mat1_ptr->indexing )
				ghost_spa_mat1.SetIndexing(1);

			spa_mat1_ptr = &ghost_spa_mat1;
		}
		if( spa_mat2.IsSymmetric() || !spa_mat2.indexing ){	//	we can only handle non-symmetric, one-based indexing matrix
			//	handle symmetric
			if( spa_mat2.IsSymmetric() ){
				if( !convertCSRSymToCSR(ghost_spa_mat2, spa_mat2) ){
					std::cout << "error:mul(csr), csr*csr failed, cannot convert spa_mat2 to non-symmetric" << std::endl;
					exit(0);
				}
			}
			else
				ghost_spa_mat2 = spa_mat2;
			//	handle zero-based indexing
			if( !spa_mat2_ptr->indexing )
				ghost_spa_mat2.SetIndexing(2);

			spa_mat2_ptr = &ghost_spa_mat2;
		}

	}

	if( mkl_useable_flag ){
		//	now, spa_mat1_ptr and spa_mat2_ptr point to the two source matric
		if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
			//	first parse, calculate size of target
			SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1_ptr->row_num, spa_mat2_ptr->col_num);
			result.SetIndexing(1);
			int info;
			mkl::c_mkl_scsrmultcsr('N', 1, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num, spa_mat2_ptr->col_num,
				(float*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0],
				(float*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				NULL, NULL, (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:mul(csr), csr*csr float failed at step 1, code: " << info << std::endl;
				exit(0);
			}
			//	second parse, calculate the mul
			int value_size = result.row_start[spa_mat1_ptr->row_num] - 1;
			assert( value_size > 0 );
			result.value.resize(value_size);
			result.col_id.resize(value_size);
			mkl::c_mkl_scsrmultcsr('N', 2, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num, spa_mat2_ptr->col_num,
				(float*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0],
				(float*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				(float*)&result.value[0], (int*)&result.col_id[0], (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:mul(csr), csr*csr float failed at step 2, code: " << info << std::endl;
				exit(0);
			}

			return result;
		}
		else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
			//	first parse, calculate size of target
			SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1_ptr->row_num, spa_mat2_ptr->col_num);
			result.SetIndexing(1);
			int info;
			mkl::c_mkl_dcsrmultcsr('N', 1, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num, spa_mat2_ptr->col_num,
				(double*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0],
				(double*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				NULL, NULL, (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:mul(csr), csr*csr double failed at step 1, code: " << info << std::endl;
				exit(0);
			}
			//	second parse, calculate the mul
			int value_size = result.row_start[spa_mat1_ptr->row_num] - 1;
			assert( value_size > 0 );
			result.value.resize(value_size);
			result.col_id.resize(value_size);
			mkl::c_mkl_dcsrmultcsr('N', 2, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num, spa_mat2_ptr->col_num,
				(double*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0],
				(double*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				(double*)&result.value[0], (int*)&result.col_id[0], (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:mul(csr), csr*csr double failed at step 2, code: " << info << std::endl;
				exit(0);
			}

			return result;
		}
	}
#endif

	std::cout << "CSR * CSR not implemented, "
		<< "try to include mkl_spblas.h and use float or double" << std::endl;
	SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1.row_num, spa_mat2.col_num);
	return result;
}

///@}

//	========================================
///@{
/**	@name Sparse Matrix + Sparse Matrix
*/
//	========================================
///@}

//	========================================
///@{
/**	@name Sparse Matrix * Sparse Matrix
*/
//	========================================
/**
	CSR Sparse Matrix + CSR Sparse Matrix

	we use mkl to calculate the result, so only float and double are allowed
*/
template<typename T1, typename T2>
inline SparseMatrixCSR<PROMOTE_T1_T2> add(SparseMatrixCSR<T1>&	spa_mat1, 
										  SparseMatrixCSR<T2>&	spa_mat2){
	if( (spa_mat1.row_num != spa_mat2.row_num) || (spa_mat1.col_num != spa_mat2.col_num) ){
		std::cout << "error:add(csr), csr+csr failed, dimension doesn't match" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_spblas_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	SparseMatrixCSR<T1> ghost_spa_mat1;
	SparseMatrixCSR<T2> ghost_spa_mat2;
	SparseMatrixCSR<T1>* spa_mat1_ptr = &spa_mat1;
	SparseMatrixCSR<T2>* spa_mat2_ptr = &spa_mat2;
	bool mkl_useable_flag = false;

	if( (Is_float<T1>::check_type && Is_float<T2>::check_type) ||	//	this problem can be solved using mkl only if both 
		(Is_double<T1>::check_type && Is_double<T2>::check_type) ){	//	T1 and T2 are float or double

		mkl_useable_flag = true;
		if( spa_mat1.IsSymmetric() || !spa_mat1.indexing ){	//	we can only handle non-symmetric, one-based indexing matrix
			//	handle symmetric
			if( spa_mat1.IsSymmetric() ){
				if( !convertCSRSymToCSR(ghost_spa_mat1, spa_mat1) ){
					std::cout << "error:add(csr), csr+csr failed, cannot convert spa_mat1 to non-symmetric" << std::endl;
					exit(0);
				}
			}
			else
				ghost_spa_mat1 = spa_mat1;
			//	handle zero-based indexing
			if( !spa_mat1_ptr->indexing )
				ghost_spa_mat1.SetIndexing(1);

			spa_mat1_ptr = &ghost_spa_mat1;
		}
		if( spa_mat2.IsSymmetric() || !spa_mat2.indexing ){	//	we can only handle non-symmetric, one-based indexing matrix
			//	handle symmetric
			if( spa_mat2.IsSymmetric() ){
				if( !convertCSRSymToCSR(ghost_spa_mat2, spa_mat2) ){
					std::cout << "error:add(csr), csr+csr failed, cannot convert spa_mat2 to non-symmetric" << std::endl;
					exit(0);
				}
			}
			else
				ghost_spa_mat2 = spa_mat2;
			//	handle zero-based indexing
			if( !spa_mat2_ptr->indexing )
				ghost_spa_mat2.SetIndexing(2);

			spa_mat2_ptr = &ghost_spa_mat2;
		}

	}

	if( mkl_useable_flag ){
		//	now, spa_mat1_ptr and spa_mat2_ptr point to the two source matric
		if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
			//	first parse, calculate size of target
			SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1_ptr->row_num, spa_mat1_ptr->col_num);
			result.SetIndexing(1);
			int info;
			mkl::c_mkl_scsradd('N', 1, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num,
				(float*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0], 1,
				(float*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				NULL, NULL, (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:add(csr), csr+csr float failed at step 1, code: " << info << std::endl;
				exit(0);
			}
			//	second parse, calculate the mul
			int value_size = result.row_start[spa_mat1_ptr->row_num] - 1;
			assert( value_size > 0 );
			result.value.resize(value_size);
			result.col_id.resize(value_size);
			mkl::c_mkl_scsradd('N', 2, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num,
				(float*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0], 1,
				(float*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				(float*)&result.value[0], (int*)&result.col_id[0], (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:add(csr), csr+csr float failed at step 2, code: " << info << std::endl;
				exit(0);
			}

			return result;
		}
		else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
			//	first parse, calculate size of target
			SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1_ptr->row_num, spa_mat1_ptr->col_num);
			result.SetIndexing(1);
			int info;
			mkl::c_mkl_dcsradd('N', 1, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num,
				(double*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0], 1,
				(double*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				NULL, NULL, (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:add(csr), csr+csr double failed at step 1, code: " << info << std::endl;
				exit(0);
			}
			//	second parse, calculate the mul
			int value_size = result.row_start[spa_mat1_ptr->row_num] - 1;
			assert( value_size > 0 );
			result.value.resize(value_size);
			result.col_id.resize(value_size);
			mkl::c_mkl_dcsradd('N', 2, 8, 
				spa_mat1_ptr->row_num, spa_mat1_ptr->col_num,
				(double*)&spa_mat1_ptr->value[0], (int*)&spa_mat1_ptr->col_id[0], (int*)&spa_mat1_ptr->row_start[0], 1,
				(double*)&spa_mat2_ptr->value[0], (int*)&spa_mat2_ptr->col_id[0], (int*)&spa_mat2_ptr->row_start[0],
				(double*)&result.value[0], (int*)&result.col_id[0], (int*)&result.row_start[0],
				0, &info);
			if( info!=0 ){
				std::cout << "error:add(csr), csr+csr double failed at step 2, code: " << info << std::endl;
				exit(0);
			}

			return result;
		}
	}
#endif

	std::cout << "CSR + CSR not implemented, "
		<< "try to include mkl_spblas.h and use float or double" << std::endl;
	SparseMatrixCSR<PROMOTE_T1_T2> result(spa_mat1.row_num, spa_mat1.col_num);
	return result;
}



///@}

//	========================================
///@{
/**	@name Dense Vector / Sparse Matrix, solve linear system
*/
//	========================================
/**
	Solve Linear System Ax = b

	If mkl_spblas.h included, then we use mkl to calculate
	if both matrix and vector are float or double	

	\param	b	vector b in Ax = b
	\param	A	CSR matrix A in Ax = b
	\return		vector x in Ax = b
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> solve_dss(DenseVector<T1>&		b,
												   SparseMatrixCSR<T2>&	A){
	//	dimension check
	if( A.row_num != A.col_num ){
		std::cout << "dense vector / sparse matrix failed, matrix is not square matrix" << std::endl;
		exit(0);
	}
	if( A.col_num != b.Dim() ){
		std::cout << "dense vector / sparse matrix failed, matrix column number doesn't equal vector dimension" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_dss_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
		DenseVector<T1> result(b.Dim());

		MKLDirectSparseSolver<T1>	solver;
		solver.SetupMatrix(A);

		solver.Solve((T1*)&result[0], (T1*)&b[0]);

		return result;
	}
	else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
		DenseVector<T1> result(b.Dim());

		MKLDirectSparseSolver<T1>	solver;
		solver.SetupMatrix(A);

		solver.Solve((T1*)&result[0], (T1*)&b[0]);

		return result;
	}
#endif	//	if we don't have mkl or the type is not legal, we use the calculation ourselves

	std::cout << "dense vector / CSR sparse matrix not implemented, "
		<< "try to include mkl_dss.h and use float or double" << std::endl;
	DenseVector<TYPE_PROMOTE(T1, T2)> result(b.Dim());
	return result;	

}

/**
	Solve Linear System Ax = b while A is symmetric structure 

	If mkl_spblas.h included, then we use mkl to calculate
	if both matrix and vector are float or double	

	\param	b	vector b in Ax = b
	\param	A	CSR symmetric matrix A in Ax = b
	\return		vector x in Ax = b
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> solve_dss(DenseVector<T1>&			b,
												   SparseMatrixCSRSym<T2>&	A){
	//	dimension check
	if( A.row_num != A.col_num ){
		std::cout << "dense vector / sparse matrix failed, matrix is not square matrix" << std::endl;
		exit(0);
	}
	if( A.col_num != b.Dim() ){
		std::cout << "dense vector / sparse matrix failed, matrix column number doesn't equal vector dimension" << std::endl;
		exit(0);
	}

#ifdef YZ_mkl_dss_h		//	if mkl is used and the type is legal to call mkl, we use mkl to perform the calculation
	if( Is_float<T1>::check_type && Is_float<T2>::check_type ){
		DenseVector<T1> result(b.Dim());

		MKLDirectSparseSolver<T1>	solver;
		solver.SetupMatrix(A);

		solver.Solve((T1*)&result[0], (T1*)&b[0]);

		return result;
	}
	else if( Is_double<T1>::check_type && Is_double<T2>::check_type ){
		DenseVector<T1> result(b.Dim());

		MKLDirectSparseSolver<T1>	solver;
		solver.SetupMatrix(A);

		solver.Solve((T1*)&result[0], (T1*)&b[0]);

		return result;
	}
#endif	//	if we don't have mkl or the type is not legal, we use the calculation ourselves

	std::cout << "dense vector / CSR sparse matrix not implemented, "
		<< "try to include mkl_dss.h and use float or double" << std::endl;
	DenseVector<TYPE_PROMOTE(T1, T2)> result(b.Dim());
	return result;	
}

/**
	Solve Linear System Ax = b, A in coordinate format

	\param	b	vector b in Ax = b
	\param	A	Coo matrix A in Ax = b
	\return		vector x in Ax = b
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> solve_dss(DenseVector<T1>&		b,
												   SparseMatrixCoo<T2>&	A){
	if( A.IsSymmetric() ){
		SparseMatrixCSRSym<T2> Acsr(A.row_num);
		convertCooToCSR(Acsr, A);
		return b / Acsr;
	}
	else{
		SparseMatrixCSR<T2> Acsr(A.row_num, A.col_num);
		convertCooToCSR(Acsr, A);
		return b / Acsr;
	}
}


///@}

//	========================================
///@{
/**	@name Operator Overloading of Sparse
*/
//	========================================

/**
	Coordinate Sparse Matrix * Dense Vector
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> operator* (SparseMatrixCoo<T1>&	spa_mat, 
													DenseVector<T2>&		den_vec){
	return mul(spa_mat, den_vec);
}

/**
	CSR Sparse Matrix * Dense Vector
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> operator* (SparseMatrixCSR<T1>&	spa_mat, 
													DenseVector<T2>&		den_vec){
	return mul(spa_mat, den_vec);
}

/**
	Solve Linear System Ax = b
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> operator / (DenseVector<T1>&		b,
													 SparseMatrixCSR<T2>&	A){
	return solve_dss(b, A);
}

/**
	Solve Linear System Ax = b
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> operator / (DenseVector<T1>&			b,
													 SparseMatrixCSRSym<T2>&	A){
	return solve_dss(b, A);
}

/**
	Solve Linear System Ax = b, A in coordinate format
*/
template<typename T1, typename T2>
inline DenseVector<TYPE_PROMOTE(T1, T2)> operator / (DenseVector<T1>&		b,
													 SparseMatrixCoo<T2>&	A){
	return solve_dss(b, A);
}


///@}

}	//	end namespace yz

#endif	//	__YZ_BLAS_SPARSE_H__