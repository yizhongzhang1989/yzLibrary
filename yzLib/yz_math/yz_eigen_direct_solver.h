/***********************************************************/
/**	\file
	\brief		direct solver using eigen
	\details	
	\author		Yizhong Zhang
	\date		11/15/2018
*/
/***********************************************************/
#ifndef __YZ_EIGEN_DIRECT_SOLVER_H__
#define __YZ_EIGEN_DIRECT_SOLVER_H__

#ifndef YZ_eigen_sparse_h
#	error	yz_eigen_direct_solver.h must be included after Eigen/Sparse
#endif

#include <iostream>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"

namespace yz{

/**
	Direct Sparse Solver using Eigen, Ax = b, given A and b, calculate x

	To use this class, you need	\n
	1,	call SetupMatrix() to setup matrix A, or call SetupMatrixStructure() and Factor() independently	\n
	2,	call Solve() to solve the linear system. It is possible to solve several x-b pairs with same A

	We use default mkl options. If you want to use your own options, change the options explicitly
	before you call the corresponding functions.
*/
template<class T>
class EigenDirectSparseSolver{
public:
	//	constructor
	EigenDirectSparseSolver(){
		if (!TypeGuard()) {
			std::cout << "error: invalid TYPE passed to EigenDirectSparseSolver\n\
						 EigenDirectSparseSolver only accept float or double" << std::endl;
			exit(0);
		}

		Reset();
	}

	~EigenDirectSparseSolver(){
		Reset();
	}

	/**
		Reset to startup status
	*/
	void Reset(){
		A.resize(0, 0);
		X.resize(0, 0);
		b.resize(0, 0);
	}

	/**
		Setup Matrix A

		This function consist of setup matrix structure and factor,
		if we have multiple matrix of same structure but different value,
		we can call them individually

		\param	rows			number of rows
		\param	cols			number of columns
		\param	nnz				number of non-zero elements
		\param	value			pointer to the matrix element array
		\param	col_id			pointer to colume id array of each element
		\param	row_start		pointer to row start id array of each row
		\param	symmetric_flag	0: non-symmetric, 1: symmetric, 2: symmetric structure\n
								symmetric structure means structure is symmetric, but values are not symmetric, \n
								the storage is the same as non-symmetric
		\param	indexing		indexing of the matrix, 0: zero-based indexing, non-zero: one-based indexing
		\return					whether success
	*/
	int SetupMatrix(
		int rows, 
		int cols, 
		int nnz, 
		T* value, 
		int* col_id, 
		int* row_start, 
		int symmetric_flag=0, 
		int indexing=0
	){
		if (indexing)
			indexing = -1;	//	offset the col index
		if (symmetric_flag)
			symmetric_flag = 1;

		A.resize(rows, cols);
		A.reserve(nnz + symmetric_flag*nnz);
		for (int row_id = 0; row_id < rows; row_id++) {
			for (int id = row_start[row_id]; id < row_start[row_id + 1]; id++) {
				A.insert(row_id, col_id[id] + indexing) = value[id];
				if (symmetric_flag && row_id != col_id[id] + indexing)	//	solve symmetric matrix
					A.insert(col_id[id] + indexing, row_id) = value[id];
			}
		}

		A.makeCompressed();

		solver.compute(A);

		return 1;
	}

	/**
		Setup Matrix A with COO input

		This function consist of setup matrix structure and factor,
		if we have multiple matrix of same structure but different value,
		we can call them individually

		\param	rows			number of rows
		\param	cols			number of columns
		\param	nnz				number of non-zero elements
		\param	value			pointer to the matrix element array
		\param	row_id			pointer to row id array of each element
		\param	col_id			pointer to colume id array of each element
		\param	symmetric_flag	0: non-symmetric, 1: symmetric, 2: symmetric structure\n
								symmetric structure means structure is symmetric, but values are not symmetric, \n
								the storage is the same as non-symmetric
		\param	indexing		indexing of the matrix, 0: zero-based indexing, non-zero: one-based indexing
		\return					whether success
	*/
	int SetupMatrixCOO(
		int rows,
		int cols,
		int nnz,
		T* value,
		int* row_id,
		int* col_id,
		int symmetric_flag = 0,
		int indexing = 0
	) {
		if (indexing)
			indexing = -1;	//	offset the col index
		if (symmetric_flag)
			symmetric_flag = 1;

		A.resize(rows, cols);
		A.reserve(nnz + symmetric_flag*nnz);
		for (int i = 0; i < nnz; i++) {
			A.insert(row_id[i] - indexing, col_id[i] - indexing) = value[i];
			if (symmetric_flag && row_id[i] != col_id[i])
				A.insert(col_id[i] - indexing, row_id[i] - indexing) = value[i];
		}

		A.makeCompressed();

		solver.compute(A);

		return 1;
	}

	/**
		Setup Matrix A from CSR sparse matrix structure directly
	*/
	int SetupMatrix(SparseMatrixCSR<T>& spa_mat){
		int ret = SetupMatrix(spa_mat.row_num, spa_mat.col_num, spa_mat.NNZ(),
			&spa_mat.value[0], &spa_mat.col_id[0], &spa_mat.row_start[0], 
			spa_mat.IsSymmetric(), spa_mat.indexing);
		return ret;
	}

	/**
		Setup Matrix A from COO sparse matrix structure directly
	*/
	int SetupMatrix(SparseMatrixCoo<T>& spa_mat) {
		int ret = SetupMatrix(spa_mat.row_num, spa_mat.col_num, spa_mat.NNZ(),
			&spa_mat.value[0], &spa_mat.row_id[0], &spa_mat.col_id[0],
			spa_mat.IsSymmetric(), spa_mat.indexing);
		return ret;
	}


	/**
		Solve the factored linear system, Ax = b

		If address of x b are identical, then the result is write back to
		the array. If the space overlap, then the result may be incorrect

		\param	x			x in Ax = b
		\param	b			b in Ax = b
		\param	vec_num		how many vectors contained in the array\n
							mkl is possible to calculate several by one call,
							default is 1
		\return				whether success
	*/
	int Solve(T* x_ptr, T* b_ptr, int vec_num = 1){
		if (A.rows() == 0)	//	A is empty
			return 0;
		if (vec_num < 1)
			return 0;

		b.resize(A.rows(), vec_num);
		for (int j = 0; j < vec_num; j++) {
			for (int i = 0; i < A.rows(); i++) {
				b(i, j) = b_ptr[j*A.rows() + i];
			}
		}

		//	solve the system
		X = solver.solve(b);

		//	write data
		for (int j = 0; j < vec_num; j++) {
			for (int i = 0; i < A.rows(); i++) {
				x_ptr[j*A.rows() + i] = X(i, j);
			}
		}

		return 1;
	}

protected:
	/**
		Check the type of this class

		MKL only accept float or double type, so if user try to pass
		other type just let the constructor fail
	*/
	inline int TypeGuard(){
		return 0;
	}

	/**
		default double precision

		single precision is achieved in template specialization
	*/
	inline int SinglePrecisionFlag(){
		return 0;	//	double precision
	}

protected:
	//	Eigen variables
	Eigen::SparseMatrix<T>								A;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	X;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	b;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<T> >		solver;
};

template<> inline int EigenDirectSparseSolver<float>::TypeGuard() {
	return 1;
}

template<> inline int EigenDirectSparseSolver<double>::TypeGuard() {
	return 1;
}


typedef EigenDirectSparseSolver<float>	EigenDirectSparseSolverf;
typedef EigenDirectSparseSolver<double>	EigenDirectSparseSolverd;

}	//	end namespace yz

#endif	//	__YZ_EIGEN_DIRECT_SOLVER_H__