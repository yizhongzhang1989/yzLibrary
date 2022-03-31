/***********************************************************/
/**	\file
	\brief		direct solver using mkl
	\details	mkl_dss.h required
	\author		Yizhong Zhang
	\date		10/6/2012
*/
/***********************************************************/
#ifndef __YZ_MKL_DIRECT_SOLVER_H__
#define __YZ_MKL_DIRECT_SOLVER_H__

#ifndef YZ_mkl_types_h
#	error	yz_mkl_direct_solver.h must be included after mkl_types.h or mkl.h
#endif
#ifndef YZ_mkl_dss_h
#	error	yz_mkl_direct_solver.h must be included after mkl_dss.h
#endif

#include <iostream>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"

namespace yz{


/**
	Direct Sparse Solver using MKL, Ax = b, given A and b, calculate x

	This class only allow float or double type due to mkl constraint, 
	other type will result in construction error

	To use this class, you need	\n
	1,	call SetupMatrix() to setup matrix A, or call SetupMatrixStructure() and Factor() independently	\n
	2,	call Solve() to solve the linear system. It is possible to solve several x-b pairs with same A

	We use default mkl options. If you want to use your own options, change the options explicitly
	before you call the corresponding functions.
*/
template<class T>
class MKLDirectSparseSolver{
public:
	//	constructor
	MKLDirectSparseSolver(){
		if( !TypeGuard() ){
			std::cout << "error: invalid TYPE passed to MKLDirectSparseSolver\n\
						 MKLDirectSparseSolver only accept float or double" << std::endl;
			exit(0);
		}

		handle = NULL;
		sym_type[0] = MKL_DSS_NON_SYMMETRIC;
		sym_type[1] = MKL_DSS_SYMMETRIC;
		sym_type[2] = MKL_DSS_SYMMETRIC_STRUCTURE;

		Reset();
	}

	~MKLDirectSparseSolver(){
		Reset();
	}

	/**
		Reset to startup status
	*/
	void Reset(){
		symmetric_flag	= 0;
		indexing		= 0;
		row_num			= 0;
		col_num			= 0;
		nnz				= 0;
		value_ptr		= NULL;
		col_id_ptr		= NULL;
		row_start_ptr	= NULL;
		x_ptr			= NULL;
		b_ptr			= NULL;
		nrhs			= 0;

		create_option	= MKL_DSS_MSG_LVL_WARNING | MKL_DSS_TERM_LVL_ERROR | 
			(indexing? 0 : MKL_DSS_ZERO_BASED_INDEXING) | SinglePrecisionFlag();
		reorder_option	= MKL_DSS_DEFAULTS;
		factor_option	= MKL_DSS_POSITIVE_DEFINITE;
		solve_option	= MKL_DSS_DEFAULTS;
		delete_option	= MKL_DSS_DEFAULTS;


		if(handle){		//	release the handle
			SafeCall( dss_delete(handle, delete_option), "delete handle in Reset" );
			handle = NULL;
		}

		stage			= DSS_NONE;
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
	int SetupMatrix(int rows, int cols, int nnz, T* value, int* col_id, int* row_start, int symmetric_flag=0, int indexing=0){
		assert( bool(this->indexing) != bool(create_option & MKL_DSS_ZERO_BASED_INDEXING) );

		//	setup structure
		SetupMatrixStructure(rows, cols, nnz, col_id, row_start, symmetric_flag, indexing);

		//	factor
		Factor(value);

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
		setup matrix A structure

		\param	rows			number of rows
		\param	cols			number of columns
		\param	nnz				number of non-zero elements
		\param	col_id			pointer to colume id array of each element
		\param	row_start		pointer to row start id array of each row
		\param	symmetric_flag	0: non-symmetric, 1: symmetric, 2: symmetric structure\n
								symmetric structure means structure is symmetric, but values are not symmetric, \n
								the storage is the same as non-symmetric
		\param	indexing		indexing of the matrix, 0: zero-based indexing, 1: one-based indexing
		\return					whether success
	*/
	int SetupMatrixStructure(int rows, int cols, int nnz, int* col_id, int* row_start, int symmetric_flag=0, int indexing=0){
		if( rows != cols ){
			#ifndef BE_QUIET
				std::cout << "error: mkl only solve square matrix" << std::endl;
			#endif
			return 0;
		}

		//	stage check
		if( stage != DSS_NONE ){	//	if the solver is not new, then reset the solver
			#ifndef BE_QUIET
				std::cout << "error: MKLDirectSparseSolver not reset, cannot setup matrix structure" << std::endl;
			#endif
			//	I don't reset the solver automatically because options may not be correct.
			return 0;
		}

		//	copy matrix
		this->symmetric_flag	= symmetric_flag;
		this->indexing			= indexing;
		row_num			= rows;
		col_num			= cols;
		this->nnz		= nnz;
		col_id_ptr		= col_id;
		row_start_ptr	= row_start;

		//	set options by setting
		create_option &= ~MKL_DSS_ZERO_BASED_INDEXING;
		create_option |= (indexing ? 0 : MKL_DSS_ZERO_BASED_INDEXING);

		//	create handle
		SafeCall( dss_create(handle, create_option), "create in SetupMatrixStructure" );
		stage = DSS_CREATE;

		//	define structure
		SafeCall( dss_define_structure(handle, sym_type[symmetric_flag], 
			row_start_ptr, row_num, col_num, col_id_ptr, nnz ), 
			"define structure in SetupMatrixStructure" );
		stage = DSS_DEFINE_STRUCTURE;

		//	reorder
		SafeCall( dss_reorder(handle, reorder_option, 0), "reorder in SetupMatrixStructure" );
		stage = DSS_REORDER;

		return 1;
	}

	/**
		Factor the matrix A after structure defined

		\param	value		the matrix that match the structure which has been defined
		\return				whether success
	*/
	int Factor(T* value){
		//	stage check
		if( stage < DSS_REORDER ){	//	the stage is at least Reorder, or refuse to factor
			#ifndef BE_QUIET
				std::cout << "error: to factor the matrix, the matrix must be reordered" << std::endl;
			#endif
			return 0;
		}

		value_ptr = value;

		//	factor the matrix
		SafeCall( dss_factor_real(handle, factor_option, value_ptr), "factor In Factor" );
		stage = DSS_FACTOR;

		return 1;
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
	int Solve(T* x, T* b, int vec_num = 1){
		//	stage check
		if( stage < DSS_FACTOR ){	//	the stage is at least Factor, or refuse to solve
			#ifndef BE_QUIET
				std::cout << "error: to solve the matrix, the matrix must be factored" << std::endl;
			#endif
			return 0;
		}

		//	if result x is write back to the array b, then we must create a temp array
		//	currently, we cannot handle situation of memory overlapping
		T*	tmp = NULL;
		if( x == b )
			tmp = new T[col_num*vec_num];

		x_ptr	= (tmp? tmp : x);
		b_ptr	= b;
		nrhs	= vec_num;

		//	solve the system
		SafeCall(dss_solve_real(handle, solve_option, b_ptr, nrhs, x_ptr), "solve in Solve");
		stage = DSS_SOLVE;

		if( tmp ){
			memcpy(x, x_ptr, sizeof(T)*col_num*nrhs);
			delete tmp;
		}

		return 1;
	}

public:
	//	options used to control mkl dss
	MKL_INT	create_option;
	MKL_INT	reorder_option;
	MKL_INT	factor_option;
	MKL_INT	solve_option;
	MKL_INT	delete_option;

	//	mkl resources
	_MKL_DSS_HANDLE_t	handle;		///<	mkl solver handle
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

	/**
		call mkl dss functions with safe guard

		\param	error	the function return value
		\param	msg		message to pring
		\return			whether success
	*/
	inline int SafeCall(int error, char* msg=NULL){
		if( error != MKL_DSS_SUCCESS ){
			std::cout << msg << ": error in MKLDirectSparseSolver, code: " << error << std::endl;
			return 0;
		}
		return 1;
	}


protected:
	//	Matrix A in CSR format
	int		symmetric_flag;			///<	0: non-symmetric, 1: symmetric, 2: symmetric structure
	int		indexing;				///<	0: zero-based indexing, non-zero: one-based indexing
	int		row_num, col_num;		///<	matrix dimension
	int		nnz;					///<	number of non-zero elements
	T*		value_ptr;				///<	non-zero elements
	int*	col_id_ptr;				///<	non-zero element column index
	int*	row_start_ptr;			///<	start position of each row

	//	vectors
	T*		x_ptr;
	T*		b_ptr;
	int		nrhs;					///<	how many vectors are solved one time

	//	symmetric type keywords
	int sym_type[3];

	//	solver stage flag
	int	stage;
	enum{ DSS_NONE, DSS_CREATE, DSS_DEFINE_STRUCTURE, DSS_REORDER, DSS_FACTOR, DSS_SOLVE };

};

template<> inline int MKLDirectSparseSolver<float>::TypeGuard(){
	return 1;
}

template<> inline int MKLDirectSparseSolver<double>::TypeGuard(){
	return 1;
}

template<> inline int MKLDirectSparseSolver<float>::SinglePrecisionFlag(){
	return MKL_DSS_SINGLE_PRECISION;
}

typedef MKLDirectSparseSolver<float>	MKLDirectSparseSolverf;
typedef MKLDirectSparseSolver<double>	MKLDirectSparseSolverd;



}	//	end namespace yz

#endif	//	__YZ_MKL_DIRECT_SOLVER_H__