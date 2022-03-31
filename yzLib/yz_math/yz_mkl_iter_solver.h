/***********************************************************/
/**	\file
	\brief		iterative solver using mkl
	\details	mkl_rci.h, mkl_spblas.h required
	\author		Yizhong Zhang
	\date		10/6/2012
*/
/***********************************************************/
#ifndef __YZ_MKL_ITER_SOLVER_H__
#define __YZ_MKL_ITER_SOLVER_H__

#ifndef YZ_mkl_rci_h
#	error	yz_mkl_iter_solver.h must be included after mkl_rci.h
#endif
#ifndef YZ_mkl_spblas_h
#	error	yz_mkl_iter_solver.h must be included after mkl_spblas.h
#endif

#include <iostream>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_math/yz_mkl_spblas_c.h"


namespace yz{

/**
	The virtual base of CG solver in mkl
*/
template<class T>
class MKLCGBase{
public:
	MKLCGBase(){
		if( !TypeGuard() ){
			std::cout << "error: invalid TYPE passed to MKLCGBase\n\
						 MKL CG only accept double" << std::endl;
			exit(0);
		}

		ipar = NULL;
		dpar = NULL;
		tmp = NULL;

		Reset();
	}

	~MKLCGBase(){
		Reset();

		//	this function may cause error, remove it if you don't care memory leak
		MKL_FreeBuffers();
	}

	/**
		Reset to startup status
	*/
	inline void Reset(){
		symmetric_flag	= 0;
		dims			= 0;
		nnz				= 0;
		value_ptr		= NULL;
		col_id_ptr		= NULL;
		row_start_ptr	= NULL;
		x_ptr			= NULL;
		b_ptr			= NULL;
		nrhs			= 0;
		x_tmp_ptr		= NULL;

		if(ipar){
			delete ipar;
			ipar = NULL;
		}
		if(dpar){
			delete dpar;
			dpar = NULL;
		}
		if(tmp){
			delete tmp;
			tmp = NULL;
		}

		stage = CG_NONE;
	}

	/**
		Setup Square Matrix A

		\param	dims			dimension of the matrix, must be square matrix
		\param	nnz				number of non-zero elements
		\param	value			pointer to the matrix element array
		\param	col_id			pointer to colume id array of each element
		\param	row_start		pointer to row start id array of each row
		\param	symmetric_flag	0: non-symmetric, 1: symmetric
		\return					whether success
	*/
	int SetupMatrix(int dims, int nnz, T* value, int* col_id, int* row_start, int symmetric_flag=0){
		//	stage check
		if( stage != CG_NONE ){
			#ifndef BE_QUIET
				std::cout << "error: MKLCGBase not reset, cannot setup matrix" << std::endl;
			#endif
			return 0;
		}

		//	copy matrix
		this->symmetric_flag	= symmetric_flag;
		this->dims		= dims;
		this->nnz		= nnz;
		value_ptr		= value;
		col_id_ptr		= col_id;
		row_start_ptr	= row_start;

		stage = CG_SETUP_MATRIX;

		return 1;
	}

	/**
		Setup Matrix A from CSR sparse matrix structure directly
	*/
	int SetupMatrix(SparseMatrixCSR<T>& spa_mat){
		if( spa_mat.row_num != spa_mat.col_num ){
			#ifndef BE_QUIET
				std::cout << "error: MKLCGBase::SetupMatrix, matrix must be square matrix" << std::endl;
			#endif
			return 0;
		}

		int ret = SetupMatrix(spa_mat.row_num, spa_mat.NNZ(),
			&spa_mat.value[0], &spa_mat.col_id[0], &spa_mat.row_start[0], 
			spa_mat.IsSymmetric());
		return ret;
	}

	/**
		Setup x and b

		this function will setup MKL resources and call initialize function
		of CG in mkl. InitParameters() is virtual function
	*/
	int SetupVectors(T* x, T* b, int vec_num = 1){
		//	stage check
		if( stage != CG_SETUP_MATRIX ){
			#ifndef BE_QUIET
				std::cout << "error: MKLCGBase::SetupVectors, matrix not set, cannot setup vector" << std::endl;
			#endif
			return 0;
		}

		//	if result x is write back to the array b, then we must create a temp array
		//	currently, we cannot handle situation of memory overlapping
		if( x == b )
			x_tmp_ptr = new T[dims*vec_num];

		x_ptr	= (x_tmp_ptr? x_tmp_ptr : x);
		b_ptr	= b;
		nrhs	= vec_num;

		//	setup resources
		if( !SetupMKLResources() )
			return 0;

		stage = CG_SETUP_VECTOR;

		//	initialize ipar and dpar
		if( !InitParameters() )
			return 0;

		return 1;
	}

	/**
		set parameters of iterations.

		you can also set ipar and dpar directly instead of calling this function

		\param	max_iterations		max iterations
		\param	tolerance			residual tolerance
		\return						whether succeed
	*/
	int SetParameters(int max_iterations=150, T tolerance=1e-6){
		//	stage check
		if( stage != CG_SETUP_VECTOR ){
			#ifndef BE_QUIET
				std::cout << "error: MKLCGBase::SetParameters, vector not set, cannot set parameters" << std::endl;
			#endif
			return 0;
		}

		//	iterations
		ipar[4] = max_iterations;
		ipar[7] = 1;

		//	tolerance
		dpar[0] = tolerance;
		ipar[8] = 1;

		//	disable user defined stopping test
		ipar[9] = 0;

		stage = CG_SET_PARAM;

		return 1;
	}

	/**
		solve the system. 

		\param	guess	initual value, 0: default value; 1: all zero; 2: same as b
		\return			whether succeed
	*/
	virtual int Solve(int guess = 1){
		//	stage check
		if( stage!=CG_SETUP_VECTOR && stage!=CG_SET_PARAM ){
			#ifndef BE_QUIET
				std::cout << "error: MKLCGBase::Solve, vector not set, cannot solve" << std::endl;
			#endif
			return 0;
		}

		//	initual guess of solution
		if (guess == 1)		//	all zero
			std::fill(&x_ptr[0], &x_ptr[dims * nrhs], 0); //memset(x_ptr, 0, sizeof(T)*dims*nrhs);
		else if( guess == 2 )	//	value of b
			memcpy(x_ptr, b_ptr, sizeof(T)*dims*nrhs);

		//	solve using mkl
		if( !MKLCGSolve() )
			return 0;

		//	write back data
		if( x_tmp_ptr ){
			memcpy(b_ptr, x_ptr, sizeof(T)*dims*nrhs);
			delete x_tmp_ptr;
			x_tmp_ptr = NULL;
		}

		return 1;
	}

	/**
		Solve the linear system, Ax = b

		combime SetupVectors() and Solve(). If you call this function, 
		you cannot set parameters youself. So if you don't want to use
		default parameters, call SetupVectors(), SetParameters() (optional),
		Solve() step by step.

		\param	x			x in Ax = b
		\param	b			b in Ax = b
		\param	vec_num		how many vectors contained in the array\n
							in CG solver, only 1 is allowed
		\param	guess		initual guess of x, 0: value in x; 1: 0; 2: b
		\return				whether success
	*/
	int Solve(T* x, T* b, int vec_num = 1, int guess = 1){

		if( !SetupVectors(x, b, vec_num) )
			return 0;

		if( !Solve(guess) )
			return 0;

		return 1;
	}

public:
	//	mkl cg variables to use
	MKL_INT		rci_request, itercount;
	MKL_INT*	ipar;
	T*			dpar;
	T*			tmp;

protected:
	/**
		Check the type of this class

		MKL only accept double type, so if user try to pass
		other type just let the constructor fail
	*/
	inline int TypeGuard(){
		return 0;
	}

	/**
		setup arrays used by MKL
	*/
	int SetupMKLResources(){
		int size = (nrhs==1 ? 128 : 128+2*nrhs);

		ipar	= new MKL_INT[size];
		dpar	= new T[size];
		tmp		= new T[dims * (3+nrhs)];

		if( !ipar || !dpar || !tmp ){
			std::cout << "error: MKLCGBase::SetupMKLResources, alloc memory failed" << std::endl;
			return 0;
		}

		return 1;
	}

	/**
		Initialize parameters of ipar and dpar

		\return		whether initialize succeed
	*/
	virtual int InitParameters() = 0;

	/**
		Solve function using mkl

		\return			whether solve succeed
	*/
	virtual int MKLCGSolve() = 0;

protected:
	//	Matrix A in CSR format
	int		symmetric_flag;			///<	0: non-symmetric, 1: symmetric. Only upper triangular is allowed
	int		dims;					///<	matrix dimension, must be square matrix
	int		nnz;					///<	number of non-zero elements
	T*		value_ptr;				///<	non-zero elements
	int*	col_id_ptr;				///<	non-zero element column index
	int*	row_start_ptr;			///<	start position of each row

	//	vectors
	T*		x_ptr;
	T*		b_ptr;
	int		nrhs;					///<	how many vectors are solved one time
	T*		x_tmp_ptr;				///<	if x, b overlap, use use this tmp array

	//	solver stage flag
	int	stage;
	enum{ CG_NONE, CG_SETUP_MATRIX, CG_SETUP_VECTOR, CG_SET_PARAM, CG_SOLVE };
};

template<> inline int MKLCGBase<double>::TypeGuard(){
	return 1;
}


/**
	mkl CG solver for symmetric matrix with single right hand side

	To use this solve, call SetupMatrix(), Solve() in sequence. In this way,
	default control parameters are used.

	If you want to set parameters directly, call SetupMatrix(),
	SetupVectors(), SetParameters() and Solve() step by step. If you want 
	to control MKL directly, replace SetParameters() by write ipar and dpar
	directly
*/
template<class T>
class MKLCGSolver : public MKLCGBase<T>{

protected:
	/**
		initialize
	*/
	int InitParameters(){
		//	initialize solver
		dcg_init(&dims, x_ptr, b_ptr, &rci_request, ipar, dpar, tmp);
		if( rci_request != 0 ){
			std::cout << "error: dcg_init in MKLCGSolver" << std::endl;
			return 0;
		}

		//	disable user specified end condition
		ipar[9] = 0;

		return 1;
	}

	/**
		the core of solving, check, cg and get
	*/
	int MKLCGSolve(){
		//	check
		dcg_check(&dims, x_ptr, b_ptr, &rci_request, ipar, dpar, tmp);
		if( rci_request != 0 ){
			std::cout << "error: dcg_check in MKLCGSolver, error code: " << rci_request << std::endl;
			return 0;
		}

		//	iterative solve
		dcg(&dims, x_ptr, b_ptr, &rci_request, ipar, dpar, tmp);

		while( rci_request ){
			if( rci_request == 1 ){
				//	multiply A
				if( symmetric_flag )
					mkl::c_mkl_cspblas_dcsrsymv('U', dims, value_ptr, row_start_ptr, col_id_ptr, tmp, &tmp[dims]);
				else
					mkl::c_mkl_cspblas_dcsrgemv('N', dims, value_ptr, row_start_ptr, col_id_ptr, tmp, &tmp[dims]);

				//	CG again
				dcg(&dims, x_ptr, b_ptr, &rci_request, ipar, dpar, tmp);

			}
			else{
				std::cout << "error: dcg in MKLCGSolver, error code: " << rci_request << std::endl;
				return 0;
			}
		}

		//	get result
		dcg_get(&dims, x_ptr, b_ptr, &rci_request, ipar, dpar, tmp, &itercount);

		return 1;
	}

};

typedef MKLCGSolver<double>	MKLCGSolverd;

/**
	mkl CG solver for symmetric matrix with multiple right hand side

	To use this solve, call SetupMatrix(), Solve() in sequence. In this way,
	default control parameters are used.

	If you want to set parameters directly, call SetupMatrix(),
	SetupVectors(), SetParameters() and Solve() step by step. If you want 
	to control MKL directly, replace SetParameters() by write ipar and dpar
	directly*/
template<class T>
class MKLCGMRHSSolver : public MKLCGBase<T>{

protected:
	/**
		initialize
	*/
	int InitParameters(){
		int method = 1;
		//	initialize solver
		dcgmrhs_init(&dims, x_ptr, &nrhs, b_ptr, &method, &rci_request, ipar, dpar, tmp);
		if( rci_request != 0 ){
			std::cout << "error: dcgmrhs_init in MKLCGMRHSSolver" << std::endl;
			return 0;
		}

		//	disable user specified end condition
		ipar[9] = 0;

		return 1;
	}

	/**
		the core of solving, check, cg and get
	*/
	int MKLCGSolve(){
		//	check
		dcgmrhs_check(&dims, x_ptr, &nrhs, b_ptr, &rci_request, ipar, dpar, tmp);
		if( rci_request != 0 ){
			std::cout << "error: dcgmrhs_check in MKLCGMRHSSolver, error code: " << rci_request << std::endl;
			return 0;
		}

		//	iterative solve
		dcgmrhs(&dims, x_ptr, &nrhs, b_ptr, &rci_request, ipar, dpar, tmp);

		while( rci_request ){
			if( rci_request == 1 ){
				//	multiply A
				if( symmetric_flag )
					mkl::c_mkl_cspblas_dcsrsymv('U', dims, value_ptr, row_start_ptr, col_id_ptr, tmp, &tmp[dims]);
				else
					mkl::c_mkl_cspblas_dcsrgemv('N', dims, value_ptr, row_start_ptr, col_id_ptr, tmp, &tmp[dims]);

				//	CG again
				dcgmrhs(&dims, x_ptr, &nrhs, b_ptr, &rci_request, ipar, dpar, tmp);

			}
			else{
				std::cout << "error: dcg in MKLCGSolver, error code: " << rci_request << std::endl;
				return 0;
			}
		}

		//	get result
		dcgmrhs_get(&dims, x_ptr, &nrhs, b_ptr, &rci_request, ipar, dpar, tmp, &itercount);

		return 1;
	}

};

typedef MKLCGMRHSSolver<double>	MKLCGMRHSSolverd;


}	//	end namespace yz

#endif	//	__YZ_MKL_ITER_SOLVER_H__