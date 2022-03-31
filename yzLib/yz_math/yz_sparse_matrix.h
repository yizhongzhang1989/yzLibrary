/***********************************************************/
/**	\file
	\brief		Sparse Matrix in 0-based indexing
	\details	Sparse matrix don't have fixed storage space,
				so we use vector to hold the data. \n
				Since sparse matrix operation is much complex
				than dense matrix, so the class itself doesn't 
				provide too many operations. You can treate it 
				as a data container. \n
				We do some check in debug mode, but no check 
				in release mode for performance. It is the
				responsibility of the user to make sure the code
				doesn't not make mistakes
	\author		Yizhong Zhang
	\date		9/6/2012
*/
/***********************************************************/
#ifndef __YZ_SPARSE_MATRIX_H__
#define __YZ_SPARSE_MATRIX_H__

#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_reorder.h"
#ifdef YZ_mkl_spblas_h	//	if mkl_spblas.h included, we can use mkl functions
#	include "yzLib/yz_math/yz_mkl_spblas_c.h"
#endif

namespace yz{

//	========================================
///@{
/**	@name compare functions used in sort
*/
//	========================================
inline bool smIsSmallerInt3xy(int3 v1, int3 v2){
	return v1.x<v2.x ? true : ( v1.x>v2.x ? false : (v1.y<v2.y) );
}

///@}

//	========================================
/**
	virtual base for sparse matrix

	All sparse matrix should be dirived from this class
*/
//	========================================
template<class T>
class SparseMatrixBase{
public:
	int row_num;	///<	rows of the matrix
	int	col_num;	///<	columns of the matrix
	int	indexing;	///<	0: zero-based indexing, 1: one-based indexing

public:

	/**
		Constructor, specify dimension of the matrix

		\param	m		rows
		\param	n		columns
		\param	idxing	indexing of the matrix, 0: zero-based indexing, 1: one-based indexing
	*/
	SparseMatrixBase(int m = 0, int n = 0, int idxing = 0){
		assert((m==0 && n==0) || (m>0 && n>0));	//	we only accept default, or legal dimension value
		row_num		= m;
		col_num		= n;
		indexing	= idxing;
	}

	//	========================================
	///@{
	/**	@name Functions must be derived
	*/
	//	========================================
	/**
		Clear the matrix, ready to get new matrix
	*/
	virtual inline void Reset() = 0;

	/**
		Set the dimension of the matrix explicitly
	*/
	virtual inline void SetDimension(int m, int n) = 0;

	/**
		Set indexing of the matrix
		
		0: zero-based indexing, non-zero: one-based indexing
	*/
	virtual inline void SetIndexing(int idxing) = 0;

	/**
		insert a single element to the matrix
	*/
	virtual inline void InsertElement(T val, int m, int n) = 0;

	/**
		whether the derived class is symmetric structure

		symmetric matrix can be stored as upper triangular 
		matrix to minimize storage. This function only check
		whether the class is symmetric structure, not whether
		the matrix is symmetric.
	*/
	virtual inline int IsSymmetric() const = 0;

	///@}

	/**
		reorder the elements

		if the storage is ordered, then no need to derive this function
	*/
	virtual inline void Reorder(){}
};


//	========================================
/**
	coordinate format sparse matrix
*/
//	========================================
template<class T>
class SparseMatrixCoo : public SparseMatrixBase<T>{
public:
	using SparseMatrixBase<T>::row_num;
	using SparseMatrixBase<T>::col_num;
	using SparseMatrixBase<T>::indexing;

	std::vector<T>		value;		//	size of value, row_id, col_id should be identical
	std::vector<int>	row_id;		//	we don't use int3, because libraries for sparse
	std::vector<int>	col_id;		//	matrix usually accept independent arrays

public:
	SparseMatrixCoo(int m=0, int n=0) : SparseMatrixBase<T>(m, n){}

	/**
		Reset the matrix, clear everything
	*/
	inline void Reset(){
		row_num		= 0;
		col_num		= 0;
		indexing	= 0;

		value.clear();
		row_id.clear();
		col_id.clear();
	}

	/**
		set the dimension of the matrix

		we can only use this function when the dimension of the 
		matrix is not set yet

		\param	m		row number
		\param	n		column number
	*/
	inline void SetDimension(int m, int n){
		assert(row_num==0 && col_num==0);	//	dimension is not set
		assert(m>0 && n>0);					//	dimension is legal
		row_num = m;
		col_num = n;
	}

	/**
		Set indexing of the matrix
		
		0: zero-based indexing, non-zero: one-based indexing
	*/
	inline void SetIndexing(int idxing){
		if( !indexing && idxing ){		//	change from zero to one
			for(int i=0; i<row_id.size(); i++)
				row_id[i] ++;
			for(int i=0; i<col_id.size(); i++)
				col_id[i] ++;
			indexing = 1;
		}
		else if( indexing && !idxing ){//	change from one to zero
			for(int i=0; i<row_id.size(); i++)
				row_id[i] --;
			for(int i=0; i<col_id.size(); i++)
				col_id[i] --;
			indexing = 0;
		}
		//	otherwise, no need to change current indexing
	}
	/**
		insert a single element to the matrix

		\param	val		value of element to insert
		\param	m		row index
		\param	n		column index
	*/
	inline void InsertElement(T val, int m, int n){
		assert(m>=0 && m<row_num && n>=0 && n<col_num);	//	legal position
		assert( !CheckExist(m, n) );					//	the element doesn't exist
		value.push_back(val);
		row_id.push_back(m);
		col_id.push_back(n);
	}


	/**
		SparseMatrixCoo is not symmetric structure
	*/
	inline int IsSymmetric() const{
		return 0;
	}

	/**
		Sort the elements in row major format.

		Some libraries require the data to be sorted in row-major,
		but we don't confirm this when insert elements. Then we can
		call this function to sort the data in row-major order.
	*/
	inline void SortRowMajor(){
		assert( value.size()==row_id.size() && value.size()==col_id.size() );
		if( value.size() < 1 )	return;

		//	create temp array to sort the elements
		std::vector<int3>	tmp;
		tmp.resize(value.size());
		for(int i=0; i<value.size(); i++){
			tmp[i].x = row_id[i];
			tmp[i].y = col_id[i];
			tmp[i].z = i;			//	z record its original position
		}

		std::sort(tmp.begin(), tmp.end(), smIsSmallerInt3xy);

		//	reset old array
		std::vector<T>	tmp_value;
		tmp_value.resize(value.size());
		for(int i=0; i<tmp_value.size(); i++){
			row_id[i]	= tmp[i].x;
			col_id[i]	= tmp[i].y;
			tmp_value[i]= value[ tmp[i].z ];
		}

		value.swap( tmp_value );
	}

	/**
		Default reorder elements function: sort row jamor
	*/
	inline void Reorder(){
		SortRowMajor();
	}

	/**
		return number of non-zero elements
	*/
	inline int NNZ() const{
		return value.size();
	}

protected:
	/**
		check whether an element already exist
	*/
	inline int CheckExist(int m, int n){
		assert(row_id.size() == col_id.size());
		for( int i=0; i<row_id.size(); i++ )
			if( row_id[i]==m && col_id[i]==n )
				return 1;
		return 0;
	}



};

//	========================================
/**
	coordinate format symmetric sparse matrix

	This class restricts coordinate format sparse matrix
	to be symmetric. It must be square matrix. Data are 
	stored in upper triangular format.
*/
//	========================================
template<class T>
class SparseMatrixCooSym : public SparseMatrixCoo<T>{
public:
	/**
		constructor, just one dimension
	*/
	SparseMatrixCooSym(int n=0) : SparseMatrixCoo<T>(n, n) {}

	/**
		set the dimension of the square matrix

		we can only use this function when the dimension of the 
		matrix is not set yet

		\param	n		dimension number
	*/
	inline void SetDimension(int n){
		SparseMatrixCoo<T>::SetDimension(n, n);
	}

	/**
		insert a single element to the matrix

		As the matrix is symmetric, we insert to upper triangular
		position only. 

		\param	val		value of element to insert
		\param	m		row index
		\param	n		column index
	*/
	inline void InsertElement(T val, int m, int n){
		if(m > n)	mySwap(m, n);
		SparseMatrixCoo<T>::InsertElement(val, m, n);
	}

	/**
		SparseMatrixCooSym is symmetric sturcture
	*/
	inline int IsSymmetric() const{
		return 1;
	}
};


//	========================================
/**
	compressed sparse row format sparse matrix (CSR)
*/
//	========================================
template<class T>
class SparseMatrixCSR : public SparseMatrixBase<T>{
public:
	using SparseMatrixBase<T>::row_num;
	using SparseMatrixBase<T>::col_num;
	using SparseMatrixBase<T>::indexing;

	std::vector<T>		value;
	std::vector<int>	col_id;
	std::vector<int>	row_start;	///< start index of each row

public:
	SparseMatrixCSR(int m=0, int n=0) : SparseMatrixBase<T>(m, n){
		row_start.resize(m+1, 0);
	}

	/**
		Reset the matrix, clear everything
	*/
	inline void Reset(){
		row_num		= 0;
		col_num		= 0;
		indexing	= 0;

		value.clear();
		col_id.clear();
		row_start.clear();
	}

	/**
		set the dimension of the matrix

		we can only use this function when the dimension of the 
		matrix is not set yet

		\param	m		row number
		\param	n		column number
	*/
	inline void SetDimension(int m, int n){
		assert(row_num==0 && col_num==0);	//	dimension is not set
		assert(m>0 && n>0);					//	dimension is legal
		row_num = m;
		col_num = n;
		row_start.resize(m+1, 0);
	}

	/**
		Set indexing of the matrix
		
		0: zero-based indexing, non-zero: one-based indexing
	*/
	inline void SetIndexing(int idxing){
		if( !indexing && idxing ){		//	change from zero to one
			for(int i=0; i<col_id.size(); i++)
				col_id[i] ++;
			for(int i=0; i<row_start.size(); i++)
				row_start[i] ++;
			indexing = 1;
		}
		else if( indexing && !idxing ){//	change from one to zero
			for(int i=0; i<col_id.size(); i++)
				col_id[i] --;
			for(int i=0; i<row_start.size(); i++)
				row_start[i] --;
			indexing = 0;
		}
		//	otherwise, no need to change current indexing
	}
	/**
		insert a single element to the matrix

		the storage of the matrix is ordered after insertion

		\param	val		value of element to insert
		\param	m		row index
		\param	n		column index
	*/
	inline void InsertElement(T val, int m, int n){
		assert(m>=0 && m<row_num && n>=0 && n<col_num);	//	legal position
		assert( !CheckExist(m, n) );					//	the element doesn't exist
		//	search insert position
		int insert_pos;
		for(insert_pos = row_start[m]; insert_pos < row_start[m+1]; insert_pos++)
			if( col_id[insert_pos] > n )	break;

		//	insert into list
		value.insert(value.begin()+insert_pos, val);
		col_id.insert(col_id.begin()+insert_pos, n);
		for( int i=m+1; i<=row_num; i++ ){
			row_start[i] ++;
		}
	}

	/**
		SparseMatrixCoo is not symmetric structure
	*/
	inline int IsSymmetric() const{
		return 0;
	}

	/**
		Sort each row increasing colume index
	*/
	inline void Reorder(){
		for(int i=0; i<row_num; i++){
			int* begin	= &col_id[0] + row_start[i];
			int* end	= &col_id[0] + row_start[i+1];
			utils::sortAndReorder(begin, end, &value[0] + row_start[i]);
		}
	}
	/**
		return number of non-zero elements
	*/
	inline int NNZ() const{
		return value.size();
	}
	/**
		Return transposed matrix, but don't change old data
	*/
	inline SparseMatrixCSR Transpose() const{
		SparseMatrixCSR result(*this);
		result.SetTranspose();
		return result;
	}
	/**
		Set the matrix to be its transpose
	*/
	inline SparseMatrixCSR& SetTranspose(){
#ifdef YZ_mkl_spblas_h
		if(Is_float<T>::check_type && row_num==col_num){		//	T is float, the matrix is square matrix
			SparseMatrixCSR<T> csc = *this;
			int job[6] = {1, bool(indexing), 0, 0, 0, 1};
			int info;
			mkl::c_mkl_scsrcsc(job, row_num, 
				(float*)&value[0], (int*)&col_id[0], (int*)&row_start[0],
				(float*)&csc.value[0], (int*)&csc.col_id[0], (int*)&csc.row_start[0],
				&info);
			return *this;
		}
		else if(Is_double<T>::check_type && row_num==col_num){	//	T is double, the matrix is square matrix
			SparseMatrixCSR<T> csc = *this;
			int job[6] = {1, bool(indexing), 0, 0, 0, 1};
			int info;
			mkl::c_mkl_dcsrcsc(job, row_num, 
				(double*)&value[0], (int*)&col_id[0], (int*)&row_start[0],
				(double*)&csc.value[0], (int*)&csc.col_id[0], (int*)&csc.row_start[0],
				&info);
			return *this;
		}
#endif
		std::cout << "SparseMatrixCSR::SetTranspose not implemented, "
			<< "try to include mkl_spblas.h and use float or double" << std::endl;
		return *this;
	}
protected:
	/**
		check whether an element already exist
	*/
	inline int CheckExist(int m, int n){
		assert(row_start.size() == row_num+1);
		assert(m>=0 && m<row_num);
		for(int i=row_start[m]; i<row_start[m+1]; i++){
			assert(i>=0 && i<col_id.size());
			if( col_id[i] == n )
				return 1;
		}
		return 0;
	}
};

//	========================================
/**
	CSR symmetric sparse matrix

	This class restricts coordinate format sparse matrix
	to be symmetric. It must be square matrix. Data are 
	stored in upper triangular format.
*/
//	========================================
template<class T>
class SparseMatrixCSRSym : public SparseMatrixCSR<T>{
public:
	/**
		constructor, just one dimension
	*/
	SparseMatrixCSRSym(int n=0) : SparseMatrixCSR<T>(n, n) {}

	/**
		set the dimension of the square matrix

		we can only use this function when the dimension of the 
		matrix is not set yet

		\param	n		dimension number
	*/
	inline void SetDimension(int n){
		SparseMatrixCSR<T>::SetDimension(n, n);
	}

	/**
		insert a single element to the matrix

		As the matrix is symmetric, we insert to upper triangular
		position only. 

		\param	val		value of element to insert
		\param	m		row index
		\param	n		column index
	*/
	inline void InsertElement(T val, int m, int n){
		if(m > n)	mySwap(m, n);
		SparseMatrixCSR<T>::InsertElement(val, m, n);
	}

	/**
		SparseMatrixCSRSym is symmetric sturcture
	*/
	inline int IsSymmetric() const{
		return 1;
	}
};



}	//	end namespace yz

#endif	//	__YZ_SPARSE_MATRIX_H__