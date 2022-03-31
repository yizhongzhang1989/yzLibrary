/***********************************************************/
/**	\file
	\brief		Sparse Vector
	\details	
	\author		Yizhong Zhang
	\date		9/6/2012
*/
/***********************************************************/
#ifndef __YZ_SPARSE_VECTOR_H__
#define __YZ_SPARSE_VECTOR_H__

#include <assert.h>
#include <iostream>
#include <vector>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
//	virtual base for sparse vector
//	========================================
template<class T>
class BaseSparseVector{
public:
	int dim_num;

public:
	BaseSparseVector(int dims = 0){
		assert( dims >= 0 );
		dim_num = dims;
	}

	virtual void SetDimension(int dims) = 0;
};


//	========================================
//	sparse vector stored by index
//	========================================
template<class T>
class SparseVector : public BaseSparseVector<T>{
public:
	using BaseSparseVector<T>::dim_num;

	std::vector<T>		value;
	std::vector<int>	index;

public:
	/**
		constructor, if dimension set, setup storage
	*/
	SparseVector(int dims = 0) : BaseSparseVector<T>(dims) {
		assert(dims >= 0);
		if( dims > 0 )
			value.resize( dims );
	}

	/**
		Set dimension number
	*/
	void SetDimension(int dims){
		assert(dims > 0);
		dim_num = dims;
		if( dim_num != value.size() )
			value.resize( dim_num );
	}

};



}	//	end namespace yz

#endif	//	__YZ_SPARSE_VECTOR_H__