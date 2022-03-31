/***********************************************************/
/**	\file
	\brief		Statistics
	\author		Yizhong Zhang
	\date		6/11/2012
*/
/***********************************************************/
#ifndef __YZ_STATISTICS_BASIC_H__
#define __YZ_STATISTICS_BASIC_H__

#include <vector>
#include "yzLib/yz_math/yz_statistics/yz_kernel_functions.h"

namespace yz{

//	========================================
///@{
/**	@name	Mean
*/
//	========================================

/**
	calculate mean value

	\param	data_ptr		data list
	\param	data_number		number of data
	\return					the mean of data, meaningless if data_number <= 0
*/
template<typename T>
inline T calculateMeanValue(const T* data_ptr, int data_number){
	if( data_number <= 0 ){
		T zero;	//	we cannot guarentee the return value is zero, we just return a meaningless number
		return zero;
	}
	T sum = data_ptr[0];
	for(int i=1; i<data_number; i++){
		sum += data_ptr[i];
	}
	return sum *= (1.0/data_number);
}

/**
	calculate mean value

	\param	data	data 
	\return			the mean of data, meaningless if data_number <= 0
*/
template<typename T>
inline T calculateMeanValue(const std::vector<T>& data){
	if( data.empty() ){
		T zero;	//	we cannot guarentee the return value is zero, we just return a meaningless number
		return zero;
	}
	return calculateMeanValue((T*)&data[0], data.size());
}

///@}

}	//	namespace yz

#endif	//	__YZ_STATISTICS_H__