/***********************************************************/
/**	\file
	\brief		Array Check Functions
	\author		Yizhong Zhang
	\date		9/28/2012
*/
/***********************************************************/
#ifndef __YZ_ARRAY_CHECK_H__
#define __YZ_ARRAY_CHECK_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>
#include <string.h>
#include "yzLib/yz_utils/yz_hash.h"

namespace yz{ namespace utils{

//	========================================
///@{
/**	@name Array Check Functions
*/
//	========================================

/**
	check whether one array is bijection with another array

	bijection means each element in data array is unique, ranging
	from first_index to size-first_index-1

	This function can be used to assert an array, so we set print_flag

	\param	data			the array to check
	\param	size			number of elements
	\param	first_index		first index of the projection array
	\param	print_flag		whether print not bijection information
	\return					whether is bijection
*/
inline int isBijection(const int* data, int size, int first_index = 0, int print_flag = 0){
	//	make sure size is legal
	if( size <= 0 ){
		std::cout << "error: isBijection, data size illegal" << std::endl;
		return 0;
	}

	//	check the mapping
	int* mapping = new int[size];
	memset(mapping, 0, sizeof(int)*size);

	for(int i=0; i<size; i++){
		//	check data in legal range
		if( data[i] < first_index || data[i] >= size+first_index ){
			if( print_flag ){
				std::cout << "data is not bijection, data[" << i << "] = " << data[i] <<
					"range is not in [" << first_index << ", " << size+first_index-1 << "]"  << std::endl;
			}
			delete[] mapping;
			return 0;
		}
		else
			mapping[data[i]] ++;
	}

	//	check whether mapping is complete
	for(int i=0; i<size; i++){
		if( mapping[i] != 1 ){
			if( mapping[i] == 0 ){
				if( print_flag )
					std::cout << "data is not bijection, " << i << " is never mapped" << std::endl;
				delete[] mapping;
				return 0;
			}
			else{
				if( print_flag )
					std::cout << "data is not bijection, " << i << " is mapped " << mapping[i] << " times" << std::endl;
				delete[] mapping;
				return 0;
			}
		}
	}

	delete[] mapping;

	return 1;
}

/**
	check whether one array is bijection with another array

	bijection means each element in data array is unique, ranging
	from first_index to size-first_index-1

	This function can be used to assert an array, so we set print_flag

	\param	data			the array to check
	\param	first_index		first index of the projection array
	\param	print_flag		whether print not bijection information
	\return					whether is bijection
*/
inline int isBijection(const std::vector<int>& data, int first_index = 0, int print_flag = 0){

	return isBijection((int*)&data[0], int(data.size()), first_index, print_flag);
}

///@}

//	========================================
///@{
/**	@name Array Util Functions
*/
//	========================================

/**
	Find the elements that appear most in data array

	This function use hash to detect duplicates. Data are treated the same only if they are the same of all bits in data array

	\param	most_element	return the element that appear most in data array
	\param	data			the data array
	\param	data_number		number of data in the array
	\return					number of appearance of the most_element
*/
template<class T>
inline int getMostAppearance(T& most_element, const T* data, int data_number) {
	if (data_number < 1)
		return 0;

	//	create hash. first: data value as key; second: duplicates of this data value
	std::unordered_map<T, unsigned int, utils::BitwiseHasher<T>> hash_map;
	hash_map.reserve(data_number);

	//	parse all data
	for (int i = 0; i < data_number; i++) {
		//std::unordered_map<T, unsigned int, utils::BitwiseHasher<T>>::iterator iter = hash_map.find(data[i]);
		auto iter = hash_map.find(data[i]);
		if (iter == hash_map.end()) {	//	data not exist, insert into hash
			hash_map.insert(std::pair<T, unsigned int>(data[i], 1));
		}
		else {		//	data already exist, increase the count
			iter->second++;
		}
	}

	//	parse tha hash table and get most appearence
	int most_amount = 0;
	//for (std::unordered_map<T, unsigned int, utils::BitwiseHasher<T>>::iterator iter = hash_map.begin(); iter != hash_map.end(); iter++) {
	for (auto iter = hash_map.begin(); iter != hash_map.end(); iter++) {
		if (iter->second > most_amount) {
			most_amount = iter->second;
			most_element = iter->first;
		}
	}

	return most_amount;
}



///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_ARRAY_CHECK_H__