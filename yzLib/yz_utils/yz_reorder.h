/***********************************************************/
/**	\file
	\brief		Reorder
	\details	This file contain reorder functions. including 
				sort, arbitrary order
	\author		Yizhong Zhang
	\date		9/28/2012
*/
/***********************************************************/
#ifndef __YZ_REORDER_H__
#define __YZ_REORDER_H__

#include <vector>
#include <unordered_set>
#include "yzLib/yz_utils/yz_array_check.h"
#include "yzLib/yz_utils/yz_hash.h"


namespace yz{ namespace utils{

template<typename T>
int reorderDataBySourceOrder(T* data, int size, const int* order);

//	========================================
///@{
/**	@name Sort Functions
*/
//	========================================


/**
	sort the index of a array by operator<

	This function is an extension of std::sort(), but instead of
	sort the array, we sort the index only. After calling this 
	function, each element in index array indicates the original 
	index in the given data array (source order).

	\param	_index		return index table
	\param	_First		start of data
	\param	_Last		end of data
	\author				Dongping Li
*/
template<class _RanIt> 
inline void sortIndex(int* _index, _RanIt _First, _RanIt _Last){
	struct SortServer{
		_RanIt	m_begin;

		SortServer(const _RanIt begin) 
			: m_begin(begin){}

		bool operator()(int a, int b) const{
			return m_begin[a] < m_begin[b];
		}
	};

	int n = _Last - _First;
	for(int i=0; i<n; i++)
		_index[i] = i;

	std::sort(_index, _index+n, SortServer(_First));
}

/**
	sort the index of a array by operator<

	This function is an extension of std::sort(), but instead of
	sort the array, we sort the index only. After calling this 
	function, each element in index array indicates the original 
	index in the given data array (source order).

	\param	index		return index table
	\param	_First		start of data
	\param	_Last		end of data
	\author				Dongping Li
*/
template<class _RanIt> 
inline void sortIndex(std::vector<int>& index, _RanIt _First, _RanIt _Last){
	int n = _Last - _First;
	if( n <= 0 )		//	no data
		return;

	index.resize(n);
	sortIndex(&index[0], _First, _Last);
}

/**
	sort the index of a array by given comparison function

	This function is an extension of std::sort(), but instead of
	sort the array, we sort the index only. After calling this 
	function, each element in index array indicates the original 
	index in the given data array (source order).

	\param	_index		return index table
	\param	_First		start of data
	\param	_Last		end of data
	\param	_Pred		comparison function
	\author				Dongping Li
*/
template<class _RanIt, class _Pr> 
inline void sortIndex(int* _index, _RanIt _First, _RanIt _Last, _Pr _Pred){
	struct SortServer{
		_RanIt	m_begin;
		_Pr		m_less;

		SortServer(const _RanIt begin, _Pr less) 
			: m_begin(begin), m_less(less){}

		bool operator()(int a, int b) const{
			return m_less(m_begin[a], m_begin[b]);
		}
	};

	int n = _Last - _First;
	for(int i=0; i<n; i++)
		_index[i] = i;

	std::sort(_index, _index+n, SortServer(_First, _Pred));
}

/**
	sort the index of a array by given comparison function

	This function is an extension of std::sort(), but instead of
	sort the array, we sort the index only. After calling this 
	function, each element in index array indicates the original 
	index in the given data array (source order).

	\param	index		return index table
	\param	_First		start of data
	\param	_Last		end of data
	\param	_Pred		comparison function
	\author				Dongping Li
*/
template<class _RanIt, class _Pr> 
inline void sortIndex(std::vector<int>& index, _RanIt _First, _RanIt _Last, _Pr _Pred){
	index.resize(_Last - _First);
	if ( index.size() < 2 )
		return;
	sortIndex(&index[0], _First, _Last, _Pred);
}

/**
	sort the array by operator<, then reorder another arrays using
	the same reorder index. 

	\param	_First		start of data
	\param	_Last		end of data
	\param	_Array2		the second array
*/
template<class T, class T2> 
inline void sortAndReorder(T*	_First, 
						   T*	_Last, 
						   T2*	_Array2){
	int n = _Last - _First;
	if( n < 2 )	return;
	int* index = new int[n];

	sortIndex(index, _First, _Last);
	reorderDataBySourceOrder(_First, n, index);
	reorderDataBySourceOrder(_Array2, n, index);

	delete[] index;
}

/**
	sort the array by given comparison function, then reorder another array
	using the same reorder index. 

	\param	_First		start of data
	\param	_Last		end of data
	\param	_Pred		the comparison function
	\param	_Array2		the second array
*/
template<class T, class _Pr, class T2> 
inline void sortAndReorder(T*	_First, 
						   T*	_Last, 
						   _Pr	_Pred,
						   T2*	_Array2){
	int n = _Last - _First;
	if( n < 2 )	return;
	int* index = new int[n];

	sortIndex(index, _First, _Last, _Pred);
	reorderDataBySourceOrder(_First, n, index);
	reorderDataBySourceOrder(_Array2, n, index);

	delete[] index;
}

///@}

//	========================================
///@{
/**	@name Reorder Functions
*/
//	========================================

/**
	invert the order from target order to source order, and 
	vice versa
*/
inline void invertOrder(int* order, int size){
	//	perform strict check to see whether order is bijection, with range 0 ~ size-1
	#ifdef STRICT_CHECK
		if( !isBijection(order, size, 0, 1) ){
			std::cout << "error: invertOrder, order is not bijection, fail to invert reorder" << std::endl;
			return 0;
		}
	#endif

	int* tmp_order = new int[size];
	for(int i=0; i<size; i++)
		tmp_order[order[i]] = i;
	memcpy(order, tmp_order, sizeof(int)*size);
}

/**
	reorder a data array to target order

	example: order[2] = 7; means data[2] should be moved to
	the position of data[7]

	order array must be bijection, starting from 0.
	If STRICT_CHECK is defined, we perform bijection check

	\param	data		the data array 
	\param	size		size of the array
	\param	order		the target order of each element
	\return				whether reorder succeed
*/
template<typename T>
int reorderDataToTargetOrder(T*			data,
							 int		size,
							 const int*	order){
	//	perform strict check to see whether order is bijection, with range 0 ~ size-1
	#ifdef STRICT_CHECK
		if( !isBijection(order, size, 0, 1) ){
			std::cout << "error: reorderDataToTargetOrder, order is not bijection, fail to reorder" << std::endl;
			return 0;
		}
	#endif

	//	find the start position that need to be reordered
	int start = 0;
	for(start=0; start<size; start++)
		if( order[start] != start )
			break;
	if( start == size )		//	if data has been ordered, just stop here
		return 1;

	//	create a temp array
	std::vector<T>	tmp_data;
	tmp_data.resize(size-start);

	//	write the ordered data in tmp_data array. Since order is bijection, so it is safe to do this
	for(int i=start; i<size; i++)
		tmp_data[order[i]-start] = data[i];

	for(int i=start; i<size; i++)
		data[i] = tmp_data[i-start];

	return 1;
}

/**
	reorder a data array to target order

	example: order[2] = 7; means data[2] should be moved to
	the position of data[7]

	order array must be bijection, or the function may
	cause error. STRICT_CHECK can be used to check whether
	the array is bijection

	\param	data	the data array 
	\param	order	new position of each element
	\return			whether reorder succeed
*/
template<typename T>
int reorderDataToTargetOrder(std::vector<T>&			data,
							 const std::vector<int>&	order){
	if( data.size() != order.size() ){
		std::cout << "error: reorderDataToTargetOrder, data size, order size don't match" << std::endl;
		return 0;
	}
	if (data.empty())
		return 1;
	return reorderDataToTargetOrder((T*)&data[0], data.size(), (int*)&order[0]);
}


/**
	reorder a data array by source order

	example: order[7] = 2; means data[7] was originally at 
	the position of data[2], so move data[2] to data[7]

	order array must be bijection, starting from 0.
	If STRICT_CHECK is defined, we perform bijection check

	\param	data		the data array 
	\param	size		size of the array
	\param	order		the source order of each element
	\return				whether reorder succeed
*/
template<typename T>
int reorderDataBySourceOrder(T*			data,
							 int		size,
							 const int*	order){
	//	perform strict check to see whether order is bijection, with range 0 ~ size-1
	#ifdef STRICT_CHECK
		if( !isBijection(order, size, 0, 1) ){
			std::cout << "error: reorderDataToTargetOrder, order is not bijection, fail to reorder" << std::endl;
			return 0;
		}
	#endif

	//	find the start position that need to be reordered
	int start = 0;
	for(start=0; start<size; start++)
		if( order[start] != start )
			break;
	if( start == size )		//	if data has been ordered, just stop here
		return 1;

	//	create a temp array
	std::vector<T>	tmp_data;
	tmp_data.resize(size-start);

	//	write the ordered data in tmp_data array. Since order is bijection, so it is safe to do this
	for(int i=start; i<size; i++)
		tmp_data[i-start] = data[order[i]];

	for(int i=start; i<size; i++)
		data[i] = tmp_data[i-start];

	return 1;
}

/**
	reorder a data array by source order

	example: order[7] = 2; means data[7] was originally at 
	the position of data[2], so move data[2] to data[7]

	order array must be bijection, starting from 0.
	If STRICT_CHECK is defined, we perform bijection check

	\param	data		the data array 
	\param	order		the source order of each element
	\return				whether reorder succeed
*/
template<typename T>
int reorderDataBySourceOrder(std::vector<T>&			data,
							 const std::vector<int>&	order){
	if( data.size() != order.size() ){
		std::cout << "error: reorderDataBySourceOrder, data size, order size don't match" << std::endl;
		return 0;
	}
	if (data.empty())
		return 1;
	return reorderDataBySourceOrder((T*)&data[0], data.size(), (int*)&order[0]);
}


///@}

//	========================================
///@{
/**	@name Remove Duplicates
*/
//	========================================

/**
	Remove duplicate in vector

	only elements identical in all bits are duplicates. This is achieved by hashing
*/
template<typename T>
inline void removeDuplicates(std::vector<T>& vec){
	std::unordered_set<T, utils::BitwiseHasher<T>> hash_vec;
	hash_vec.reserve(vec.size());

	//	parse all data
	int count = 0;
	for (int i = 0; i < vec.size(); i++) {
		//std::unordered_set<T, utils::BitwiseHasher<T>>::const_iterator iter = hash_vec.find(vec[i]);
		auto iter = hash_vec.find(vec[i]);
		if (iter == hash_vec.end()) {
			vec[count++] = vec[i];
			hash_vec.insert(vec[i]);
		}
	}

	vec.resize(count);
}
///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_REORDER_H__