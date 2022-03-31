/***********************************************************/
/**	\file
	\brief		Set memory to a given pattern
	\details	Functions in this file aims to set memory of an image,
				so width and height are need to give the size of memory.
	\author		Yizhong Zhang
	\date		5/26/2012
*/
/***********************************************************/
#ifndef __YZ_MEM_PATTERN_H__
#define __YZ_MEM_PATTERN_H__

#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{ namespace utils{
//	========================================
///@{
/**	@name Set Memory to A Given Pattern
*/
//	========================================

/**
	Set Memory to the Same Value

	\param	ptr		memory ptr
	\param	width	width of memory image
	\param	height	height of memory image
	\param	value	value to be set
*/
template<typename T, typename TYPE>
inline void setMemSolid(T* ptr, int width, int height, TYPE value){
	int size = width*height;
	for( int i=0; i<size; i++ )
		ptr[i] = value;
}

/**
	Set Memory to Vertical Sine Wave

	\param	ptr			memory ptr
	\param	width		width of memory image
	\param	height		height of memory image
	\param	min_value	minimal value of sin function
	\param	max_value	maximal value of sin function
	\param	period		period of the wave, in pixels
*/
template<typename T, typename T1, typename T2>
inline void setMemVerticalWave(T* ptr, int width, int height, T1 min_value, T2 max_value, int period){
	TYPE_PROMOTE(T, float) amp = ((TYPE_PROMOTE(T, float))max_value - min_value) / 2;
	TYPE_PROMOTE(T, float) mid = ((TYPE_PROMOTE(T, float))max_value + min_value) / 2;
	for( int i=0; i<width; i++ )
		ptr[i] = amp * sin(2*YZ_PI*i/period) + mid;
	for( int j=1; j<height; j++ )
		memcpy(ptr+j*width, ptr, sizeof(T)*width);
}

/**
	Set Memory to Horizontal Sine Wave

	\param	ptr			memory ptr
	\param	width		width of memory image
	\param	height		height of memory image
	\param	min_value	minimal value of sin function
	\param	max_value	maximal value of sin function
	\param	period		period of the wave, in pixels
*/
template<typename T, typename T1, typename T2>
inline void setMemHorizontalWave(T* ptr, int width, int height, T1 min_value, T2 max_value, int period){
	TYPE_PROMOTE(T, float) amp = ((TYPE_PROMOTE(T, float))max_value - min_value) / 2;
	TYPE_PROMOTE(T, float) mid = ((TYPE_PROMOTE(T, float))max_value + min_value) / 2;
	//	calculate first period
	for( int j=0; j<(myMin(height, period)); j++ ){
		T value = amp * sin(2*YZ_PI*j/period) + mid;
		for( int i=0; i<width; i++ )
			ptr[j*width+i] = value;
	}
	//	fill other period
	for( int j=period; j<height; j++ )
		memcpy(ptr+j*width, ptr+(j-period)*width, sizeof(T)*width);
}

///@}

//	========================================
///@{
/**	@name Check Data in Memory
*/
//	========================================

/**
	Convert data to array

	To print the data in binary format, call \n
	std::cout << MemoryBinary<DataType>(data) << std::endl;
*/
template<class T>
class MemoryBinary {
public:
	std::vector<unsigned char> data;

	MemoryBinary(T n) {
		data.resize(sizeof(T));
		*(T*)&data[0] = n;
	}

	void Print() const {
		for (int i = 0; i < data.size(); i++) {
			unsigned char tmp = data[i];
			for (int j = 0; j < 8; j++) {
				printf("%d", (tmp & 0x80 ? 1 : 0));
				tmp <<= 1;
			}
			if (i + 1 != data.size())
				printf("|");
		}
	}
};

template<class T>
inline std::ostream& operator << (
	std::ostream&			stream,
	const MemoryBinary<T>&	mb)
{
	mb.Print();
	return stream;
}

/**
	Convert data to array

	To print the data in binary format, call \n
	std::cout << memoryBinary(data) << std::endl;
*/
template<typename T>
MemoryBinary<T> memoryBinary(T n) {
	return MemoryBinary<T>(n);
}


/**
	Convert data to array

	To print the data in hex format, call \n
	std::cout << MemoryHex<DataType>(data) << std::endl;
*/
template<class T>
class MemoryHex : public MemoryBinary<T> {
public:
	MemoryHex(T n) : MemoryBinary<T>(n) {	}

	void Print() const {
		for (int i = 0; i < this->data.size(); i++) {
			printf("%02X", this->data[i]);
		}
	}
};

template<class T>
inline std::ostream& operator << (
	std::ostream&		stream,
	const MemoryHex<T>&	mh)
{
	mh.Print();
	return stream;
}

/**
	Convert data to array

	To print the data in hex format, call \n
	std::cout << memoryHex(data) << std::endl;
*/
template<typename T>
MemoryHex<T> memoryHex(T n) {
	return MemoryHex<T>(n);
}


///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_MEM_PATTERN_H__