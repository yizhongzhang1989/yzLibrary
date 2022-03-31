/***********************************************************/
/**	\file
	\brief		Data Filter
	\details	This file provide filters to make data more smooth.
	\author		Yizhong Zhang
	\date		6/9/2012
*/
/***********************************************************/
#ifndef __YZ_FILTER_H__
#define __YZ_FILTER_H__

#include <stdlib.h>
#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_numerical_utils.h" 

namespace yz{

//	========================================
///@{
/**	@name 1D array filter
*/
//	========================================

/**
	Filter by average neighbor data, including front and end.

	Non-basic type is also acceptable, if operator +=, -= are defined

	the old array is used to store filtered data

	\param	data		the array to be filtered
	\param	size		length of the array
	\param	diameter	diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void setFilterByAverage(T* data, int size, int diameter=3){
	if( size <= 1 )
		return;
	if( diameter <= 1 )
		return;

	//	create tmp array
	T* tmp = new T[size+diameter];
	T* ptr = tmp + (diameter-1)/2;
	memcpy(ptr, data, sizeof(T)*size);
	for( T* p=tmp; p!=ptr; p++ )
		*p = *ptr;
	for( T *p=ptr+size, *end = tmp+size+diameter; p!=end; p++ )
		*p = *(ptr+size-1);

	//	create initial average
	T avg	= *tmp;
	for( int i=1; i<diameter; i++ )
		avg += tmp[i];

	//	parse the array
	T* data_ptr = data;
	T* data_end = data + size;
	ptr = tmp;
	double inv_diameter = 1.0 / diameter;
	while( data_ptr != data_end ){
		*data_ptr = avg * inv_diameter;
		avg -= *ptr;
		avg += *(ptr + diameter);
		ptr ++;
		data_ptr ++;
	}

	delete[] tmp;
}

/**
	Filter by average neighbor data stored in vector, including front and end.

	Non-basic type is also acceptable, if operator +=, -= are defined

	\param	data		data stored in std::vector
	\param	diameter	diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void setFilterByAverage(std::vector<T>& data, int diameter=3){
	if( !data.empty() )
		setFilterByAverage(&data[0], data.size(), diameter);
}

/**
	Filter by average neighbor data, including front and end.

	Non-basic type is also acceptable, if operator +=, -=, *double are defined.

	Filtering Method:

	1,	expand the array front and end, insert same data of front and end to expanded space

	2,	calculate the average of neighbor data and write to destination

	if des == src, we filter on the old array. If they are not the same,
	des and src cannot overlap, or the result may not be correct

	\param	des			filter destination
	\param	src			filter source
	\param	size		length of the array
	\param	diameter	diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void filterByAverage(T* des, T* src, int size, int diameter=3){
	if( size <= 0 )
		return;

	if( des != src )
		memcpy(des, src, sizeof(T)*size);

	setFilterByAverage(des, size, diameter);
}

/**
	Filter by average neighbor data, keep two ends fixed.

	Non-basic type is also acceptable, if operator +=, -= are defined

	the old array is used to store filtered data

	\param	data			the array to be filtered
	\param	size			length of the array
	\param	fixed_length	length of fixed area, default 1: just fixed the ends, no smooth
	\param	diameter		diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void setFilterByAverageFixEnds(T* data, int size, int fixed_length=1, int diameter=3){
	if( fixed_length*2 > size )		//	fixed area cannot overlap
		fixed_length = size / 2;

	if( fixed_length <= 0 ){
		setFilterByAverage(data, size, diameter);
		return;
	}

	T front = data[0];
	T end	= data[size-1];

	setFilterByAverage(data, size, diameter);

	//	smooth fixed area
	front -= data[0];
	end -= data[size-1];
	for(int i=0; i<fixed_length; i++ ){
		double coef = double(fixed_length-i) / fixed_length;
		data[i] += front * coef;
		data[size-i] += end * coef;
	}

}
/**
	Filter by average neighbor data stored in vector, keep two ends fixed.

	Non-basic type is also acceptable, if operator +=, -= are defined

	\param	data		data stored in std::vector
	\param	fixed_length	length of fixed area, default 1: just fixed the ends, no smooth
	\param	diameter	diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void setFilterByAverageFixEnds(std::vector<T>& data, int fixed_length=1, int diameter=3){
	if( !data.empty() )
		setFilterByAverageFixEnds(&data[0], data.size(), fixed_length, diameter);
}
/**
	Filter by average neighbor data, keep two ends fixed.

	Non-basic type is also acceptable, if operator +=, -=, *double are defined.

	Filtering Method:

	1,	record front and end data

	2,	filter data using non-fixed method

	3,	recover front and end, smooth that area

	if des == src, we filter on the old array. If they are not the same,
	des and src cannot overlap, or the result may not be correct

	\param	des			filter destination
	\param	src			filter source
	\param	size		length of the array
	\param	fixed_length	length of fixed area, default 1: just fixed the ends, no smooth
	\param	diameter	diameter of average(how many numbers are used to average)
*/
template<typename T>
inline void filterByAverageFixEnds(T* des, const T* src, int size, int fixed_length=1, int diameter=3){
	if( size <= 0 )
		return;

	if( des != src )
		memcpy(des, src, sizeof(T)*size);

	setFilterByAverageFixEnds(des, size, fixed_length, diameter);
}


/**
	bilateral filter

	only basic data type could be used, because we need substraction of data

	destination array could be the same as source array, but they cannot overlap if not the same

	\param	des			filter destination
	\param	src			filter source
	\param	size		length of the array
	\param	radius		radius of the kernel
	\param	euclidean	standard variation of quantity
	\param	gaussian	standard variation of distance
*/
template<typename T1, typename T2, typename T3>
inline void filterBilateral(T1* des, T1* src, int size, int radius, T2 euclidean, T3 gaussian){
	//	if the size of radius is small, we use existing kernel space, or we have to alloc a new one
	double kernel_space[16];
	double* kernel = kernel_space;
	if( radius > 16 ){
		kernel = new double[radius];
		assert(kernel);
	}

	//	create tmp array if necessary
	T1* tmp_array = NULL;
	if( des == src ){
		tmp_array = new T1[size];
		assert(tmp_array);
		memcpy(tmp_array, src, sizeof(T1)*size);
		src = tmp_array;
	}

	//	calculate gaussian kernel
	double coef1 = - 1.0 / ( 2.0 * gaussian * gaussian );
	for(int i=0; i<radius; i++){
		kernel[i] = exp( coef1 * i * i );
	}

	//	filter the data
	double coef2 = - 1.0 / ( 2.0 * euclidean * euclidean );
	for(int i=0; i<size; i++){		//	for each element
		double center_data = src[i];
		double val = 0.0;
		double sum = 0.0;
		for(int j=myMax(i-radius+1, 0); j<myMin(i+radius, size); j++){	//	for each neighbor in kernel
			int dist = i>j ? i-j : j-i;
			double curr_data = src[j];
			double euc_kernel = exp( coef2 * (curr_data - center_data) * (curr_data - center_data) );
			double factor = kernel[dist] * euc_kernel;
			val += factor * curr_data;
			sum += factor;
		}

		//	write data
		des[i] = val / sum;
	}

	//	delete extra space
	if( tmp_array )
		delete tmp_array;
	if( kernel != kernel_space )
		delete kernel;
}


/**
	gaussian filter

	only basic data type could be used, because we need substraction of data

	destination array could be the same as source array, but they cannot overlap if not the same

	\param	des			filter destination
	\param	src			filter source
	\param	size		length of the array
	\param	radius		radius of the kernel
	\param	sigma		standard variation
*/
template<typename T1, typename T2>
inline void filterGaussian(T1* des, T1* src, int size, int radius, T2 sigma){
	//	if the size of radius is small, we use existing kernel space, or we have to alloc a new one
	double kernel_space[16];
	double* kernel = kernel_space;
	if( radius > 16 ){
		kernel = new double[radius];
		assert(kernel);
	}

	//	create tmp array if necessary
	T1* tmp_array = NULL;
	if( des == src ){
		tmp_array = new T1[size];
		assert(tmp_array);
		memcpy(tmp_array, src, sizeof(T1)*size);
		src = tmp_array;
	}

	//	calculate gaussian kernel
	double coef1 = - 1.0 / ( 2.0 * sigma * sigma );
	for(int i=0; i<radius; i++)
		kernel[i] = exp( coef1 * i * i );

	//	normalize the kernel, note we only store half the kernel
	double sum = kernel[0];		//	this value should be 1
	for(int i=1; i<radius; i++)
		sum += kernel[i] * 2.0;	//	each value appear twice
	sum = 1.0 / sum;
	for(int i=0; i<radius; i++)
		kernel[i] *= sum;

	//	filter the data
	for(int i=0; i<size; i++){		//	for each element
		double center_data = src[i];
		double val = 0.0;
		sum = 0.0;
		for(int j=yz::myMax(i-radius+1, 0); j<yz::myMin(i+radius, size); j++){	//	for each neighbor in kernel
			int dist = i>j ? i-j : j-i;
			double factor = kernel[dist];
			val += factor * src[j];
			sum += factor;
		}

		//	write data
		des[i] = val / sum;
	}

	//	delete extra space
	if( tmp_array )
		delete tmp_array;
	if( kernel != kernel_space )
		delete kernel;
}
///@}

//	========================================
///@{
/**	@name Soft Trajectory
*/
//	========================================

/**
	Soft start function, defined in 0 - 1

	f(0) = 0, f(1) = 1, f'(0) = 0, f'(1) = 0;

	\param	x	x
	\return		f(x)
*/
inline double softStartFunc(double x){
	x = clamp(x, 0.0, 1.0);
	if(x <= 0.5)
		return 4.0*x*x*x;
	else
		return 1.0 - 4.0*(1-x)*(1-x)*(1-x);
}

/**
	Set the trajectory to start(f'(0) = 0) using softStartFunc()

	Non-basic type is also acceptable, if operator +, -, *double are defined

	the old array is used to store filtered data

	\param	data			trajectory to be filtered
	\param	size			length of the array
	\param	soft_length		length of soft area
*/
template<typename T>
inline void setSoftStart(T* data, int size, int soft_length){
	if( size <= 1 )
		return;
	if( soft_length <= 1 )
		return;
	if( soft_length > size )
		soft_length = size;

	T d0 = data[0];
	for( int i=1; i<soft_length; i++ ){
		double x = double(i)/soft_length;
		data[i] = d0 + (data[i]-d0) * softStartFunc(x);
	}
}
/**
	Set the trajectory in vector to start(f'(0) = 0) using softStartFunc()

	Non-basic type is also acceptable, if operator +, -, *double are defined

	the old array is used to store filtered data

	\param	data			trajectory in vector
	\param	soft_length		length of soft area
*/
template<typename T>
inline void setSoftStart(std::vector<T>& data, int soft_length){
	if(	!data.empty() )
		setSoftStart(&data[0], data.size(), soft_length);
}



/**
	Set the trajectory to stop(f'(0) = 0) using softStartFunc()

	Non-basic type is also acceptable, if operator +, -, *double are defined

	the old array is used to store filtered data

	\param	data			trajectory to be filtered
	\param	size			length of the array
	\param	soft_length		length of soft area
*/
template<typename T>
inline void setSoftStop(T* data, int size, int soft_length){
	if( size <= 1 )
		return;
	if( soft_length <= 1 )
		return;
	if( soft_length > size )
		soft_length = size;

	T d0 = data[size-1];
	for( int i=1; i<soft_length; i++ ){
		double x = double(i)/soft_length;
		data[size-1-i] = d0 + (data[size-1-i]-d0) * softStartFunc(x);
	}
}
/**
	Set the trajectory in vector to stop(f'(0) = 0) using softStartFunc()

	Non-basic type is also acceptable, if operator +, -, *double are defined

	the old array is used to store filtered data

	\param	data			trajectory in vector
	\param	soft_length		length of soft area
*/
template<typename T>
inline void setSoftStop(std::vector<T>& data, int soft_length){
	if(	!data.empty() )
		setSoftStop(&data[0], data.size(), soft_length);
}

///@}

}	//	namespace yz

#endif	//	__YZ_FILTER_H__