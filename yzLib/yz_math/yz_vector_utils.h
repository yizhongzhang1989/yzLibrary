/***********************************************************/
/**	\file
	\brief		Vector Utilities
	\details	functions for general vector and yz::Vec
	\author		Yizhong Zhang
	\date		6/5/2012
*/
/***********************************************************/
#ifndef __YZ_VECTOR_UTILS_H__
#define __YZ_VECTOR_UTILS_H__

#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_numerical_utils.h"

namespace yz{

//	========================================
///@{
/**	@name General Vector Utilities
*/
//	========================================
/**
	set all the vector component to the same value
*/
template<typename T1, typename T2>
inline void setVectorValue(T1* vec, T2 value, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] = value;
}
/**
	calculate square sum of each element of a vector
*/
template<typename T>
inline T calculateVectorSquareSum(const T* vec, int dim){
	if( dim < 1 )
		return 0;
	T square_sum = vec[0] * vec[0];
	for( int i=1; i<dim; i++ )
		square_sum += vec[i] * vec[i];
	return square_sum;
}

/**
	calculate length of the vector
*/
template<typename T>
inline T calculateVectorLength(const T* vec, int dim){
	return sqrt( calculateVectorSquareSum(vec, dim) );
}

/**
	calculate sum of each element of a vector

	this function is also able to calculate non-basic 
	data type vectors if dimension is legal and + operator
	is defined
*/
template<typename T>
inline T calculateVectorSum(const T* vec, int dim){
	if( dim < 1 )
		return 0;
	T sum = vec[0];
	for( int i=1; i<dim; i++ )
		sum += vec[i];
	return sum;
}

/**
	calculate mean of each element of a vector

	this function is also able to calculate non-basic 
	data type vectors if dimension is legal and + * operators
	are defined
*/
template<typename T>
inline T calculateVectorMean(const T* vec, int dim){
	if( dim < 1 )
		return 0;
	T sum = vec[0];
	for( int i=1; i<dim; i++ )
		sum += vec[i];
	return sum * (1.0/dim);
}

/**
	calculate dot product of two vectors of the same dimension
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2 calculateVectorDot(const T1* vec1, const T2* vec2, int dim){
	if( dim < 1 )
		return 0;
	PROMOTE_T1_T2 sum = vec1[0] * vec2[0];
	for( int i=1; i<dim; i++ )
		sum += vec1[i] * vec2[i];
	return sum;
}

/**
	scale a vector, new vector overwrite the old vector
*/
template<typename T1, typename T2>
inline void setVectorScale(T1* vec, T2 scale, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] *= scale;
}

/**
	add the second vector to the original

	\param	vec			the original vector, value will be modified
	\param	vec_add		the vector to add
	\param	dim			dimension of the two vectors
*/
template<typename T1, typename T2>
inline void setVectorAddVector(T1* vec, const T2* vec_add, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] += vec_add[i];
}

/**
	add value to each component of the vector

	\param	vec			the original vector, value will be modified
	\param	value_add	the value to add
	\param	dim			dimension of the vector
*/
template<typename T1, typename T2>
inline void setVectorAddNumber(T1* vec, T2 value_add, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] += value_add;
}

/**
	sub the second vector from the original

	\param	vec			the original vector, value will be modified
	\param	vec_sub		the vector to subtract
	\param	dim			dimension of the two vectors
*/
template<typename T1, typename T2>
inline void setVectorSubVector(T1* vec, const T2* vec_sub, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] -= vec_sub[i];
}

/**
	sub value to each component of the vector

	\param	vec			the original vector, value will be modified
	\param	value_sub	the value to sub
	\param	dim			dimension of the two vectors
*/
template<typename T1, typename T2>
inline void setVectorSubNumber(T1* vec, T2 value_sub, int dim){
	for( int i=0; i<dim; i++ )
		vec[i] -= value_sub;
}

/**
	set the vector to normalize

	\param	vec			the vector to be normalized
	\param	dim			dimension of the vector
*/
template<typename T>
inline void setVectorNormalize(T* vec, int dim){
	T inv_length = 1.0 / calculateVectorLength(vec, dim);
	setVectorScale(vec, inv_length, dim);
}
/**
	get the min value from the vector
*/
template<typename T>
inline T getVectorMin(T* vec, int dim){
	assert(dim>0);
	T min_val = vec[0];
	for( int i=1; i<dim; i++ )
		if( vec[i] < min_val )
			min_val = vec[i];
	return min_val;
}

/**
	get the max value from the vector
*/
template<typename T>
inline T getVectorMax(T* vec, int dim){
	assert(dim>0);
	T max_val = vec[0];
	for( int i=1; i<dim; i++ )
		if( vec[i] > max_val )
			max_val = vec[i];
	return max_val;
}




///@}

//	========================================
///@{
/**	@name Geometry functions of angle between vectors
*/
//	========================================
/**
	angle in degrees between two vectors, -180 ~ 180 degree, rotate r1 to r2
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT angleDegBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	return rad2deg( angleRadBetweenVectors(r1, r2) );
}

/**
	angle in degrees between two vectors, 0 ~ 180 degree
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT angleDegBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	return rad2deg( angleRadBetweenVectors(r1, r2) );
}

/**
	angle in rad between two vectors, -PI ~ PI, rotate r1 to r2
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT angleRadBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	r1.SetNormalize();
	r2.SetNormalize();
	TYPE_PROMOTE(T1, T2) cross_val	= cross(r1, r2);
	TYPE_PROMOTE(T1, T2) dot_val	= clamp(dot(r1, r2), -1.0f, 1.0f);

	if( abs(cross_val) < 1e-5 ){
		if( dot_val > 0 )
			return 0;
		else
			return YZ_PI;	
	}
	else{
		if( cross_val > 0 )
			return acos(dot_val);
		else
			return -acos(dot_val);
	}
}

/**
	angle in rad between two vectors, 0 ~ PI
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT angleRadBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	r1.SetNormalize();
	r2.SetNormalize();

	return acos( clamp(dot(r1, r2), -1.0f, 1.0f) );
}

/**
	sin value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT sinBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	return cross(r1, r2) / sqrt( r1.SquareLength() * r2.SquareLength() );
}
/**
	sin value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT sinBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	return cross(r1, r2).Length() / sqrt( r1.SquareLength() * r2.SquareLength() );
}

/**
	cos value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT cosBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	return dot(r1, r2) / sqrt( r1.SquareLength() * r2.SquareLength() );
}
/**
	cos value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT cosBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	return dot(r1, r2) / sqrt( r1.SquareLength() * r2.SquareLength() );
}

/**
	tan value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT tanBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	return cross(r1, r2) / dot(r1, r2);
}
/**
	tan value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT tanBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	return cross(r1, r2).Length() / dot(r1, r2);
}

/**
	cot value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT cotBetweenVectors(Vec2<T1> r1, Vec2<T2> r2){
	return dot(r1, r2) / cross(r1, r2);
}
/**
	cot value of angle between two vectors
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT cotBetweenVectors(Vec3<T1> r1, Vec3<T2> r2){
	return dot(r1, r2) / cross(r1, r2).Length();
}

///@}


}	//	namespace yz

#endif	//	__YZ_VECTOR_UTILS_H__