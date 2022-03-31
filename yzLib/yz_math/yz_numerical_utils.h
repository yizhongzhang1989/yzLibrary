/***********************************************************/
/**	\file
	\brief		Commonly Used Numerical Utilities
	\details	Random Number, deg-rad conversion, clamp.
				Functions in this file should be functions 
				that used very often. Specific functions 
				should not be placed here. 
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_NUMERICAL_UTILS_H__
#define __YZ_NUMERICAL_UTILS_H__

#include <iostream>
#include <math.h>
#include <time.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
///@{
/**	@name Random Number
*/
//	========================================

/**
	Return random floating point number between min_val and max_val. 

	return type is promoted to the biggest one among float, T1 and T2
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT randFloatingPointNumber(T1 min_val, T2 max_val){
	return TYPE_PROMOTE(TYPE_PROMOTE(T1,T2), float)(rand()%RAND_MAX)/(RAND_MAX-1)*(max_val-min_val) + min_val;
}

/**
	Return random number between min_val and max_val. 

	return type is promoted to bigger one of min_val, max_val
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2 randNumber(T1 min_val, T2 max_val){
	return randFloatingPointNumber(min_val, max_val);
}

/**
	Return random float number between 0 and 1
*/
inline float rand0to1f(){
	return rand()%RAND_MAX /  float(RAND_MAX-1);
}

/**
	Return random double number between 0 and 1
*/
inline double rand0to1d(){
	return rand()%RAND_MAX /  double(RAND_MAX-1);
}

/**
	Return random float number between 0 and 2*PI
*/
inline float rand0to2PIf(){
	return rand0to1f() * 2.0f * float(YZ_PI);
}

/**
	Return random double number between 0 and 2*PI
*/
inline float rand0to2PId(){
	return float(rand0to1d() * 2.0f * float(YZ_PI));
}

/**
	Set current time to seed, so random number sequence 
	will be different each time of execution
*/
inline void setRandWithTime(){
	srand( (unsigned int)time(NULL) );
}
///@}

//	========================================
///@{
/**	@name Deg-Rad Conversion
*/
//	========================================
template<typename T>
inline PROMOTE_T_TO_FLOAT deg2rad(T deg){
	return deg * YZ_PI / 180;
}

template<typename T>
inline PROMOTE_T_TO_FLOAT rad2deg(T rad){
	return rad * 180 / YZ_PI;
}

///@}

//	========================================
///@{
/**	@name Clamp and Round
*/
//	========================================
/**
	Clamp value to min_val and max_val
*/
template<typename T, typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT clamp(T value, T1 min_val, T2 max_val){
	if( value > max_val )
		return max_val;
	if( value < min_val )
		return min_val;
	return value;
}

/**
	round a floating point number to closest integer
*/
template<typename T>
inline T roundToClosestInteger(T num){
	return floor(num+T(0.5));
}

///@}

//	========================================
///@{
/**	@name	Numerical Functions
*/
//	========================================

/**
	calculate the power of integer number

	if exponent < 0, just retun 0	\n
*/
inline int power(int base, int exponent){
	if( exponent < 0 )
		return 0;
	else if( exponent == 0 )
		return 1;

	int acc = base;
	while(--exponent)
		acc *= base;
	return acc;
}


///@}

//	========================================
///@{
/**	@name Two Number Operation
*/
//	========================================

/**
	swap two numbers

	swap conflicts with c++, so I use the name mySwap
*/
template<typename T>
inline void mySwap(T& a, T& b){
	T tmp = a;
	a = b;
	b = tmp;
}

/**
	return smaller number
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2 myMin(T1 a, T2 b){
	return a<b ? a : b;
}

/**
	return bigger number
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2 myMax(T1 a, T2 b){
	return a>b ? a : b;
}


///@}

}	//	namespace yz

#endif	//	__YZ_NUMERICAL_UTILS_H__