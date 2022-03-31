/***********************************************************/
/**	\file
	\brief		Kernel Functions
	\details	a collection of kernel functions. 
				Most kernel functions only open window in -1 to 1.

				Detailed function description can be found on wiki:
				http://en.wikipedia.org/wiki/Kernel_(statistics)
	\author		Yizhong Zhang
	\date		6/11/2012
*/
/***********************************************************/
#ifndef __YZ_KERNEL_FUNCTIONS_H__
#define __YZ_KERNEL_FUNCTIONS_H__

#include <stdlib.h>
#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
///@{
/**	@name Kernel Functions used in Statistics
*/
//	========================================

/**
	Uniform Kernel, -1<x<1

	K(x) = 1/2
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelUniform(T x){
	if( x<-1 || x>1 )
		return 0;

	return 0.5;
}

/**
	Triangular Kernel, -1<x<1

	K(x) = 1-abs(x)
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelTriangular(T x){
	if( x<-1 || x>1 )
		return 0;
	
	if( x<0 )
		return 1+x;
	else
		return 1-x;
}

/**
	Epanechnikov Kernel, -1<x<1

	K(x) = 3/4 * (1-x^2)
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelEpanechnikov(T x){
	if( x<-1 || x>1 )
		return 0;
	
	return (3.0/4.0)*(1-x*x);
}

/**
	Quartic Kernel, -1<x<1

	K(x) = 15/16 * (1-x^2)^2
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelQuartic(T x){
	if( x<-1 || x>1 )
		return 0;

	return (15.0/16.0)*(1-x*x)*(1-x*x);
}

/**
	Biweight Kernel, same as Quartic, -1<x<1

	K(x) = 15/16 * (1-x^2)^2
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelBiweight(T x){
	return kernelQuartic(x);
}

/**
	Triweight Kernel, -1<x<1

	K(x) = 35/32 * (1-x^2)^3
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelTriweight(T x){
	if( x<-1 || x>1 )
		return 0;

	return (35.0/32.0)*(1-x*x)*(1-x*x)*(1-x*x);
}

/**
	Tricube Kernel, -1<x<1

	K(x) = 70/81 * (1-abs(x)^3)^3
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelTricube(T x){
	if( x<-1 || x>1 )
		return 0;

	x = (x<0 ? -x : x);
	return (70.0/81.0)*(1-x*x*x)*(1-x*x*x)*(1-x*x*x);
}

/**
	Gaussian Kernel, x belong to Real

	K(x) = 1/sqrt(2*PI) * e^(-1/2 * x^2)
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelGaussian(T x){
	return (1.0/sqrt(2*YZ_PI))*exp(-0.5*x*x);
}

/**
	Cosine Kernel, -1<x<1

	K(x) = PI/4 * cos(PI/2 * x)
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelCosine(T x){
	if( x<-1 || x>1 )
		return 0;

	return (YZ_PI/4.0)*cos(YZ_PI/2.0 * x);
}

/**
	SharpQuadratic Kernel, -1<x<1

	K(x) = 3/2 * (x-1)^2
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelSharpQuadratic(T x){
	if( x<-1 || x>1 )
		return 0;

	x = (x<0 ? -x : x);
	return 3.0/2.0 * (x-1)*(x-1);
}

/**
	SharpCubic Kernel, -1<x<1

	K(x) = 2 * (x-1)^3
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT kernelSharpCubic(T x){
	if( x<-1 || x>1 )
		return 0;

	x = (x<0 ? -x : x);
	return 2 * (1-x)*(1-x)*(1-x);
}

///@}

}	//	namespace yz

#endif	//	__YZ_KERNEL_FUNCTIONS_H__