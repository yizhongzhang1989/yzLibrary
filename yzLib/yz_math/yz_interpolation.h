/***********************************************************/
/**	\file
	\brief		Some Interpolation functions
	\details	Interpolation on different dimensions, with 
				different method.

				Type promote is used in implementing the functions.
				However, these promotion only work on basic data 
				types. If user defined types are passed as function
				arguments, original type is returned, since type 
				promotion only work for basic data type.
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_INTERPOLATION_H__
#define __YZ_INTERPOLATION_H__

#include <assert.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
///@{
/**	@name 1D Interpolation
*/
//	========================================
/**	
	Nearest Interpolation on Unit Length

	\param	val0	value at position 0
	\param	val1	value at position 1
	\param	t		interpolation position, should be 0 ~ 1
	\return			interplated value
*/
template <typename T1, typename T2>
inline PROMOTE_T1_T2 interpNearest(T1 val0, T2 val1, double t){
	if( t < 0.5 )
		return val0;
	else
		return val1;
}

/**
	Linear Interpolation on Unit Length

	\param	val0	value at position 0
	\param	val1	value at position 1
	\param	t		interpolation position, should be 0 ~ 1
	\return			interplated value
*/
template <typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT interpLinear(T1 val0, T2 val1, double t){
	return val0 * (1.0 - t) + val1 * t;
}

/**
	Interpolation using cosine function

	\param	val0	value at position 0
	\param	val1	value at position 1
	\param	t		interpolation position, should be 0 ~ 1
	\return			interplated value
*/
template<typename T1, typename T2>
inline PROMOTE_T1_T2_TO_FLOAT interpCosine(T1 val0, T2 val1, double t){
	t = (1 - cos(t*YZ_PI)) / 2;
	return val0 * (1.0 - t) + val1 * t;
}

/**
	Cubic Interpolation, 4 points are needed. p0-p1-(interpolation area)-p2-p3

	\param	val0	value at position 0
	\param	val1	value at position 1
	\param	val2	value at position 2
	\param	val3	value at position 3	
	\param	t		interpolation position between point 1 and 2, should be 0 ~ 1
	\return			interplated value
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT interpCubic(T val0, T val1, T val2, T val3, double t){
	double t2 = t * t;

	return
		(val3 - val2 - val0 + val1) * t2*t +
		(val0*2 - val1*2 - val3 + val2) * t2 +
		(val2 - val0) * t + val1;
}

/**
	Catmull-Rom Spline Interpolation, 4 points are needed. p0-p1-(interpolation area)-p2-p3

	Catmull-Rom is an improved version of Cubic interpolation with some parameters changed. 
	It can provide better result in some cases. 

	\param	val0	value at position 0
	\param	val1	value at position 1
	\param	val2	value at position 2
	\param	val3	value at position 3	
	\param	t		interpolation position between point 1 and 2, should be 0 ~ 1
	\return			interplated value
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT interpCatmullRom(T val0, T val1, T val2, T val3, double t){
	double t2 = t * t;

	return 
		(val0*(-0.5)	+ val1*1.5		+ val2*(-1.5)	+ val3*(0.5)) * t2*t +
		(val0			+ val1*(-2.5)	+ val2*2		+ val3*(-0.5)) * t2 +
		(val0*(-0.5)	+ val2*(0.5)) * t + val1;
}

/**
	Cubic Bezier Interpolation

	Can be used to interpolate user defined data, make sure T * double is defined

	If you want to generate a curve, just pass position Vec2-4 of control points as val0-3. 
	Vec2-4 have defined * double operator

	\param	val0	value of point 0
	\param	val1	value of point 1
	\param	val2	value of point 2
	\param	val3	value of point 3
	\param	t		interpolation position, , should be 0 ~ 1
	\return			interpolated value
*/
template <typename T>
inline PROMOTE_T_TO_FLOAT interpCubicBezier(T val0, T val1, T val2, T val3, double t){
	return	val0 * (1-t)*(1-t)*(1-t) +
			val1 * 3*(1-t)*(1-t)*t + 
			val2 * 3*(1-t)*t*t +
			val3 * t*t*t;
}

///@}

//	========================================
///@{
/**	@name 2D Interpolation
*/
//	========================================
/**
	Nearest Interpolation on Unit Square

	^ y
	|   01  11
	|   00  10
	o-------->  x

	\param	val00	value at position (0, 0)
	\param	val10	value at position (1, 0)
	\param	val01	value at position (0, 1)
	\param	val11	value at position (1, 1)
	\param	tx		interpolation position on x axis, should be 0 ~ 1
	\param	ty		interpolation position on y axis, should be 0 ~ 1
	\return			interpolated value
*/
template <typename T>
inline T interpNearest(T val00, T val10, T val01, T val11, double tx, double ty){
	if( tx<0.5 && ty<0.5 )
		return val00;
	else if( tx<0.5 )
		return val01;
	else if( ty<0.5 )
		return val10;
	else
		return val11;
}

/**
	Bilinear Interpolation on Unit Square

	^ y
	|   01  11
	|   00  10
	o--------> x

	\param	val00	value at position (0, 0)
	\param	val10	value at position (1, 0)
	\param	val01	value at position (0, 1)
	\param	val11	value at position (1, 1)
	\param	tx		interpolation position on x axis, should be 0 ~ 1
	\param	ty		interpolation position on y axis, should be 0 ~ 1
	\return			interpolated value
*/
template <typename T>
inline PROMOTE_T_TO_FLOAT interpBilinear(T val00, T val10, T val01, T val11, double tx, double ty){
	return (1.0-tx)*(1.0-ty)*val00 + tx*(1.0-ty)*val10 + (1.0-tx)*ty*val01 + tx*ty*val11;
}

/**
	Bicubic Interpolation on Unit Square, square(11-22)

	^ y
	|	03	13	23	33
	|	02	12	22	32
	|   01  11	21	31
	|   00  10	20	30
	o------------------> x

	\param	val00	value at position (0, 0)
	\param	val10	value at position (1, 0)
	\param	val20	value at position (2, 0)
	\param	val30	value at position (3, 0)
	\param	val01	value at position (0, 1)
	\param	val11	value at position (1, 1)
	\param	val21	value at position (2, 1)
	\param	val31	value at position (3, 1)
	\param	val02	value at position (0, 2)
	\param	val12	value at position (1, 2)
	\param	val22	value at position (2, 2)
	\param	val32	value at position (3, 2)
	\param	val03	value at position (0, 3)
	\param	val13	value at position (1, 3)
	\param	val23	value at position (2, 3)
	\param	val33	value at position (3, 3)
	\param	tx		interpolation position on x axis, should be 0 ~ 1
	\param	ty		interpolation position on y axis, should be 0 ~ 1
	\return			interpolated value
*/
template <typename T>
inline PROMOTE_T_TO_FLOAT interpBicubic(T val00, T val10, T val20, T val30, 
										T val01, T val11, T val21, T val31, 
										T val02, T val12, T val22, T val32, 
										T val03, T val13, T val23, T val33, 
										double tx, double ty){
	return interpCubic(
		interpCubic(val00, val10, val20, val30, tx),
		interpCubic(val01, val11, val21, val31, tx),
		interpCubic(val02, val12, val22, val32, tx),
		interpCubic(val03, val13, val23, val33, tx),
		ty );
}

///@}

//	========================================
///@{
/**	@name Triangle Interpolation
*/
//	========================================

/**
	Cubic Bezier Interpolation on Triangle

	given 10 sampling points, interpolate the triangle at given barcentric coordinate

	sequence of control points on the triangle are dipicted in the following

	b003	\n
	b102	b012	\n
	b201	b111	b021	\n
	b300	b210	b120	b030	\n

	\param	b300	value in corner 0
	\param	b030	value in corner 1
	\param	b003	value in corner 2
	\param	b210	value on edge 0-1, closer to 0
	\param	b120	value on edge 0-1, closer to 1
	\param	b021	value on edge 1-2, closer to 1
	\param	b012	value on edge 1-2, closer to 2
	\param	b102	value on edge 0-2, closer to 2
	\param	b201	value on edge 0-2, closer to 0
	\param	b111	value at center
	\param	u		barcentric coordinate, coef on b300-b030
	\param	v		barcentric coordinate, coef on b300-b003
	\return			interpolated value
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT interpCubicBezierTriangle(T b300,	T b030,	T b003,
													T b210,	T b120, T b021, T b012, T b102, T b201,
													T b111,
													double u, double v){
	double	t = 1 - u - v;
	double	t2 = t*t;
	double	u2 = u*u;
	double	v2 = v*v;
	return b300*t2*t + b030*u2*u + b003*v2*v +
		3*(b210*t2*u + b120*t*u2 + b021*u2*v + b012*u*v2 + b102*v2*t + b201*v*t2) +
		6*b111*t*u*v;
}

/**
	Cubic Bezier Interpolation on Triangle

	given 10 sampling points, interpolate the triangle at given barcentric coordinate

	sequence of control points on the triangle are dipicted in the following

	b003	\n
	b102	b012	\n
	b201	b111	b021	\n
	b300	b210	b120	b030	\n

	\param	data_on_control_points	data on control points arranged in sequence of 
					b300, b030, b003, b210, b120, b021, b012, b102, b201, b111
	\param	u		barcentric coordinate, coef on b300-b030
	\param	v		barcentric coordinate, coef on b300-b003
	\return			interpolated value
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT interpCubicBezierTriangle(T* data_on_control_points, double u, double v){
	return interpCubicBezierTriangle(
		data_on_control_points[0], data_on_control_points[1], data_on_control_points[2],
		data_on_control_points[3], data_on_control_points[4], data_on_control_points[5],
		data_on_control_points[6], data_on_control_points[7], data_on_control_points[8],
		data_on_control_points[9],
		u, v);
}

///@}

}	//	namespace yz

#endif	//	__YZ_INTERPOLATION_H__