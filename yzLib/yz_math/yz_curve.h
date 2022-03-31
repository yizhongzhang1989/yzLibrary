/***********************************************************/
/**	\file
	\brief		Functions of Curves
	\author		Yizhong Zhang
	\date		11/3/2012
*/
/***********************************************************/
#ifndef __YZ_CURVE_H__
#define __YZ_CURVE_H__

#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_interpolation.h"

namespace yz{

//	========================================
///@{
/**	@name	Length of Curves
*/
//	========================================

/**
	calculate length of 2D curve

	the curve is represented by a list of points, 
	each line segment is treated as straight line

	\param	curve			key points of the curve
	\param	point_number	number of key points
	\return					the length of the whole curve
*/
template<typename T>
inline T calculateCurveLength2D(const T* curve, int point_number){
	if( point_number < 2 )	//	too short, cannot form a curve
		return 0;
	T length = 0;
	for( int i=0; i<point_number-1; i++ ){
		Vec2<T> r = Vec2<T>(curve+i*2) - Vec2<T>(curve+i*2+2);
		length += r.Length();
	}
	return length;
}

/**
	calculate length of 3D curve

	the curve is represented by a list of points, 
	each line segment is treated as straight line

	\param	curve			key points of the curve
	\param	point_number	number of key points
	\return					the length of the whole curve
*/
template<typename T>
inline T calculateCurveLength3D(const T* curve, int point_number){
	if( point_number < 2 )	//	too short, cannot form a curve
		return 0;
	T length = 0;
	for( int i=0; i<point_number-1; i++ ){
		Vec3<T> r = Vec3<T>(curve+i*3) - Vec3<T>(curve+i*3+3);
		length += r.Length();
	}
	return length;
}

/**
	calculate length of curve

	the curve is represented by a list of points, 
	each line segment is treated as straight line

	\param	curve			key points of the curve
	\return					the length of the whole curve
*/
template<typename T>
inline T calculateCurveLength(const std::vector<Vec2<T>>& curve){
	if( curve.empty() )
		return 0;
	return calculateCurveLength2D((T*)&curve[0], curve.size());
}

/**
	calculate length of curve

	the curve is represented by a list of points, 
	each line segment is treated as straight line

	\param	curve			key points of the curve
	\return					the length of the whole curve
*/
template<typename T>
inline T calculateCurveLength(const std::vector<Vec3<T>>& curve){
	if( curve.empty() )
		return 0;
	return calculateCurveLength3D((T*)&curve[0], curve.size());
}

/**
	calculate length of cubic bezier

	\param	v0		position of point 0
	\param	v1		position of point 1
	\param	v2		position of point 2
	\param	v3		position of point 3
	\param	slices	how many slices is the curve to be discretized
	\return			the length of the whole curve
*/
template<typename T>
inline T calculateCubicBezierLength(Vec2<T> v0, Vec2<T> v1, Vec2<T> v2, Vec2<T> v3, int slices = 16){
	std::vector<Vec2<T>> point;
	point.resize(slices+1);
	point[0] = v0;
	point[slices] = v3;
	for(int i=1; i<slices; i++){
		point[i] = interpCubicBezier(v0, v1, v2, v3, double(i)/slices);
	}
	return calculateCurveLength(point);
}

/**
	calculate length of cubic bezier

	\param	v0		position of point 0
	\param	v1		position of point 1
	\param	v2		position of point 2
	\param	v3		position of point 3
	\param	slices	how many slices is the curve to be discretized
	\return			the length of the whole curve
*/
template<typename T>
inline T calculateCubicBezierLength(Vec3<T> v0, Vec3<T> v1, Vec3<T> v2, Vec3<T> v3, int slices = 16){
	std::vector<Vec3<T>> point;
	point.resize(slices+1);
	point[0] = v0;
	point[slices] = v3;
	for(int i=1; i<slices; i++){
		point[i] = interpCubicBezier(v0, v1, v2, v3, double(i)/slices);
	}
	return calculateCurveLength(point);
}


///@}

//	========================================
///@{
/**	@name	Area of Polygon
*/
//	========================================

/**
	calculate area of 2D polygon represented by surrounding points

	the polygon is expected to be self-intersection free

	\param	point_ptr		point list in xy_xy format
	\param	point_number	number of points
	\return					the area calculated, positive: counter-clockwise, negetive: clockwise
*/
template<typename T>
inline T calculatePolygonArea2D(const T* point_ptr, int point_number){
	T	area = 0;
	for(int i=0, j = point_number-1; i<point_number; i++) {
		area += (point_ptr[j*2]+point_ptr[i*2]) * (point_ptr[j*2+1]-point_ptr[i*2+1]); 
		j = i; 
	}
	return area * -0.5; 
}

/**
	calculate area of 2D polygon represented by surrounding points

	the polygon is expected to be self-intersection free

	\param	point		point list
	\return				the area calculated, positive: counter-clockwise, negetive: clockwise
*/
template<typename T>
inline T calculatePolygonArea(const std::vector<Vec2<T>>& point){
	if( point.empty() )
		return 0;
	return calculatePolygonArea2D((T*)&point[0], point.size());
}


///@}

}	//	namespace yz

#endif	//	__YZ_CURVE_H__