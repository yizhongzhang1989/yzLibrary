/***********************************************************/
/**	\file
	\brief		Intersection Test of Geometry Elements
	\author		Yizhong Zhang
	\date		6/12/2012
*/
/***********************************************************/
#ifndef __YZ_INTERSECTION_TEST_H__
#define __YZ_INTERSECTION_TEST_H__

#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_matrix.h"
#include "yzLib/yz_geometry/yz_clipping.h"

namespace yz{	namespace geometry{
//	========================================
///@{
/**	@name Line-Point Projection
*/
//	========================================

/**
	Get Projection Coefficient of Line and Point

	\param	coef_prj	coefficient of projection point with v0-v1
	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
*/
template<typename T, typename TV>
inline void getLinePointProjectionCoef(T& coef_prj,
									   Vec2<TV> linev0,
									   Vec2<TV> linev1,
									   Vec2<TV> point ){
	Vec2<TV> r = linev1 - linev0;
	Vec2<TV> rp = point - linev0;

	T len = r.Length();
	coef_prj = dot(rp, r) / (len*len);
}

/**
	Get Projection Point of Line and Point

	return value doesn't affect write result to prj_point

	\param	prj_point	the projection point of the point on the line
	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
	\return				whether the projection point is on the line segment
*/
template<typename T>
inline int getLinePointProjectionPoint(Vec2<T>& prj_point,
									   Vec2<T> linev0,
									   Vec2<T> linev1,
									   Vec2<T> point ){
	T prj_coef;
	getLinePointProjectionCoef(prj_coef, linev0, linev1, point);

	prj_point = (linev1-linev0) * prj_coef + linev0;

	if( prj_coef>=0 && prj_coef<=1 )
		return 1;
	else
		return 0;
}

/**
	get the distance of point and line

	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
	\return				the distance between point and line
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT getLinePointDistance(Vec2<T> linev0,
											   Vec2<T> linev1,
											   Vec2<T> point )
{
	float line_len = (linev1 - linev0).Length();

	if( line_len < 1e-5f ){
		//	the line is so short that it reduced to a point
		return (point - linev0).Length();
	}
	else{
		return fabs( (linev1.x-linev0.x)*(linev0.y-point.y) -
			(linev0.x-point.x)*(linev1.y-linev0.y) ) / line_len;
	}
}

/**
	get the distance of point and a line segment. 

	If the projection point on the line is not inside the line segment,
	then the distance is the point to the nearest end point of the line segment

	\param	segv0		coordinate of v0 of line segment
	\param	segv1		coordinate of v1 of line segment
	\param	point		coordinate of the point
	\return				the distance between point and line segment
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT getSegmentPointDistance(Vec2<T> segv0,
												  Vec2<T> segv1,
												  Vec2<T> point ){
	T lamda;
	getLinePointProjectionCoef(lamda, segv0, segv1, point);
	if( lamda <= 0 )
		return (point - segv0).Length();
	else if( lamda >= 1 )
		return (point - segv1).Length();
	else{
		Vec2<T> prj_point = (segv1 - segv0) * lamda + segv0;
		return (prj_point-point).Length();
	}
}

/**
	Get Projection Coefficient of Line and Point in 3D space

	\param	coef_prj	coefficient of projection point with v0-v1
	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
*/
template<typename T, typename TV>
inline void getLinePointProjectionCoef(T& coef_prj,
									   Vec3<TV> linev0,
									   Vec3<TV> linev1,
									   Vec3<TV> point ){
	Vec3<TV> r = linev1 - linev0;
	Vec3<TV> rp = point - linev0;

	T len = r.Length();
	coef_prj = dot(rp, r) / (len*len);
}

/**
	Get Projection Point of Line and Point in 3D space

	return value doesn't affect write result to prj_point

	\param	prj_point	the projection point of the point on the line
	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
	\return				whether the projection point is on the line segment
*/
template<typename T>
inline int getLinePointProjectionPoint(Vec3<T>& prj_point,
									   Vec3<T> linev0,
									   Vec3<T> linev1,
									   Vec3<T> point ){
	T prj_coef;
	getLinePointProjectionCoef(prj_coef, linev0, linev1, point);

	prj_point = (linev1-linev0) * prj_coef + linev0;

	if( prj_coef>=0 && prj_coef<=1 )
		return 1;
	else
		return 0;
}



/**
	get the distance of point and line

	If the projection point on the line is not inside the line segment,
	then the distance is the point to the nearest end point of the line segment

	\param	linev0		coordinate of v0 of line
	\param	linev1		coordinate of v1 of line
	\param	point		coordinate of the point
	\return				the distance between point and line
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT getLinePointDistance(Vec3<T> linev0,
											   Vec3<T> linev1,
											   Vec3<T> point ){
	Vec3<T> prj_point;
	getLinePointProjectionPoint(prj_point, linev0, linev1, point);
	return (prj_point-point).Length();
}


/**
	get the distance of point and line segment

	\param	segv0		coordinate of v0 of line segment
	\param	segv1		coordinate of v1 of line segment
	\param	point		coordinate of the point
	\return				the distance between point and line segment
*/
template<typename T>
inline PROMOTE_T_TO_FLOAT getSegmentPointDistance(Vec3<T> segv0,
												  Vec3<T> segv1,
												  Vec3<T> point ){
	Vec3<T> nearest_point;
	getNearestPointOnSegment(nearest_point, point, segv0, segv1);
	return (nearest_point - point).Length();
}

/**
	get the nearest point on a line segment to the point

	\param	nearest_point	the nearest point on the line segment from the point
	\param	segv0			coordinate of v0 of line segment
	\param	segv1			coordinate of v1 of line segment
	\param	point			coordinate of the point
	\return					1:	the nearest point is on the line segment, not including the end points
							0:	the nearest point is the end point
*/
template<typename T>
inline int getNearestPointOnSegment(Vec2<T>&	nearest_point,
									Vec2<T>		point,
									Vec2<T>		segv0,
									Vec2<T>		segv1 ){
	T lamda;
	getLinePointProjectionCoef(lamda, segv0, segv1, point);
	if( lamda <= 0 ){
		nearest_point = segv0;
		return 0;
	}
	else if( lamda >= 1 ){
		nearest_point = segv1;
		return 0;
	}
	else{
		nearest_point = (segv1 - segv0) * lamda + segv0;
		return 1;
	}
}
/**
	get the nearest point on a line segment to the point

	\param	nearest_point	the nearest point on the line segment from the point
	\param	segv0			coordinate of v0 of line segment
	\param	segv1			coordinate of v1 of line segment
	\param	point			coordinate of the point
	\return					1:	the nearest point is on the line segment, not including the end points
							0:	the nearest point is the end point
*/
template<typename T>
inline int getNearestPointOnSegment(Vec3<T>&	nearest_point,
									Vec3<T>		point,
									Vec3<T>		segv0,
									Vec3<T>		segv1 ){
	T lamda;
	getLinePointProjectionCoef(lamda, segv0, segv1, point);
	if( lamda <= 0 ){
		nearest_point = segv0;
		return 0;
	}
	else if( lamda >= 1 ){
		nearest_point = segv1;
		return 0;
	}
	else{
		nearest_point = (segv1 - segv0) * lamda + segv0;
		return 1;
	}
}
///@}

//	========================================
///@{
/**	@name Line-Line Intersection and Projection
*/
//	========================================

/**
	Get intersection coefficient of two lines 2D

	P0 = coef0 * (P2 - P1) = coef1 * (P4 - P3)

	If two lines are not parallel, they must intersect.

	Lines are represented by two points on that line. The line must be
	legal (distinct end points), or error may accour calling this function.

	\param	coef0			return value, P0 = coef0 * (P2 - P1)
	\param	coef1			return value, P0 = coef1 * (P4 - P3)
	\param	line0v0			P1, end point 0 of line0
	\param	line0v1			P2, end point 1 of line0
	\param	line1v0			P3, end point 0 of line1
	\param	line1v1			P4, end point 1 of line1
	\return					whether the two lines intersect, 0: don't intersect, 1: intersect
*/
template<typename T, typename TV>
inline int getLineLineIntersectionCoef(T& coef0, 
									   T& coef1, 
									   Vec2<TV> line0v0, 
									   Vec2<TV> line0v1, 
									   Vec2<TV> line1v0, 
									   Vec2<TV> line1v1){
	Vec2<TV>	r0 = line0v1 - line0v0;
	Vec2<TV>	r1 = line1v1 - line1v0;

	PROMOTE_T_TO_FLOAT denominator = cross(r0, r1);
	if( denominator < 1e-6 && denominator > -1e-6 )	//	the two line segments are parellel
		return 0;

	coef0 = ((line1v1[0]-line1v0[0])*(line0v0[1]-line1v0[1])-(line1v1[1]-line1v0[1])*(line0v0[0]-line1v0[0])) / denominator;
	coef1 = ((line0v1[0]-line0v0[0])*(line0v0[1]-line1v0[1])-(line0v1[1]-line0v0[1])*(line0v0[0]-line1v0[0])) / denominator;

	return 1;
}

/**
	Get intersection point of two lines 2D

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	line0v0					end point 0 of line0
	\param	line0v1					end point 1 of line0
	\param	line1v0					end point 0 of line1
	\param	line1v1					end point 1 of line1
	\return							whether the two lines intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getLineLineIntersectionPoint(Vec2<T>& intersection_point, 
										Vec2<T> line0v0, 
										Vec2<T> line0v1, 
										Vec2<T> line1v0, 
										Vec2<T> line1v1){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, line0v0, line0v1, line1v0, line1v1);

	if( intersect_flag ){
		intersection_point = line0v0 + coef0 * (line0v1-line0v0);
		return 1;
	}

	return 0;
}

/**
	Get intersection point of line and ray

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	linev0					end point 0 of line
	\param	linev1					end point 1 of line
	\param	ray_origin				origin of ray
	\param	ray_next				next point of ray, next = origin + direction
	\return							whether they intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getLineRayIntersectionPoint(Vec2<T>& intersection_point, 
									   Vec2<T> linev0, 
									   Vec2<T> linev1, 
									   Vec2<T> ray_origin, 
									   Vec2<T> ray_next){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, linev0, linev1, ray_origin, ray_next);

	if( intersect_flag && coef1>=0 ){
		intersection_point = ray_origin + coef1 * (ray_next - ray_origin);
		return 1;
	}

	return 0;
}

/**
	Get intersection point of a line and a line segment 2D

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	linev0					end point 0 of line
	\param	linev1					end point 1 of line
	\param	segv0					end point 0 of line segment
	\param	segv1					end point 1 of line segment
	\return							whether they intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getLineSegmentIntersectionPoint(Vec2<T>& intersection_point, 
										   Vec2<T> linev0, 
										   Vec2<T> linev1, 
										   Vec2<T> segv0, 
										   Vec2<T> segv1){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, linev0, linev1, segv0, segv1);

	if( intersect_flag && coef1>=0 && coef1<=1 ){
		intersection_point = segv0 + coef1 * (segv1-segv0);
		return 1;
	}

	return 0;
}

/**
	Get intersection point of two rays 2D

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	ray0_origin				start point of ray0
	\param	ray0_next				next point of ray0, next = origin + direction
	\param	ray1_origin				start point of ray1
	\param	ray1_next				next point of ray1, next = origin + direction
	\return							whether they intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getRayRayIntersectionPoint(Vec2<T>& intersection_point, 
									  Vec2<T> ray0_origin, 
									  Vec2<T> ray0_next, 
									  Vec2<T> ray1_origin, 
									  Vec2<T> ray1_next){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, 
		ray0_origin, ray0_next, ray1_origin, ray1_next);

	if( intersect_flag && coef0>=0 && coef1>=0 ){
		intersection_point = ray1_origin + coef1 * (ray1_next - ray1_origin);
		return 1;
	}

	return 0;
}

/**
	Get intersection point of ray and a line segment 2D

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	ray_origin				start point of the ray
	\param	ray_next				next point of ray, next = origin + direction
	\param	segv0					end point 0 of line segment
	\param	segv1					end point 1 of line segment
	\return							whether they intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getRaySegmentIntersectionPoint(Vec2<T>& intersection_point, 
										  Vec2<T> ray_origin, 
										  Vec2<T> ray_next,
										  Vec2<T> segv0, 
										  Vec2<T> segv1){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, ray_origin, ray_next, segv0, segv1);

	if( intersect_flag && coef0>=0 && coef1>=0 && coef1<=1 ){
		intersection_point = ray_origin + coef0 * (ray_next - ray_origin);
		return 1;
	}

	return 0;
}

/**
	Get intersection point of two line segments 2D

	\param	intersection_point		return value, if don't intersect, don't change its value
	\param	seg0v0					end point 0 of line segment0
	\param	seg0v1					end point 1 of line segment0
	\param	seg1v0					end point 0 of line segment1
	\param	seg1v1					end point 1 of line segment1
	\return							whether the two segments intersect, 0: don't intersect, 1: intersect
*/
template<typename T>
inline int getSegmentSegmentIntersectionPoint(Vec2<T>& intersection_point, 
											  Vec2<T> seg0v0, 
											  Vec2<T> seg0v1, 
											  Vec2<T> seg1v0, 
											  Vec2<T> seg1v1){
	PROMOTE_T_TO_FLOAT coef0, coef1;
	int intersect_flag = getLineLineIntersectionCoef(coef0, coef1, seg0v0, seg0v1, seg1v0, seg1v1);

	if( intersect_flag && coef0>=0 && coef0<=1 && coef1>=0 && coef1<=1 ){
		intersection_point = seg0v0 + coef0 * (seg0v1-seg0v0);
		return 1;
	}

	return 0;
}


/**
	Get projection coefficient of two lines in 3D space

	If two lines don't intersect and parallel to each other,
	they must have common perpendicular line segment. We set 
	P0-P1 to be common perpendicular line segment of the two lines,
	and P0 is on line0, P1 is on line1.

	\param	coef0			return value, P0 = coef0 * (P2 - P1)
	\param	coef1			return value, P1 = coef1 * (P4 - P3)
	\param	line0v0			P1, end point 0 of line0
	\param	line0v1			P2, end point 1 of line0
	\param	line1v0			P3, end point 0 of line1
	\param	line1v1			P4, end point 1 of line1
	\return					whether the two lines are parallel, 0: parallel, 1: not parallel
*/
template<typename T, typename TV>
inline int getLineLineProjectionCoef(T& coef0, 
									 T& coef1, 
									 Vec3<TV> line0v0, 
									 Vec3<TV> line0v1, 
									 Vec3<TV> line1v0, 
									 Vec3<TV> line1v1){
	Vec3<TV>	r0 = line0v1 - line0v0;
	Vec3<TV>	r1 = line1v1 - line1v0;

	Vec3<TV> com_p = cross(r0, r1);
	double len = com_p.Length();
	if( len < 1e-6 && len > -1e-6 )	//	the two lines are parallel
		return 0;

	Matrix3x3<T> mat(-r0, r1, -com_p);
	mat.SetInverse();

	Vec3<T>	coef_vec = mat * (line0v0 - line1v0);

	coef0 = coef_vec[0];
	coef1 = coef_vec[1];

	return 1;
}

/**
	Get Common Perpendicular of two lines in 3D space

	If two lines in 3D space are not parallel, they must intersect 
	or have common perpendicular line segment. This function return
	the two end points of the common perpendicular line segment.

	\param	p0				common perpendicular end on line0
	\param	p1				common perpendicular end on line1
	\param	line0v0			P1, end point 0 of line0
	\param	line0v1			P2, end point 1 of line0
	\param	line1v0			P3, end point 0 of line1
	\param	line1v1			P4, end point 1 of line1
	\return					whether the two lines are parallel, 0: parallel, 1: not parallel
*/
template<typename T>
inline int getLineLineCommonPerpendicular(Vec3<T>& p0,
										  Vec3<T>& p1,
										  Vec3<T> line0v0, 
										  Vec3<T> line0v1, 
										  Vec3<T> line1v0, 
										  Vec3<T> line1v1){
	T coef0, coef1;
	int ret_val = getLineLineProjectionCoef(coef0, coef1, line0v0, line0v1, line1v0, line1v1);
	if( ret_val ){
		p0 = line0v0 + (line0v1 - line0v0) * coef0;
		p1 = line1v0 + (line1v1 - line1v0) * coef1;
	}
	return ret_val;
}
///@}

//	========================================
///@{
/**	@name Point-Triangle Intersection Projection
*/
//	========================================

/**
	Get point coordinate in triangle coordinate system.

	Given coordinate of three vertices of a triangle and a point in
	global coordinate system, get the coordinate of the point in
	local coordinate system of the triangle. Let v0 be the original 
	point, v1-v0 be x-axis, v2-v0 be y-axis.

	If coef0 > 0 && coef1 > 0 && coef0+coef1 < 1, the point is strictly inside the triangle.

	\param	coef1		return local coordinate of point on v1-v0
	\param	coef2		return local coordinate of point on v2-v0
	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
*/
template<typename T, typename TV>
inline void getPointTriangleIntersectionCoef(T& coef1,
											 T& coef2,
											 Vec2<TV> point,
											 Vec2<TV> v0,
											 Vec2<TV> v1,
											 Vec2<TV> v2 ){
	Vec2<T> r1	= v1 - v0;
	Vec2<T> r2	= v2 - v0;
	Vec2<T> r	= point - v0;

	Matrix2x2<T> mat(r1, r2);

	T det = mat.Det();
	if (det != 0){
		//	the triangle is regular
		mySwap(mat[0][0], mat[1][1]);
		mat[0][0] /= det;
		mat[0][1] /= -det;
		mat[1][0] /= -det;
		mat[1][1] /= det;

		Vec2<T> coef_vec = mat * r;
		coef1 = coef_vec[0];
		coef2 = coef_vec[1];
	}
	else{
		//	the triangle is irregular, return NaN to coefficient
		coef1 = std::numeric_limits<T>::quiet_NaN();
		coef2 = std::numeric_limits<T>::quiet_NaN();
	}

}

/**
	Is a point inside a triangle, including on boundary

	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				test result, 1: is inside triangle; 0: not inside triangle
*/
template<typename T>
inline int isPointInsideTriangle(Vec2<T> point,
								 Vec2<T> v0,
								 Vec2<T> v1,
								 Vec2<T> v2	){
	T coef1, coef2;
	getPointTriangleIntersectionCoef(coef1, coef2, point, v0, v1, v2);
	if( coef1>=0 && coef2>=0 && coef1+coef2<=1 )
		return 1;
	else
		return 0;
}

/**
	Is a point inside a triangle, not including on boundary

	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				test result, 1: is inside triangle; 0: not inside triangle
*/
template<typename T>
inline int isPointStrictlyInsideTriangle(Vec2<T> point,
										 Vec2<T> v0,
										 Vec2<T> v1,
										 Vec2<T> v2	){
	T coef1, coef2;
	getPointTriangleIntersectionCoef(coef1, coef2, point, v0, v1, v2);
	if( coef1>0 && coef2>0 && coef1+coef2<1 )
		return 1;
	else
		return 0;
}

/**
	Project the point on the triangle plane, and return coordinate
	of the projected point in triangle local coordinate system.

	Given coordinate of three vertices of a triangle and a point in
	global coordinate system, get the coordinate of the point in
	local coordinate system of the triangle. Let v0 be the original 
	point, v1-v0 be x-axis, v2-v0 be y-axis.

	If coef0 > 0 && coef1 > 0 && coef0+coef1 < 1, the point is strictly inside the triangle.

	\param	coef1		return local coordinate of point on v1-v0
	\param	coef2		return local coordinate of point on v2-v0
	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				return which side does the point lie compared with the triangle.
						1: on the normal size; -1: behind the normal side
*/
template<typename T, typename TV>
inline int getPointTriangleProjectionCoef(T& coef1,
										  T& coef2,
										  Vec3<TV> point,
										  Vec3<TV> v0,
										  Vec3<TV> v1,
										  Vec3<TV> v2 ){
	Vec3<T> r1	= v1 - v0;
	Vec3<T> r2	= v2 - v0;
	Vec3<T> nor = cross(r1, r2).Normalize();
	Vec3<T> r	= point - v0;

	Matrix3x3<T> mat(r1, r2, nor);
	Vec3<T> coef_vec = mat.Inverse() * r;
	coef1 = coef_vec[0];
	coef2 = coef_vec[1];
	if( coef_vec[2] >= 0 )
		return 1;
	else
		return -1;
}

/**
	Project the 3D point onto the plane of the triangle

	\param	prj_point	return the projected point
	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				whether the projected point is inside the triangle(including boundary)
						1: inside the triangle; 0: not inside the triangle
*/
template<typename T>
inline int getPointTriangleProjectionPoint(Vec3<T>& prj_point,
										   Vec3<T> point,
										   Vec3<T> v0,
										   Vec3<T> v1,
										   Vec3<T> v2){
	T coef1, coef2;
	getPointTriangleProjectionCoef(coef1, coef2, point, v0, v1, v2);
	prj_point = v0 + (v1-v0) * coef1 + (v2-v0) * coef2;

	if( coef1>=0 && coef2>=0 && coef1+coef2<=1 )
		return 1;
	else
		return 0;
}

/**
	Get the nearest point of a point to a triangle. If the projection 
	point is inside the triangle, then the nearest point is the projection
	point. Otherwise, the nearest point is on the boundary of the triangle

	\param	nearest_point	return the projected point
	\param	point			coordinate of the point in global coordinate system
	\param	v0				coordinate of triangle v0 in global coordinate system
	\param	v1				coordinate of triangle v1 in global coordinate system
	\param	v2				coordinate of triangle v2 in global coordinate system
	\return					whether the projected point is inside the triangle(including boundary)
							1: inside the triangle; 0: not inside the triangle
*/
template<typename T>
inline int getNearestPointOnTriangle(Vec3<T>&	nearest_point,
									 Vec3<T>	point,
									 Vec3<T>	v0,
									 Vec3<T>	v1,
									 Vec3<T>	v2){
	T coef1, coef2;
	getPointTriangleProjectionCoef(coef1, coef2, point, v0, v1, v2);

	//	projection point is inside of triangle
	if( coef1>=0 && coef2>=0 && coef1+coef2<=1 ){
		nearest_point = v0 + (v1-v0) * coef1 + (v2-v0) * coef2;
		return 1;
	}

	//	projection point is outside of triangle
	Vec3<T>	tmp[3];
	getNearestPointOnSegment(tmp[0], point, v0, v1);
	getNearestPointOnSegment(tmp[1], point, v0, v2);
	getNearestPointOnSegment(tmp[2], point, v1, v2);

	T squ_dist[3];
	for(int i=0; i<3; i++)
		squ_dist[i] = (point - tmp[i]).SquareLength();

	if( squ_dist[0]<=squ_dist[1] && squ_dist[0]<=squ_dist[2] )
		nearest_point = tmp[0];
	else if( squ_dist[1]<=squ_dist[0] && squ_dist[1]<=squ_dist[2] )
		nearest_point = tmp[1];
	else
		nearest_point = tmp[2];

	return 0;
}

/**
	Whether the projection of the point onto the triangle plane is inside the 
	triangle, including the boundary

	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				test result, 1: is inside triangle; 0: not inside triangle
*/
template<typename T>
inline int isPointProjectionInsideTriangle(Vec3<T> point,
										   Vec3<T> v0,
										   Vec3<T> v1,
										   Vec3<T> v2	){
	T coef1, coef2;
	getPointTriangleProjectionCoef(coef1, coef2, point, v0, v1, v2);
	if( coef1>=0 && coef2>=0 && coef1+coef2<=1 )
		return 1;
	else
		return 0;
}

/**
	Whether the projection of the point onto the triangle plane is strictly 
	inside the triangle, including the boundary

	\param	point		coordinate of the point in global coordinate system
	\param	v0			coordinate of triangle v0 in global coordinate system
	\param	v1			coordinate of triangle v1 in global coordinate system
	\param	v2			coordinate of triangle v2 in global coordinate system
	\return				test result, 1: is inside triangle; 0: not inside triangle
*/
template<typename T>
inline int isPointProjectionStrictlyInsideTriangle(Vec3<T> point,
												   Vec3<T> v0,
												   Vec3<T> v1,
												   Vec3<T> v2	){
	T coef1, coef2;
	getPointTriangleProjectionCoef(coef1, coef2, point, v0, v1, v2);
	if( coef1>0 && coef2>0 && coef1+coef2<1 )
		return 1;
	else
		return 0;
}

///@}

//	========================================
///@{
/**	@name Line-Triangle Intersection
*/
//	========================================

/**
	Get intersection point of a line and the plane, the plane is 
	represented using a point on the plane and the normal of the plane

	\param	intersection_point	return the intersection point
	\param	linev0				vertex 0 on the line
	\param	linev1				vertex 1 on the line
	\param	plane_v				vertex on the plane
	\param	plane_normal		normal of the plane
	\return						whether the line and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getLinePlaneIntersectionPoint(Vec3<T>& intersection_point,
										 Vec3<T> linev0,
										 Vec3<T> linev1,
										 Vec3<T> plane_v,
										 Vec3<T> plane_normal ){
	Vec3<T> u1 = linev0 - plane_v;
	Vec3<T> u2 = linev1 - plane_v;
	T c1 = dot(u1, plane_normal);
	T c2 = dot(u2, plane_normal);
	if( abs(c1-c2) < 1e-6 )	//	the line is almost parallel to the plane
		return 0;			//	there is no intersection point

	intersection_point = (u2 - u1)*(c1/(c1-c2)) + u1 +  plane_v;
	return 1;
}


/**
	Get intersection point of a line and the plane, the plane 
	is represented using three non-colinear vertices on the plane

	\param	intersection_point	return the intersection point
	\param	linev0				vertex 0 on the line
	\param	linev1				vertex 1 on the line
	\param	planev0				vertex 0 on the plane
	\param	planev1				vertex 1 on the plane
	\param	planev2				vertex 2 on the plane
	\return						whether the line and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getLinePlaneIntersectionPoint(Vec3<T>& intersection_point,
										 Vec3<T> linev0,
										 Vec3<T> linev1,
										 Vec3<T> planev0,
										 Vec3<T> planev1,
										 Vec3<T> planev2 ){
	Vec3<T> r1 = planev1 - planev0;
	Vec3<T> r2 = planev2 - planev0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	return getLinePlaneIntersectionPoint(intersection_point, linev0, linev1, planev0, nor);
}

/**
	Get intersection point of a ray and the plane, the plane is represented 
	using a point on the plane and the normal of the plane

	\param	intersection_point	return the intersection point
	\param	ray_origin			start point of a ray
	\param	ray_next			next point of a ray
	\param	plane_v				vertex on the plane
	\param	plane_normal		normal of the plane
	\return						whether the ray and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getRayPlaneIntersectionPoint(Vec3<T>& intersection_point,
										Vec3<T> ray_origin,
										Vec3<T> ray_next,
										Vec3<T> plane_v,
										Vec3<T> plane_normal ){
	Vec3<T> u1 = ray_origin - plane_v;
	Vec3<T> u2 = ray_next - plane_v;
	T c1 = dot(u1, plane_normal);
	T c2 = dot(u2, plane_normal);
	if( (c1>0 && c2>0 && c2>c1) || (c1<0 && c2<0 && c2<c1) || abs(c1-c2) < 1e-6 )
		return 0;	//	the ray is not pointing towards the plane, or the ray is parallel to the plane

	intersection_point = (u2 - u1)*(c1/(c1-c2)) + u1 +  plane_v;
	return 1;
}

/**
	Get intersection point of a ray and the plane, the plane 
	is represented using three non-coliner vertices on the plane

	\param	intersection_point	return the intersection point
	\param	ray_origin			start point of a ray
	\param	ray_next			next point of a ray
	\param	planev0				vertex 0 on the plane
	\param	planev1				vertex 1 on the plane
	\param	planev2				vertex 2 on the plane
	\return						whether the ray and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getRayPlaneIntersectionPoint(Vec3<T>& intersection_point,
										Vec3<T> ray_origin,
										Vec3<T> ray_next,
										Vec3<T> planev0,
										Vec3<T> planev1,
										Vec3<T> planev2 ){
	Vec3<T> r1 = planev1 - planev0;
	Vec3<T> r2 = planev2 - planev0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	return getRayPlaneIntersectionPoint(intersection_point, ray_origin, ray_next, planev0, nor);
}

/**
	Get intersection point of a line segment and the plane, the plane 
	is represented using a point on the plane and the normal of the plane

	\param	intersection_point	return the intersection point
	\param	segv0				end point 0 of line segment
	\param	segv1				end point 1 of line segment
	\param	plane_v				vertex on the plane
	\param	plane_normal		normal of the plane
	\return						whether the line segment and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getSegmentPlaneIntersectionPoint(Vec3<T>& intersection_point,
											Vec3<T> segv0,
											Vec3<T> segv1,
											Vec3<T> plane_v,
											Vec3<T> plane_normal ){
	Vec3<T> u1 = segv0 - plane_v;
	Vec3<T> u2 = segv1 - plane_v;
	T c1 = dot(u1, plane_normal);
	T c2 = dot(u2, plane_normal);
	if( (c1>0 && c2>0) || (c1<0 && c2<0) )
		return 0;	//	two end points are on the same side of the plane

	intersection_point = (u2 - u1)*(c1/(c1-c2)) + u1 +  plane_v;
	return 1;
}

/**
	Get intersection point of a line segment and the plane, the plane 
	is represented using three non-coliner vertices on the plane

	\param	intersection_point	return the intersection point
	\param	segv0				end point 0 of line segment
	\param	segv1				end point 1 of line segment
	\param	planev0				vertex 0 on the plane
	\param	planev1				vertex 1 on the plane
	\param	planev2				vertex 2 on the plane
	\return						whether the line segment and the plane intersect
								1: intersect, 0: the line is parallel to the plane
*/
template<typename T>
inline int getSegmentPlaneIntersectionPoint(Vec3<T>& intersection_point,
											Vec3<T> segv0,
											Vec3<T> segv1,
											Vec3<T> planev0,
											Vec3<T> planev1,
											Vec3<T> planev2 ){
	Vec3<T> r1 = planev1 - planev0;
	Vec3<T> r2 = planev2 - planev0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	return getSegmentPlaneIntersectionPoint(intersection_point, segv0, segv1, planev0, nor);
}

/**
	Get whether a line intersect a triangle in 3D space

	\param	intersection_point	the intersection point
	\param	linev0				end vertex 0 of the line
	\param	linev1				end vertex 1 of the line
	\param	triv0				triangle vertex 0
	\param	triv1				triangle vertex 1
	\param	triv2				triangle vertex 2
	\return						whether the line has intersected the triangle, including boundary
								1: intersected; 0: not intersected
*/
template<typename T>
inline int getLineTriangleIntersectionPoint(Vec3<T>& intersection_point,
											Vec3<T> linev0, 
											Vec3<T> linev1,
											Vec3<T> triv0,
											Vec3<T> triv1,
											Vec3<T> triv2 ){
	Vec3<T> r1 = triv1 - triv0;
	Vec3<T> r2 = triv2 - triv0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	int flag = getLinePlaneIntersectionPoint(intersection_point, linev0, linev1, triv0, nor);
	if( flag )
		return isPointProjectionInsideTriangle(intersection_point, triv0, triv1, triv2);
	else
		return 0;
}

/**
	Get whether a ray intersect a triangle in 3D space

	\param	intersection_point	the intersection point
	\param	ray_origin			start point of the ray
	\param	ray_next			next point of the ray
	\param	triv0				triangle vertex 0
	\param	triv1				triangle vertex 1
	\param	triv2				triangle vertex 2
	\return						whether the ray has intersected the triangle, including boundary
								1: intersected; 0: not intersected
*/
template<typename T>
inline int getRayTriangleIntersectionPoint(Vec3<T>& intersection_point,
										   Vec3<T> ray_origin, 
										   Vec3<T> ray_next,
										   Vec3<T> triv0,
										   Vec3<T> triv1,
										   Vec3<T> triv2 ){
	Vec3<T> r1 = triv1 - triv0;
	Vec3<T> r2 = triv2 - triv0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	int flag = getRayPlaneIntersectionPoint(intersection_point, ray_origin, ray_next, triv0, nor);
	if( flag )
		return isPointProjectionInsideTriangle(intersection_point, triv0, triv1, triv2);
	else
		return 0;
}

/**
	Get whether a X ray intersect a triangle in 3D space

	\param	intersection_point	the intersection point
	\param	ray_origin			start point of the ray
	\param	triv0				triangle vertex 0
	\param	triv1				triangle vertex 1
	\param	triv2				triangle vertex 2
	\return						whether the ray has intersected the triangle, including boundary
	1: intersected; 0: not intersected
*/
template<typename T>
inline int getXRayTriangleIntersectionPoint(
	Vec3<T>&	intersection_point,
	Vec3<T>		ray_origin,
	Vec3<T>		triv0,
	Vec3<T>		triv1,
	Vec3<T>		triv2)
{
	//	reduce to 2D
	Vec2<T> p(ray_origin.y, ray_origin.z);
	Vec2<T> v0(triv0.y, triv0.z);
	Vec2<T> v1(triv1.y, triv1.z);
	Vec2<T> v2(triv2.y, triv2.z);

	T coef1, coef2;
	getPointTriangleIntersectionCoef(coef1, coef2, p, v0, v1, v2);

	if (coef1 > 0 && coef2 > 0 && coef1 + coef2 < 1){
		//	the line intersect this triangle, now check whether intersection point is on the right of ray origin
		intersection_point = triv0 *(1 - coef1 - coef2) + triv1*coef1 + triv2*coef2;
		if (intersection_point.x > ray_origin.x)
			return 1;
		else
			return 0;
	}
	else
		return 0;
}

/**
	Get whether a line segment intersect a triangle in 3D space

	\param	intersection_point	the intersection point
	\param	segv0				end vertex 0 of line segment
	\param	segv1				end vertex 1 of line segment
	\param	triv0				triangle vertex 0
	\param	triv1				triangle vertex 1
	\param	triv2				triangle vertex 2
	\return						whether the ray has intersected the triangle, including boundary
								1: intersected; 0: not intersected
*/
template<typename T>
inline int getSegmentTriangleIntersectionPoint(Vec3<T>& intersection_point,
											   Vec3<T> segv0, 
											   Vec3<T> segv1,
											   Vec3<T> triv0,
											   Vec3<T> triv1,
											   Vec3<T> triv2 ){
	Vec3<T> r1 = triv1 - triv0;
	Vec3<T> r2 = triv2 - triv0;
	Vec3<T> nor = cross(r1, r2).Normalize();

	int flag = getSegmentPlaneIntersectionPoint(intersection_point, segv0, segv1, triv0, nor);
	if( flag )
		return isPointProjectionInsideTriangle(intersection_point, triv0, triv1, triv2);
	else
		return 0;
}

///@}

//	========================================
///@{
/**	@name Bonding Box Intersection Test
*/
//	========================================
/**
	Check whether point is inside bonding box 2D

	\param		point		the point
	\param		aabb_min	AABB min
	\param		aabb_max	AABB max
	\return					whether point is inside
*/
template<typename T>
inline int isPointInsideAABB(Vec2<T> point,
							 Vec2<T> aabb_min,
							 Vec2<T> aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1]);

	if( point[0]>=aabb_min[0] && point[0]<=aabb_max[0] && 
		point[1]>=aabb_min[1] && point[1]<=aabb_max[1] )
		return 1;
	else
		return 0;
}

/**
	Get nearest point of a point to AABB

	\param		nearest_point	return the nearest point on AABB
	\param		point			the point
	\param		aabb_min		AABB min
	\param		aabb_max		AABB max
	\return						whether point is inside the AABB
*/
template<typename T>
inline int getNearestPointOnAABB(Vec2<T>&	nearest_point,
								 Vec2<T>	point,
								 Vec2<T>	aabb_min,
								 Vec2<T>	aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1]);

	int outside_flag = 0;
	int space[2] = {0, 0};

	//	scan each dimension
	for(int dim=0; dim<2; dim++){
		if( point[dim] <= aabb_min[dim] ){
			nearest_point[dim]	= aabb_min[dim];
			outside_flag = 1;
		}
		else if( point[dim] >= aabb_max[dim] ){
			nearest_point[dim]	= aabb_max[dim];
			outside_flag = 1;
		}
		else{
			if( point[dim]-aabb_min[dim] < aabb_max[dim]-point[dim] ){
				nearest_point[dim]	= point[dim];
				space[dim] = 1;
			}
			else{
				nearest_point[dim]	= point[dim];
				space[dim] = 2;
			}
		}
	}

	if( outside_flag )		//	if the vertex is outside of AABB, we have calculated nearest point
		return 0;

	//	if the point is inside AABB, move nearest point to surface of AABB
	int moving_dim = 0;
	T	min_dist = (space[0]==1 ? point[0]-aabb_min[0] : aabb_max[0]-point[0]);
	for(int dim=1; dim<2; dim++){
		T	this_dist = (space[dim]==1 ? point[dim]-aabb_min[dim] : aabb_max[dim]-point[dim]);
		if( this_dist < min_dist ){
			min_dist	= this_dist;
			moving_dim	= dim;
		}
	}

	nearest_point[moving_dim] = (space[moving_dim]==1 ? aabb_min[moving_dim] : aabb_max[moving_dim]);

	return 1;
}

/**
	Get farthest point of a point to AABB

	\param		farthest_point	return the farest point on AABB
	\param		point			the point
	\param		aabb_min		AABB min
	\param		aabb_max		AABB max
	\return						whether point is inside the AABB
*/
template<typename T>
inline int getFarthestPointOnAABB(Vec2<T>&	farthest_point,
								  Vec2<T>	point,
								  Vec2<T>	aabb_min,
								  Vec2<T>	aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1]);

	int outside_flag = 0;

	for(int dim=0; dim<2; dim++){
		if( point[dim] <= aabb_min[dim] ){
			farthest_point[dim]	= aabb_max[dim];
			outside_flag = 1;
		}
		else if( point[dim] >= aabb_max[dim] ){
			farthest_point[dim]	= aabb_min[dim];
			outside_flag = 1;
		}
		else{
			if( point[dim]-aabb_min[dim] < aabb_max[dim]-point[dim] )
				farthest_point[dim]	= aabb_max[dim];
			else
				farthest_point[dim]	= aabb_min[dim];
		}
	}

	return !outside_flag;
}

/**
	Check whether point is inside bonding box

	\param		point		the point
	\param		aabb_min	AABB min
	\param		aabb_max	AABB max
	\return					whether point is inside
*/
template<typename T>
inline int isPointInsideAABB(Vec3<T> point,
							 Vec3<T> aabb_min,
							 Vec3<T> aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1] && aabb_min[2]<=aabb_max[2] );

	if( point[0]>=aabb_min[0] && point[0]<=aabb_max[0] && 
		point[1]>=aabb_min[1] && point[1]<=aabb_max[1] &&
		point[2]>=aabb_min[2] && point[2]<=aabb_max[2] )
		return 1;
	else
		return 0;
}

/**
	Get nearest point of a point to AABB

	\param		nearest_point	return the nearest point on AABB
	\param		point			the point
	\param		aabb_min		AABB min
	\param		aabb_max		AABB max
	\return						whether point is inside the AABB
*/
template<typename T>
inline int getNearestPointOnAABB(Vec3<T>&	nearest_point,
								 Vec3<T>	point,
								 Vec3<T>	aabb_min,
								 Vec3<T>	aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1] && aabb_min[2]<=aabb_max[2]);

	int outside_flag = 0;
	int space[3] = {0, 0, 0};

	//	scan each dimension
	for(int dim=0; dim<3; dim++){
		if( point[dim] <= aabb_min[dim] ){
			nearest_point[dim]	= aabb_min[dim];
			outside_flag = 1;
		}
		else if( point[dim] >= aabb_max[dim] ){
			nearest_point[dim]	= aabb_max[dim];
			outside_flag = 1;
		}
		else{
			if( point[dim]-aabb_min[dim] < aabb_max[dim]-point[dim] ){
				nearest_point[dim]	= point[dim];
				space[dim] = 1;
			}
			else{
				nearest_point[dim]	= point[dim];
				space[dim] = 2;
			}
		}
	}

	if( outside_flag )		//	if the vertex is outside of AABB, we have calculated nearest point
		return 0;

	//	if the point is inside AABB, move nearest point to surface of AABB
	int moving_dim = 0;
	T	min_dist = (space[0]==1 ? point[0]-aabb_min[0] : aabb_max[0]-point[0]);
	for(int dim=1; dim<3; dim++){
		T	this_dist = (space[dim]==1 ? point[dim]-aabb_min[dim] : aabb_max[dim]-point[dim]);
		if( this_dist < min_dist ){
			min_dist	= this_dist;
			moving_dim	= dim;
		}
	}

	nearest_point[moving_dim] = (space[moving_dim]==1 ? aabb_min[moving_dim] : aabb_max[moving_dim]);

	return 1;
}

/**
	Get farthest point of a point to AABB

	\param		farthest_point	return the farest point on AABB
	\param		point			the point
	\param		aabb_min		AABB min
	\param		aabb_max		AABB max
	\return						whether point is inside the AABB
*/
template<typename T>
inline int getFarthestPointOnAABB(Vec3<T>&	farthest_point,
								  Vec3<T>	point,
								  Vec3<T>	aabb_min,
								  Vec3<T>	aabb_max ){
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1] && aabb_min[2]<=aabb_max[2]);

	int outside_flag = 0;

	for(int dim=0; dim<3; dim++){
		if( point[dim] <= aabb_min[dim] ){
			farthest_point[dim]	= aabb_max[dim];
			outside_flag = 1;
		}
		else if( point[dim] >= aabb_max[dim] ){
			farthest_point[dim]	= aabb_min[dim];
			outside_flag = 1;
		}
		else{
			if( point[dim]-aabb_min[dim] < aabb_max[dim]-point[dim] )
				farthest_point[dim]	= aabb_max[dim];
			else
				farthest_point[dim]	= aabb_min[dim];
		}
	}

	return !outside_flag;
}

/**
	Check whether line segment intersect bonding box 2D

	\param		seg_v0		end point 0 of line segment
	\param		seg_v1		end point 1 of line segment
	\param		aabb_min	AABB min
	\param		aabb_max	AABB max
	\return					whether line segment intersect bonding box
*/
template<typename T>
inline int isSegmentAABBIntersect(Vec2<T> seg_v0,
								  Vec2<T> seg_v1,
								  Vec2<T> aabb_min,
								  Vec2<T> aabb_max){
	assert(aabb_min[0]<aabb_max[0] && aabb_min[1]<aabb_max[1]);

	return (clipCohenSutherland(seg_v0, seg_v1, aabb_min, aabb_max) ? 1 : 0);	//	if intersect, return 1, else return 0
}

/**
	Check whether line segment intersect bonding box

	\param		seg_v0		end point 0 of line segment
	\param		seg_v1		end point 1 of line segment
	\param		aabb_min	AABB min
	\param		aabb_max	AABB max
	\return					whether line segment intersect bonding box
*/
template<typename T>
inline int isSegmentAABBIntersect(Vec3<T> seg_v0,
								  Vec3<T> seg_v1,
								  Vec3<T> aabb_min,
								  Vec3<T> aabb_max){
	assert(aabb_min[0]<aabb_max[0] && aabb_min[1]<aabb_max[1] && aabb_min[2]<aabb_max[2] );

	return (clipCohenSutherland(seg_v0, seg_v1, aabb_min, aabb_max) ? 1 : 0);	//	if intersect, return 1, else return 0
}

/**
	Check whether ray and AABB intersect 2D

	\param	ray_origin			origin point of the ray
	\param	ray_next			next point of the ray
	\param	aabb_min			AABB min
	\param	aabb_max			AABB max
	\return						whether they intersect
*/
template<typename T>
inline int isRayAABBIntersect(Vec2<T> ray_origin,
							  Vec2<T> ray_next,
							  Vec2<T> aabb_min,
							  Vec2<T> aabb_max){
	Vec2<T> intp;
	return getRayAABBIntersectPoint(intp, ray_origin, ray_next, aabb_min, aabb_max);
}
/**
	Check whether ray and AABB intersect

	\param	ray_origin			origin point of the ray
	\param	ray_next			next point of the ray
	\param	aabb_min			AABB min
	\param	aabb_max			AABB max
	\return						whether they intersect
*/
template<typename T>
inline int isRayAABBIntersect(Vec3<T> ray_origin,
							  Vec3<T> ray_next,
							  Vec3<T> aabb_min,
							  Vec3<T> aabb_max){
	Vec3<T> intp;
	return getRayAABBIntersectPoint(intp, ray_origin, ray_next, aabb_min, aabb_max);
}
/**
	Get Ray AABB Intersection Point 2D

	If ray origin is inside the aabb, return the origin

	Fast Ray-Box Intersection by Andrew Woo \n
	from "Graphics Gems", Academic Press, 1990 \n
	http://tog.acm.org/resources/GraphicsGems/gems/RayBox.c

	\param	intersection_point	intersection point
	\param	ray_origin			origin point of the ray
	\param	ray_next			next point of the ray
	\param	aabb_min			AABB min
	\param	aabb_max			AABB max
	\return						whether they intersect
*/
template<typename T>
inline int getRayAABBIntersectPoint(Vec2<T>& intersection_point,
									Vec2<T> ray_origin,
									Vec2<T> ray_next,
									Vec2<T> aabb_min,
									Vec2<T> aabb_max){
	//	bonding box must be legal, or the algorithm may fail
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1] );
	Vec2<T>	dir = (ray_next - ray_origin).Normalize();

	bool inside = true;
	char quadrant[2];
	register int i;
	int whichPlane;
	T maxT[2];
	T candidatePlane[2];

	/* Find candidate planes; this loop can be avoided if
	rays cast all from the eye(assume perpsective view) */
	for (i=0; i<2; i++){
		if(ray_origin[i] < aabb_min[i]) {
			quadrant[i] = 1;
			candidatePlane[i] = aabb_min[i];
			inside = false;
		}
		else if (ray_origin[i] > aabb_max[i]){
			quadrant[i] = 0;
			candidatePlane[i] = aabb_max[i];
			inside = false;
		}
		else{
			quadrant[i] = 2;
		}
	}

	/* Ray origin inside bounding box */
	if(inside)	{
		intersection_point = ray_origin;
		return 1;
	}

	/* Calculate T distances to candidate planes */
	for (i = 0; i < 2; i++){
		if (quadrant[i] != 2 && fabs(dir[i])>1e-5)
			maxT[i] = (candidatePlane[i]-ray_origin[i]) / dir[i];
		else
			maxT[i] = -1;
	}

	/* Get largest of the maxT's for final choice of intersection */
	whichPlane = 0;
	for (i = 1; i < 2; i++)
		if (maxT[whichPlane] < maxT[i])
			whichPlane = i;

	/* Check final candidate actually inside box */
	if (maxT[whichPlane] < 0) 
		return 0;
	for (i = 0; i < 2; i++){
		if (whichPlane != i) {
			intersection_point[i] = ray_origin[i] + maxT[whichPlane] *dir[i];
			if (intersection_point[i] < aabb_min[i] || intersection_point[i] > aabb_max[i])
				return 0;
		} 
		else {
			intersection_point[i] = candidatePlane[i];
		}
	}
	return 1;				/* ray hits box */
}

/**
	Get Ray AABB Intersection Point

	If ray origin is inside the aabb, return the origin

	Fast Ray-Box Intersection by Andrew Woo \n
	from "Graphics Gems", Academic Press, 1990 \n
	http://tog.acm.org/resources/GraphicsGems/gems/RayBox.c

	\param	intersection_point	intersection point
	\param	ray_origin			origin point of the ray
	\param	ray_next			next point of the ray
	\param	aabb_min			AABB min
	\param	aabb_max			AABB max
	\return						whether they intersect
*/
template<typename T>
inline int getRayAABBIntersectPoint(Vec3<T>& intersection_point,
									Vec3<T> ray_origin,
									Vec3<T> ray_next,
									Vec3<T> aabb_min,
									Vec3<T> aabb_max){
	//	bonding box must be legal, or the algorithm may fail
	assert(aabb_min[0]<=aabb_max[0] && aabb_min[1]<=aabb_max[1] && aabb_min[2]<=aabb_max[2] );
	Vec3<T>	dir = (ray_next - ray_origin).Normalize();

	bool inside = true;
	char quadrant[3];
	register int i;
	int whichPlane;
	T maxT[3];
	T candidatePlane[3];

	/* Find candidate planes; this loop can be avoided if
	rays cast all from the eye(assume perpsective view) */
	for (i=0; i<3; i++){
		if(ray_origin[i] < aabb_min[i]) {
			quadrant[i] = 1;
			candidatePlane[i] = aabb_min[i];
			inside = false;
		}
		else if (ray_origin[i] > aabb_max[i]){
			quadrant[i] = 0;
			candidatePlane[i] = aabb_max[i];
			inside = false;
		}
		else{
			quadrant[i] = 2;
		}
	}

	/* Ray origin inside bounding box */
	if(inside)	{
		intersection_point = ray_origin;
		return 1;
	}

	/* Calculate T distances to candidate planes */
	for (i = 0; i < 3; i++){
		if (quadrant[i] != 2 && fabs(dir[i])>1e-5)
			maxT[i] = (candidatePlane[i]-ray_origin[i]) / dir[i];
		else
			maxT[i] = -1;
	}

	/* Get largest of the maxT's for final choice of intersection */
	whichPlane = 0;
	for (i = 1; i < 3; i++)
		if (maxT[whichPlane] < maxT[i])
			whichPlane = i;

	/* Check final candidate actually inside box */
	if (maxT[whichPlane] < 0) 
		return 0;
	for (i = 0; i < 3; i++){
		if (whichPlane != i) {
			intersection_point[i] = ray_origin[i] + maxT[whichPlane] *dir[i];
			if (intersection_point[i] < aabb_min[i] || intersection_point[i] > aabb_max[i])
				return 0;
		} 
		else {
			intersection_point[i] = candidatePlane[i];
		}
	}
	return 1;				/* ray hits box */
}

/**
	check whether axis aligned ray intersect AABB 2D

	\param	ray_origin	origin of the ray
	\param	ray_dir		direction of the ray, 1 x, 2 y, -1 -x, -2 -y
	\param	aabb_min	aabb min
	\param	aabb_max	aabb max
	\return				whether the ray intersect bonding box
*/
template<typename T>
inline int isAARayAABBIntersect(Vec2<T>	ray_origin,
								int		ray_dir,
								Vec2<T> aabb_min,
								Vec2<T> aabb_max){
	assert( ray_dir>=-2 && ray_dir<=2 && ray_dir!=0 );

	int sgn = ray_dir > 0;
	int idx = sgn ? ray_dir-1 : -1-ray_dir;
	int cmp_idx = !idx;

	if( sgn && ray_origin[idx] > aabb_max[idx] )
		return 0;
	if( !sgn && ray_origin[idx] < aabb_min[idx] )
		return 0;

	if(	ray_origin[cmp_idx]>=aabb_min[cmp_idx] && ray_origin[cmp_idx]<=aabb_max[cmp_idx] )
		return 1;

	return 0;
}

/**
	check whether axis aligned ray intersect AABB 3D

	\param	ray_origin	origin of the ray
	\param	ray_dir		direction of the ray, 1 x, 2 y, 3 z, -1 -x, -2 -y, -3 -z
	\param	aabb_min	aabb min
	\param	aabb_max	aabb max
	\return				whether the ray intersect bonding box
*/
template<typename T>
inline int isAARayAABBIntersect(Vec3<T>	ray_origin,
								int		ray_dir,
								Vec3<T> aabb_min,
								Vec3<T> aabb_max){
	assert( ray_dir>=-3 && ray_dir<=3 && ray_dir!=0 );

	int pos_dir = ray_dir > 0;
	int idx = pos_dir ? ray_dir-1 : -1-ray_dir;
	int cmp_idx1 = (idx + 1) % 3;
	int cmp_idx2 = (idx + 2) % 3;

	if( pos_dir && ray_origin[idx] > aabb_max[idx] )
		return 0;
	if( !pos_dir && ray_origin[idx] < aabb_min[idx] )
		return 0;

	if(	ray_origin[cmp_idx1]>=aabb_min[cmp_idx1] && ray_origin[cmp_idx1]<=aabb_max[cmp_idx1] && 
		ray_origin[cmp_idx2]>=aabb_min[cmp_idx2] && ray_origin[cmp_idx2]<=aabb_max[cmp_idx2] )
		return 1;

	return 0;
}


/**
	check whether +x direction ray intersect AABB 2D

	\param	xray_origin	origin of the +x direction ray
	\param	aabb_min	aabb min
	\param	aabb_max	aabb max
	\return				whether the ray intersect bonding box
*/
template<typename T>
inline int isXRayAABBIntersect(Vec2<T> xray_origin,
							   Vec2<T> aabb_min,
							   Vec2<T> aabb_max){
	if( xray_origin[0] > aabb_max[0] )
		return 0;

	if(	xray_origin[1]>=aabb_min[1] && xray_origin[1]<=aabb_max[1] )
		return 1;
	else
		return 0;
}

/**
	check whether +x direction ray intersect AABB

	\param	xray_origin	origin of the +x direction ray
	\param	aabb_min	aabb min
	\param	aabb_max	aabb max
	\return				whether the ray intersect bonding box
*/
template<typename T>
inline int isXRayAABBIntersect(Vec3<T> xray_origin,
							   Vec3<T> aabb_min,
							   Vec3<T> aabb_max){
	if( xray_origin[0] > aabb_max[0] )
		return 0;

	if(	xray_origin[1]>=aabb_min[1] && xray_origin[1]<=aabb_max[1] && 
		xray_origin[2]>=aabb_min[2] && xray_origin[2]<=aabb_max[2] )
		return 1;
	else
		return 0;
}

/**
	check whether -x direction ray intersect AABB 2D

	\param	nxray_origin	origin of the -x direction ray
	\param	aabb_min		aabb min
	\param	aabb_max		aabb max
	\return					whether the ray intersect bonding box
*/
template<typename T>
inline int isNXRayAABBIntersect(Vec2<T> nxray_origin,
								Vec2<T> aabb_min,
								Vec2<T> aabb_max){
	if( nxray_origin[0] < aabb_min[0] )
		return 0;

	if(	nxray_origin[1]>=aabb_min[1] && nxray_origin[1]<=aabb_max[1] )
		return 1;
	else
		return 0;
}

/**
	check whether -x direction ray intersect AABB

	\param	nxray_origin	origin of the -x direction ray
	\param	aabb_min		aabb min
	\param	aabb_max		aabb max
	\return					whether the ray intersect bonding box
*/
template<typename T>
inline int isNXRayAABBIntersect(Vec3<T> nxray_origin,
								Vec3<T> aabb_min,
								Vec3<T> aabb_max){
	if( nxray_origin[0] < aabb_min[0] )
		return 0;

	if(	nxray_origin[1]>=aabb_min[1] && nxray_origin[1]<=aabb_max[1] && 
		nxray_origin[2]>=aabb_min[2] && nxray_origin[2]<=aabb_max[2] )
		return 1;
	else
		return 0;
}

/**
	Check whether two axis aligned bonding boxes intersect

	\param		aabb1_min	AABB 1 min
	\param		aabb1_max	AABB 1 max
	\param		aabb2_min	AABB 2 min
	\param		aabb2_max	AABB 2 max
	\return					whether the two bonding boxes intersect
*/
template<typename T>
inline int isAABBsIntersect(Vec2<T> aabb1_min,
							Vec2<T> aabb1_max,
							Vec2<T> aabb2_min,
							Vec2<T> aabb2_max){
	//	bonding box must be legal, or the algorithm may fail
	assert(aabb1_min[0]<=aabb1_max[0] && aabb1_min[1]<=aabb1_max[1]);
	assert(aabb2_min[0]<=aabb2_max[0] && aabb2_min[1]<=aabb2_max[1]);

	Vec2<T> t = ((aabb1_max - aabb1_min) + (aabb2_max - aabb2_min)) * 0.5;
	Vec2<T> r = ((aabb1_min + aabb1_max) - (aabb2_min + aabb2_max)) * 0.5;
	if( r[0] < 0 )	r[0] = - r[0];
	if( r[1] < 0 )	r[1] = - r[1];
	if( r[0] <= t[0] && r[1] <= t[1] )
		return 1;
	else
		return 0;
}

/**
	Check whether two axis aligned bonding boxes intersect

	\param		aabb1_min	AABB 1 min
	\param		aabb1_max	AABB 1 max
	\param		aabb2_min	AABB 2 min
	\param		aabb2_max	AABB 2 max
	\return					whether the two bonding boxes intersect
*/
template<typename T>
inline int isAABBsIntersect(Vec3<T> aabb1_min,
							Vec3<T> aabb1_max,
							Vec3<T> aabb2_min,
							Vec3<T> aabb2_max){
	//	bonding box must be legal, or the algorithm may fail
	assert(aabb1_min[0]<=aabb1_max[0] && aabb1_min[1]<=aabb1_max[1] && aabb1_min[2]<=aabb1_max[2] );
	assert(aabb2_min[0]<=aabb2_max[0] && aabb2_min[1]<=aabb2_max[1] && aabb2_min[2]<=aabb2_max[2] );

	Vec3<T> t = ((aabb1_max - aabb1_min) + (aabb2_max - aabb2_min)) * 0.5;
	Vec3<T> r = ((aabb1_min + aabb1_max) - (aabb2_min + aabb2_max)) * 0.5;
	if( r[0] < 0 )	r[0] = - r[0];
	if( r[1] < 0 )	r[1] = - r[1];
	if( r[2] < 0 )	r[2] = - r[2];
	if( r[0] <= t[0] && r[1] <= t[1] && t[2] <= t[2] )
		return 1;
	else
		return 0;
}

///@}



}}	//	namespace yz::geometry

#endif	//	__YZ_INTERSECTION_TEST_H__