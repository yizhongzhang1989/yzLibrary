/***********************************************************/
/**	\file
	\brief		regression of 3d data
	\author		Yizhong Zhang
	\date		8/27/2013
*/
/***********************************************************/
#ifndef __YZ_REGRESSION_3D_H__
#define __YZ_REGRESSION_3D_H__

#include <vector>
#include "yzLib/yz_math/yz_matrix.h"
#include "yzlib/yz_math/yz_vector_utils.h"

namespace yz{	namespace statistics{

//	========================================
///@{
/**	@name	Fitting Plane in 3D
*/
//	========================================

/**
	fit a plane with 3d points, the plane is z = Ax + By + C

	the fitting is better if the plane is perpendicular to z axis

	\param	A			return coef A
	\param	B			return coef B
	\param	C			return coef C
	\param	xyz_ptr		pointer to 3d point cloud
	\param	point_num	how many points
	\return				whether fit succeed
*/
template<typename T>
inline int fitPlaneLeastSquare(T& A, T& B, T& C, const T* xyz_ptr, int point_num){
	if( point_num < 3 )
		return 0;

	Matrix3x3<PROMOTE_T_TO_FLOAT>	mat(0, 0, 0, 0, 0, 0, 0, 0, 0);
	Vec3<PROMOTE_T_TO_FLOAT>		rhs(0, 0, 0);

	//	fill the matrix
	for(int i=0; i<point_num; i++){
		T x = xyz_ptr[i*3];
		T y = xyz_ptr[i*3+1];
		T z = xyz_ptr[i*3+2];
		mat.data[0][0] += x * x;
		mat.data[0][1] += x * y;
		mat.data[0][2] += x;
		mat.data[1][1] += y * y;
		mat.data[1][2] += y;
		rhs.x += x * z;
		rhs.y += y * z;
		rhs.z += z;
	}
	mat.data[1][0] = mat.data[0][1];
	mat.data[2][0] = mat.data[0][2];
	mat.data[2][1] = mat.data[1][2];
	mat.data[2][2] = point_num;

	Vec3<PROMOTE_T_TO_FLOAT>	res = mat.Inverse() * rhs;
	A = res[0];
	B = res[1];
	C = res[2];
	return 1;
}

/**
	fit a plane with 3d points, the plane is z = Ax + By + C

	the fitting is better if the plane is perpendicular to z axis

	\param	A		return coef A
	\param	B		return coef B
	\param	C		return coef C
	\param	xyz		list of the point
	\return			whether fit succeed
*/
template<typename T>
inline int fitPlaneLeastSquare(T& A, T& B, T& C, const std::vector<Vec3<T>>& xyz){
	if( xyz.empty() )
		return 0;

	return fitPlaneLeastSquare( A, B, C, (T*)&xyz[0], xyz.size() );
}

/**
	fit a plane with 3d points which minimize the distance in normal direction

	\param	center_x	return center coordinate of the plane
	\param	center_y	return center coordinate of the plane
	\param	center_z	return center coordinate of the plane
	\param	normal_x	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	normal_y	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	normal_z	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	xyz_ptr		pointer to 3d point cloud
	\param	point_num	how many points
	\return				whether fit succeed
*/
template<typename T>
inline int fitPlaneLeastSquareNormalDirection(T& center_x, T& center_y, T& center_z, 
											  T& normal_x, T& normal_y, T& normal_z, 
											  const T* xyz_ptr, int point_num){
	if( point_num < 3 )
		return 0;

	//	calculate rotation matrix
	Matrix3x3<PROMOTE_T_TO_FLOAT>	rot;
	Vec3<PROMOTE_T_TO_FLOAT>		normal(normal_x, normal_y, normal_z);
	Vec3<PROMOTE_T_TO_FLOAT>		old_normal(normal);
	Vec3<PROMOTE_T_TO_FLOAT>		z(0, 0, 1);
	if( normal.Length() < 1e-4 )		//	if invalid guess normal is given, we just treate it to be +z
		normal = yz::Vec3d(0, 0, 1);
	else if( normal.z < 0 )
		normal *= -1;
	if( normal.z > 0.995 )	//	very close to z direction, don't need to rotate
		rot.SetIdentity();
	else					//	need to rotate in order to improve fitting result
		rot.SetRotationRad( yz::cross(normal, z), yz::angleRadBetweenVectors(normal, z) );

	Matrix3x3<PROMOTE_T_TO_FLOAT>	mat(0, 0, 0, 0, 0, 0, 0, 0, 0);
	Vec3<PROMOTE_T_TO_FLOAT>		rhs(0, 0, 0);
	Vec3<T>							center(0, 0, 0);

	//	fill the matrix
	for(int i=0; i<point_num; i++){
		Vec3<T> p(xyz_ptr+i*3);
		p = rot * p;
		mat.data[0][0] += p.x * p.x;
		mat.data[0][1] += p.x * p.y;
		mat.data[0][2] += p.x;
		mat.data[1][1] += p.y * p.y;
		mat.data[1][2] += p.y;
		rhs.x += p.x * p.z;
		rhs.y += p.y * p.z;
		rhs.z += p.z;

		center.x += p.x;
		center.y += p.y;
	}
	mat.data[1][0] = mat.data[0][1];
	mat.data[2][0] = mat.data[0][2];
	mat.data[2][1] = mat.data[1][2];
	mat.data[2][2] = point_num;
	center.z = rhs.z;

	center /= point_num;

	//	calculate plane function
	Vec3<PROMOTE_T_TO_FLOAT>	res = mat.Inverse() * rhs;
	normal.x = -res[0];
	normal.y = -res[1];
	normal.z = 1;
	normal.SetNormalize();

	//	rotate back
	rot.SetInverse();
	center = rot * center;
	normal = rot * normal;
	if( yz::dot(old_normal, normal) < 0 )
		normal *= -1;

	center_x = center.x;
	center_y = center.y;
	center_z = center.z;
	normal_x = normal.x;
	normal_y = normal.y;
	normal_z = normal.z;

	return 1;
}

/**
	fit a plane with 3d points which minimize the distance in normal direction

	\param	center	return center coordinate of the plane
	\param	normal	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	xyz		list of the point
	\return			whether fit succeed
*/
template<typename T>
inline int fitPlaneLeastSquareNormalDirection(Vec3<T>& center, Vec3<T>& normal, const std::vector<Vec3<T>>& xyz){
	if( xyz.empty() )
		return 0;

	return fitPlaneLeastSquareNormalDirection(center.x, center.y, center.z,
		normal.x, normal.y, normal.z, (T*)&xyz[0], xyz.size() );
}
/**
	fit a plane with 3d points which minimize the distance in normal direction

	\param	center_x	return center coordinate of the plane
	\param	center_y	return center coordinate of the plane
	\param	center_z	return center coordinate of the plane
	\param	normal_x	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	normal_y	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	normal_z	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	xyz_ptr		pointer to 3d point cloud
	\param	w_ptr		pointer to weight of each point
	\param	point_num	how many points
	\return				whether fit succeed
*/
template<typename T>
inline int fitPlaneWeightedLeastSquareNormalDirection(T& center_x, T& center_y, T& center_z, 
													  T& normal_x, T& normal_y, T& normal_z, 
													  const T* xyz_ptr, const T* w_ptr, int point_num){
	if( point_num < 3 )
		return 0;

	//	calculate rotation matrix
	Matrix3x3<PROMOTE_T_TO_FLOAT>	rot;
	Vec3<PROMOTE_T_TO_FLOAT>		normal(normal_x, normal_y, normal_z);
	Vec3<PROMOTE_T_TO_FLOAT>		old_normal(normal);
	Vec3<PROMOTE_T_TO_FLOAT>		z(0, 0, 1);
	if( normal.Length() < 1e-4 )		//	if invalid guess normal is given, we just treate it to be +z
		normal = yz::Vec3d(0, 0, 1);
	else if( normal.z < 0 )
		normal *= -1;
	if( normal.z > 0.995 )	//	very close to z direction, don't need to rotate
		rot.SetIdentity();
	else					//	need to rotate in order to improve fitting result
		rot.SetRotationRad( yz::cross(normal, z), yz::angleRadBetweenVectors(normal, z) );

	Matrix3x3<PROMOTE_T_TO_FLOAT>	mat(0, 0, 0, 0, 0, 0, 0, 0, 0);
	Vec3<PROMOTE_T_TO_FLOAT>		rhs(0, 0, 0);
	Vec3<T>							center(0, 0, 0);

	//	fill the matrix
	for(int i=0; i<point_num; i++){
		Vec3<T> p(xyz_ptr+i*3);
		T		w = w_ptr[i];
		p = rot * p;
		mat.data[0][0] += p.x * p.x * w;
		mat.data[0][1] += p.x * p.y * w;
		mat.data[0][2] += p.x * w;
		mat.data[1][1] += p.y * p.y * w;
		mat.data[1][2] += p.y * w;
		mat.data[2][2] += w;
		rhs.x += p.x * p.z * w;
		rhs.y += p.y * p.z * w;
		rhs.z += p.z * w;

		center.x += p.x * w;
		center.y += p.y * w;
	}
	mat.data[1][0] = mat.data[0][1];
	mat.data[2][0] = mat.data[0][2];
	mat.data[2][1] = mat.data[1][2];
	center.z = rhs.z;

	center /= mat.data[2][2];

	//	calculate plane function
	Vec3<PROMOTE_T_TO_FLOAT>	res = mat.Inverse() * rhs;
	normal.x = -res[0];
	normal.y = -res[1];
	normal.z = 1;
	normal.SetNormalize();

	//	rotate back
	rot.SetInverse();
	center = rot * center;
	normal = rot * normal;
	if( yz::dot(old_normal, normal) < 0 )
		normal *= -1;

	center_x = center.x;
	center_y = center.y;
	center_z = center.z;
	normal_x = normal.x;
	normal_y = normal.y;
	normal_z = normal.z;

	return 1;
}

/**
	fit a plane with 3d points which minimize the distance in normal direction

	\param	center	return center coordinate of the plane
	\param	normal	i/o, given the approximate of the plane, return the calculated normal of the plane
	\param	xyz		list of the point
	\param	w		list of weight of each point
	\return			whether fit succeed
*/
template<typename T>
inline int fitPlaneWeightedLeastSquareNormalDirection(Vec3<T>&						center, 
													  Vec3<T>&						normal, 
													  const std::vector<Vec3<T>>&	xyz,
													  const std::vector<T>&			w){
	if( xyz.empty() || xyz.size() != w.size() )
		return 0;

	return fitPlaneWeightedLeastSquareNormalDirection(center.x, center.y, center.z,
		normal.x, normal.y, normal.z, (T*)&xyz[0], (T*)&w[0], xyz.size() );
}


///@}


}}	//	namespace yz::statistics

#endif	//	__YZ_STATISTICS_H__