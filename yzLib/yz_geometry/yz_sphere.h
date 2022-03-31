/***********************************************************/
/**	\file
	\brief		Sphere
	\author		Yizhong Zhang
	\date		10/23/2012
*/
/***********************************************************/
#ifndef __YZ_SPHERE_H__
#define __YZ_SPHERE_H__

#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

/**
	A sphere in 3D space
*/
template<class T>
class Sphere{
public:
	Vec3<T>	center;
	T		radius;

public:
	/**
		default constructor, at original point, radius 1
	*/
	Sphere() : radius(1){};

	template<typename T1, typename T2> Sphere(Vec3<T1> sphere_center, T2 sphere_radius){
		center	= sphere_center;
		radius	= sphere_radius;
	}
};


}}	//	namespace yz::geometry

#endif	//	__YZ_SPHERE_H__