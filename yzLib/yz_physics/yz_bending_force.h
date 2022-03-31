/***********************************************************/
/**	\file
	\brief		Bending Force
	\details	
	\author		Yizhong Zhang
	\date		10/18/2012
*/
/***********************************************************/
#ifndef __YZ_BENDING_FORCE_H__
#define __YZ_BENDING_FORCE_H__

#include <iostream>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace physics{

/**
	calculate the bending force of the four vertices

	v1, v2 are wing vertices; v3, v4 are common vertices. 
	refer to Figure 1. in paper : 
	Simulation of Clothing with Folds and Wrinkles

		  v4			\n
		 /|\			\n
	v1	/ | \ v2		\n
		\ | /			\n
		 \|/			\n
		  v3

	\param	force_v1	return force on vertex 1
	\param	force_v2	return force on vertex 2
	\param	force_v3	return force on vertex 3
	\param	force_v4	return force on vertex 4
	\param	v1			vertex 1
	\param	v2			vertex 2
	\param	v3			vertex 3, on common edge of the two faces
	\param	v4			vertex 4, on common edge of the two faces
	\param	k_bending	bending stiffness
*/
template<typename T>
inline void calculateBendingForces(Vec3<T>&	force_v1, 
								   Vec3<T>&	force_v2, 
								   Vec3<T>&	force_v3, 
								   Vec3<T>&	force_v4,
								   Vec3<T>	v1,
								   Vec3<T>	v2,
								   Vec3<T>	v3,
								   Vec3<T>	v4,
								   T		k_bending = 1){
	Vec3<T>	n1	= cross(v1-v3, v1-v4);	//	weighted normal, don't normalize
	Vec3<T> n2	= cross(v2-v4, v2-v3);	//	weighted normal, don't normalize
	Vec3<T> e	= v4 - v3;

	Vec3<T> u1	= e.Length() / n1.SquareLength() * n1;
	Vec3<T> u2	= e.Length() / n2.SquareLength() * n2;
	Vec3<T> u3	= dot(v1-v4, e) / (e.Length() * n1.SquareLength()) * n1 
		+ dot(v2-v4, e) / (e.Length() * n2.SquareLength()) * n2;
	Vec3<T> u4	= -dot(v1-v3, e) / (e.Length() * n1.SquareLength()) * n1 
		- dot(v2-v3, e) / (e.Length() * n2.SquareLength()) * n2;

	T	coef = k_bending * e.SquareLength() / (n1.Length() + n2.Length());
	n1.SetNormalize();
	n2.SetNormalize();
	T	sin_theta		= dot( cross(n1, n2), e );
	T	value = (1 - dot(n1, n2)) * 0.5;
	if( value < 0 ){
		value = 0;
	}
	T	sin_half_theta	= sqrt( value );
	if( sin_theta<0 )	sin_half_theta = -sin_half_theta;

	force_v1 = coef * sin_half_theta * u1;
	force_v2 = coef * sin_half_theta * u2;
	force_v3 = coef * sin_half_theta * u3;
	force_v4 = coef * sin_half_theta * u4;

	if( force_v1.x != force_v1.x ){
		std::cout << "error" << std::endl;
	}
}



}}	//	end namespace yz::physics

#endif	//	__YZ_BENDING_FORCE_H__