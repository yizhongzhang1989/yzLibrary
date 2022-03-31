/***********************************************************/
/**	\file
	\brief		Bonding Box
	\author		Yizhong Zhang
	\date		6/17/2012
*/
/***********************************************************/
#ifndef __YZ_BONDING_BOX_H__
#define __YZ_BONDING_BOX_H__

#include <vector>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Calculate AABB
*/
//	========================================

/**
	calculate bonding box from a bundle of 2D vertices

	\param	bb_min			return min coordinate of bonding box
	\param	bb_max			return max coordinate of bonding box
	\param	vertex_2d		vertex list, arranged in xy_xy_
	\param	vertex_number	number of vertices
	\return					whether the bonding box coef is valid, 0: not valid, 1: valid
*/
template<typename T>
inline int getAABBCoef(Vec2<T>& bb_min, Vec2<T>& bb_max, const T* vertex_2d, int vertex_number){
	if( vertex_number <= 0 )
		return 0;

	bb_min.x = vertex_2d[0];
	bb_min.y = vertex_2d[1];
	bb_max = bb_min;

	for( int i=1; i<vertex_number; i++ ){
		if( vertex_2d[i*2  ]<bb_min.x )	bb_min.x = vertex_2d[i*2  ];
		if( vertex_2d[i*2  ]>bb_max.x )	bb_max.x = vertex_2d[i*2  ];
		if( vertex_2d[i*2+1]<bb_min.y )	bb_min.y = vertex_2d[i*2+1];
		if( vertex_2d[i*2+1]>bb_max.y )	bb_max.y = vertex_2d[i*2+1];
	}

	return 1;
}

/**
	calculate bonding box from a bundle of 2D vertices stored in vector

	\param	bb_min			return min coordinate of bonding box
	\param	bb_max			return max coordinate of bonding box
	\param	vertex_2d		vector to hold vertices
	\return					whether the bonding box coef is valid, 0: not valid, 1: valid
*/
template<typename T>
inline int getAABBCoef(Vec2<T>& bb_min, Vec2<T>& bb_max, const std::vector<Vec2<T>>& vertex_2d){
	if( vertex_2d.empty() )
		return 0;

	return getAABBCoef( bb_min, bb_max, (T*)&vertex_2d[0], vertex_2d.size() );
}

/**
	calculate bonding box from a bundle of 3D vertices

	\param	bb_min			return min coordinate of bonding box
	\param	bb_max			return max coordinate of bonding box
	\param	vertex_3d		vertex list, arranged in xyz_xyz_
	\param	vertex_number	number of vertices
	\return					whether the bonding box coef is valid, 0: not valid, 1: valid
*/
template<typename T>
inline int getAABBCoef(Vec3<T>& bb_min, Vec3<T>& bb_max, const T* vertex_3d, int vertex_number){
	if( vertex_number <= 0 )
		return 0;

	bb_min.x = vertex_3d[0];
	bb_min.y = vertex_3d[1];
	bb_min.z = vertex_3d[2];
	bb_max = bb_min;

	for( int i=1; i<vertex_number; i++ ){
		if( vertex_3d[i*3  ]<bb_min.x )	bb_min.x = vertex_3d[i*3  ];
		if( vertex_3d[i*3  ]>bb_max.x )	bb_max.x = vertex_3d[i*3  ];
		if( vertex_3d[i*3+1]<bb_min.y )	bb_min.y = vertex_3d[i*3+1];
		if( vertex_3d[i*3+1]>bb_max.y )	bb_max.y = vertex_3d[i*3+1];
		if( vertex_3d[i*3+2]<bb_min.z )	bb_min.z = vertex_3d[i*3+2];
		if( vertex_3d[i*3+2]>bb_max.z )	bb_max.z = vertex_3d[i*3+2];
	}

	return 1;
}

/**
	calculate bonding box from a bundle of 3D vertices stored in vector

	\param	bb_min			return min coordinate of bonding box
	\param	bb_max			return max coordinate of bonding box
	\param	vertex_3d		vector to hold vertices
	\return					whether the bonding box coef is valid, 0: not valid, 1: valid
*/
template<typename T>
inline int getAABBCoef(Vec3<T>& bb_min, Vec3<T>& bb_max, const std::vector<Vec3<T>>& vertex_3d){
	if( vertex_3d.empty() )
		return 0;

	return getAABBCoef( bb_min, bb_max, (T*)&vertex_3d[0], vertex_3d.size() );
}

///@}


//	========================================
///@{
/**	@name AABB Parameters
*/
//	========================================

/**
	get the 12 edges of an AABB
*/
template<typename T>
void getAABBEdges(
	yz::Vec3<T>					edge_v0[12],
	yz::Vec3<T>					edge_v1[12],
	yz::Vec3<T>					bb_min,
	yz::Vec3<T>					bb_max)
{
	edge_v0[0] = yz::Vec3f(bb_min.x, bb_min.y, bb_min.z);
	edge_v0[1] = yz::Vec3f(bb_min.x, bb_max.y, bb_min.z);
	edge_v0[2] = yz::Vec3f(bb_min.x, bb_min.y, bb_max.z);
	edge_v0[3] = yz::Vec3f(bb_min.x, bb_max.y, bb_max.z);
	edge_v0[4] = yz::Vec3f(bb_min.x, bb_min.y, bb_min.z);
	edge_v0[5] = yz::Vec3f(bb_max.x, bb_min.y, bb_min.z);
	edge_v0[6] = yz::Vec3f(bb_min.x, bb_min.y, bb_max.z);
	edge_v0[7] = yz::Vec3f(bb_max.x, bb_min.y, bb_max.z);
	edge_v0[8] = yz::Vec3f(bb_min.x, bb_min.y, bb_min.z);
	edge_v0[9] = yz::Vec3f(bb_max.x, bb_min.y, bb_min.z);
	edge_v0[10] = yz::Vec3f(bb_min.x, bb_max.y, bb_min.z);
	edge_v0[11] = yz::Vec3f(bb_max.x, bb_max.y, bb_min.z);

	edge_v1[0] = yz::Vec3f(bb_max.x, bb_min.y, bb_min.z);
	edge_v1[1] = yz::Vec3f(bb_max.x, bb_max.y, bb_min.z);
	edge_v1[2] = yz::Vec3f(bb_max.x, bb_min.y, bb_max.z);
	edge_v1[3] = yz::Vec3f(bb_max.x, bb_max.y, bb_max.z);
	edge_v1[4] = yz::Vec3f(bb_min.x, bb_max.y, bb_min.z);
	edge_v1[5] = yz::Vec3f(bb_max.x, bb_max.y, bb_min.z);
	edge_v1[6] = yz::Vec3f(bb_min.x, bb_max.y, bb_max.z);
	edge_v1[7] = yz::Vec3f(bb_max.x, bb_max.y, bb_max.z);
	edge_v1[8] = yz::Vec3f(bb_min.x, bb_min.y, bb_max.z);
	edge_v1[9] = yz::Vec3f(bb_max.x, bb_min.y, bb_max.z);
	edge_v1[10] = yz::Vec3f(bb_min.x, bb_max.y, bb_max.z);
	edge_v1[11] = yz::Vec3f(bb_max.x, bb_max.y, bb_max.z);
};

///@}

/**
Axis aligned bonding box

Template is required to specify vector representation, such as Vec3<double>, Vec2<float>
*/
template<class VEC_NT>
class AABB {
public:
	VEC_NT	bb_min, bb_max;

public:
	inline void GetAABBCoef(const typename VEC_NT::Type* vertex_nd, int vertex_number) {
		getAABBCoef(bb_min, bb_max, vertex_nd, vertex_number);
	}

	inline void GetAABBCoef(const std::vector<VEC_NT>& vertex_nd) {
		getAABBCoef(bb_min, bb_max, vertex_nd);
	}

	inline void Expand(AABB<VEC_NT> aabb) {
		for (int i = 0; i < VEC_NT::dims; i++) {
			if (aabb.bb_min[i] < bb_min[i])		bb_min[i] = aabb.bb_min[i];
			if (aabb.bb_max[i] > bb_max[i])		bb_max[i] = aabb.bb_max[i];
		}
	}

	inline void Expand(VEC_NT p) {
		for (int i = 0; i < VEC_NT::dims; i++) {
			if (p[i] < bb_min[i])	bb_min[i] = p[i];
			if (p[i] > bb_max[i])	bb_max[i] = p[i];
		}
	}
};


/**
2D bonding box
*/
template<class T>
class AABB2D : public AABB<Vec2<T>> {
};


/**
2D bonding box
*/
template<class T>
class AABB3D : public AABB<Vec3<T>> {
};



}}	//	namespace yz::geometry

#endif	//	__YZ_BONDING_BOX_H__