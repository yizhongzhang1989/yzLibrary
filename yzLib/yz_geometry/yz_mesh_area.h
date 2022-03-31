/***********************************************************/
/**	\file
	\brief		calculate areas of mesh
	\author		Yizhong Zhang
	\date		12/31/2015
*/
/***********************************************************/
#ifndef __YZ_MESH_AREA_H__
#define __YZ_MESH_AREA_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_vector_utils.h"


namespace yz{	namespace geometry{

/**
	Calculate mixed area of each vertex

	This function is copied from PtrCurvTriMesh::CalculateEdgeWeightMixedArea(),
	if any bug detected in this function, check that function as well.

	\param	mixed_area		return the mixed area of each vertex
	\param	vertex			input vertex of the mesh
	\param	face			input face of the mesh
	\return					sum of all mixed area values
*/
template<typename T>
inline T calculateMixedArea(
	std::vector<T>&				mixed_area,
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face )
{
	mixed_area.clear();
	mixed_area.resize(vertex.size(), 0);

	T area_sum = 0;

	//	for each face, calculate new data
	for (int f = 0; f<face.size(); f++){
		int		v_idx[3] = { face[f].x, face[f].y, face[f].z };
		Vec3<T> v[3] = { vertex[v_idx[0]], vertex[v_idx[1]], vertex[v_idx[2]] };

		//	edge weight
		T cot_val[3];
		for (int i = 0; i<3; i++){
			cot_val[i] = cotBetweenVectors(v[(i + 1) % 3] - v[i], v[(i + 2) % 3] - v[i]);
		}

		//	mixed area
		if (cot_val[0] >= 0 && cot_val[1] >= 0 && cot_val[2] >= 0){	//	non-obtuse triangle
			for (int i = 0; i < 3; i++){
				T add1 = cot_val[(i + 1) % 3] * (v[(i + 2) % 3] - v[i]).SquareLength() / 8;
				T add2 = cot_val[(i + 2) % 3] * (v[(i + 1) % 3] - v[i]).SquareLength() / 8;
				mixed_area[v_idx[i]] += add1;
				mixed_area[v_idx[i]] += add2;
				area_sum += (add1 + add2);
			}
		}
		else{													//	obtuse triangle
			T area = cross(v[1] - v[0], v[2] - v[0]).Length() / 2;
			for (int i = 0; i < 3; i++){
				T add = area / (cot_val[i] < 0 ? 2 : 4);	//	obtuse angle take half the area
				mixed_area[v_idx[i]] += add;
				area_sum += add;
			}
		}
	}

	return area_sum;
}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_AREA_H__