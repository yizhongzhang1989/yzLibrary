/***********************************************************/
/**	\file
	\brief		polygon util functions
	\author		Yizhong Zhang
	\date		4/22/2014
*/
/***********************************************************/
#ifndef __YZ_POLYGON_UTILS_H__
#define __YZ_POLYGON_UTILS_H__

#include <iostream>
#include <vector>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_intersection_test.h"

namespace yz{  namespace geometry{  namespace polygon{

//	========================================
///@{
/**	@name Polygon Area
*/
//	========================================
/**
	calculate the area of the polygon with sign

	\param	vertex			vertex list of polygon
	\param	edge			edge list of polygon
	\return					the area of the polygon, positive for counter-clockwise, negitive otherwise
*/
template<typename T>
T calculatePolygonAreaSigned(
	const std::vector<Vec2<T>>&	vertex, 
	const std::vector<int2>&	edge )
{
	T	area = 0;
	for(int i=0; i<edge.size(); i++){
		Vec2<T> v1 = vertex[edge[i].x];
		Vec2<T> v2 = vertex[edge[i].y];
		area += v1.x * v2.y - v1.y * v2.x;
	}
	return area * 0.5; 
}

///@}

}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_POLYGON_UTILS_H__