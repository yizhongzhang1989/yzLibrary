/***********************************************************/
/**	\file
	\brief		distance field of polygon
	\author		Yizhong Zhang
	\date		1/17/2013
*/
/***********************************************************/
#ifndef __YZ_POLYGON_DISTANCE_FIELD_H__
#define __YZ_POLYGON_DISTANCE_FIELD_H__

#include <iostream>
#include <vector>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_intersection_test.h"

namespace yz{  namespace geometry{  namespace polygon{

//	========================================
///@{
/**	@name Nearest Point On Polygon
*/
//	========================================
/**
	get the nearest point of a given point on the polygon

	\param	nearest_point	return the nearest point on polygon
	\param	point			the testing point
	\param	vertex			vertex list of polygon
	\param	edge			edge list of polygon
	\return					the edge id containning the nearest point
*/
template<typename T>
int getNearestPointOnPolygon(
	Vec2<T>&						nearest_point,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge
) {
	int nearest_edge_id;
	return getNearestPointOnPolygon(nearest_point, nearest_edge_id, point, vertex, edge);
}

/**
	get the nearest point and edge of a given point on the polygon

	\param	nearest_point	return the nearest point on mesh
	\param	nearest_edge_id	return id of the edge that the nearest point lies in
	\param	point			the testing point
	\param	vertex			vertex list of polygon
	\param	edge			edge list of polygon
	\return					the edge id containning the nearest point
*/
template<typename T>
int getNearestPointOnPolygon(
	Vec2<T>&						nearest_point,
	int&							nearest_edge_id,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge
) {
	//	create face aabb tree
	int i;
	AABBTree2D<T>	edge_aabb;
	edge_aabb.BuildEdgeAABBTree(vertex, edge);

	//	calculate nearest point using AABB Tree
	return getNearestPointOnPolygon(nearest_point, nearest_edge_id,
		point, vertex, edge, edge_aabb);
}

/**
recursively get the nearest point on polygon using AABB tree

\param	nearest_point		return the nearest point position
\param	nearest_edge_id		the id of the edge that contain the nearest point
\param	nearest_squ_dist	square distance of point to nearest point
\param	point				the testing point
\param	vertex				vertex list of polygon
\param	edge				edge list of polygon
\param	edge_aabb			edge aabb tree of the polygon
\param	node_id				current checking node id
*/
template<typename T>
void recursiveGetNearestPointOnPolygon(
	Vec2<T>&						nearest_point,
	int&							nearest_edge_id,
	T&								nearest_squ_dist,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb,
	int								node_id
) {
	//	leaf has reached
	if (node_id < edge_aabb.leaf_number) {
		Vec2<T>	this_nearest_point;
		getNearestPointOnSegment(this_nearest_point, point,
			vertex[edge[node_id].x], vertex[edge[node_id].y]);

		T this_squ_dist = (this_nearest_point - point).SquareLength();
		if (this_squ_dist < nearest_squ_dist) {
			nearest_squ_dist = this_squ_dist;
			nearest_point = this_nearest_point;
			nearest_edge_id = node_id;
		}

		return;
	}

	//	check two children, if point is inside the bonding box or 
	//	closer than current nearest distance, then go deeper
	//	in order to accelerate this process, we first pick out the node 
	//	that is closer to the point, in this way we have a chance that 
	//	we don't need to parse the second node any more
	int child[2] = { edge_aabb.node[node_id].left, edge_aabb.node[node_id].right };
	int inside_flag[2];
	Vec2<T>	this_farthest_point[2];
	for (int i = 0; i<2; i++) {
		inside_flag[i] = getFarthestPointOnAABB(this_farthest_point[i], point,
			edge_aabb.node[child[i]].bb_min, edge_aabb.node[child[i]].bb_max);
	}
	if (inside_flag[1] && !inside_flag[0])
		mySwap(child[0], child[1]);
	else if ((this_farthest_point[1] - point).SquareLength() < (this_farthest_point[0] - point).SquareLength())
		mySwap(child[0], child[1]);

	//	we have stored the more likely to be closer node in child[0]
	for (int i = 0; i<2; i++) {
		Vec2<T>	this_nearest_point;
		inside_flag[i] = getNearestPointOnAABB(this_nearest_point, point,
			edge_aabb.node[child[i]].bb_min, edge_aabb.node[child[i]].bb_max);
		if (inside_flag[i] || (this_nearest_point - point).SquareLength() < nearest_squ_dist) {
			recursiveGetNearestPointOnPolygon(nearest_point, nearest_edge_id, nearest_squ_dist, point,
				vertex, edge, edge_aabb, child[i]);
		}
	}
}

/**
	get the nearest point of a given point on the polygon using AABB tree

	\param	nearest_point	return the nearest point on polygon
	\param	nearest_edge_id	return id of the edge that the nearest point lies in
	\param	point			the testing point
	\param	vertex			vertex list of polygon
	\param	edge			edge list of polygon
	\param	edge_aabb		edge aabb tree of polygon
	\return					the edge id containning the nearest point
*/
template<typename T>
int getNearestPointOnPolygon(
	Vec2<T>&						nearest_point,
	int&							nearest_edge_id,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb
) {
	assert( edge.size() == edge_aabb.leaf_number );

	//	if the mesh contain just one triangle, calculate directly
	if( edge_aabb.leaf_number == 1 ){
		getNearestPointOnSegment(nearest_point, point,
			vertex[edge[0].x], vertex[edge[0].y]);
		return 0;
	}

	//	we start from the root with farthest distance as starting point
	int node_id = edge_aabb.leaf_number;	//	index of root
	getFarthestPointOnAABB(nearest_point, point, 
		edge_aabb.node[node_id].bb_min, edge_aabb.node[node_id].bb_max);

	T nearest_squ_dist = (nearest_point - point).SquareLength();

	recursiveGetNearestPointOnPolygon(nearest_point, nearest_edge_id, nearest_squ_dist, point,
		vertex, edge, edge_aabb, node_id);

	return nearest_edge_id;
}

/**
	get the nearest point of a given point on the polygon using AABB tree

	\param	nearest_point	return the nearest point on polygon
	\param	point			the testing point
	\param	vertex			vertex list of polygon
	\param	edge			edge list of polygon
	\param	edge_aabb		edge aabb tree of polygon
	\return					the edge id containning the nearest point
*/
template<typename T>
int getNearestPointOnPolygon(
	Vec2<T>&						nearest_point,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb) {
	int nearest_edge_id;
	return getNearestPointOnPolygon(nearest_point, nearest_edge_id, point, vertex, edge, edge_aabb);
}

///@}

}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_POLYGON_DISTANCE_FIELD_H__