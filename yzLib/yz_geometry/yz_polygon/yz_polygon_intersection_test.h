/***********************************************************/
/**	\file
	\brief		Intersection Test of Polygon
	\author		Yizhong Zhang
	\date		1/17/2013
*/
/***********************************************************/
#ifndef __YZ_POLYGON_INTERSECTION_TEST_H__
#define __YZ_POLYGON_INTERSECTION_TEST_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_intersection_test.h" 
#include "yzLib/yz_geometry/yz_clipping.h"

namespace yz{  namespace geometry{  namespace polygon{

//	========================================
///@{
/**	@name Line Segment Polygon Intersection Test
*/
//	========================================
/**
	get intersection points of polygon and line segment given aabb tree

	if the polygon is degenerate (contain horizontal or vertical lines), 
	edge_aabb must be expanded, or the intersection test may fail

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\param	edge_aabb			edge aabb tree
	\return						the number of intersection points
*/
template<typename T>
int getSegmentPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							seg_v0,
	Vec2<T>							seg_v1,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb
) {
	//	clear old data
	intersection_points.clear();
	edge_id.clear();

	//	check aabb tree
	if( edge_aabb.leaf_number*2-1 != edge_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}	
	if( edge_aabb.leaf_number != edge.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree is not built by this mesh, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}

	//	clip the line segment to root bonding box
	int inters_flag = clipCohenSutherland(seg_v0, seg_v1,
		edge_aabb.node[edge_aabb.leaf_number].bb_min, edge_aabb.node[edge_aabb.leaf_number].bb_max);	//	node[leaf_number] is root

	if(inters_flag){
		recursiveGetSegmentPolygonIntersectionPoints(intersection_points, edge_id,
			seg_v0, seg_v1, vertex, edge, edge_aabb, edge_aabb.leaf_number);	//	scan from the root, whose id is leaf_number
	}


	return intersection_points.size();
}

/**
	recursively get intersection points of polygon and line segment with AABB tree acceleration

	if the polygon is degenerate (contain horizontal or vertical lines), 
	edge_aabb must be expanded, or the intersection test may fail

	called by getSegmentMeshIntersectionPoints()

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\param	edge_aabb			aabb tree created for polygon edge
	\param	node_id				current scan node id
	\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
void recursiveGetSegmentPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							seg_v0,
	Vec2<T>							seg_v1,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb,
	int								node_id
){
	//	if the node is leaf, then check line segment - triangle intersection
	if( edge_aabb.node[node_id].left == -1 || edge_aabb.node[node_id].left == -1 ){	
		Vec2<T>	v0(vertex[edge[node_id].x]);
		Vec2<T>	v1(vertex[edge[node_id].y]);

		Vec2<T>	inter_p;
		int flag = getSegmentSegmentIntersectionPoint(inter_p, 
			seg_v0, seg_v1, v0, v1);

		if( flag ){
			intersection_points.push_back(inter_p);
			edge_id.push_back(node_id);
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id		= edge_aabb.node[node_id].left;
	Vec2<T> new_seg_v0(seg_v0), new_seg_v1(seg_v1);
	int inters_flag = clipCohenSutherland(new_seg_v0, new_seg_v1, 
		edge_aabb.node[left_id].bb_min, edge_aabb.node[left_id].bb_max);	//	line segment is cliped to the bonding box
	if(inters_flag){
		recursiveGetSegmentPolygonIntersectionPoints(intersection_points, edge_id,
			new_seg_v0, new_seg_v1, vertex, edge, edge_aabb, left_id);
	}

	//	right sub tree
	int right_id	= edge_aabb.node[node_id].right;
	inters_flag = clipCohenSutherland(seg_v0, seg_v1, 
		edge_aabb.node[right_id].bb_min, edge_aabb.node[right_id].bb_max);
	if(inters_flag){
		recursiveGetSegmentPolygonIntersectionPoints(intersection_points, edge_id,
			seg_v0, seg_v1, vertex, edge, edge_aabb, right_id);
	}
}

/**
	get intersection points of polygon and line segment

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\return						the number of intersection points
*/
template<typename T>
int getSegmentMeshIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							seg_v0,
	Vec2<T>							seg_v1,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge
) {
	//	clear old data
	intersection_points.clear();
	edge_id.clear();

	//	scan each edge
	for( int i=0; i<edge.size(); i++ ){
		Vec2<T>	v0(vertex[edge[i].x]);
		Vec2<T>	v1(vertex[edge[i].y]);

		Vec2<T>	inter_p;
		int flag = getSegmentSegmentIntersectionPoint(inter_p, 
			seg_v0, seg_v1, v0, v1);

		if( flag ){
			intersection_points.push_back(inter_p);
			edge_id.push_back(i);
		}
	}

	return intersection_points.size();
}

///@}

//	========================================
///@{
/**	@name Ray Polygon Intersection Test
*/
//	========================================

/**
	get intersection points of polygon and ray with help of AABB tree

	if the polygon is degenerate (contain horizontal or vertical lines), 
	edge_aabb must be expanded, or the intersection test may fail

	we just clip the ray to line segment to do intersection test

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	ray_origin			origin point of ray
	\param	ray_next			next point of ray
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\param	edge_aabb			aabb tree of edge
	\return						the number of intersection points
*/
template<typename T>
int getRayPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							ray_origin,
	Vec2<T>							ray_next,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb
) {
	//	clear old data
	intersection_points.clear();
	edge_id.clear();

	//	calculate line segmet
	Vec2<T>	point_enter;
	int flag = getRayAABBIntersectPoint(point_enter, ray_origin, ray_next,
		edge_aabb.node[edge_aabb.leaf_number].bb_min, edge_aabb.node[edge_aabb.leaf_number].bb_max);
	if (!flag)	return 0;		//	ray bonding box don't intersect at all

	//	ray intersect the bonding box, now calculate the exit point
	Vec2<T>	ray_dir = (ray_next - ray_origin).Normalize();
	Vec2<T>	bbmin = edge_aabb.node[edge_aabb.leaf_number].bb_min;
	Vec2<T>	bbmax = edge_aabb.node[edge_aabb.leaf_number].bb_max;
	T lamda[2], distance[2];
	for( int i=0; i<2; i++ ){
		distance[i] = (ray_dir[i] > 0 ? bbmax[i]-point_enter[i] : bbmin[i]-point_enter[i]);	//	distance to bondary
		lamda[i] = ( fabs(ray_dir[i])<1e-5 ? 1e6 : distance[i]/ray_dir[i] );				//	steps go to bondary
	}
	lamda[0] = myMin(lamda[0], lamda[1]);	//	which direction cost minimal steps

	Vec2<T>	point_exit = point_enter + ray_dir * lamda[0];

	//	intersection test using line segment and mesh
	return getSegmentPolygonIntersectionPoints(intersection_points, edge_id,
		point_enter, point_exit, vertex, edge, edge_aabb);
}

/**
recursively get intersection points of polygon and +x ray with AABB tree acceleration

if the polygon is degenerate (contain horizontal or vertical lines),
edge_aabb must be expanded, or the intersection test may fail

\param	intersection_points	all intersection points
\param	edge_id				penetrated edge id, correspond to intersection_points
\param	xray_origin			origin of +x direction ray
\param	vertex				vertex list of the polygon
\param	edge				edge list of the polygon
\param	edge_aabb			aabb tree created for polygon edge
\param	node_id				current scan node id
\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
void recursiveGetXRayPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							xray_origin,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb,
	int								node_id
) {
	//	if the node is leaf, then check ray - triangle intersection
	if (edge_aabb.node[node_id].left == -1 || edge_aabb.node[node_id].left == -1) {
		Vec2<T>	v0(vertex[edge[node_id].x]);
		Vec2<T>	v1(vertex[edge[node_id].y]);

		Vec2<T>	inter_p;
		int flag = getRaySegmentIntersectionPoint(inter_p,
			xray_origin, xray_origin + Vec2<T>(1, 0), v0, v1);

		if (flag) {
			intersection_points.push_back(inter_p);
			edge_id.push_back(node_id);
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id = edge_aabb.node[node_id].left;
	int inters_flag = isXRayAABBIntersect(xray_origin,
		edge_aabb.node[left_id].bb_min, edge_aabb.node[left_id].bb_max);
	if (inters_flag) {
		recursiveGetXRayPolygonIntersectionPoints(intersection_points, edge_id,
			xray_origin, vertex, edge, edge_aabb, left_id);
	}

	//	right sub tree
	int right_id = edge_aabb.node[node_id].right;
	inters_flag = isXRayAABBIntersect(xray_origin,
		edge_aabb.node[right_id].bb_min, edge_aabb.node[right_id].bb_max);
	if (inters_flag) {
		recursiveGetXRayPolygonIntersectionPoints(intersection_points, edge_id,
			xray_origin, vertex, edge, edge_aabb, right_id);
	}
}

/**
	get intersection points of mesh and +x direction ray with help of AABB tree

	if the polygon is degenerate (contain horizontal or vertical lines), 
	edge_aabb must be expanded, or the intersection test may fail

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	xray_origin			origin point of +x ray
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\param	edge_aabb			aabb tree of edge
	\return						the number of intersection points
*/
template<typename T>
int getXRayPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							xray_origin,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb
) {
	//	clear old data
	intersection_points.clear();
	edge_id.clear();

	//	check aabb tree
	if( edge_aabb.leaf_number*2-1 != edge_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}
	if( edge_aabb.leaf_number != edge.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree is not built by this mesh, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}

	//	recursively get intersection points
	recursiveGetXRayPolygonIntersectionPoints(intersection_points, edge_id,
		xray_origin, vertex, edge, edge_aabb, edge_aabb.leaf_number);

	return intersection_points.size();
}

/**
	get intersection points of polygon and ray

	\param	intersection_points	all intersection points
	\param	edge_id				penetrated edge id, correspond to intersection_points
	\param	ray_origin			origin point of ray
	\param	ray_next			next point of ray
	\param	vertex				vertex list of the polygon
	\param	edge				edge list of the polygon
	\return						the number of intersection points
*/
template<typename T>
int getRayPolygonIntersectionPoints(
	std::vector<Vec2<T>>&			intersection_points,
	std::vector<int>&				edge_id,
	Vec2<T>							ray_origin,
	Vec2<T>							ray_next,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge
) {
	//	clear old data
	intersection_points.clear();
	edge_id.clear();

	//	scan each edge
	for (int i = 0; i < edge.size(); i++) {
		Vec2<T>	v0(vertex[edge[i].x]);
		Vec2<T>	v1(vertex[edge[i].y]);

		Vec2<T>	inter_p;
		int flag = getRaySegmentIntersectionPoint(inter_p, 
			ray_origin, ray_next, v0, v1);

		if( flag ){
			intersection_points.push_back(inter_p);
			edge_id.push_back(i);
		}
	}

	return intersection_points.size();
}


///@}

//	========================================
///@{
/**	@name Point Polygon Intersection Test
*/
//	========================================
/**
check whether a point is inside of a closed polygon with aabb tree acceleration

use +x ray to intersect with the mesh and count the intersection point number,
odd number means the point is inside mesh while even number indicates outside.

This method is only valid when the polygon is closed, and all the vertices are in
general position (such as aligned to a plane). If degenerate case appear, add
random noise to each vertex.

\param	point				the point to check
\param	vertex				vertex list of the polygon
\param	edge				edge list of the polygon
\param	edge_aabb			aabb tree created for polygon edge
\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
int isPointInsidePolygon(
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int2>&		edge,
	const AABBTree2D<T>&			edge_aabb
) {
	std::vector<Vec2<T>>	inter_p;
	std::vector<int>		edge_id;
	int inter_num = getXRayPolygonIntersectionPoints(inter_p, edge_id,
		point, vertex, edge, edge_aabb);
	if (inter_num % 2)
		return 1;
	else
		return 0;
}

///@}

}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_POLYGON_INTERSECTION_TEST_H__