/***********************************************************/
/**	\file
	\brief		distance field of mesh
	\author		Yizhong Zhang
	\date		10/7/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_NN_H__
#define __YZ_MESH_NN_H__

#include <iostream>
#include <vector>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_intersection_test.h"


namespace yz{	namespace geometry{

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
inline int getNearestPointOnTriangle(
	Vec2<T>&	nearest_point,
	Vec2<T>		point,
	Vec2<T>		v0,
	Vec2<T>		v1,
	Vec2<T>		v2
) {
	T coef1, coef2;
	getPointTriangleIntersectionCoef(coef1, coef2, point, v0, v1, v2);

	//	the point is inside of triangle
	if (coef1 >= 0 && coef2 >= 0 && coef1 + coef2 <= 1){
		nearest_point = point;
		return 1;
	}

	//	projection point is outside of triangle
	Vec2<T>	tmp[3];
	getNearestPointOnSegment(tmp[0], point, v0, v1);
	getNearestPointOnSegment(tmp[1], point, v0, v2);
	getNearestPointOnSegment(tmp[2], point, v1, v2);

	T squ_dist[3];
	for (int i = 0; i<3; i++)
		squ_dist[i] = (point - tmp[i]).SquareLength();

	if (squ_dist[0] <= squ_dist[1] && squ_dist[0] <= squ_dist[2])
		nearest_point = tmp[0];
	else if (squ_dist[1] <= squ_dist[0] && squ_dist[1] <= squ_dist[2])
		nearest_point = tmp[1];
	else
		nearest_point = tmp[2];

	return 0;
}

//	========================================
///@{
/**	@name Nearest Point On 2D Mesh
*/
//	========================================

/**
	get the nearest point of a given point on the mesh using AABB tree

	\param	nearest_point	return the nearest point on mesh
	\param	nearest_face_id	return id of the face that the nearest point lies in
	\param	point			the testing point
	\param	vertex			vertex list of mesh
	\param	face			face list of mesh
	\param	face_aabb		face aabb tree of mesh
	\return					the face id containning the nearest point
*/
template<typename T>
int getNearestPointOnMesh(
	Vec2<T>&						nearest_point,
	int&							nearest_face_id,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int3>&		face,
	const AABBTree2D<T>&			face_aabb
) {
	assert(face.size() == face_aabb.leaf_number);

	//	if the mesh contain just one triangle, calculate directly
	if (face_aabb.leaf_number == 1){
		getNearestPointOnTriangle(nearest_point, point,
			vertex[face[0].x], vertex[face[0].y], vertex[face[0].z]);
		return 0;
	}

	//	we start from the root with farthest distance as starting point
	int node_id = face_aabb.leaf_number;	//	index of root
	getFarthestPointOnAABB(nearest_point, point,
		face_aabb.node[node_id].bb_min, face_aabb.node[node_id].bb_max);

	T nearest_squ_dist = (nearest_point - point).SquareLength();

	recursiveGetNearestPointOnMesh(nearest_point, nearest_face_id, nearest_squ_dist, point,
		vertex, face, face_aabb, node_id);

	return nearest_face_id;
}

/**
	recursively get the nearest point on mesh using AABB tree

	\param	nearest_point		return the nearest point position
	\param	nearest_face_id		the id of the face that contain the nearest point
	\param	nearest_squ_dist	square distance of point to nearest point
	\param	point				the testing point
	\param	vertex				vertex list of mesh
	\param	face				face list of mesh
	\param	face_aabb			face aabb tree of the mesh
	\param	node_id				current checking node id
*/
template<typename T>
void recursiveGetNearestPointOnMesh(
	Vec2<T>&						nearest_point,
	int&							nearest_face_id,
	T&								nearest_squ_dist,
	Vec2<T>							point,
	const std::vector<Vec2<T>>&		vertex,
	const std::vector<int3>&		face,
	const AABBTree2D<T>&			face_aabb,
	int								node_id
) {
	//	leaf has reached
	if (node_id < face_aabb.leaf_number){
		Vec2<T>	this_nearest_point;
		getNearestPointOnTriangle(this_nearest_point, point,
			vertex[face[node_id].x], vertex[face[node_id].y], vertex[face[node_id].z]);

		T this_squ_dist = (this_nearest_point - point).SquareLength();
		if (this_squ_dist < nearest_squ_dist){
			nearest_squ_dist = this_squ_dist;
			nearest_point = this_nearest_point;
			nearest_face_id = node_id;
		}

		return;
	}

	//	check two children, if point is inside the bonding box or 
	//	closer than current nearest distance, then go deeper
	//	in order to accelerate this process, we first pick out the node 
	//	that is closer to the point, in this way we have a chance that 
	//	we don't need to parse the second node any more
	int child[2] = { face_aabb.node[node_id].left, face_aabb.node[node_id].right };
	int inside_flag[2];
	Vec2<T>	this_farthest_point[2];
	for (int i = 0; i<2; i++){
		inside_flag[i] = getFarthestPointOnAABB(this_farthest_point[i], point,
			face_aabb.node[child[i]].bb_min, face_aabb.node[child[i]].bb_max);
	}
	if (inside_flag[1] && !inside_flag[0])
		mySwap(child[0], child[1]);
	else if ((this_farthest_point[1] - point).SquareLength() < (this_farthest_point[0] - point).SquareLength())
		mySwap(child[0], child[1]);

	//	we have stored the more likely to be closer node in child[0]
	for (int i = 0; i<2; i++){
		Vec2<T>	this_nearest_point;
		inside_flag[i] = getNearestPointOnAABB(this_nearest_point, point,
			face_aabb.node[child[i]].bb_min, face_aabb.node[child[i]].bb_max);
		if (inside_flag[i] || (this_nearest_point - point).SquareLength() < nearest_squ_dist){
			recursiveGetNearestPointOnMesh(nearest_point, nearest_face_id, nearest_squ_dist, point,
				vertex, face, face_aabb, child[i]);
		}
	}
}

///@}

//	========================================
///@{
/**	@name Nearest Point On 3D Mesh
*/
//	========================================
/**
	get the nearest point of a given point on the mesh using AABB tree

	\param	nearest_point	return the nearest point on mesh
	\param	point			the testing point
	\param	vertex			vertex list of mesh
	\param	face			face list of mesh
	\return					the face id containning the nearest point
*/
template<typename T>
int getNearestPointOnMesh(
	Vec3<T>&						nearest_point,
	Vec3<T>							point,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face
){
	int nearest_face_id;
	return getNearestPointOnMesh(nearest_point, nearest_face_id, point, vertex, face);
}

/**
	get the nearest point of a given point on the mesh using AABB tree

	\param	nearest_point	return the nearest point on mesh
	\param	nearest_face_id	return id of the face that the nearest point lies in
	\param	point			the testing point
	\param	vertex			vertex list of mesh
	\param	face			face list of mesh
	\return					the face id containning the nearest point
*/
template<typename T>
int getNearestPointOnMesh(
	Vec3<T>&						nearest_point,
	int&							nearest_face_id,
	Vec3<T>							point,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face
) {
	//	create face aabb tree
	int i;
	AABBTree3D<T>	face_aabb;
	face_aabb.BuildTriangleAABBTree(vertex, face);

	//	calculate nearest point using AABB Tree
	return getNearestPointOnMesh(nearest_point, nearest_face_id,
		point, vertex, face, face_aabb);
}

/**
	get the nearest point of a given point on the mesh using AABB tree

	\param	nearest_point	return the nearest point on mesh
	\param	point			the testing point
	\param	vertex			vertex list of mesh
	\param	face			face list of mesh
	\param	face_aabb		face aabb tree of mesh
	\return					the face id containning the nearest point
*/
template<typename T>
int getNearestPointOnMesh(
	Vec3<T>&						nearest_point,
	Vec3<T>							point,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face,
	const AABBTree3D<T>&			face_aabb
) {
	int nearest_face_id;
	return getNearestPointOnMesh(nearest_point, nearest_face_id, point, vertex, face, face_aabb);
}

/**
	get the nearest point of a given point on the mesh using AABB tree

	\param	nearest_point	return the nearest point on mesh
	\param	nearest_face_id	return id of the face that the nearest point lies in
	\param	point			the testing point
	\param	vertex			vertex list of mesh
	\param	face			face list of mesh
	\param	face_aabb		face aabb tree of mesh
	\return					the face id containning the nearest point
*/
template<typename T>
int getNearestPointOnMesh(
	Vec3<T>&						nearest_point,
	int&							nearest_face_id,
	Vec3<T>							point,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face,
	const AABBTree3D<T>&			face_aabb
) {
	assert( face.size() == face_aabb.leaf_number );

	//	if the mesh contain just one triangle, calculate directly
	if( face_aabb.leaf_number == 1 ){
		getNearestPointOnTriangle(nearest_point, point,
			vertex[face[0].x], vertex[face[0].y], vertex[face[0].z]);
		return 0;
	}

	//	we start from the root with farthest distance as starting point
	int node_id = face_aabb.leaf_number;	//	index of root
	getFarthestPointOnAABB(nearest_point, point, 
		face_aabb.node[node_id].bb_min, face_aabb.node[node_id].bb_max);

	T nearest_squ_dist = (nearest_point - point).SquareLength();

	recursiveGetNearestPointOnMesh(nearest_point, nearest_face_id, nearest_squ_dist, point,
		vertex, face, face_aabb, node_id);

	return nearest_face_id;
}

/**
	recursively get the nearest point on mesh using AABB tree

	\param	nearest_point		return the nearest point position
	\param	nearest_face_id		the id of the face that contain the nearest point
	\param	nearest_squ_dist	square distance of point to nearest point
	\param	point				the testing point
	\param	vertex				vertex list of mesh
	\param	face				face list of mesh
	\param	face_aabb			face aabb tree of the mesh
	\param	node_id				current checking node id
*/
template<typename T>
void recursiveGetNearestPointOnMesh(
	Vec3<T>&						nearest_point,
	int&							nearest_face_id,
	T&								nearest_squ_dist,
	Vec3<T>							point,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face,
	const AABBTree3D<T>&			face_aabb,
	int								node_id
) {
	//	leaf has reached
	if( node_id < face_aabb.leaf_number ){
		Vec3<T>	this_nearest_point;
		getNearestPointOnTriangle(this_nearest_point, point, 
			vertex[face[node_id].x], vertex[face[node_id].y], vertex[face[node_id].z]);

		T this_squ_dist = (this_nearest_point - point).SquareLength();
		if( this_squ_dist < nearest_squ_dist ){
			nearest_squ_dist	= this_squ_dist;
			nearest_point		= this_nearest_point;
			nearest_face_id		= node_id;
		}

		return;
	}

	//	check two children, if point is inside the bonding box or 
	//	closer than current nearest distance, then go deeper
	//	in order to accelerate this process, we first pick out the node 
	//	that is closer to the point, in this way we have a chance that 
	//	we don't need to parse the second node any more
	int child[2] = {face_aabb.node[node_id].left, face_aabb.node[node_id].right};
	int inside_flag[2];
	Vec3<T>	this_farthest_point[2];
	for(int i=0; i<2; i++){
		inside_flag[i] = getFarthestPointOnAABB(this_farthest_point[i], point,
			face_aabb.node[child[i]].bb_min, face_aabb.node[child[i]].bb_max);
	}
	if( inside_flag[1] && !inside_flag[0] )
		mySwap(child[0], child[1]);
	else if( (this_farthest_point[1]-point).SquareLength() < (this_farthest_point[0]-point).SquareLength() )
		mySwap(child[0], child[1]);

	//	we have stored the more likely to be closer node in child[0]
	for(int i=0; i<2; i++){
		Vec3<T>	this_nearest_point;
		inside_flag[i] = getNearestPointOnAABB(this_nearest_point, point,
			face_aabb.node[child[i]].bb_min, face_aabb.node[child[i]].bb_max);
		if( inside_flag[i] || (this_nearest_point-point).SquareLength() < nearest_squ_dist ){
			recursiveGetNearestPointOnMesh( nearest_point, nearest_face_id, nearest_squ_dist, point,
				vertex, face, face_aabb, child[i]);
		}
	}
}

///@}

//	========================================
///@{
/**	@name Nearest Point On Point Cloud
*/
//	========================================

/**
get the nearest point of a given point in the point cloud using AABB tree

\param	nearest_point		return the nearest point
\param	nearest_point_id	return id of the nearest vertex
\param	point				the testing point
\param	vertex				vertex list of target point cloud
\param	vertex_aabb			vertex aabb tree of point cloud
\return						the vertex id of the nearest point
*/
template<typename VEC_NT>
int getNearestPointOnVertices(
	VEC_NT&							nearest_point,
	int&							nearest_point_id,
	VEC_NT							point,
	const std::vector<VEC_NT>&		vertex,
	const AABBTree<VEC_NT>&			vertex_aabb
) {
	assert(vertex.size() == vertex_aabb.leaf_number);

	//	if the mesh contain just one triangle, calculate directly
	if (vertex_aabb.leaf_number == 1) {
		nearest_point = vertex[0];
		nearest_point_id = 0;
		return 0;
	}

	//	we start from the root with farthest distance as starting point
	int node_id = vertex_aabb.leaf_number;	//	index of root
	getFarthestPointOnAABB(nearest_point, point,
		vertex_aabb.node[node_id].bb_min, vertex_aabb.node[node_id].bb_max);

	typename VEC_NT::Type nearest_squ_dist = (nearest_point - point).SquareLength();

	recursiveGetNearestPointOnVertices(
		nearest_point, nearest_point_id, nearest_squ_dist, point,
		vertex, vertex_aabb, node_id);

	return nearest_point_id;
}

/**
recursively get the nearest point on mesh using AABB tree

\param	nearest_point		return the nearest point position
\param	nearest_face_id		the id of the face that contain the nearest point
\param	nearest_squ_dist	square distance of point to nearest point
\param	point				the testing point
\param	vertex				vertex list
\param	vertex_aabb			vertex aabb tree of the point cloud
\param	node_id				current checking node id
*/
template<typename VEC_NT>
void recursiveGetNearestPointOnVertices(
	VEC_NT&							nearest_point,
	int&							nearest_point_id,
	typename VEC_NT::Type&			nearest_squ_dist,
	VEC_NT							point,
	const std::vector<VEC_NT>&		vertex,
	const AABBTree<VEC_NT>&			vertex_aabb,
	int								node_id
) {
	//	leaf has reached
	if (node_id < vertex_aabb.leaf_number) {
		typename VEC_NT::Type this_squ_dist = (vertex[node_id] - point).SquareLength();
		if (this_squ_dist < nearest_squ_dist) {
			nearest_point = vertex[node_id];
			nearest_point_id = node_id;
			nearest_squ_dist = this_squ_dist;
		}

		return;
	}

	//	check two children, if point is inside the bonding box or 
	//	closer than current nearest distance, then go deeper
	//	in order to accelerate this process, we first pick out the node 
	//	that is closer to the point, in this way we have a chance that 
	//	we don't need to parse the second node any more
	int child[2] = { vertex_aabb.node[node_id].left, vertex_aabb.node[node_id].right };
	int inside_flag[2];
	VEC_NT	this_farthest_point[2];
	for (int i = 0; i<2; i++) {
		inside_flag[i] = getFarthestPointOnAABB(this_farthest_point[i], point,
			vertex_aabb.node[child[i]].bb_min, vertex_aabb.node[child[i]].bb_max);
	}
	if (inside_flag[1] && !inside_flag[0])
		mySwap(child[0], child[1]);
	else if ((this_farthest_point[1] - point).SquareLength() < (this_farthest_point[0] - point).SquareLength())
		mySwap(child[0], child[1]);

	//	we have stored the more likely to be closer node in child[0]
	for (int i = 0; i<2; i++) {
		VEC_NT	this_nearest_point;
		inside_flag[i] = getNearestPointOnAABB(this_nearest_point, point,
			vertex_aabb.node[child[i]].bb_min, vertex_aabb.node[child[i]].bb_max);
		if (inside_flag[i] || (this_nearest_point - point).SquareLength() < nearest_squ_dist) {
			recursiveGetNearestPointOnVertices(
				nearest_point, nearest_point_id, nearest_squ_dist, point,
				vertex, vertex_aabb, child[i]);
		}
	}
}

/**
get the k nearest points of a given point in the point cloud using AABB tree

\param	k_nearest_point		return the k nearest points
\param	k_nearest_point_id	return id of the nearest vertices
\param	k					number of nearest points
\param	point				the testing point
\param	vertex				vertex list of target point cloud
\param	vertex_aabb			vertex aabb tree of point cloud
\return						distance of the k-th vertex (negative means invalid)
*/
template<typename VEC_NT>
typename VEC_NT::Type getKNearestPointsOnVertices(
	std::vector<VEC_NT>&			k_nearest_point,
	std::vector<int>&				k_nearest_point_id,
	int								k,
	VEC_NT							point,
	const std::vector<VEC_NT>&		vertex,
	const AABBTree<VEC_NT>&			vertex_aabb
) {
	assert(vertex.size() == vertex_aabb.leaf_number);

	if (k <= 0) {
		k_nearest_point.clear();
		k_nearest_point_id.clear();
		return -1;
	}
	if (k == 1) {
		VEC_NT	np;
		int		npid;
		getNearestPointOnVertices(np, npid, point, vertex, vertex_aabb);
		k_nearest_point.clear();
		k_nearest_point_id.clear();
		k_nearest_point.push_back(np);
		k_nearest_point_id.push_back(npid);
		return (np - point).Length();
	}

	//	create initial list of the first k points
	if (vertex.size() < k)
		k = vertex.size();

	std::vector<typename VEC_NT::Type>	k_nearest_point_squ_dist;
	k_nearest_point.clear();
	k_nearest_point_id.clear();
	for (int i = 0; i < k; i++) {
		k_nearest_point.push_back(vertex[i]);
		k_nearest_point_id.push_back(i);
		k_nearest_point_squ_dist.push_back((point - vertex[i]).SquareLength());
	}

	//	order from smallest to largest
	std::vector<int> order;
	for (int i = 0; i < k; i++)
		order.push_back(i);
	utils::sortIndex(order, k_nearest_point_squ_dist.begin(), k_nearest_point_squ_dist.end());
	utils::reorderDataBySourceOrder(k_nearest_point, order);
	utils::reorderDataBySourceOrder(k_nearest_point_id, order);
	utils::reorderDataBySourceOrder(k_nearest_point_squ_dist, order);

	if (vertex.size() == k)
		return sqrt(k_nearest_point_squ_dist.back());

	//	we start from the root with farthest distance as starting point
	int node_id = vertex_aabb.leaf_number;	//	index of root
	recursiveGetKNearestPointsOnVertices(
		k_nearest_point, k_nearest_point_id, k_nearest_point_squ_dist, k, point,
		vertex, vertex_aabb, node_id);

	return sqrt(k_nearest_point_squ_dist.back());
}


/**
recursively get the nearest point on mesh using AABB tree

\param	k_nearest_point				return the current k nearest point position
\param	k_nearest_point_id			the id of the face that contain the nearest point
\param	k_nearest_point_squ_dist	square distance of point to nearest point
\param	k							number of nearest points
\param	point						the testing point
\param	vertex						vertex list
\param	vertex_aabb					vertex aabb tree of the point cloud
\param	node_id						current checking node id
*/
template<typename VEC_NT>
void recursiveGetKNearestPointsOnVertices(
	std::vector<VEC_NT>&			k_nearest_point,
	std::vector<int>&				k_nearest_point_id,
	std::vector<typename VEC_NT::Type>&	k_nearest_point_squ_dist,
	int								k,
	VEC_NT							point,
	const std::vector<VEC_NT>&		vertex,
	const AABBTree<VEC_NT>&			vertex_aabb,
	int								node_id
) {
	//	leaf has reached
	if (node_id < vertex_aabb.leaf_number) {
		if (node_id < k)	//	we already processed the first k nodes, skip
			return;

		//	if this vertex is closer to the current k-th vertex, insert it
		typename VEC_NT::Type this_squ_dist = (vertex[node_id] - point).SquareLength();
		if (this_squ_dist < k_nearest_point_squ_dist.back()) {
			//	insert this point into the increasing array
			int idx = k - 1;
			while (idx > 0 && k_nearest_point_squ_dist[idx - 1] > this_squ_dist) {
				k_nearest_point[idx] = k_nearest_point[idx - 1];
				k_nearest_point_id[idx] = k_nearest_point_id[idx - 1];
				k_nearest_point_squ_dist[idx] = k_nearest_point_squ_dist[idx - 1];
				idx--;
			}

			k_nearest_point[idx] = vertex[node_id];
			k_nearest_point_id[idx] = node_id;
			k_nearest_point_squ_dist[idx] = this_squ_dist;
		}
		return;
	}

	//	check two children, if point is inside the bonding box or 
	//	closer than current nearest distance, then go deeper
	//	in order to accelerate this process, we first pick out the node 
	//	that is closer to the point, in this way we have a chance that 
	//	we don't need to parse the second node any more
	int child[2] = { vertex_aabb.node[node_id].left, vertex_aabb.node[node_id].right };
	int inside_flag[2];
	VEC_NT	this_farthest_point[2];
	for (int i = 0; i<2; i++) {
		inside_flag[i] = getFarthestPointOnAABB(this_farthest_point[i], point,
			vertex_aabb.node[child[i]].bb_min, vertex_aabb.node[child[i]].bb_max);
	}
	if (inside_flag[1] && !inside_flag[0])
		mySwap(child[0], child[1]);
	else if ((this_farthest_point[1] - point).SquareLength() < (this_farthest_point[0] - point).SquareLength())
		mySwap(child[0], child[1]);

	//	we have stored the more likely to be closer node in child[0]
	for (int i = 0; i<2; i++) {
		VEC_NT	this_nearest_point;
		inside_flag[i] = getNearestPointOnAABB(this_nearest_point, point,
			vertex_aabb.node[child[i]].bb_min, vertex_aabb.node[child[i]].bb_max);
		if (inside_flag[i] || (this_nearest_point - point).SquareLength() < k_nearest_point_squ_dist.back()) {
			recursiveGetKNearestPointsOnVertices(
				k_nearest_point, k_nearest_point_id, k_nearest_point_squ_dist, k, point,
				vertex, vertex_aabb, child[i]);
		}
	}
}

///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_NN_H__