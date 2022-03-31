/***********************************************************/
/**	\file
	\brief		Intersection Test of Mesh
	\author		Yizhong Zhang
	\date		9/17/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_INTERSECTION_TEST_H__
#define __YZ_MESH_INTERSECTION_TEST_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_intersection_test.h" 
#include "yzLib/yz_geometry/yz_clipping.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Point Mesh Intersection Test
*/
//	========================================
/**
	check whether a point is inside of a closed mesh with aabb tree acceleration

	use an axis aligned ray to intersect with the mesh and count the intersection point number,
	odd number means the point is inside mesh while even number indicates outside.

	This method is only valid when the mesh is closed, and all the vertices are in
	general position (such as aligned to a plane). If degenerate case appear, add
	random noise to each vertex.

	\param	point				the point to check
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\param	ray_dir				the axis and direction used to project ray
								'x' negative x direction; 'X' positive x direction;
								so for 'y' 'Y' 'z' 'Z', other value treate as 'X'
	\return						whether the point is inside
*/
template<typename T>
int isPointInsideClosedMesh(
	Vec3<T>						point,
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const AABBTree3D<T>&		face_aabb,
	char						ray_dir = 'X')
{
	std::vector<Vec3<T>>	inter_p;
	std::vector<int>		face_id;
	
	int dir = 1;
	if (ray_dir == 'Y')			dir = 2;
	else if (ray_dir == 'Z')	dir = 3;
	else if (ray_dir == 'x')	dir = -1;
	else if (ray_dir == 'y')	dir = -2;
	else if (ray_dir == 'z')	dir = -3;

	int inter_num = getAARayMeshIntersectionPoints(
		inter_p, face_id, point, dir, vertex, face, face_aabb);
	if( inter_num % 2 )
		return 1;
	else
		return 0;
}

/**
	check whether a point is inside the open mesh by axis aligned check

	We check 3 axis, xyz. A legal direction has even number of intersection points

	the point is treated to be inside if directions that got even number is more than that of odd number

	\param	point				the point to check
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\return						whether the point is inside
*/
template<typename T>
int isPointInsideOpenMeshByAxisAlignedCheck(Vec3<T>						point,
											const std::vector<Vec3<T>>&	vertex,
											const std::vector<int3>&	face,
											const AABBTree3D<T>&		face_aabb){
	int count[3] = {0, 0, 0};	//	inside, outside, error
	for( int dir = 1; dir <= 3; dir ++ ){	//	check each axis
		std::vector<Vec3<T>>	inter_p;
		std::vector<int>		face_id;

		//	get intersection point number of both positive and negtive directions
		int inter_num_p = getAARayMeshIntersectionPoints(inter_p, face_id,
			point, dir, vertex, face, face_aabb);
		int inter_num_n = getAARayMeshIntersectionPoints(inter_p, face_id,
			point, -dir, vertex, face, face_aabb);

		//	count number
		if( (inter_num_p + inter_num_n) % 2 )	//	number of intersection point is odd number
			count[2] ++;
		else if( inter_num_p%2 )				//	ray shooting each direction with odd number, inside
			count[0] ++;
		else
			count[1] ++;
	}

	if( count[0] > count[1] )
		return 1;
	else
		return 0;
}

/**
	check whether a 2d point is on 2d mesh with aabb tree acceleration

	\param	point				the point to check
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\return						whether the point is on the mesh
*/
template<typename T>
int isPointOn2DMesh(Vec2<T>						point,
					const std::vector<Vec2<T>>&	vertex,
					const std::vector<int3>&	face,
					const AABBTree2D<T>&		face_aabb ){
	std::vector<int> face_id;
	getFacesContainingPoint(face_id, point, vertex, face, face_aabb);
	return !face_id.empty();
}

/**
	get face id on 2d mesh that contains the point

	it is possible that the 2D mesh overlap, so multi-faces may exist

	\param	face_id				return face id that contains the point
	\param	point				the point to check
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\return						the number of faces that contain the point
*/
template<typename T>
int getFacesContainingPoint(std::vector<int>&			face_id,
							Vec2<T>						point,
							const std::vector<Vec2<T>>&	vertex,
							const std::vector<int3>&	face,
							const AABBTree2D<T>&		face_aabb ){
	//	clear old data
	face_id.clear();

	//	check aabb tree
	if( face_aabb.leaf_number*2-1 != face_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot get face containing point" << std::endl;
		#endif
		return -1;
	}
	if( face_aabb.leaf_number != face.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot get face containing point" << std::endl;
		#endif
		return -1;
	}

	//	clip the line segment to root bonding box
	int inters_flag = isPointInsideAABB(point, 
		face_aabb.node[face_aabb.leaf_number].bb_min, face_aabb.node[face_aabb.leaf_number].bb_max);	//	node[leaf_number] is root

	if(inters_flag){
		recursiveGetFacesContainingPoint(face_id, point, vertex, face, face_aabb, face_aabb.leaf_number);	//	scan from the root, whose id is leaf_number
	}

	return face_id.size();	
}

/**
	recursively get face id on 2d mesh that contains the point

	it is possible that the 2D mesh overlap, so multi-faces may exist

	\param	face_id				return face id that contains the point
	\param	point				the point to check
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\param	node_id				current scan node id
*/
template<typename T>
void recursiveGetFacesContainingPoint(std::vector<int>&				face_id,
									  Vec2<T>						point,
									  const std::vector<Vec2<T>>&	vertex,
									  const std::vector<int3>&		face,
									  const AABBTree2D<T>&			face_aabb,
									  int							node_id){
	//	if the node is leaf, then check whether point in triangle
	if( face_aabb.node[node_id].left == -1 || face_aabb.node[node_id].left == -1 ){	
		Vec2<T>	v0(vertex[face[node_id].x]);
		Vec2<T>	v1(vertex[face[node_id].y]);
		Vec2<T>	v2(vertex[face[node_id].z]);

		int flag = isPointInsideTriangle(point, v0, v1, v2);

		if( flag ){
			face_id.push_back( node_id );
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id		= face_aabb.node[node_id].left;
	int inters_flag = isPointInsideAABB(point, face_aabb.node[left_id].bb_min, face_aabb.node[left_id].bb_max);
	if( inters_flag )
		recursiveGetFacesContainingPoint(face_id, point, vertex, face, face_aabb, left_id);

	//	right sub tree
	int right_id	= face_aabb.node[node_id].right;
	inters_flag = isPointInsideAABB(point, face_aabb.node[right_id].bb_min, face_aabb.node[right_id].bb_max);
	if( inters_flag )
		recursiveGetFacesContainingPoint(face_id, point, vertex, face, face_aabb, right_id);
}

///@}

//	========================================
///@{
/**	@name Line Segment Mesh Intersection Test
*/
//	========================================
/**
	get intersection points of mesh and line segment given aabb tree

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			face aabb tree
	\return						the number of intersection points
*/
template<typename T>
int getSegmentMeshIntersectionPoints(std::vector<Vec3<T>>&			intersection_points,
									 std::vector<int>&				face_id,
									 Vec3<T>						seg_v0,
									 Vec3<T>						seg_v1,
									 const std::vector<Vec3<T>>&	vertex,
									 const std::vector<int3>&		face,
									 const AABBTree3D<T>&			face_aabb ){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	check aabb tree
	if( face_aabb.leaf_number*2-1 != face_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}
	if( face_aabb.leaf_number != face.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree is not built by this mesh, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}

	//	clip the line segment to root bonding box
	int inters_flag = clipCohenSutherland(seg_v0, seg_v1, 
		face_aabb.node[face_aabb.leaf_number].bb_min, face_aabb.node[face_aabb.leaf_number].bb_max);	//	node[leaf_number] is root

	if(inters_flag){
		recursiveGetSegmentMeshIntersectionPoints(intersection_points, face_id,
			seg_v0, seg_v1, vertex, face, face_aabb, face_aabb.leaf_number);	//	scan from the root, whose id is leaf_number
	}


	return intersection_points.size();
}


/**
	recursively get intersection points of mesh and line segment with AABB tree acceleration

	called by getSegmentMeshIntersectionPoints()

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\param	node_id				current scan node id
	\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
void recursiveGetSegmentMeshIntersectionPoints(std::vector<Vec3<T>>&		intersection_points,
											   std::vector<int>&			face_id,
											   Vec3<T>						seg_v0,
											   Vec3<T>						seg_v1,
											   const std::vector<Vec3<T>>&	vertex,
											   const std::vector<int3>&		face,
											   const AABBTree3D<T>&			face_aabb,
											   int							node_id ){
	//	if the node is leaf, then check line segment - triangle intersection
	if( face_aabb.node[node_id].left == -1 || face_aabb.node[node_id].left == -1 ){	
		Vec3<T>	v0(vertex[face[node_id].x]);
		Vec3<T>	v1(vertex[face[node_id].y]);
		Vec3<T>	v2(vertex[face[node_id].z]);

		Vec3<T>	inter_p;
		int flag = getSegmentTriangleIntersectionPoint(inter_p, 
			seg_v0, seg_v1, v0, v1, v2);

		if( flag ){
			intersection_points.push_back(inter_p);
			face_id.push_back(node_id);
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id		= face_aabb.node[node_id].left;
	Vec3<T> new_seg_v0(seg_v0), new_seg_v1(seg_v1);
	int inters_flag = clipCohenSutherland(new_seg_v0, new_seg_v1, 
		face_aabb.node[left_id].bb_min, face_aabb.node[left_id].bb_max);	//	line segment is cliped to the bonding box
	if(inters_flag){
		recursiveGetSegmentMeshIntersectionPoints(intersection_points, face_id,
			new_seg_v0, new_seg_v1, vertex, face, face_aabb, left_id);
	}

	//	right sub tree
	int right_id	= face_aabb.node[node_id].right;
	inters_flag = clipCohenSutherland(seg_v0, seg_v1, 
		face_aabb.node[right_id].bb_min, face_aabb.node[right_id].bb_max);
	if(inters_flag){
		recursiveGetSegmentMeshIntersectionPoints(intersection_points, face_id,
			seg_v0, seg_v1, vertex, face, face_aabb, right_id);
	}
}

/**
	get intersection points of mesh and line segment

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	seg_v0				end point 0 of line segment
	\param	seg_v1				end point 1 of line segment
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\return						the number of intersection points
*/
template<typename T>
int getSegmentMeshIntersectionPoints(std::vector<Vec3<T>>&			intersection_points,
									 std::vector<int>&				face_id,
									 Vec3<T>						seg_v0,
									 Vec3<T>						seg_v1,
									 const std::vector<Vec3<T>>&	vertex,
									 const std::vector<int3>&		face){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	scan each face
	for( int i=0; i<face.size(); i++ ){
		Vec3<T>	v0(vertex[face[i].x]);
		Vec3<T>	v1(vertex[face[i].y]);
		Vec3<T>	v2(vertex[face[i].z]);

		Vec3<T>	inter_p;
		int flag = getSegmentTriangleIntersectionPoint(inter_p, 
			seg_v0, seg_v1, v0, v1, v2);

		if( flag ){
			intersection_points.push_back(inter_p);
			face_id.push_back(i);
		}
	}

	return intersection_points.size();
}

///@}

//	========================================
///@{
/**	@name Ray Mesh Intersection Test
*/
//	========================================


/**
	get intersection points of mesh and ray with help of AABB tree

	we just clip the ray to line segment to do intersection test

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	ray_origin			origin point of ray
	\param	ray_next			next point of ray
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree of face
	\return						the number of intersection points
*/
template<typename T>
int getRayMeshIntersectionPoints(std::vector<Vec3<T>>&			intersection_points,
								 std::vector<int>&				face_id,
								 Vec3<T>						ray_origin,
								 Vec3<T>						ray_next,
								 const std::vector<Vec3<T>>&	vertex,
								 const std::vector<int3>&		face,
								 const AABBTree3D<T>&			face_aabb ){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	calculate line segmet
	Vec3<T>	point_enter;
	int flag = getRayAABBIntersectPoint(point_enter, ray_origin, ray_next, 
		face_aabb.node[face_aabb.leaf_number].bb_min, face_aabb.node[face_aabb.leaf_number].bb_max);
	if( !flag )	return 0;		//	ray bonding box don't intersect at all

	//	ray intersect the bonding box, now calculate the exit point
	Vec3<T>	ray_dir = (ray_next - ray_origin).Normalize();
	Vec3<T>	bbmin = face_aabb.node[face_aabb.leaf_number].bb_min;
	Vec3<T>	bbmax = face_aabb.node[face_aabb.leaf_number].bb_max;
	T lamda[3], distance[3];
	for( int i=0; i<3; i++ ){
		distance[i] = (ray_dir[i] > 0 ? bbmax[i]-point_enter[i] : bbmin[i]-point_enter[i]);	//	distance to bondary
		lamda[i] = ( fabs(ray_dir[i])<1e-5 ? 1e6 : distance[i]/ray_dir[i] );				//	steps go to bondary
	}
	lamda[0] = myMin(lamda[0], myMin(lamda[1], lamda[2]));	//	which direction cost minimal steps

	Vec3<T>	point_exit = point_enter + ray_dir * lamda[0];

	//	intersection test using line segment and mesh
	return getSegmentMeshIntersectionPoints(intersection_points, face_id,
		point_enter, point_exit, vertex, face, face_aabb);
}

/**
	get intersection points of mesh and axis aligned ray with help of AABB tree

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	ray_origin			origin point of ray
	\param	ray_dir				direction of the ray, 1 x, 2 y, 3 z, -1 -x, -2 -y, -3 -z
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree of face
	\return						the number of intersection points
*/
template<typename T>
int getAARayMeshIntersectionPoints(std::vector<Vec3<T>>&		intersection_points,
								   std::vector<int>&			face_id,
								   Vec3<T>						ray_origin,
								   int							ray_dir,
								   const std::vector<Vec3<T>>&	vertex,
								   const std::vector<int3>&		face,
								   const AABBTree3D<T>&			face_aabb ){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	check aabb tree
	if( face_aabb.leaf_number*2-1 != face_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}
	if( face_aabb.leaf_number != face.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree is not built by this mesh, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}

	//	recursively get intersection points
	recursiveGetAARayMeshIntersectionPoints(intersection_points, face_id,
		ray_origin, ray_dir, vertex, face, face_aabb, face_aabb.leaf_number);

	return intersection_points.size();
}

/**
	recursively get intersection points of mesh and axis aligned ray with AABB tree acceleration

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	ray_origin			origin of axis aligned ray
	\param	ray_dir				direction of the ray, 1 x, 2 y, 3 z, -1 -x, -2 -y, -3 -z
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\param	node_id				current scan node id
	\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
void recursiveGetAARayMeshIntersectionPoints(std::vector<Vec3<T>>&		intersection_points,
											std::vector<int>&			face_id,
											Vec3<T>						ray_origin,
											int							ray_dir,
											const std::vector<Vec3<T>>&	vertex,
											const std::vector<int3>&	face,
											const AABBTree3D<T>&		face_aabb,
											int							node_id ){
	//	if the node is leaf, then check ray - triangle intersection
	if( face_aabb.node[node_id].left == -1 || face_aabb.node[node_id].left == -1 ){	
		Vec3<T>	v0(vertex[face[node_id].x]);
		Vec3<T>	v1(vertex[face[node_id].y]);
		Vec3<T>	v2(vertex[face[node_id].z]);

		int sgn = ray_dir>0 ? 1 : -1;
		int idx = sgn>0 ? ray_dir-1 : -1-ray_dir;
		Vec3<T>	dir_vec(0, 0, 0);
		dir_vec[idx] = sgn;

		Vec3<T>	inter_p;
		int flag = getRayTriangleIntersectionPoint(inter_p,
			ray_origin, ray_origin+dir_vec, v0, v1, v2);

		if( flag ){
			intersection_points.push_back(inter_p);
			face_id.push_back(node_id);
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id		= face_aabb.node[node_id].left;
	int inters_flag = isAARayAABBIntersect(ray_origin, ray_dir,
		face_aabb.node[left_id].bb_min, face_aabb.node[left_id].bb_max);
	if(inters_flag){
		recursiveGetAARayMeshIntersectionPoints(intersection_points, face_id,
			ray_origin, ray_dir, vertex, face, face_aabb, left_id);
	}

	//	right sub tree
	int right_id	= face_aabb.node[node_id].right;
	inters_flag = isAARayAABBIntersect(ray_origin, ray_dir,
		face_aabb.node[right_id].bb_min, face_aabb.node[right_id].bb_max);
	if(inters_flag){
		recursiveGetAARayMeshIntersectionPoints(intersection_points, face_id,
			ray_origin, ray_dir, vertex, face, face_aabb, right_id);
	}
}

/**
	get intersection points of mesh and +x direction ray with help of AABB tree

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	xray_origin			origin point of +x ray
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree of face
	\return						the number of intersection points
*/
template<typename T>
int getXRayMeshIntersectionPoints(std::vector<Vec3<T>>&			intersection_points,
								  std::vector<int>&				face_id,
								  Vec3<T>						xray_origin,
								  const std::vector<Vec3<T>>&	vertex,
								  const std::vector<int3>&		face,
								  const AABBTree3D<T>&			face_aabb ){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	check aabb tree
	if( face_aabb.leaf_number*2-1 != face_aabb.node.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree node number doesn't match, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}
	if( face_aabb.leaf_number != face.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: AABB Tree is not built by this mesh, cannot calculate intersection points" << std::endl;
		#endif
		return -1;
	}

	//	recursively get intersection points
	recursiveGetXRayMeshIntersectionPoints(intersection_points, face_id,
		xray_origin, vertex, face, face_aabb, face_aabb.leaf_number);

	return intersection_points.size();
}

/**
	recursively get intersection points of mesh and +x ray with AABB tree acceleration

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	xray_origin			origin of +x direction ray
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	face_aabb			aabb tree created for mesh face
	\param	node_id				current scan node id
	\return						the number of intersection points, -1 if aabb tree error
*/
template<typename T>
void recursiveGetXRayMeshIntersectionPoints(std::vector<Vec3<T>>&		intersection_points,
											std::vector<int>&			face_id,
											Vec3<T>						xray_origin,
											const std::vector<Vec3<T>>&	vertex,
											const std::vector<int3>&	face,
											const AABBTree3D<T>&		face_aabb,
											int							node_id ){
	//	if the node is leaf, then check ray - triangle intersection
	if( face_aabb.node[node_id].left == -1 || face_aabb.node[node_id].left == -1 ){	
		Vec3<T>	v0(vertex[face[node_id].x]);
		Vec3<T>	v1(vertex[face[node_id].y]);
		Vec3<T>	v2(vertex[face[node_id].z]);

		Vec3<T>	inter_p;
		int flag = getXRayTriangleIntersectionPoint(inter_p, xray_origin, v0, v1, v2);

		if( flag ){
			intersection_points.push_back(inter_p);
			face_id.push_back(node_id);
		}

		return;
	}

	//	the node is not leaf, then check left and right sub trees
	//	left sub tree
	int left_id		= face_aabb.node[node_id].left;
	int inters_flag = isXRayAABBIntersect(xray_origin, 
		face_aabb.node[left_id].bb_min, face_aabb.node[left_id].bb_max);
	if(inters_flag){
		recursiveGetXRayMeshIntersectionPoints(intersection_points, face_id,
			xray_origin, vertex, face, face_aabb, left_id);
	}

	//	right sub tree
	int right_id	= face_aabb.node[node_id].right;
	inters_flag = isXRayAABBIntersect(xray_origin, 
		face_aabb.node[right_id].bb_min, face_aabb.node[right_id].bb_max);
	if(inters_flag){
		recursiveGetXRayMeshIntersectionPoints(intersection_points, face_id,
			xray_origin, vertex, face, face_aabb, right_id);
	}
}

/**
	get intersection points of mesh and ray

	\param	intersection_points	all intersection points
	\param	face_id				penetrated face id, correspond to intersection_points
	\param	ray_origin			origin point of ray
	\param	ray_next			next point of ray
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\return						the number of intersection points
*/
template<typename T>
int getRayMeshIntersectionPoints(std::vector<Vec3<T>>&			intersection_points,
								 std::vector<int>&				face_id,
								 Vec3<T>						ray_origin,
								 Vec3<T>						ray_next,
								 const std::vector<Vec3<T>>&	vertex,
								 const std::vector<int3>&		face){
	//	clear old data
	intersection_points.clear();
	face_id.clear();

	//	scan each face
	for( int i=0; i<face.size(); i++ ){
		Vec3<T>	v0(vertex[face[i].x]);
		Vec3<T>	v1(vertex[face[i].y]);
		Vec3<T>	v2(vertex[face[i].z]);

		Vec3<T>	inter_p;
		int flag = getRayTriangleIntersectionPoint(inter_p, 
			ray_origin, ray_next, v0, v1, v2);

		if( flag ){
			intersection_points.push_back(inter_p);
			face_id.push_back(i);
		}
	}

	return intersection_points.size();
}


///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_INTERSECTION_TEST_H__