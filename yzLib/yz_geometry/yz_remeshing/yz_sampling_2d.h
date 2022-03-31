/***********************************************************/
/**	\file
	\brief		Sampling on 2D mesh
	\author		Yizhong Zhang
	\date		1/7/2013
*/
/***********************************************************/
#ifndef __YZ_SAMPLING_2D_H__
#define __YZ_SAMPLING_2D_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_mesh_intersection_test.h"

namespace yz{  namespace geometry{  namespace remeshing{

//	========================================
///@{
/**	@name Uniform Sampling on 2D mesh
*/
//	========================================
/**
	spawn sampling points on 2d mesh with staggered grid (equilateral triangle)

	\param	sample_point		return sample points
	\param	vertex				vertex of the 2D mesh
	\param	face				face of the 2D mesh
	\param	distance			edge length of the equilateral triangle pattern
	\return						the number of sampling points
*/
template<typename T>
int samplingPointsOn2DMeshStaggered(std::vector<Vec2<T>>&		sample_point, 
									const std::vector<Vec2<T>>&	vertex,
									const std::vector<int3>&	face,
									T							distance ){
	//	check
	if( vertex.empty() || face.empty() )
		return 0;
	if( distance < 1e-6 ){
		std::cout << "Error: samplingPointsOn2DMeshStaggered, distance too small: " << distance << std::endl;
		return 0;
	}

	//	clear old data
	sample_point.clear();

	//	create aabb tree
	AABBTree2D<T>	face_aabb;
	face_aabb.BuildTriangleAABBTree(vertex, face);

	//	spawn points
	T x0 = face_aabb.node[face_aabb.leaf_number].bb_min.x + distance/2;
	T y0 = face_aabb.node[face_aabb.leaf_number].bb_min.y + distance/2;
	T dx = distance;
	T dy = distance * sqrtf(3.0) / 2.0;
	int x_num = (face_aabb.node[face_aabb.leaf_number].bb_max.x - x0) / dx + 1;
	int y_num = (face_aabb.node[face_aabb.leaf_number].bb_max.y - y0) / dy + 1;

	for(int j=0; j<y_num; j++){
		T y = y0 + j * dy;
		T x0_new = x0 + dx/2 * (j%2);
		for(int i=0; i<x_num; i++){
			T x = x0_new + i * dx;
			yz::Vec2<T> p(x, y);
			if( isPointOn2DMesh(p, vertex, face, face_aabb) )
				sample_point.push_back(p);
		}
	}

	return sample_point.size();
}

///@}

//	========================================
///@{
/**	@name Poisson Disk Sampling
*/
//	========================================

/**
	poisson disk sampling on 2d mesh with different methods

	\param	sample_point	return sample points
	\param	vertex			vertex of the 2D mesh
	\param	face			face of the 2D mesh
	\param	distance		edge length of the equilateral triangle pattern
	\param	method			character string used to choose method of poisson disk sampling, we use author of the method(case insensitive):	\n
							RobertBridson: method in paper Fast Poisson Disk Sampling in Arbitrary Dimensions	\n
	\return					the number of sampling points
*/
template<typename T>
int samplingPoissonDiskOn2DMesh(std::vector<Vec2<T>>&		sample_point, 
								const std::vector<Vec2<T>>&	vertex,
								const std::vector<int3>&	face,
								T							distance,
								const char*					method = NULL){
	//	check
	if( vertex.empty() || face.empty() )
		return 0;
	if( distance < 1e-6 ){
		std::cout << "Error: samplingPointsOn2DMeshStaggered, distance too small: " << distance << std::endl;
		return 0;
	}

	//	clear old data
	sample_point.clear();

	//	create aabb tree
	AABBTree2D<T>	face_aabb;
	face_aabb.BuildTriangleAABBTree(vertex, face);
	AABB2D<T>		bonding_box = face_aabb.node[face_aabb.leaf_number];

	//	get potential sample points
	std::vector<Vec2<T>> potential_sample_point;
	samplingPoissonDisk(potential_sample_point, distance, bonding_box.bb_min, bonding_box.bb_max, method);

	//	add points on the mesh to sampling point list
	for(int i=0; i<potential_sample_point.size(); i++){
		if( isPointOn2DMesh(potential_sample_point[i], vertex, face, face_aabb) ){
			sample_point.push_back(potential_sample_point[i]);
		}
	}

	return sample_point.size();
}

/**
	Poisson Disk Sampling on 2d rectangle

	\param	sample_point	return the sample points generated on rectangle
	\param	distance		minimal distance between any points
	\param	bb_min			min corner of the rectangle
	\param	bb_max			max corner of the rectangle
	\param	method			character string used to choose method of poisson disk sampling, we use author of the method(case insensitive):	\n
							RobertBridson: method in paper Fast Poisson Disk Sampling in Arbitrary Dimensions	\n
	\return	
*/
template<typename T>
int samplingPoissonDisk(std::vector<Vec2<T>>&	sample_point,
						T						distance,
						Vec2<T>					bb_min,
						Vec2<T>					bb_max,
						const char*				method = NULL){
	//	get lower case method string
	std::string method_str;
	if( method ){
		method_str = method;
		std::transform(method_str.begin(), method_str.end(), method_str.begin(), ::tolower);
	}

	//	choose method
	if( method_str == "robertbridson" ){	//	Robert Bridson
		return samplingPoissonDiskRobertBridson(sample_point, distance, bb_min, bb_max);
	}
	else{								//	default method, RobertBridson
		#ifndef	BE_QUIET
		std::cout << "warning: samplingPoissonDisk, unrecognized method: " << method_str << ", use default method: RobertBridson" << std::endl;
		#endif
		return samplingPoissonDiskRobertBridson(sample_point, distance, bb_min, bb_max);
	}
}


/**
	Poisson Disk Sampling on 2d rectangle

	method comes from the paper:	\n
	Fast Poisson Disk Sampling in Arbitrary Dimensions	\n
	SIGGRAPH 2007 sketch, Robert Bridson

	\param	sample_point	return the sample points generated on rectangle
	\param	distance		minimal distance between any points
	\param	bb_min			min corner of the rectangle
	\param	bb_max			max corner of the rectangle
	\return	
*/
template<typename T>
int samplingPoissonDiskRobertBridson(std::vector<Vec2<T>>&	sample_point,
									 T						distance,
									 Vec2<T>				bb_min,
									 Vec2<T>				bb_max){
	const int k = 30;	//	number of samples each point, provided in the paper
	const T square_distance = distance * distance;

	//	check size
	if( distance < 1e-6 ){
		std::cout << "Error: samplingPoissonDiskRobertBridson, distance too small: " << distance << std::endl;
		return 0;
	}
	if( (bb_max.x - bb_min.x) < 1e-6 || (bb_max.y - bb_min.y) < 1e-6 ){
		std::cout << "Error: samplingPoissonDiskRobertBridson, illegal rectangle" << std::endl;
		return 0;
	}

	//	initialize background grid
	T	cell_size = distance / sqrtf(2);
	int dim_x = (bb_max - bb_min).x / cell_size + 1;
	int dim_y = (bb_max - bb_min).y / cell_size + 1;
	std::vector<int> grid;
	grid.resize(dim_x*dim_y, -1);

	//	find first sampling point
	sample_point.clear();
	Vec2<T>	point(bb_min.x + (bb_max.x - bb_min.x)*rand0to1d(), 
		bb_min.y + (bb_max.y - bb_min.y)*rand0to1d());	//	this point must be within the rectangle
	sample_point.push_back(point);	//	save point
	int x = clamp((point.x-bb_min.x)/cell_size, 0, dim_x-1);
	int y = clamp((point.y-bb_min.y)/cell_size, 0, dim_y-1);
	int idx = x + y * dim_x;
	grid[idx] = 0;					//	add to grid
	std::vector<int> active_list;
	active_list.push_back(0);		//	add to active list

	//	loop until the whole space is filled
	while( !active_list.empty() ){
		//	randomly choose a seed point
		int seed_id = rand() % active_list.size();
		int seed_point_id = active_list[seed_id];
		Vec2<T>	seed_point = sample_point[seed_point_id];

		//	generate k points
		bool find_flag = false;
		for(int i=0; i<k; i++){
			//	calculate a random point around the seed point
			T angle	= YZ_PI * 2 * rand0to1d();
			T radius = distance * (1 + rand0to1d());
			Vec2<T>	r0 = Vec2d(1, 0).SetRotateRad(angle);
			point = seed_point + r0 * radius;
			if( !isPointInsideAABB(point, bb_min, bb_max) )		//	point not inside rectangle, just skip
				continue;

			//	check whether potential point is legal
			x = clamp((point.x-bb_min.x)/cell_size, 0, dim_x-1);
			y = clamp((point.y-bb_min.y)/cell_size, 0, dim_y-1);
			idx = x + y * dim_x;
			if( grid[idx] != -1 )	//	this grid is already filled with one point
				continue;

			//	check neighbor grids, we have to check 25 neighbor grids
			bool near_flag = false;
			for(int v=myMax(0, y-2); v<=myMin(dim_y-1, y+2); v++){
				for(int u=myMax(0, x-2); u<=myMin(dim_x-1, x+2); u++){
					int this_idx = u + v * dim_x;
					if( grid[this_idx] != -1 ){	//	find a existing point
						assert( grid[this_idx] < sample_point.size() );
						Vec2<T> this_point = sample_point[grid[this_idx]];
						if( (this_point - point).SquareLength() < square_distance )	//	find a existing point that is close to this point
							near_flag = true;
					}
				}
			}
			if( !near_flag ){	//	all existing points are at least distance away from point
				grid[idx] = sample_point.size();
				active_list.push_back(sample_point.size());
				sample_point.push_back(point);
				find_flag = true;
			}
		}

		//	no point found, remove the point from active list
		if( !find_flag ){
			active_list[seed_id] = active_list.back();
			active_list.pop_back();
		}
	}

	return sample_point.size();
}


///@}

}}}	//	namespace yz::geometry::remeshing

#endif	//	__YZ_SAMPLING_2D_H__