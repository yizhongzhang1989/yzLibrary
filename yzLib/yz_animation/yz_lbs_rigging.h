/***********************************************************/
/**	\file
	\brief		Rigging by Linear Blend Skinning
	\author		Yizhong Zhang
	\date		9/17/2012
*/
/***********************************************************/
#ifndef __YZ_LBS_RIGGING_H__
#define __YZ_LBS_RIGGING_H__

#include <vector>
#include <algorithm>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_geometry/yz_mesh_matrix.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_animation/yz_skeleton.h"
#include "yzLib/yz_animation/yz_skin.h"
#include "yzLib/yz_animation/yz_skeleton_utils.h"
#include "yzLib/yz_utils/yz_color_bar.h"

namespace yz{	namespace animation{

//	========================================
///@{
/**	@name Linear Blend Rigging
*/
//	========================================

//	some functions that should only be used by rigging functions
namespace rigging{
	
/**
	calculate diagomal matrix H and nearest bone id of each vertex

	H is described in Ilya Baran's paper.

	\param	H					return the diagonal of the matrix as vector
	\param	nearest_bone_id		id of the nearest bone of each vertex
	\param	skeleton			the skeleton at rigging pose
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	distance_threshold	minimal tolerance distance of vertex to bone.
								Trancate to this value is distance is smaller.
	\param	heating_coef		coefficient of heating of each bone. typically 1
*/
template<typename T>
void calculateH(std::vector<T>&				H,
				std::vector<int>&			nearest_bone_id,
				const Skeleton<T>&			skeleton,
				const std::vector<Vec3<T>>&	vertex,
				const std::vector<int3>&	face,
				T							distance_threshold,
				T							heating_coef){
	//	prepare space for H and bone
	H.resize(vertex.size());
	nearest_bone_id.resize(vertex.size());

	//	create BVH for the mesh, currently we use AABB Tree
	geometry::AABBTree3D<T>	aabb_tree;
	aabb_tree.BuildTriangleAABBTree(vertex, face);

	//	create VF connectivity,
	geometry::TriMeshVF	vf;
	vf.CreateVF(vertex.size(), face);

	//	calculate start point of each bone
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	check whether line to nearest bone totally inside the bone
#ifdef ENABLE_OPENMP
#	pragma omp parallel for
#endif
	for( int i=0; i<vertex.size(); i++ ){
		//	get the nearest point on the bone
		int		bone_id;
		Vec3<T>	prj_on_bone;
		T		distance;
		getNearestPointOnSkeleton(bone_id, prj_on_bone, distance, skeleton, bone_start, vertex[i]);
		nearest_bone_id[i] = bone_id;

		//	check intersection 
		std::vector<Vec3<T>>	inter_p;
		std::vector<int>		face_id;
		int inter_p_num = geometry::getSegmentMeshIntersectionPoints(inter_p, face_id, 
			vertex[i], prj_on_bone, vertex, face, aabb_tree);

		//	ignore vertices that is vertex[i] itself
		for(int j=0; j<face_id.size(); j++){
			for(int k=vf.vf_start[i]; k<vf.vf_start[i+1]; k++){
				if( vf.vf[k] == face_id[j] ){	//	the intersection point is just vertex[i]
					inter_p_num --;				//	ignore this intersection point
					break;						//	check next intersection face
				}
			}
		}

		//	no intersection means the line segment must be totally inside or outside the mesh
		int inside_flag = 0;
		if( !inter_p_num ){	//	check whether projection point is inside the mesh. If so, the line is inside
			inside_flag = geometry::isPointInsideClosedMesh(prj_on_bone, vertex, face, aabb_tree);
		}

		if( inside_flag ){
			if( distance < distance_threshold )
				distance = distance_threshold;
			H[i] = heating_coef / (distance*distance);
		}
		else{
			H[i] = 0;
		}
	}

}


/**
	calculate diagomal matrix H and nearest bone id of each vertex on open mesh

	H is described in Ilya Baran's paper.

	\param	H					return the diagonal of the matrix as vector
	\param	nearest_bone_id		id of the nearest bone of each vertex
	\param	skeleton			the skeleton at rigging pose
	\param	vertex				vertex list of the mesh
	\param	face				face list of the mesh
	\param	distance_threshold	minimal tolerance distance of vertex to bone.
								Trancate to this value is distance is smaller.
	\param	heating_coef		coefficient of heating of each bone. typically 1
*/
template<typename T>
void calculateHOpenMesh(std::vector<T>&				H,
						std::vector<int>&			nearest_bone_id,
						const Skeleton<T>&			skeleton,
						const std::vector<Vec3<T>>&	vertex,
						const std::vector<int3>&	face,
						T							distance_threshold,
						T							heating_coef){
	//	prepare space for H and bone
	H.resize(vertex.size());
	nearest_bone_id.resize(vertex.size());

	//	calculate vertex normal of the mesh
	std::vector<Vec3<T>>	vertex_normal;
	geometry::calculateVertexNormal(vertex_normal, vertex, face);

	//	create BVH for the mesh, currently we use AABB Tree
	geometry::AABBTree3D<T>	aabb_tree;
	aabb_tree.BuildTriangleAABBTree(vertex, face);

	//	create VF connectivity
	geometry::TriMeshVF	vf;
	vf.CreateVF(vertex.size(), face);

	//	calculate start point of each bone
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	check whether line to nearest bone totally inside the bone
#ifdef ENABLE_OPENMP
#	pragma omp parallel for
#endif
	for( int i=0; i<vertex.size(); i++ ){
		//	get the nearest point on the bone
		int		bone_id;
		Vec3<T>	prj_on_bone;
		T		distance;
		getNearestPointOnSkeleton(bone_id, prj_on_bone, distance, skeleton, bone_start, vertex[i]);
		nearest_bone_id[i] = bone_id;

		//	check intersection 
		std::vector<Vec3<T>>	inter_p;
		std::vector<int>		face_id;
		int inter_p_num = geometry::getSegmentMeshIntersectionPoints(inter_p, face_id, 
			vertex[i], prj_on_bone, vertex, face, aabb_tree);

		//	ignore vertices that is vertex[i] itself
		for(int j=0; j<face_id.size(); j++){
			for(int k=vf.vf_start[i]; k<vf.vf_start[i+1]; k++){
				if( vf.vf[k] == face_id[j] ){	//	the intersection point is just vertex[i]
					inter_p_num --;				//	ignore this intersection point
					break;						//	check next intersection face
				}
			}
		}

		//	no intersection means the line segment must be totally inside or outside the mesh
		int inside_flag = 0;
		T	normal_coef = 1;
		if( !inter_p_num ){	//	check whether projection point is inside the mesh. If so, the line is inside
			Vec3<T> r = prj_on_bone - vertex[i];
			Vec3<T> n = - vertex_normal[i];
			normal_coef = dot(r, n);
			if( normal_coef > 0 )
				inside_flag = 1;
		}

		if( inside_flag ){
			if( distance < distance_threshold )
				distance = distance_threshold;
			H[i] = heating_coef * normal_coef / (distance*distance);
		}
		else{
			H[i] = 0;
		}
	}

}

}	//	namespace yz::animation::rigging

/**
	Rigging by bone heat

	Method from: Automatic Rigging and Animation of 3D Characters \n
	Ilya Baran, Jovan Popovic, SIGGRAPH 2007

	\param	rigging				rigging weight
	\param	skeleton			skeleton at the rest pose
	\param	transformer			rigging pose transformer
	\param	vertex				vertex of the skin mesh
	\param	face				face of the skin mesh
	\param	weight_threshold	the minimal weight we accept, below this, we don't take into account
	\param	c					coefficient c described in the paper, default 1
	\param	open_mesh_flag		if the mesh is open mesh, then set this flag and we will calculate parameters in different ways
	\return						success flag
*/
template<typename T>
int riggingBoneHeat(RiggingLBS<T>&				rigging,
					Skeleton<T>					skeleton,
					SkeletonTransformer<T>		transformer, 
					const std::vector<Vec3<T>>&	vertex,
					const std::vector<int3>&	face,
					T							weight_threshold = 0.001,
					T							c = 1,
					int							open_mesh_flag = 0){
	if( skeleton.bone.size() != transformer.transformer.size() )	//	assure skeleton match transformer
		return 0;

	//	calculate global transformer, transform skeleton to rigging pose,
	//	and calculate start position of each bone at rigging pose
	setTransformerToGlobal(transformer, skeleton);
	transformSkeleton(skeleton, transformer);
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	inverse all rotation matrix in transformer
	for(int i=0; i<transformer.transformer.size(); i++){
		transformer.transformer[i].rotation.SetInverse();
	}

	//	transform as rest pose
	riggingBoneHeat(rigging, skeleton, vertex, face, weight_threshold, c, open_mesh_flag);

	//	rotate offset vector
	for(int i=0; i<rigging.bone_weight.size(); i++){
		int		bone_id = rigging.bone_weight[i].bone_id;
		rigging.bone_weight[i].offset	= transformer.transformer[bone_id].rotation * rigging.bone_weight[i].offset;
	}

	return 1;
}

/**
	Rigging by bone heat

	Method from: Automatic Rigging and Animation of 3D Characters \n
	Ilya Baran, Jovan Popovic, SIGGRAPH 2007

	\param	rigging				rigging weight
	\param	skeleton			skeleton at the rest pose, also rigging pose
	\param	vertex_of_mesh		vertex of the skin mesh
	\param	face				face of the skin mesh
	\param	weight_threshold	the minimal weight we accept, below this, we don't take into account
	\param	c					coefficient c described in the paper, default 1
	\param	open_mesh_flag		if the mesh is open mesh, then set this flag and we will calculate parameters in different ways
	\return						success flag
*/
template<typename T>
int riggingBoneHeat(RiggingLBS<T>&				rigging,
					const Skeleton<T>&			skeleton,
					const std::vector<Vec3<T>>&	vertex_of_mesh,
					const std::vector<int3>&	face,
					T							weight_threshold = 0.001,
					T							c = 1,
					int							open_mesh_flag = 0){

	//	we copy the vertex, and apply noise to each vertex
	std::vector<Vec3<T>> vertex = vertex_of_mesh;
	geometry::applyNoiseToVertices(vertex, 0.0001);

	//	step 1, calculate laplacian matrix
	#ifndef	BE_QUIET
		std::cout << "creating laplacian matrix" << std::endl;
	#endif
	SparseMatrixCSR<T>	lap_mat(vertex.size(), vertex.size());
	geometry::createLaplacianMatrixForMesh(lap_mat, vertex, face);

	//	step 2, calculate diagonal matrix H, we just represent as vector
	#ifndef	BE_QUIET
		std::cout << "creating matrix H" << std::endl;
	#endif
	std::vector<T>		H;
	std::vector<int>	nearest_bone_id;	//	this vector used in step 4
	if( open_mesh_flag )
		rigging::calculateHOpenMesh(H, nearest_bone_id, skeleton, vertex, face, T(0.002), c);
	else
		rigging::calculateH(H, nearest_bone_id, skeleton, vertex, face, T(0.002), c);

	//	step 3, calculate lap_mat + H
	#ifndef	BE_QUIET
		std::cout << "creating matrix A" << std::endl;
	#endif
	for(int i=0; i<lap_mat.row_num; i++){
		for(int k=lap_mat.row_start[i]; k<lap_mat.row_start[i+1]; k++){
			int j = lap_mat.col_id[k];
			if( i == j )	//	diagonal, add H[i]
				lap_mat.value[k] += H[i];
		}
	}

	//	step 4, solve bone heat equation for every bone
#if defined(YZ_eigen_sparse_h)
	EigenDirectSparseSolver<T>	solver;
	solver.SetupMatrix(lap_mat.row_num, lap_mat.col_num, lap_mat.NNZ(),		//	laplacian matrix must be symmetric structure
		(T*)&lap_mat.value[0], (int*)&lap_mat.col_id[0], (int*)&lap_mat.row_start[0], 2);

#elif defined(YZ_mkl_dss_h)	//	mkl dss can be used to solve
	MKLDirectSparseSolver<T>	solver;
	solver.SetupMatrix(lap_mat.row_num, lap_mat.col_num, lap_mat.NNZ(),		//	laplacian matrix must be symmetric structure
		(T*)&lap_mat.value[0], (int*)&lap_mat.col_id[0], (int*)&lap_mat.row_start[0], 2);
#endif

	std::vector<std::pair<int, int>>	v_p_list;	//	vertex index and corresponding position in b_w_list
	std::vector<std::pair<int, T>>		b_w_list;	//	bone weight list
	for(int b_id=0; b_id<skeleton.bone.size(); b_id++){
		#ifndef	BE_QUIET
			std::cout << "\rrigging bone " << b_id+1 << "/" << skeleton.bone.size();
			if( b_id == skeleton.bone.size()-1 )
				std::cout << std::endl;
		#endif
		DenseVector<T>	Hp(vertex.size());
		for(int i=0; i<Hp.Dim(); i++){
			Hp[i] = (nearest_bone_id[i] == b_id) ? H[i] : 0;
		}

		DenseVector<T>	w(vertex.size());
		#if defined(YZ_eigen_sparse_h)
			solver.Solve((T*)&w[0], (T*)&Hp[0], 1);
		#elif defined(YZ_mkl_dss_h)	//	use MKLDirectSparseSolver
			solver.Solve((T*)&w[0], (T*)&Hp[0], 1);
		#else				//	use default solver
			w = Hp / lap_mat;
		#endif

		//	add weight of this bone to potential list
		for(int i=0; i<w.Dim(); i++){
			if(w[i] > weight_threshold ){
				v_p_list.push_back( std::pair<int, int>(i, b_w_list.size()) );
				b_w_list.push_back( std::pair<int, T>(b_id, w[i]) );
			}
		}
	}

	//	copy bone weight from potential list to weight
	std::vector<Vec3<T>> node;
	getSkeletonNodePosition(node, skeleton);

	rigging.Reset();
	rigging.bone_weight.resize( b_w_list.size() );
	rigging.start.resize( vertex.size()+1 );

	std::sort(v_p_list.begin(), v_p_list.end());
	for(int i=0; i<v_p_list.size(); i++){
		int id = v_p_list[i].second;
		int v_id = v_p_list[i].first;
		int b_id = b_w_list[id].first;

		rigging.bone_weight[i].offset	= vertex[v_id] - node[b_id];

		rigging.bone_weight[i].bone_id	= b_id;	//	copy weight
		rigging.bone_weight[i].weight	= b_w_list[id].second;
	}
	rigging.start[0] = 0;
	int end_pos = 0;
	for(int v_id=0; v_id<vertex.size(); v_id++){
		while( end_pos != v_p_list.size() && v_p_list[end_pos].first <= v_id )
			end_pos ++;
		rigging.start[v_id+1] = end_pos;
	}

	return 1;
}


///@}

//	========================================
///@{
/**	@name Rigid Rigging
*/
//	========================================
/**
	Rigging each vertex to its nearest bone on skeleton

	This function binds every vertex in vertex list to the skeleton,
	including stand alone vertices. 

	\param	rigging				rigging information, each vertex is rigged to one bone with weight = 1
	\param	vertex				vertex of the skin mesh
	\param	skeleton			skeleton at the rest pose, it is very important to pass rest pose as parameter
	\param	transformer			rigging pose transformer
	\return						success flag
*/
template<typename T>
int riggingNearest(RiggingLBS<T>&					rigging,
				   const std::vector<Vec3<T>>&		vertex,
				   Skeleton<T>						skeleton,
				   SkeletonTransformer<T>			transformer ){
	if( skeleton.bone.size() != transformer.transformer.size() )	//	assure skeleton match transformer
		return 0;

	//	calculate global transformer, transform skeleton to rigging pose,
	//	and calculate start position of each bone at rigging pose
	setTransformerToGlobal(transformer, skeleton);
	transformSkeleton(skeleton, transformer);
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	inverse all rotation matrix in transformer
	for(int i=0; i<transformer.transformer.size(); i++){
		transformer.transformer[i].rotation.SetInverse();
	}

	//	transform as rest pose
	riggingNearest(rigging, vertex, skeleton);

	//	rotate offset vector
	for(int i=0; i<vertex.size(); i++){
		int		bone_id = rigging.bone_weight[i].bone_id;
		rigging.bone_weight[i].offset	= transformer.transformer[bone_id].rotation * rigging.bone_weight[i].offset;
	}

	return 1;
}

/**
	Rigging each vertex to its nearest bone on skeleton

	This function assumes that rest pose of skeleton is the rigging
	pose. If the skeleton is transformed skeleton, the rigging will
	not be correct

	This function binds every vertex in vertex list to the skeleton,
	including stand alone vertices. 

	\param	rigging				rigging information, each vertex is rigged to one bone with weight = 1
	\param	vertex				vertex of the skin mesh
	\param	skeleton			skeleton at the rest pose, it is very important to pass rest pose as parameter
	\return						success flag
*/
template<typename T>
int riggingNearest(RiggingLBS<T>&					rigging,
				   const std::vector<Vec3<T>>&		vertex,
				   const Skeleton<T>&				skeleton ){

	rigging.bone_weight.resize(vertex.size());
	rigging.start.resize(vertex.size()+1);

	//	calculate start position of each bone at rest pose
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	calculate each vertex binding position
#ifdef ENABLE_OPENMP
#	pragma omp parallel for
#endif
	for(int i=0; i<vertex.size(); i++){
		//	find the nearest bone of the skeleton, bind the vertex to that bone with weight = 1
		int		bone_id;
		Vec3<T>	prj_point;
		T		distance;
		getNearestPointOnSkeleton(bone_id, prj_point, distance, skeleton, bone_start, vertex[i]);

		//	copy information
		rigging.bone_weight[i].bone_id	= bone_id;
		rigging.bone_weight[i].offset	= vertex[i] - bone_start[bone_id];
		rigging.bone_weight[i].weight	= 1;
		rigging.start[i] = i;
	}
	rigging.start[vertex.size()] = vertex.size();

	return 1;
}

/**
	Rigging each vertex to the given bone on skeleton

	This function binds every vertex in vertex list to the skeleton,
	including stand alone vertices. 

	\param	rigging				rigging information, each vertex is rigged to one bone with weight = 1
	\param	vertex				vertex of the skin mesh
	\param	skeleton			skeleton at the rest pose, it is very important to pass rest pose as parameter
	\param	bone_id				the bone index that is rigged to
	\param	transformer			rigging pose transformer
	\return						success flag
*/
template<typename T>
int riggingToGivenBone(RiggingLBS<T>&					rigging,
					   const std::vector<Vec3<T>>&		vertex,
					   Skeleton<T>						skeleton,
					   SkeletonTransformer<T>			transformer,
					   int								bone_id){
	if( skeleton.bone.size() != transformer.transformer.size() )	//	assure skeleton match transformer
		return 0;

	//	calculate global transformer, transform skeleton to rigging pose,
	//	and calculate start position of each bone at rigging pose
	setTransformerToGlobal(transformer, skeleton);
	transformSkeleton(skeleton, transformer);
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	inverse all rotation matrix in transformer
	for(int i=0; i<transformer.transformer.size(); i++){
		transformer.transformer[i].rotation.SetInverse();
	}

	//	transform as rest pose
	riggingToGivenBone(rigging, vertex, skeleton, bone_id);

	//	rotate offset vector
	for(int i=0; i<vertex.size(); i++){
		int		bone_id = rigging.bone_weight[i].bone_id;
		rigging.bone_weight[i].offset	= transformer.transformer[bone_id].rotation * rigging.bone_weight[i].offset;
	}

	return 1;
}

/**
	Rigging each vertex to the given bone on skeleton

	This function assumes that rest pose of skeleton is the rigging
	pose. If the skeleton is transformed skeleton, the rigging will
	not be correct

	This function binds every vertex in vertex list to the skeleton,
	including stand alone vertices. 

	\param	rigging				rigging information, each vertex is rigged to one bone with weight = 1
	\param	vertex				vertex of the skin mesh
	\param	skeleton			skeleton at the rest pose, it is very important to pass rest pose as parameter
	\param	bone_id				the bone index that is rigged to
	\return						success flag
*/
template<typename T>
int riggingToGivenBone(RiggingLBS<T>&					rigging,
					   const std::vector<Vec3<T>>&		vertex,
					   const Skeleton<T>&				skeleton,
					   int								bone_id){

	rigging.bone_weight.resize(vertex.size());
	rigging.start.resize(vertex.size()+1);

	//	calculate start position of each bone at rest pose
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	calculate each vertex binding position
	for(int i=0; i<vertex.size(); i++){
		rigging.bone_weight[i].bone_id	= bone_id;
		rigging.bone_weight[i].offset	= vertex[i] - bone_start[bone_id];
		rigging.bone_weight[i].weight	= 1;
		rigging.start[i] = i;
	}
	rigging.start[vertex.size()] = vertex.size();

	return 1;
}

///@}

//	========================================
///@{
/**	@name Update Rigging
*/
//	========================================

/**
	update rigging offset

	This function is used to calculate offset of each vertex to
	the bone that it has rigged to. Because rigging is only meaningful
	at rest pose of skeleton, so we require the input skeleton to
	be at rest pose.

	This function is useful especially when we have a new rest pose,
	and we want to use weight of the old rigging. So what we need is to
	update offset only.
	
	\param	rigging			return the rigging, offset will be updated according to vertex and skeleton
	\param	vertex			the vertex
	\param	rest_skeleton	skeleton at rest pose
	\return					1: succeed, 0: failed
*/
template<typename T>
int updateRiggingOffset(RiggingLBS<T>&				rigging,
						const std::vector<Vec3<T>>&	vertex,
						const Skeleton<T>&			rest_skeleton){
	if( rigging.start.size()-1 != vertex.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: updateRiggingOffset, rigging vertex size don't match" << std::endl;
		#endif
		return 0;
	}

	//	calculate start position of each bone at rest pose
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, rest_skeleton);

	//	update offset
	for(int i=0; i<vertex.size(); i++){
		for(int j=rigging.start[i]; j<rigging.start[i+1]; j++){
			int bone_id = rigging.bone_weight[j].bone_id;
			rigging.bone_weight[j].offset = vertex[i] - bone_start[bone_id];
		}
	}

	return 1;
}


///@}

//	========================================
///@{
/**	@name Get Rigging and Color
*/
//	========================================
/**
	get weight of each vertex to a single bone

	number of vertex can be achieved in RiggingLBS

	\param	weight		return weight of each vertex
	\param	rigging		LBS rigging information
	\param	bone_id		the id of the bone, illegal bone_id will result in returning all zero
	\return				whether succeed
*/
template<typename T>
int getWeightOfBoneFromRigging(std::vector<T>&		weight,
							   const RiggingLBS<T>&	rigging,
							   int					bone_id){
	int vertex_number = rigging.start.size() - 1;
	weight.resize( vertex_number );

	//	scan rigging, extract weight of the given bone
	for(int i=0; i<vertex_number; i++){
		weight[i] = 0;	//	if the vertex is not attached to this bone, its weight is 0
		for(int j=rigging.start[i]; j<rigging.start[i+1]; j++){
			if( bone_id == rigging.bone_weight[j].bone_id ){
				weight[i] = rigging.bone_weight[j].weight;
				break;
			}
		}
	}

	return 1;
}

/**
	get weight color of each vertex to a single bone

	number of vertex can be achieved in RiggingLBS

	\param	weight_color	return weight color of each vertex
	\param	rigging			LBS rigging information
	\param	bone_id			the id of the bone, illegal bone_id will result in returning all zero
	\param	color_func		function of color calculation. refer to yz_color_bar.h
							Default: convertToColorJet
	\return					whether succeed
*/
template<typename C, typename T>
int getWeightColorOfBoneFromRigging(std::vector<yz::Com3<C>>&	weight_color,
									const RiggingLBS<T>&		rigging,
									int							bone_id,
									void (*color_func) (C*, float, float, float, float) = utils::convertToColorJet<C> ){
	int vertex_number = rigging.start.size() - 1;
	weight_color.resize( vertex_number );

	//	scan rigging, extract weight of the given bone
	for(int i=0; i<vertex_number; i++){
		T weight = 0;	//	if the vertex is not attached to this bone, its weight is 0
		for(int j=rigging.start[i]; j<rigging.start[i+1]; j++){
			if( bone_id == rigging.bone_weight[j].bone_id ){
				weight = rigging.bone_weight[j].weight;
				break;
			}
		}

		//	calculate color
		color_func((C*)&weight_color[i], weight, 0, 1, 0);
	}

	return 1;
}


///@}

}}	//	namespace yz::animation

#endif	//	__YZ_LBS_RIGGING_H__