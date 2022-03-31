/***********************************************************/
/**	\file
	\brief		Mesh Optimization
	\author		Yizhong Zhang
	\date		12/27/2017
*/
/***********************************************************/
#ifndef __YZ_MESH_OPTIMIZATION_H__
#define __YZ_MESH_OPTIMIZATION_H__

#include <iostream>
#include <vector>
#include <functional>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "yzLib/yz_setting.h"


namespace yz{  namespace geometry{  namespace remeshing{

//	========================================
///@{
/**	@name Mesh Optimization
*/
//	========================================


/**
	Perform mesh optimization by edge collapse, edge swap, edge split. 

	We also provide smoothing functions. Non-manifold mesh is also appliable. 
	The refinement doesn't change structure of the mesh. For example, we don't seperate mesh that form a thin tube.

	This class provide a variaty of interfaces for mesh optimization. Interface functions can be called 
	in a sequence, or called individually. Interfaces can be used in the following methods:	\n
	----- Method 1:	\n
	1, call SetMesh(v, f), or construct the class with (v, f), the class will point to the mesh	\n
	2, call Refine(collapse_thre, split_thre, 1) to perform edge collapse, edge swap, edge split one step	\n
	3, call SmoothLaplacian(coef) or SmoothTaubin(iterations, coef_forward, coef_backward) to smooth the mesh	\n
	4, perform step 2, step 3 iteratively, until the mesh quality is good enough	\n
	----- Method 2:	\n
	1, call SetMesh(v, f), or construct the class with (v, f), the class will point to the mesh	\n
	2, call Refine(collapse_thre, T split_thre, 10000) to perform edge collapse, edge swap, edge split until no more edge will change	\n
	3, call SmoothLaplacian(coef) or SmoothTaubin(iterations, coef_forward, coef_backward) to smooth the mesh	\n

	You can also call EdgeCollapse, EdgeSwap or EdgeSplit to perform individual operation. 

	This class include a virtual function: ReIndexMesh(). This function performs post processing after edge operations. 
	If you have more data associated with vertex and face, you can overload this function, so that the mesh refinement
	can update these data as well. We provide ReIndexFaceData(std::vector<DATA>& f_data) of handle face data.
*/
template <class T>
class TriMeshRefiner {
public:
	std::vector<Vec3<T>>*	vertex_ptr;		///<	pointer to vertex array
	std::vector<int3>*		face_ptr;		///<	pointer to face array

public:
	TriMeshRefiner() {
		vertex_ptr = NULL;
		face_ptr = NULL;
	}

	TriMeshRefiner(std::vector<Vec3<T>>& v, std::vector<int3>& f) {
		vertex_ptr = &v;
		face_ptr = &f;

		SetMesh(v, f);
	}

	/**
		Setup mesh to refine

		\param	v	mesh vertex 
		\param	f	mesh face
	*/
	void SetMesh(std::vector<Vec3<T>>& v, std::vector<int3>& f) {
		vertex_ptr = &v;
		face_ptr = &f;

		ClearDataArray();
		SetupVVF();
	}

	/**
		Refine the mesh by performing edge collapse, edge swap, edge split iteratively

		\param	collapse_thre	min edge length
		\param	split_thre		max edge length
		\param	max_iterations	max edge collapse, swap, split iterations. By default, just one iteration
		\return					(uncollapsed_edge_count, unswapped_edge_count, unsplitted_edge_count)
	*/
	int3 Refine(T collapse_thre, T split_thre, int max_iterations = 1) {
		if (!vertex_ptr || !face_ptr) {
			std::cout << "error: TriMeshRefiner::Refine, mesh not set, cannot refine" << std::endl;
			return int3(0, 0, 0);
		}

		Vec3i remain_count(0, 0, 0);

		//	perform optimization iteratively, until no more edges will change, or reached max iterations
		for (int i = 0; i < max_iterations; i++) {
			CreateEdgeToOptimize(collapse_thre*collapse_thre, split_thre*split_thre);

			remain_count.x = SolveCollapse();
			remain_count.y = SolveSwap();
			remain_count.z = SolveSplit();

			ReIndexMesh();

			ClearDataArray();
			SetupVVF();

			//	if no strange edges exist, refine can stop
			if (remain_count == Vec3i(0, 0, 0))
				break;
			//	if this iteration cannot improve any more, refine can stop
			if (edge_to_collapse.size() == remain_count.x
				&& edge_to_swap.size() == remain_count.y
				&& edge_to_split.size() == remain_count.z)
				break;
		}

		return remain_count;
	}

	/**
		Perform edge collapse, until no more edges can collapse

		\param	collapse_thre	min edge length
		\param	max_iterations	max edge collapse iterations. By default, just one iteration
		\return					number of edges uncollapsed
	*/
	int EdgeCollapse(T collapse_thre, int max_iterations = 1) {
		if (!vertex_ptr || !face_ptr)
			return 0;

		int remain_count = 0;

		//	perform optimization iteratively, until no more edges will change, or reached max iterations
		for (int i = 0; i < max_iterations; i++) {
			CreateEdgeToOptimize(collapse_thre*collapse_thre, 1e12);

			remain_count = SolveCollapse();

			ReIndexMesh();

			ClearDataArray();
			SetupVVF();

			//	if no strange edges exist, refine can stop
			if (remain_count == 0)
				break;
			//	if this iteration cannot improve any more, refine can stop
			if (edge_to_collapse.size() == remain_count)
				break;
		}

		return remain_count;
	}

	/**
		Perform edge swap, until no more edges can swap any more

		Swap the edge if the edge angle sum is bigger than 180 degrees

		\param	max_iterations	max edge swap iterations. By default, just one iteration
		\return					number of edges unswapped
	*/
	int EdgeSwap(int max_iterations = 1) {
		if (!vertex_ptr || !face_ptr)
			return 0;

		int remain_count = 0;

		for (int i = 0; i < max_iterations; i++) {
			CreateEdgeToOptimize();

			remain_count = SolveSwap();

			ReIndexMesh();

			ClearDataArray();
			SetupVVF();

			//	if no strange edges exist, refine can stop
			if (remain_count == 0)
				break;
			//	if this iteration cannot improve any more, refine can stop
			if (edge_to_swap.size() == remain_count)
				break;
		}

		return remain_count;
	}

	/**
		Perform edge split, until no more edges can split

		\param	split_thre		max edge length
		\param	max_iterations	max edge split iterations. By default, just one iteration
		\return					number of edges split
	*/
	int EdgeSplit(T split_thre, int max_iterations = 1) {
		if (!vertex_ptr || !face_ptr)
			return -1;

		int remain_count = 0;

		for (int i = 0; i < max_iterations; i++) {
			CreateEdgeToOptimize(0, split_thre*split_thre);

			remain_count = SolveSplit();

			ReIndexMesh();

			ClearDataArray();
			SetupVVF();

			//	if no strange edges exist, refine can stop
			if (remain_count == 0)
				break;
			//	if this iteration cannot improve any more, refine can stop
			if (edge_to_split.size() == remain_count)
				break;
		}

		return remain_count;
	}

	/**
		Perform Laplacian smooth on the mesh

		\param	coef			linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
	*/
	void SmoothLaplacian(T coef = 1.0) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);

		laplacian_vertex.resize(vertex.size());

		//	perform laplacian smooth
		for (unsigned int vid = 0; vid != vertex.size(); vid++) {
			//	for each vertex, calculate the average vertex position of its neighbor
			auto its = vvf_hash.equal_range(vid);
			Vec3<T> v_sum(0, 0, 0);
			unsigned int count = 0;
			for (auto iter = its.first; iter != its.second; iter++) {
				v_sum += vertex[iter->second.vid];
				count++;
			}
			if (count) {	//	normal vertex, interpolate between original and Laplacian
				v_sum /= count;
				laplacian_vertex[vid] = (1.0 - coef) * vertex[vid] + coef * v_sum;
			}
			else {			//	isolated vertex, just original 
				laplacian_vertex[vid] = vertex[vid];
			}
		}

		vertex.swap(laplacian_vertex);
	}

	/**
		Perform Taubin smooth on the mesh

		Taubin smooth is to laplacian forward, then backward, so that the volume is almost preserved.
		suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

		\param	iterations			number of smoothing iterations
		\param	coef_forward		forward coefficient of Taubin smooth
		\param	coef_backward		backward coefficient of Taubin smooth
	*/
	void SmoothTaubin(
		int		iterations		= 10,
		T		coef_forward	= 0.4507499669,
		T		coef_backward	= -0.4720265626
	) {
		for (int i = 0; i < iterations; i++) {
			SmoothLaplacian(coef_forward);
			SmoothLaplacian(coef_backward);
		}
	}

protected:
	/**
		neighbor edge of a vertex
	*/
	struct VVF {
		int		vid;	///<	neighbor vertex index
		int		fid;	///<	neighbor face index
		VVF(int nv, int nf) : vid(nv), fid(nf) {}
	};

	//	mesh topology
	std::unordered_multimap<int, VVF>		vvf_hash;				///< vertex-vertex-face hash table, (vertex index, neighbor vertex, neighbor face)

	//	refinement parameters. Edges are represented by Vec2i, x < y
	typedef std::unordered_set<Vec2i, utils::BitwiseHasher<Vec2i>>	UNORDERED_SET_VEC2I;
	UNORDERED_SET_VEC2I						locked_edge_set;		///< locked edges, cannot perform any operation
	std::vector<std::pair<T, Vec2i>>		edge_to_collapse;		///< all edges shorter than threshold, recorded as (edge length, edge)
	std::vector<std::pair<T, Vec2i>>		edge_to_swap;			///< all edges to swap, recorded as (-angle sum, edge). add negtive to perform std::sort
	std::vector<std::pair<T, Vec2i>>		edge_to_split;			///< all edges longer than threshold, recorded as (-edge length, edge). add negtive to perform std::sort

	//	post processing parameters
	std::vector<int2>						vertex_to_merge_vector;	///< all vertices that merge together. int2(vertex to remove, merge target)
	std::unordered_map<int, int>			vertex_to_merge_map;	///< all vertices that merge together. std::pair(vertex to remove, merge target)
	std::vector<Vec3<T>>					vertex_to_add_vector;	///< vertices that need to add to the back of vertex array, should not change the order
	std::vector<int>						face_to_remove_vector;	///< all faces that to be removed
	std::vector<int3>						face_to_update_vector;	///< all faces that update vertex index. int3(face idx, original vertex idx, target vertex idx)
	std::vector<int4>						face_to_add_vector;		///< faces that need to add to the back of face array, should not change the order. (face x, y, z, original face idx)

	//	temporary data
	std::vector<Vec3<T>>					laplacian_vertex;		///< temporary vertex list

protected:
	/**
		ReIndex the mesh

		We make this function virtual, so that Refine() will call overload function automatically

		If overload class have original data associated with face, ReIndexFaceData() can be called after ReIndexMesh().
	*/
	virtual void ReIndexMesh() {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		//	add new detected vertices and faces
		for (unsigned int i = 0; i != vertex_to_add_vector.size(); i++)
			vertex.push_back(vertex_to_add_vector[i]);
		for (unsigned int i = 0; i != face_to_add_vector.size(); i++)
			face.push_back(int3(face_to_add_vector[i].x, face_to_add_vector[i].y, face_to_add_vector[i].z));

		//	if a vertex will be removed, update all faces around it
		for (unsigned int i = 0; i != vertex_to_merge_vector.size(); i++) {
			int vid_source = vertex_to_merge_vector[i].x;
			int vid_target = vertex_to_merge_vector[i].y;
			auto its = vvf_hash.equal_range(vid_source);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int fid = iter->second.fid;
				if (face[fid].x == vid_source)		face[fid].x = vid_target;
				else if (face[fid].y == vid_source)	face[fid].y = vid_target;
				else if (face[fid].z == vid_source)	face[fid].z = vid_target;
			}
		}

		//	if a face is to be updated, do accordingly
		for (unsigned int i = 0; i != face_to_update_vector.size(); i++) {
			int fid = face_to_update_vector[i].x;
			for (int j = 0; j < 3; j++) {
				if (face[fid][j] == face_to_update_vector[i].y) {
					face[fid][j] = face_to_update_vector[i].z;
					break;
				}
			}
			//	because when creating face_to_update_vector, we don't take vertex_to_merge_vector into account, so we check it explicitly
			for (int j = 0; j < 3; j++) {
				std::unordered_map<int, int>::iterator iter = vertex_to_merge_map.find(face[fid][j]);
				if (iter != vertex_to_merge_map.end())
					face[fid][j] = iter->second;
			}
		}

		//	remove faces that are labeled to be removed
		std::sort(face_to_remove_vector.begin(), face_to_remove_vector.end(), std::greater<int>());
		for (unsigned int i = 0; i != face_to_remove_vector.size(); i++) {
			face[face_to_remove_vector[i]] = face.back();
			face.pop_back();
		}

		//	it is possible that isolated vertices exist
	}

	/**
		re-index data associated with face, used for overload functions of ReIndexMesh

		This function should be called after ReIndexMesh
	*/
	template <typename DATA>
	void ReIndexFaceData(std::vector<DATA>& f_data) {
		for (unsigned int i = 0; i != face_to_add_vector.size(); i++) {
			f_data.push_back(f_data[face_to_add_vector[i].w]);
		}
		for (unsigned int i = 0; i != face_to_remove_vector.size(); i++) {
			f_data[face_to_remove_vector[i]] = f_data.back();
			f_data.pop_back();
		}
	}

	/**
		Clear all data array

		call this function every time before performing edge operations
	*/
	void ClearDataArray() {
		vvf_hash.clear();
		if (face_ptr)
			vvf_hash.reserve((*face_ptr).size() * 6);

		locked_edge_set.clear();
		edge_to_collapse.clear();
		edge_to_swap.clear();
		edge_to_split.clear();

		vertex_to_merge_vector.clear();
		vertex_to_merge_map.clear();
		vertex_to_add_vector.clear();
		face_to_remove_vector.clear();
		face_to_update_vector.clear();
		face_to_add_vector.clear();
	}

	/**
		Setup the neighborhood hash table: vvf_hash, for fast VVF lookup
	*/
	void SetupVVF() {
		if (!vertex_ptr || !face_ptr)
			return;

		std::vector<int3>& face = (*face_ptr);

		//	prepare data array
		vvf_hash.clear();
		vvf_hash.reserve(face.size() * 6);

		//	create the VVF hash table for fast lookup
		for (unsigned int fid = 0; fid != face.size(); fid++) {
			vvf_hash.insert(std::pair<int, VVF>(face[fid].x, VVF(face[fid].y, fid)));
			vvf_hash.insert(std::pair<int, VVF>(face[fid].x, VVF(face[fid].z, fid)));
			vvf_hash.insert(std::pair<int, VVF>(face[fid].y, VVF(face[fid].x, fid)));
			vvf_hash.insert(std::pair<int, VVF>(face[fid].y, VVF(face[fid].z, fid)));
			vvf_hash.insert(std::pair<int, VVF>(face[fid].z, VVF(face[fid].x, fid)));
			vvf_hash.insert(std::pair<int, VVF>(face[fid].z, VVF(face[fid].y, fid)));
		}
	}

	/**
		Create Edges to Collapse, Swap and Split

		\param	collapse_thre_square	collapse_thre * collapse_thre
		\param	split_thre_square		split_thre * split_thre
	*/
	virtual void CreateEdgeToOptimize(
		T	collapse_thre_square	= 0.0,
		T	split_thre_square		= 1e12
	) {
		if (!vertex_ptr || !face_ptr)
			return;

		std::vector<Vec3<T>>&	vertex	= (*vertex_ptr);
		std::vector<int3>&		face	= (*face_ptr);

		//	we access each edge only once, and record accessed edge with locked_edge_set
		locked_edge_set.clear();
		edge_to_collapse.clear();
		edge_to_swap.clear();
		edge_to_split.clear();

		for (unsigned int fid = 0; fid != face.size(); fid++) {
			for (int i = 0; i < 3; i++) {
				//	get the edge that is not accessed yet
				int2 edge(face[fid][i], face[fid][(i + 1) % 3]);
				if (edge.x > edge.y)
					mySwap(edge.x, edge.y);
				if (locked_edge_set.find(edge) != locked_edge_set.end())
					continue;
				locked_edge_set.insert(edge);

				//	calculate edge length, and check whether it should be collapsed or split
				T squ_len = (vertex[edge.x] - vertex[edge.y]).SquareLength();
				if (squ_len < collapse_thre_square)
					edge_to_collapse.push_back(std::pair<T, Vec2i>(squ_len, edge));
				else if (squ_len > split_thre_square)
					edge_to_split.push_back(std::pair<T, Vec2i>(-squ_len, edge));	//	we use -squ_len to perform sort

				//	calculate angle sum of this edge
				int co_edge_face_count = 0;
				T angle_sum = 0;
				auto its = vvf_hash.equal_range(edge.x);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					if (iter->second.vid != edge.y)
						continue;

					int fid = iter->second.fid;
					int third_vid = getThirdVertexOfFace(face[fid], edge.x, edge.y);

					//	count number of faces sharing this edge
					co_edge_face_count++;

					//	calculate angle sum of this edge
					angle_sum += angleRadBetweenVectors(
						vertex[edge.x] - vertex[third_vid],
						vertex[edge.y] - vertex[third_vid]
					);
				}
				//	only normal edges (nor boundary or non-manifold edge) can be swapped
				if (co_edge_face_count == 2 && angle_sum > 3.16) {
					edge_to_swap.push_back(std::pair<T, Vec2i>(-angle_sum, edge));	//	we use -angle_sum to perform sort
				}
			}
		}

		locked_edge_set.clear();
	}

	/**
		Collapse edges from edge_to_collapse list, from shortest edge to longest

		We only collapse edges that are collapsable

		\return		number of short edges remain uncollapsed
	*/
	virtual int SolveCollapse() {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		int collapsed_edge_count = 0;
		std::vector<int2> edge_to_lock;
		std::vector<int> neighbor_face;

		//	collapse from shortest to longest
		std::sort(edge_to_collapse.begin(), edge_to_collapse.end());

		for (unsigned int i = 0; i != edge_to_collapse.size(); i++) {
			int2 edge = edge_to_collapse[i].second;

			//	check whether able to collapse
			if (locked_edge_set.find(edge) != locked_edge_set.end())		//	if the edge is locked, cannot collapse
				continue;
			if (!CollapseLegalCheck(edge, neighbor_face, edge_to_lock))		//	if the edge will cause non-manifold issue due to collapse, cannot collapse
				continue;

			//	check collapse direction, if one of the vertex is on non-manifold boundary, collapse the other one
			int nonmanifold_x_flag = 0, nonmanifold_y_flag = 0;
			auto its = vvf_hash.equal_range(edge.x);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int2 nei_edge = int2(edge.x, iter->second.vid);
				std::vector<int> nei_edge_nei_face;
				GetNeighborFacesOfEdge(nei_edge_nei_face, nei_edge);
				if (nei_edge_nei_face.size() != 2) {
					nonmanifold_x_flag = 1;
					break;
				}
			}
			its = vvf_hash.equal_range(edge.y);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int2 nei_edge = int2(edge.y, iter->second.vid);
				std::vector<int> nei_edge_nei_face;
				GetNeighborFacesOfEdge(nei_edge_nei_face, nei_edge);
				if (nei_edge_nei_face.size() != 2) {
					nonmanifold_y_flag = 1;
					break;
				}
			}
			if (!nonmanifold_x_flag && nonmanifold_y_flag) {		//	x is normal, y is non-manifold, we collapse x
			}
			else if (nonmanifold_x_flag && !nonmanifold_y_flag) {	//	x is non-manifold, y is normal, we collapse y
				mySwap(edge.x, edge.y);	
			}
			else {	//	other cases, we collapse the vertex that cause less volume change
				T swap_volume_x = 0, swap_volume_y = 0;
				auto its = vvf_hash.equal_range(edge.x);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					int fid = iter->second.fid;
					if (face[fid][0] == edge.y || face[fid][1] == edge.y || face[fid][2] == edge.y)
						continue;
					Vec3<T> r1 = vertex[edge.y] - vertex[face[fid][0]];
					Vec3<T> r2 = vertex[face[fid][1]] - vertex[face[fid][0]];
					Vec3<T> r3 = vertex[face[fid][2]] - vertex[face[fid][0]];
					swap_volume_x += fabs(dot(r1, cross(r2, r3)));
				}
				its = vvf_hash.equal_range(edge.y);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					int fid = iter->second.fid;
					if (face[fid][0] == edge.x || face[fid][1] == edge.x || face[fid][2] == edge.x)
						continue;
					Vec3<T> r1 = vertex[edge.x] - vertex[face[fid][0]];
					Vec3<T> r2 = vertex[face[fid][1]] - vertex[face[fid][0]];
					Vec3<T> r3 = vertex[face[fid][2]] - vertex[face[fid][0]];
					swap_volume_y += fabs(dot(r1, cross(r2, r3)));
				}
				if (swap_volume_x > swap_volume_y)
					mySwap(edge.x, edge.y);
			}
			
			//	collapse this edge
			vertex_to_merge_vector.push_back(edge);
			vertex_to_merge_map.insert(std::pair<int, int>(edge.x, edge.y));

			for (unsigned int j = 0; j != neighbor_face.size(); j++)
				face_to_remove_vector.push_back(neighbor_face[j]);

			//	lock neighboring edges. Don't touch these edges until the vvf_hash is updated
			for (unsigned int j = 0; j != edge_to_lock.size(); j++)
				LockEdge(edge_to_lock[j]);

			collapsed_edge_count++;
		}

		return edge_to_collapse.size() - collapsed_edge_count;
	}

	/**
		Swap edges from edge_to_swap list, from most sharp edges

		We only swap edges that are swappable

		\return		number of sharp edges remain unswapped
	*/
	virtual int SolveSwap() {
		int swap_edge_count = 0;
		std::vector<int2> edge_to_lock;
		std::vector<int2> neighbor_face_vertex;

		std::sort(edge_to_swap.begin(), edge_to_swap.end());

		for (unsigned int i = 0; i != edge_to_swap.size(); i++) {
			int2 edge = edge_to_swap[i].second;

			//	check whether able to swap
			if (locked_edge_set.find(edge) != locked_edge_set.end())		//	if the edge is locked, cannot collapse
				continue;
			if (!SwapLegalCheck(edge, neighbor_face_vertex, edge_to_lock))	//	if the edge will cause non-manifold issue due to swap, cannot swap
				continue;

			//	swap this edge
			face_to_update_vector.push_back(int3(neighbor_face_vertex[0].x, edge.x, neighbor_face_vertex[1].y));
			face_to_update_vector.push_back(int3(neighbor_face_vertex[1].x, edge.y, neighbor_face_vertex[0].y));

			//	lock edge
			for (unsigned int j = 0; j != edge_to_lock.size(); j++)
				LockEdge(edge_to_lock[j]);

			swap_edge_count++;
		}

		return edge_to_swap.size() - swap_edge_count;
	}

	/**
		Collapse edges from edge_to_split list, from longest edge to shortest

		We only split edges that are splitable

		\return		number of long edges remain unsplitted
	*/
	virtual int SolveSplit() {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		int split_edge_count = 0;

		std::sort(edge_to_split.begin(), edge_to_split.end());

		for (unsigned int i = 0; i != edge_to_split.size(); i++) {
			int2 edge = edge_to_split[i].second;

			//	check whether able to split. Edge split doesn't change topology, so we don't perform topology check
			if (locked_edge_set.find(edge) != locked_edge_set.end())		//	if the edge is locked, cannot collapse
				continue;

			//	insert new vertex
			int new_vid = vertex.size() + vertex_to_add_vector.size();
			Vec3<T> new_v = (vertex[edge.x] + vertex[edge.y]) * 0.5;
			vertex_to_add_vector.push_back(new_v);

			//	this edge doesn't exist any more, lock it
			LockEdge(edge);

			//	split all faces sharing this edge
			auto its = vvf_hash.equal_range(edge.x);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				if (iter->second.vid == edge.y) {
					//	split the face, by change one to two
					int fid = iter->second.fid;
					int vid = geometry::getThirdVertexOfFace(face[fid], edge.x, edge.y);

					face_to_update_vector.push_back(int3(fid, edge.x, new_vid));

					int new_fid = face.size() + face_to_add_vector.size();
					face_to_add_vector.push_back(int4(face[fid], fid));
					face_to_update_vector.push_back(int3(new_fid, edge.y, new_vid));

					LockEdge(int2(vid, edge.x));
					LockEdge(int2(vid, edge.y));
				}
			}

			split_edge_count++;
		}

		return edge_to_split.size() - split_edge_count;
	}

	/**
		If an edge is locked, it cannot be collapsed, swapped, or split
	*/
	inline void LockEdge(int2 edge) {
		if (edge.x > edge.y)
			mySwap(edge.x, edge.y);
		locked_edge_set.insert(edge);
	}

	/**
		Check whether an edge is legal collapse. This collapse cannot cause topology issues

		If a thin tube is formed, this edge cannot collapse

		\param	edge			the edge to check
		\param	neighbor_face	all faces that share the edge as boundary
		\param	edge_to_lock	if the edge is collapsed, all edges that need to lock
		\return					0: illegal;  1: legal
	*/
	inline int CollapseLegalCheck(
		int2				edge, 
		std::vector<int>&	neighbor_face, 
		std::vector<int2>&	edge_to_lock
	) {
		neighbor_face.clear();
		edge_to_lock.clear();

		//	find all vertices that form a face with the edge
		std::vector<int>	co_face_vertex;
		auto its = vvf_hash.equal_range(edge.x);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			if (iter->second.vid == edge.y) {
				int fid = iter->second.fid;
				int vid = getThirdVertexOfFace((*face_ptr)[fid], edge.x, edge.y);
				neighbor_face.push_back(fid);
				co_face_vertex.push_back(vid);
			}
		}

		//	If co_face_vertex are connected, it forms a tetrahedron (as a single tet, or a tet on the surface), cannot collapse
		for (unsigned int i = 0; i != neighbor_face.size(); i++) {
			std::vector<int>	co_face_vertex_x, co_face_vertex_y;
			GetCoFaceVerticesOfEdge(co_face_vertex_x, int2(co_face_vertex[i], edge.x));
			GetCoFaceVerticesOfEdge(co_face_vertex_y, int2(co_face_vertex[i], edge.y));
			//	check whether co_face_vertex_x, co_face_vertex_y has overlap
			for (unsigned int x = 0; x != co_face_vertex_x.size(); x++)
				for (unsigned int y = 0; y != co_face_vertex_y.size(); y++)
					if (co_face_vertex_x[x] == co_face_vertex_y[y])
						return 0;
		}

		//	find all vertices neighboring edge.x, exclude co_face_vertex detected previously
		std::vector<int>	one_ring_vertex_x;
		its = vvf_hash.equal_range(edge.x);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			int one_ring_vid = iter->second.vid;
			if (one_ring_vid == edge.y)
				continue;
			if (std::find(co_face_vertex.begin(), co_face_vertex.end(), one_ring_vid) != co_face_vertex.end())
				continue;
			if (std::find(one_ring_vertex_x.begin(), one_ring_vertex_x.end(), one_ring_vid) == one_ring_vertex_x.end())
				one_ring_vertex_x.push_back(one_ring_vid);
		}

		//	find all vertices neighboring edge.y, exclude co_face_vertex detected previously
		std::vector<int>	one_ring_vertex_y;
		its = vvf_hash.equal_range(edge.y);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			int one_ring_vid = iter->second.vid;
			if (one_ring_vid == edge.x)
				continue;
			if (std::find(co_face_vertex.begin(), co_face_vertex.end(), one_ring_vid) != co_face_vertex.end())
				continue;
			if (std::find(one_ring_vertex_y.begin(), one_ring_vertex_y.end(), one_ring_vid) == one_ring_vertex_y.end())
				one_ring_vertex_y.push_back(one_ring_vid);
		}

		//	if the mesh forms a thin tube, cannot collapse
		for (unsigned int i = 0; i != one_ring_vertex_x.size(); i++) {
			for (unsigned int j = 0; j != one_ring_vertex_y.size(); j++) {
				if (one_ring_vertex_x[i] == one_ring_vertex_y[j])
					return 0;
				//	if this edge is collapsed, we must lock edges that may cause thin tube issues
				//	because edge will not collapse until the end of the refinement, when detecting
				//	those edges, this edge is not collapsed yet and cannot be detected
				if (IsEdgeExist(one_ring_vertex_x[i], one_ring_vertex_y[j])) {
					edge_to_lock.push_back(int2(one_ring_vertex_x[i], one_ring_vertex_y[j]));
				}
			}
		}

		//	lock the whole neighborhood of edge
		its = vvf_hash.equal_range(edge.x);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			int fid = iter->second.fid;
			edge_to_lock.push_back(int2((*face_ptr)[fid].x, (*face_ptr)[fid].y));
			edge_to_lock.push_back(int2((*face_ptr)[fid].x, (*face_ptr)[fid].z));
			edge_to_lock.push_back(int2((*face_ptr)[fid].y, (*face_ptr)[fid].z));
		}
		its = vvf_hash.equal_range(edge.y);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			int fid = iter->second.fid;
			edge_to_lock.push_back(int2((*face_ptr)[fid].x, (*face_ptr)[fid].y));
			edge_to_lock.push_back(int2((*face_ptr)[fid].x, (*face_ptr)[fid].z));
			edge_to_lock.push_back(int2((*face_ptr)[fid].y, (*face_ptr)[fid].z));
		}

		return 1;
	}

	/**
		Check whether an edge is legal swap. This swap cannot cause topology issues

		\param	edge					the edge to check
		\param	neighbor_face_vertex	all faces that share the edge as boundary, and their corresponding vertices
		\param	edge_to_lock			if the edge is swapped, all edges that need to lock
		\return							0: illegal;  1: legal
	*/
	inline int SwapLegalCheck(
		int2				edge, 
		std::vector<int2>&	neighbor_face_vertex,	
		std::vector<int2>&	edge_to_lock
	) {
		edge_to_lock.clear();
		neighbor_face_vertex.clear();

		//	create neighbor_face_vertex
		auto its = vvf_hash.equal_range(edge.x);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			if (iter->second.vid == edge.y) {
				int fid = iter->second.fid;
				int vid = geometry::getThirdVertexOfFace((*face_ptr)[fid], edge.x, edge.y);
				neighbor_face_vertex.push_back(int2(fid, vid));
			}
		}
		if (neighbor_face_vertex.size() != 2) {
			std::cout << "error: TriMeshRefiner::SwapLegalCheck, this case may not happen" << std::endl;
			return 0;
		}

		int nv[2] = { neighbor_face_vertex[0].y, neighbor_face_vertex[1].y };

		//	if thin tube exist, cannot swap
		if (IsEdgeExist(nv[0], nv[1]))
			return 0;

		//	lock neighboring edges
		edge_to_lock.push_back(edge);
		edge_to_lock.push_back(int2(nv[0], edge.x));
		edge_to_lock.push_back(int2(nv[0], edge.y));
		edge_to_lock.push_back(int2(nv[1], edge.x));
		edge_to_lock.push_back(int2(nv[1], edge.y));

		//	when this edge is swapped, we must lock edges that may cause thin tube issues
		std::unordered_set<int2, utils::BitwiseHasher<int2>>	nei_edge;
		for (int i = 0; i < 2; i++) {		
			its = vvf_hash.equal_range(nv[i]);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int fid = iter->second.fid;
				for (int j = 0; j < 3; j++) {
					if ((*face_ptr)[fid][j] == nv[i]) {
						//	create edge of the face that is opposite to nv
						int ev0 = (*face_ptr)[fid][(j + 1) % 3];
						int ev1 = (*face_ptr)[fid][(j + 2) % 3];
						if (ev0 > ev1)
							mySwap(ev0, ev1);

						int2 e(ev0, ev1);
						if (i == 0)		//	for nv[0], we record all nei_edge 
							nei_edge.insert(e);
						else if (nei_edge.find(e) != nei_edge.end())	//	for nv[1], if it share nei_edge, clock it
							edge_to_lock.push_back(e);
					}
				}
			}
		}

		return 1;
	}

protected:		//	helper functions
	/**
		Check whether the edge exist
	*/
	inline int IsEdgeExist(int v1, int v2) {
		auto its = vvf_hash.equal_range(v1);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			if (iter->second.vid == v2) {
				return 1;
			}
		}
		return 0;
	}

	/**
		Get all faces connecting the vertex

		VVF must be properly set before calling this function

		\param	nei_face	output faces
		\param	vid			input vertex index
	*/
	inline void GetNeighborFacesOfVertex(std::vector<int>& nei_face, int vid) {
		nei_face.clear();
		auto its = vvf_hash.equal_range(vid);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			nei_face.push_back(iter->second.fid);
		}
		if (nei_face.size() < 2)
			return;

		//	sort and remove duplicates
		std::sort(nei_face.begin(), nei_face.end());

		int count = 0;
		for (unsigned int i = 1; i != nei_face.size(); i++) {
			nei_face[count++] = nei_face[i];
		}

		nei_face.resize(count);
	}

	/**
		Get all faces connecting the vertex

		VVF must be properly set before calling this function

		\param	nei_face	output faces
		\param	vid			input vertex index
	*/
	inline void GetNeighborFacesOfVertex(std::set<int>& nei_face, int vid) {
		nei_face.clear();
		auto its = vvf_hash.equal_range(vid);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			nei_face.insert(iter->second.fid);
		}
	}

	/**
		Get all faces sharing the edge

		VVF must be properly set before calling this function

		\param	nei_face	output faces
		\param	edge		input edge
	*/
	inline void GetNeighborFacesOfEdge(std::vector<int>& nei_face, int2 edge) {
		nei_face.clear();
		auto its = vvf_hash.equal_range(edge.x);
		//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
		for (auto iter = its.first; iter != its.second; iter++) {
			if (iter->second.vid == edge.y) {
				nei_face.push_back(iter->second.fid);
			}
		}
	}

	/**
		Get all vertices that share the same face with edge

		\param	co_face_vertex		output vertices
		\param	edge				input edge
	*/
	inline void GetCoFaceVerticesOfEdge(std::vector<int>& co_face_vertex, int2 edge) {
		GetNeighborFacesOfEdge(co_face_vertex, edge);	//	this step, co_face_vertex stores neighbor faces
		for (unsigned int i = 0; i < co_face_vertex.size(); i++) {
			co_face_vertex[i] = getThirdVertexOfFace((*face_ptr)[co_face_vertex[i]], edge.x, edge.y);
		}
	}

};


/**
	Perform mesh optimization on non-manifold mesh

	Method to use:
	1, call SetMesh(v, f), or construct the class with (v, f), the class will point to the mesh	\n
	2, call RefineNonManifold(collapse_thre, split_thre, area_perc_thre, segment_len_thre, 10000) to perform mesh refinement \n
	3, call SplitSingularVertex() to split the mesh at singular vertices	\n
	4, call RemoveIsolatedVertices() to remove isolated vertices	\n
	5, call RemoveSmallClusters() to remove vertices that shrink to a single cluster or ForceRemoveSmallClustersStep() to force remove (with topology change)
*/
template <class T>
class NonManifoldTriMeshRefiner : protected TriMeshRefiner<T> {
public:
	using TriMeshRefiner<T>::vertex_ptr;
	using TriMeshRefiner<T>::face_ptr;
	using TriMeshRefiner<T>::vvf_hash;
	using TriMeshRefiner<T>::UNORDERED_SET_VEC2I;
	using TriMeshRefiner<T>::locked_edge_set;
	using TriMeshRefiner<T>::edge_to_collapse;
	using TriMeshRefiner<T>::edge_to_swap;
	using TriMeshRefiner<T>::edge_to_split;
	using TriMeshRefiner<T>::vertex_to_merge_vector;
	using TriMeshRefiner<T>::vertex_to_merge_map;
	using TriMeshRefiner<T>::vertex_to_add_vector;
	using TriMeshRefiner<T>::face_to_remove_vector;
	using TriMeshRefiner<T>::face_to_update_vector;
	using TriMeshRefiner<T>::face_to_add_vector;
	using TriMeshRefiner<T>::laplacian_vertex;


public:
	NonManifoldTriMeshRefiner() : TriMeshRefiner<T>() {
	}

	NonManifoldTriMeshRefiner(std::vector<Vec3<T>>& v, std::vector<int3>& f) : TriMeshRefiner<T>(v, f) {
		CreatePatch();
		CalculatePatchArea();
		CreateNonManifold();
	}

	/**
		Setup the non-manifold mesh to refine

		\param	v	mesh vertex
		\param	f	mesh face
	*/
	void SetMesh(std::vector<Vec3<T>>& v, std::vector<int3>& f) {
		TriMeshRefiner<T>::SetMesh(v, f);

		CreatePatch();
		CalculatePatchArea();
		CreateNonManifold();
	}

	/**
		Refine non-manifold mesh

		We perform: 1, shrink small patch; 2, shrink small boundary segment; 3, smooth Taubin; 4, edge operation; 5, update topology.
		then iterate these steps.

		\param	collapse_thre					min edge length
		\param	split_thre						max edge length
		\param	patch_area_percentage_thre		if patch area is smaller than this threshould with total area, shrink it
		\param	boundary_segment_length_thre	if boundary segment is shorter than this threshould, shrink it
		\param	max_iterations					max refinement iterations. By default, just one iteration
		\param	iterations_boundary				number of smoothing iterations on non-manifold boundary
		\param	iterations_face					number of smoothing iterations on manifold face
		\param	coef_forward					forward coefficient of Taubin smooth
		\param	coef_backward					backward coefficient of Taubin smooth
		\return									remain small patch number, remain short boundary segment number, remain strange edge number
	*/
	int3 RefineNonManifold(
		T		collapse_thre,
		T		split_thre,
		T		patch_area_percentage_thre = 0.0,
		T		boundary_segment_length_thre = 1e6,
		int		max_iterations = 1,
		int		iterations_boundary = 1,
		int		iterations_face = 10,
		T		coef_forward = 0.4507499669,
		T		coef_backward = -0.4720265626
	) {
		if (!vertex_ptr || !face_ptr) {
			std::cout << "error: NonManifoldTriMeshRefiner::RefineNonManifold, mesh not set, cannot refine" << std::endl;
			return int3(0, 0, 0);
		}

		Vec3i ret_num;

		for (int i = 0; i < max_iterations; i++) {
			ret_num.x = ShrinkSmallPatchStep(patch_area_percentage_thre);

			ret_num.y = ShrinkSmallBoundarySegmentStep(boundary_segment_length_thre);

			//	edge collapse, swap and split
			Vec3i remain_count = Refine(collapse_thre, split_thre, 1);
			ret_num.z = remain_count.x + remain_count.y + remain_count.z;

			//	rebuild topology
			this->ClearDataArray();
			this->SetupVVF();
			UpdatePatch();
			CalculatePatchArea();
			CreateNonManifold();

			//	if no more to do, stop
			if (ret_num == int3(0, 0, 0))
				break;
			//	if this iteration cannot improve any more, refine can stop
			if (ret_num.x == 0
				&& ret_num.y == 0
				&& edge_to_collapse.size() == remain_count.x
				&& edge_to_swap.size() == remain_count.y
				&& edge_to_split.size() == remain_count.z)
				break;
		}

		return ret_num;
	}

	/**
		Perform Taubin smooth on the mesh

		Taubin smooth is to laplacian forward, then backward, so that the volume is almost preserved.
		For non-manifold mesh, we first smooth non-manifold boundary, then manifold faces.
		Suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

		\param	iterations_boundary		number of smoothing iterations on non-manifold boundary
		\param	iterations_face			number of smoothing iterations on manifold face
		\param	coef_forward			forward coefficient of Taubin smooth
		\param	coef_backward			backward coefficient of Taubin smooth
	*/
	void SmoothTaubin(
		int		iterations_boundary = 1,
		int		iterations_face = 10,
		T		coef_forward = 0.4507499669,
		T		coef_backward = -0.4720265626
	) {
		if (!vertex_ptr || !face_ptr)
			return;

		SmoothTaubinBoundary(iterations_boundary, coef_forward, coef_backward);
		SmoothTaubinFace(iterations_face, coef_forward, coef_backward);
	}

	/**
		Perform Taubin smooth on the non-manifold boundary

		Taubin smooth is to laplacian forward, then backward, so that the overall shape of the boundary is almost preserved.
		Suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

		\param	iterations				number of smoothing iterations on non-manifold boundary
		\param	coef_forward			forward coefficient of Taubin smooth
		\param	coef_backward			backward coefficient of Taubin smooth
	*/
	void SmoothTaubinBoundary(
		int		iterations = 1,
		T		coef_forward = 0.4507499669,
		T		coef_backward = -0.4720265626
	) {
		if (!vertex_ptr || !face_ptr)
			return;

		for (int i = 0; i < iterations; i++) {
			SmoothLaplacianBoundary(coef_forward);
			SmoothLaplacianBoundary(coef_backward);
		}
	}

	/**
		Perform Taubin smooth on the mesh manifold faces

		Taubin smooth is to laplacian forward, then backward, so that the volume is almost preserved.
		This function does not move non-manifold boundary edges.
		Suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

		\param	iterations				number of smoothing iterations on manifold face
		\param	coef_forward			forward coefficient of Taubin smooth
		\param	coef_backward			backward coefficient of Taubin smooth
	*/
	void SmoothTaubinFace(
		int		iterations = 10,
		T		coef_forward = 0.4507499669,
		T		coef_backward = -0.4720265626
	) {
		if (!vertex_ptr || !face_ptr)
			return;

		for (int i = 0; i < iterations; i++) {
			SmoothLaplacianFace(coef_forward);
			SmoothLaplacianFace(coef_backward);
		}
	}

	/**
		Remove isolated vertices inside the mesh, and re-index faces

		\return		number of vertices removed
	*/
	int RemoveIsolatedVertices() {
		if (!vertex_ptr || !face_ptr) {
			return 0;
		}

		int num = removeIsolatedVerticesFromMesh(*vertex_ptr, *face_ptr);

		this->ClearDataArray();
		this->SetupVVF();
		UpdatePatch();
		CalculatePatchArea();
		CreateNonManifold();

		return num;
	}

	/**
		Remove small vertex cluster

		A vertex cluster should contain at least 3 vertices and 1 face.
		The distance between vertices should be smaller than dist_thre.
		We don't remove clusters that may cause topology change.

		\param	dist_thre	distance threshold of cluster
		\return				number of cluster removed
	*/
	int RemoveSmallClusters(T dist_thre) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		//	We label each vertex, 0: not referenced, 1: cluster, should remove, 2: one-ring neighbor of cluster
		//	clusters should not overlap or too close, so if a cluster is to be removed, we don't touch its one-ring neighbor
		std::vector<int> vertex_flag;
		vertex_flag.resize(vertex.size(), 0);

		T dist_thre_square = dist_thre * dist_thre;
		int removed_cluster_count = 0;
		std::deque<int>		cluster_queue;	// queue to perform floodfill
		std::set<int>		cluster_vertex;	// vertices of the cluster
		std::set<int>		cluster_face;	// faces connected with cluster_vertex

											//	parse all vertices and detect whether a cluster exist around it
		for (unsigned int i = 0; i != vertex.size(); i++) {
			if (vertex_flag[i])
				continue;

			//	check whether the vertex is clost to its neighbor
			auto its = vvf_hash.equal_range(i);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int nvid = iter->second.vid;
				if (vertex_flag[nvid])
					continue;
				Vec3<T> r = vertex[nvid] - vertex[i];
				//	if we found close vertex pair, insert them into the queue
				if (r.SquareLength() < dist_thre_square) {
					cluster_queue.push_back(i);
					cluster_vertex.insert(i);
					vertex_flag[i] = 1;
					cluster_queue.push_back(nvid);
					cluster_vertex.insert(nvid);
					vertex_flag[nvid] = 1;
					break;
				}
			}

			if (cluster_queue.empty())
				continue;

			//	now we got a seed of the cluster, extend it by floodfill
			int recover_flag = 0;	//	this flag records whether this cluster is legal

									//	create the cluster by floodfill
			while (!cluster_queue.empty()) {
				int vid = cluster_queue.front();
				cluster_queue.pop_front();

				//	for each vertex, check whether its neighbor can add to this cluster
				auto its = vvf_hash.equal_range(vid);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					int nvid = iter->second.vid;
					if (vertex_flag[nvid] == 1)		// this vertex is already in this cluster
						continue;
					Vec3<T> r = vertex[nvid] - vertex[i];
					if (r.SquareLength() < dist_thre_square) {	//	the neighbor vertex is very close to the seed vertex
						if (vertex_flag[nvid] == 2) {	// this cluster will touch another cluster's neighbor, this cluster is illegal
							recover_flag = 1;
							break;
						}
						cluster_queue.push_back(nvid);
						cluster_vertex.insert(nvid);
						vertex_flag[nvid] = 1;
					}
				}

				if (recover_flag) {
					cluster_queue.clear();
					break;
				}
			}

			// if this cluster is just a short edge, ignore it
			if (cluster_vertex.size() == 2) {
				recover_flag = 2;
			}

			std::vector<int> face_to_remove;	// at least 2 vertices in the cluster
			std::vector<int> face_to_update;	// only 1 vertex in the cluster

			if (!recover_flag) {
				//	create cluster faces: all faces that contain cluster vertices
				for (std::set<int>::iterator iter = cluster_vertex.begin(); iter != cluster_vertex.end(); iter++) {
					auto its = vvf_hash.equal_range(*iter);
					//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
					for (auto iter = its.first; iter != its.second; iter++) {
						int nfid = iter->second.fid;
						cluster_face.insert(nfid);
					}
				}

				//	get all face_to_remove and face_to_update
				int three_vertex_face_flag = 0; // a legal cluster should contain at least 1 face that has 3 cluster vertices
				for (std::set<int>::iterator iter = cluster_face.begin(); iter != cluster_face.end(); iter++) {
					int fid = *iter;

					//	count how many vertices of this face are in cluster_vertex
					int count = 0;
					for (int j = 0; j < 3; j++) {
						if (cluster_vertex.find(face[fid][j]) != cluster_vertex.end()) {
							count++;
						}
					}

					if (count == 3) {		//	3 vertices, this face should remove
						face_to_remove.push_back(fid);
						three_vertex_face_flag = 1;
					}
					else if (count == 2)	//	2 vertices, this face should remove
						face_to_remove.push_back(fid);
					else if (count == 1)	//	1 vertex, just need to change vertex index
						face_to_update.push_back(fid);
				}

				if (!three_vertex_face_flag)
					recover_flag = 3;
			}

			//	check whether removing this cluster may cause duplicate face
			if (!recover_flag) {
				//	we assume that duplicate face doesn't exist originally, so if duplicate face appear, they should come from face_to_update
				std::set<Vec2i> duplicate_face_check;
				for (unsigned int f = 0; f != face_to_update.size(); f++) {
					//	since all these faces will merge at a single vertex, we only check the other two vertices
					for (int j = 0; j < 3; j++) {
						int vid = face[face_to_update[f]][j];
						if (cluster_vertex.find(vid) != cluster_vertex.end()) {
							Vec2i check(face[face_to_update[f]][(j + 1) % 3], face[face_to_update[f]][(j + 2) % 3]);
							if (check.x > check.y)
								mySwap(check.x, check.y);
								
							if (duplicate_face_check.find(check) == duplicate_face_check.end())	//	this face(edge) not exist, record it
								duplicate_face_check.insert(check);
							else	// this face(edge) already exist, it means duplicate face will appear if we remove this cluster
								recover_flag = 4;

							break;
						}
					}

					if (recover_flag)
						break;
				}
			}

			//	if this cluster cannot remove, recover the labels
			if (recover_flag) {
				for (std::set<int>::iterator iter = cluster_vertex.begin(); iter != cluster_vertex.end(); iter++) {
					vertex_flag[*iter] = 0;
				}

				cluster_vertex.clear();
				cluster_face.clear();
				continue;
			}

			//	now we can remove this cluster safely
			//	merge all cluster vertices to a single vertex
			for (std::set<int>::iterator iter = cluster_vertex.begin(); iter != cluster_vertex.end(); iter++) {
				if (*iter == i)
					continue;
				vertex_to_merge_vector.push_back(int2(*iter, i));
				vertex_to_merge_map.insert(std::pair<int, int>(*iter, i));
			}

			//	remove all faces that no longer exist
			for (unsigned int j = 0; j != face_to_remove.size(); j++) {
				face_to_remove_vector.push_back(face_to_remove[j]);
			}

			//	re-label faces that only 1 vertex in cluster_vertex, and label neighborhood
			for (unsigned int j = 0; j != face_to_update.size(); j++) {
				int cluster_v;
				for (int k = 0; k < 3; k++) {
					int vid = face[face_to_update[j]][k];
					if (cluster_vertex.find(vid) != cluster_vertex.end())
						cluster_v = vid;
					else
						vertex_flag[vid] = 2;
				}

				face_to_update_vector.push_back(int3(face_to_update[j], cluster_v, i));
			}

			removed_cluster_count++;

			cluster_vertex.clear();
			cluster_face.clear();
		}

		//	update the mesh
		if (removed_cluster_count) {
			ReIndexMesh();

			this->ClearDataArray();
			this->SetupVVF();
			UpdatePatch();
			CalculatePatchArea();
			CreateNonManifold();
		}

		return removed_cluster_count;
	}

	/**
		Force remove small vertex cluster

		If removing this cluster may cause topology issue, we shrink a bigger
		area towards this cluster, making it safe to remove. It is possible that 
		removing a cluster could cause singular vertex, we don't solve this issue.

		\param	dist_thre	distance threshold of cluster
		\return				number of cluster removed
	*/
	int ForceRemoveSmallClustersStep(T dist_thre) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		int cluster_num = CreateSmallClusters(dist_thre);
		if (!cluster_num)	//	no clusters
			return 0;

		int remain_cluster_count = 0;

		//	move clusters that are adjacent towards their center
		std::vector<int> cluster_isolate_flag;
		cluster_isolate_flag.resize(cluster_num, 1);
		std::vector<Vec3<T>> cluster_offset;
		cluster_offset.resize(cluster_num);

		//	calculate displacement of adjacent clusters
		for (std::set<Vec2i>::iterator iter = neighbor_cluster_pair.begin(); iter != neighbor_cluster_pair.end(); iter++) {
			int c0 = iter->x, c1 = iter->y;
			cluster_isolate_flag[c0] = 0;
			cluster_isolate_flag[c1] = 0;

			int c0_vid = cluster_vertex[cluster_vertex_start[c0]];
			int c1_vid = cluster_vertex[cluster_vertex_start[c1]];

			Vec3<T> r = vertex[c1_vid] - vertex[c0_vid];
			cluster_offset[c0] += r * 0.45;
			cluster_offset[c1] -= r * 0.45;
		}

		//	perform offset on adjacent clusters
		for (unsigned int c = 0; c != cluster_num; c++) {
			if (cluster_isolate_flag[c])
				continue;
			for (unsigned int j = cluster_vertex_start[c]; j != cluster_vertex_start[c + 1]; j++) {
				int vid = cluster_vertex[j];
				vertex[vid] += cluster_offset[c];
			}
			remain_cluster_count++;
		}

		//	for isolated clusters, check whether it can be removed directly
		for (unsigned int c = 0; c != cluster_num; c++) {
			if (!cluster_isolate_flag[c])
				continue;

			//	create cluster faces: all faces that contain cluster vertices
			std::set<int>		cluster_face;
			for (unsigned int j = cluster_vertex_start[c]; j != cluster_vertex_start[c + 1]; j++) {
				int vid = cluster_vertex[j];
				auto its = vvf_hash.equal_range(vid);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					int nfid = iter->second.fid;
					cluster_face.insert(nfid);
				}
			}

			//	seperate faces inside this cluster
			std::vector<int> face_to_remove;	// at least 2 vertices in the cluster, will be collapsed
			std::vector<int> face_to_update;	// only 1 vertex in the cluster, one ring neighbor of this cluster

			for (std::set<int>::iterator iter = cluster_face.begin(); iter != cluster_face.end(); iter++) {
				int fid = *iter;

				//	count how many vertices of this face are in cluster_vertex
				int count = 0;
				for (int j = 0; j < 3; j++) {
					int vid = face[fid][j];
					if (vertex_cluster_index[vid] == c)
						count++;
				}

				if(count == 1)
					face_to_update.push_back(fid);
				else
					face_to_remove.push_back(fid);
			}

			//	check whether this cluster may cause duplicate faces
			std::set<Vec2i> cluster_boundary_edge;
			std::set<int> cluster_singular_vertex;
			for (unsigned int f = 0; f != face_to_update.size(); f++) {
				//	since all these faces will merge at a single vertex, we only check the other two vertices
				for (int j = 0; j < 3; j++) {
					int vid = face[face_to_update[f]][j];
					if (vertex_cluster_index[vid] == c) {
						Vec2i cluster_edge(face[face_to_update[f]][(j + 1) % 3], face[face_to_update[f]][(j + 2) % 3]);
						if (cluster_edge.x > cluster_edge.y)
							mySwap(cluster_edge.x, cluster_edge.y);

						if (cluster_boundary_edge.find(cluster_edge) == cluster_boundary_edge.end())	//	this face(edge) not exist, record it
							cluster_boundary_edge.insert(cluster_edge);
						else {	// this face(edge) already exist, it means duplicate face will appear if we remove this cluster
							cluster_singular_vertex.insert(cluster_edge.x);
							cluster_singular_vertex.insert(cluster_edge.y);
						}

						break;
					}
				}
			}

			//	duplicate face appear, so we shrink it to the center
			if (!cluster_singular_vertex.empty()) {
				int cluster_first_vid = cluster_vertex[cluster_vertex_start[c]];
				Vec3<T> center = vertex[cluster_first_vid];
				for (std::set<int>::iterator iter = cluster_singular_vertex.begin(); iter != cluster_singular_vertex.end(); iter++)
					center += vertex[*iter];
				center /= (cluster_singular_vertex.size() + 1);

				//	offset cluster
				Vec3<T> offset = (center - vertex[cluster_first_vid]) * 0.9;
				for (unsigned int j = cluster_vertex_start[c]; j != cluster_vertex_start[c + 1]; j++) {
					vertex[cluster_vertex[j]] += offset;
				}
				//	offset each singular vertex
				for (std::set<int>::iterator iter = cluster_singular_vertex.begin(); iter != cluster_singular_vertex.end(); iter++) {
					vertex[*iter] += (center - vertex[*iter]) * 0.9;
				}

				remain_cluster_count++;
			}	//	else, we can remove this cluster
			else {
				int merge_target = cluster_vertex[cluster_vertex_start[c]];

				//	now we can remove this cluster safely, merge all cluster vertices to a single vertex
				for (unsigned int j = cluster_vertex_start[c] + 1; j != cluster_vertex_start[c + 1]; j++) {
					vertex_to_merge_vector.push_back(int2(cluster_vertex[j], merge_target));
					vertex_to_merge_map.insert(std::pair<int, int>(cluster_vertex[j], merge_target));
				}

				//	remove all faces that no longer exist
				for (unsigned int j = 0; j != face_to_remove.size(); j++) {
					face_to_remove_vector.push_back(face_to_remove[j]);
				}

				//	re-label faces that only 1 vertex in cluster_vertex, and label neighborhood
				for (unsigned int j = 0; j != face_to_update.size(); j++) {
					int cluster_v;
					for (int k = 0; k < 3; k++) {
						int vid = face[face_to_update[j]][k];
						if (vertex_cluster_index[vid] == c) {
							cluster_v = vid;
							break;
						}
					}

					face_to_update_vector.push_back(int3(face_to_update[j], cluster_v, merge_target));
				}
			}
		}

		//	update the mesh
		if (remain_cluster_count != cluster_num) {
			ReIndexMesh();

			this->ClearDataArray();
			this->SetupVVF();
			UpdatePatch();
			CalculatePatchArea();
			CreateNonManifold();
		}

		return remain_cluster_count;
	}

	/**
		Split the mesh at singular vertices

		\return		number of new added vertices
	*/
	int SplitSingularVertex() {
		if (!vertex_ptr || !face_ptr)
			return 0;

		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		std::set<int> nei_face_of_v;
		std::deque<int> face_queue;
		std::vector<int> fan_face;
		std::vector<int> nei_face_of_e;

		int split_vertex_count = 0;

		//	parse all vertices, check whether it is singular vertex
		for (unsigned int vid = 0; vid != vertex.size(); vid++) {
			//	get all neighbor faces of this vertex, then detect fans
			this->GetNeighborFacesOfVertex(nei_face_of_v, vid);

			while (!nei_face_of_v.empty()) {
				//	pick one face as seed
				fan_face.clear();
				face_queue.push_back(*nei_face_of_v.begin());
				fan_face.push_back(*nei_face_of_v.begin());
				nei_face_of_v.erase(nei_face_of_v.begin());

				//	floodfill with face-edge-face connectivity to get a fan
				while (!face_queue.empty()) {
					int fid = face_queue.front();
					face_queue.pop_front();

					for (int i = 0; i < 3; i++) {
						int v0 = face[fid][i], v1 = face[fid][(i + 1) % 3];
						if (v0 != vid && v1 != vid)
							continue;

						//	get neighbor face of existing fan faces
						this->GetNeighborFacesOfEdge(nei_face_of_e, int2(v0, v1));
						for (unsigned int j = 0; j != nei_face_of_e.size(); j++) {
							int nfid = nei_face_of_e[j];
							std::set<int>::iterator iter = nei_face_of_v.find(nfid);
							if (iter != nei_face_of_v.end()) {
								//	a new face is detected in this fan, record it
								face_queue.push_back(nfid);
								fan_face.push_back(nfid);
								nei_face_of_v.erase(iter);
							}
						}
					}
				}

				//	a fan is detected, if there are non-labeled faces, we need to split the fan here
				if (!nei_face_of_v.empty()) {
					int split_vid = vertex.size() + split_vertex_count;
					split_vertex_count++;
					vertex_to_add_vector.push_back(vertex[vid]);
					for (unsigned int i = 0; i != fan_face.size(); i++) {
						face_to_update_vector.push_back(int3(fan_face[i], vid, split_vid));
					}
				}
			}
		}

		if (split_vertex_count) {
			ReIndexMesh();

			this->ClearDataArray();
			this->SetupVVF();
		}

		return split_vertex_count;
	}

	/**
		Draw the non-manifold mesh

		If openGL is not included, this function will not work properly
	*/
	void Display() {
#ifdef YZ_gl_h
		if (!vertex_ptr || !face_ptr)
			return;

		std::vector<Vec3<T>>& vertex = (*vertex_ptr);
		std::vector<int3>& face = (*face_ptr);

		//	calculate average edge length of the mesh
		T avg_edge_length = 0;
		for (unsigned int i = 0; i != face.size(); i++) {
			for (int j = 0; j < 3; j++) {
				avg_edge_length += (vertex[face[i][j]] - vertex[face[i][(j + 1) % 3]]).Length();
			}
		}
		avg_edge_length /= (face.size() * 3);

		//	calculate face normal
		std::vector<Vec3<T>>	face_normal;
		calculateFaceNormal(face_normal, *vertex_ptr, *face_ptr);

		//	draw the mesh as wireframe, with different patch color
		std::vector<uchar3>	face_color;
		face_color.resize(face.size());
		for (int i = 0; i < face_color.size(); i++) {
			opengl::getSequentialDisplayColor(&face_color[i].x, face_mark[i]);
		}

		opengl::drawFlatColorTriMeshEdgeFromFace(vertex, face, face_color, face_normal);

		//	draw non-manifold boundary edge
		glColor3f(1, 1, 1);
		for (std::map<Vec2i, int>::iterator iter = non_manifold_edge_nei_face_count.begin();
			iter != non_manifold_edge_nei_face_count.end();
			iter++)
		{
			Vec2i edge = iter->first;
			opengl::drawCylinder(vertex[edge.x], vertex[edge.y], avg_edge_length*0.1, 8);
		}

		//	draw non-manifold boundary joint vertex
		glColor3f(1, 0, 0);
		for (int i = 0; i < non_manifold_boundary_joint_vertex.size(); i++) {
			int vid = non_manifold_boundary_joint_vertex[i];
			opengl::drawPointAsCube(vertex[vid], avg_edge_length*0.5);
		}

		//	draw non-manifold boundary end vertex
		glColor3f(0, 1, 0);
		for (int i = 0; i < non_manifold_boundary_end_vertex.size(); i++) {
			int vid = non_manifold_boundary_end_vertex[i];
			opengl::drawPointAsCube(vertex[vid], avg_edge_length*0.5);
		}
#else
		std::cout << "gl.h has to be included in order to use Display() in NonManifoldTriMeshRefiner" << std::endl;
#endif

	}

protected:
	//	manifold patches, seperated by non-manifold boundary edges
	std::vector<int>		face_mark;								///< mark (indicating patch index) of each face
	std::vector<int>		patch_face_count;						///< face number of each patch
	std::vector<T>			patch_area;								///< face area sum of each patch

	//	non-manifold boundary
	std::map<Vec2i, int>	non_manifold_edge_nei_face_count;		///< all non-manifold edge and corresponding neighbor faces
	std::multimap<int, int>	non_manifold_boundary_vv;				///< vertex-vertex on non-manifold boundary
	std::vector<int>		non_manifold_boundary_joint_vertex;		///< vertex connecting different boundary segments
	std::vector<int>		non_manifold_boundary_end_vertex;		///< vertex that is the end of a boundary segment, exclude joint vertex
	std::vector<int>		non_manifold_boundary_normal_vertex;	///< boundary vertex exclude joint vertex and end vertex

	//	cluster
	std::vector<int>		vertex_cluster_index;					///< cluster index of each vertex. -1: doesn't belong to a cluster
	std::vector<int>		cluster_vertex;							///< vertices inside clusters, arranged from index 0 to the last
	std::vector<int>		cluster_vertex_start;					///< vertices of each cluster in cluster_vertex, i.e. [cluster_vertex_start[0], cluster_vertex_start[1]) is the first cluster
	std::set<Vec2i>			neighbor_cluster_pair;					///< all cluster pairs that are adjacent (cannot remove without merging)

protected:
	/**
		Overload ReIndexMesh(), add face_mark operation
	*/
	void ReIndexMesh() {
		//	original re-index 
		TriMeshRefiner<T>::ReIndexMesh();

		//	perform re-index on face_mark
		this->ReIndexFaceData(face_mark);
	}

	/**
		Create patches seperated by non-manifold edges using floodfill algorithm

		\return		number of patches
	*/
	int CreatePatch() {
		std::vector<int3>& face = (*face_ptr);

		face_mark.clear();
		face_mark.resize(face.size(), -1);

		std::deque<int>	face_queue;
		std::vector<int> nei_face;
		int seed_value = 0;

		for (int seed_fid = 0; seed_fid < face.size(); seed_fid++) {
			if (face_mark[seed_fid] != -1)
				continue;

			//	get a seed face
			face_mark[seed_fid] = seed_value;
			face_queue.push_back(seed_fid);

			//	floodfill until no more faces can be added to the patch
			while (!face_queue.empty()) {
				int fid = face_queue.front();
				face_queue.pop_front();

				//	get neighbor face of the seed face, extend the patch if possible
				for (int i = 0; i < 3; i++) {
					int2 nei_e(face[fid][i], face[fid][(i + 1) % 3]);

					this->GetNeighborFacesOfEdge(nei_face, nei_e);

					if (nei_face.size() != 2)
						continue;

					int nei_fid = nei_face[0] == fid ? nei_face[1] : nei_face[0];

					if (face_mark[nei_fid] == -1) {
						//	check normal direction of this face and insure face normal on the same patch are aligned
						for (int j = 0; j < 3; j++) {
							if (face[nei_fid][j] == nei_e.x && face[nei_fid][(j + 1) % 3] == nei_e.y) {
								mySwap(face[nei_fid][1], face[nei_fid][2]);
								break;
							}
						}

						face_mark[nei_fid] = seed_value;
						face_queue.push_back(nei_fid);
					}
				}
			}

			seed_value++;
		}

		return seed_value;
	}

	/**
		Update patch after edge optimization

		patch may appear and disappear during edge optimization, this function make face mark as continuous as possible
	*/
	int UpdatePatch() {
		//	record old face_mark
		std::vector<int>	old_face_mark;
		old_face_mark.swap(face_mark);

		//	create patch
		int curr_patch_num = CreatePatch();

		if (old_face_mark.size() != face_mark.size())
			return 0;

		//	check correspondence between current and old
		std::set<Vec2i> mapping;
		for (int i = 0; i < face_mark.size(); i++) {
			mapping.insert(Vec2i(old_face_mark[i], face_mark[i]));
		}

		std::vector<int> relabel, old_mark_used_flag;
		relabel.resize(curr_patch_num, -1);
		old_mark_used_flag.resize(patch_face_count.size(), 0);
		for (std::set<Vec2i>::iterator iter = mapping.begin(); iter != mapping.end(); iter++) {
			int old_label = iter->x;
			int curr_label = iter->y;
			if (relabel[curr_label] < 0) {
				if (old_mark_used_flag[old_label]) {
					relabel[curr_label] = old_mark_used_flag.size();
					old_mark_used_flag.push_back(1);
				}
				else {
					old_mark_used_flag[old_label] = 1;
					relabel[curr_label] = old_label;
				}
			}
		}

		for (int i = 0; i < relabel.size(); i++) {
			if (relabel[i] == -1) {
				std::cout << "error: NonManifoldTriMeshRefiner::UpdatePatch, mapping didn't found" << std::endl;
			}
		}

		//	relabel current face_mark
		for (int i = 0; i < face_mark.size(); i++) {
			face_mark[i] = relabel[face_mark[i]];
		}

		return curr_patch_num;
	}

	/**
		Calculate face number and face area sum of each patch
	*/
	void CalculatePatchArea() {
		patch_face_count.clear();
		patch_area.clear();

		for (unsigned int fid = 0; fid != face_ptr->size(); fid++) {
			int mark_id = face_mark[fid];

			if (patch_face_count.size() <= mark_id) {
				patch_face_count.resize(mark_id + 1, 0);
				patch_area.resize(mark_id + 1, 0);
			}

			patch_face_count[mark_id] ++;
			patch_area[mark_id] += FaceArea(fid);
		}
	}

	/**
		Create non-manifold structure
	*/
	void CreateNonManifold() {
		//	clear existing data
		non_manifold_edge_nei_face_count.clear();
		non_manifold_boundary_vv.clear();
		non_manifold_boundary_joint_vertex.clear();
		non_manifold_boundary_end_vertex.clear();
		non_manifold_boundary_normal_vertex.clear();

		//	get all non-manifold edges, and corresponding neighbor face number
		std::vector<int3>& face = (*face_ptr);
		std::vector<int> nei_face;
		for (int fid = 0; fid < face.size(); fid++) {
			for (int i = 0; i < 3; i++) {
				int2 edge(face[fid][i], face[fid][(i + 1) % 3]);
				if (edge.x > edge.y)
					mySwap(edge.x, edge.y);

				this->GetNeighborFacesOfEdge(nei_face, edge);

				if (nei_face.size() != 2)	//	non-manifold edge has 1 or more than 2 neighbor faces
					non_manifold_edge_nei_face_count.insert(std::pair<Vec2i, int>(edge, nei_face.size()));
			}
		}

		//	create non-manifold vertex-vertex
		for (std::map<Vec2i, int>::iterator iter = non_manifold_edge_nei_face_count.begin();
			iter != non_manifold_edge_nei_face_count.end();
			iter++)
		{
			Vec2i edge = iter->first;
			non_manifold_boundary_vv.insert(std::pair<int, int>(edge.x, edge.y));
			non_manifold_boundary_vv.insert(std::pair<int, int>(edge.y, edge.x));
		}

		//	parse vv
		std::map<int, int> vv_count;
		for (std::multimap<int, int>::iterator iter = non_manifold_boundary_vv.begin();
			iter != non_manifold_boundary_vv.end();
			iter++)
		{
			std::map<int, int>::iterator vv_iter = vv_count.find(iter->first);
			if (vv_iter == vv_count.end())
				vv_count.insert(std::pair<int, int>(iter->first, 1));
			else
				vv_iter->second++;
		}

		//	create non-manifold vertices of different type
		for (std::map<int, int>::iterator iter = vv_count.begin(); iter != vv_count.end(); iter++) {
			if (iter->second == 1)		// only 1 neighbor boundary vertex
				non_manifold_boundary_end_vertex.push_back(iter->first);
			else if (iter->second > 2)	// 3 and more neighbor boundary vertices
				non_manifold_boundary_joint_vertex.push_back(iter->first);
			else {	// 2 neighbor boundary vertices
				std::multimap<int, int>::iterator vv_iter = non_manifold_boundary_vv.find(iter->first);
				int2 edge_0(iter->first, vv_iter->second);
				vv_iter++;
				int2 edge_1(iter->first, vv_iter->second);
				if (edge_0.x > edge_0.y)
					mySwap(edge_0.x, edge_0.y);
				if (edge_1.x > edge_1.y)
					mySwap(edge_1.x, edge_1.y);
				std::map<Vec2i, int>::iterator iter0 = non_manifold_edge_nei_face_count.find(edge_0);
				std::map<Vec2i, int>::iterator iter1 = non_manifold_edge_nei_face_count.find(edge_1);
				if (iter0->second == iter1->second)
					non_manifold_boundary_normal_vertex.push_back(iter->first);
				else
					non_manifold_boundary_joint_vertex.push_back(iter->first);
			}				
		}
	}

	/**
		Shrink small patches (with non-manifold boundary) toward its center

		We don't move boundary. For adjacent small patches, we shrink the smaller patch first

		\param	patch_area_percentage_thre		if patch area is smaller than this threshould with total area, shrink it
		\return									number of small patches
	*/
	int ShrinkSmallPatchStep(T patch_area_percentage_thre = 0.0) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		//	calculate patch area threshold
		T patch_area_sum = 0;
		for (unsigned int i = 0; i != patch_area.size(); i++) {
			patch_area_sum += patch_area[i];
		}
		T patch_area_thre = patch_area_sum * patch_area_percentage_thre;

		//	get small patch faces
		std::vector<std::vector<int>> small_patch_fid;
		small_patch_fid.resize(patch_area.size());
		for (unsigned int fid = 0; fid != face.size(); fid++) {
			int mark_id = face_mark[fid];
			if (patch_area[mark_id] < patch_area_thre)
				small_patch_fid[mark_id].push_back(fid);
		}

		//	sort patch, and shrink from smallest to biggest
		std::vector<int> index;
		utils::sortIndex(index, patch_area.begin(), patch_area.end());
		std::map<int, std::pair<int, Vec3<T>>> vertex_patch_displace;	//	vertex index, (patch index, displace vector)
		int small_patch_count = 0;

		for (unsigned int idx = 0; idx != index.size(); idx++) {
			int patch_idx = index[idx];
			if (patch_area[patch_idx] > patch_area_thre)	//	if the patch area is big enough, there is no more small patches
				break;
			if (!patch_face_count[patch_idx])				//	if this patch no longer exist, skip it
				continue;

			small_patch_count++;
			std::vector<int> nei_f;

			//	this is a small patch, shrink it
			for (unsigned int j = 0; j != small_patch_fid[patch_idx].size(); j++) {
				int fid = small_patch_fid[patch_idx][j];
				//	for each edge of the face
				for (int i = 0; i < 3; i++) {
					int vid0 = face[fid][i], vid1 = face[fid][(i + 1) % 3];

					this->GetNeighborFacesOfEdge(nei_f, int2(vid0, vid1));	//	check whether this edge is boundary
					if (nei_f.size() == 1)								//	don't move boundary edge
						continue;

					Vec3<T> r = (vertex[vid1] - vertex[vid0]) * (T)0.1;	//	we use a small coefficient

					//std::map<int, std::pair<int, Vec3<T>>>::iterator iter0 = vertex_patch_displace.find(vid0);
					auto iter0 = vertex_patch_displace.find(vid0);
					if (iter0 == vertex_patch_displace.end())	//	this vertex is never displaced
						vertex_patch_displace.insert(std::pair<int, std::pair<int, Vec3<T>>>(vid0, std::pair<int, Vec3<T>>(patch_idx, r)));
					else if (iter0->second.first == patch_idx)	//	this vertex is already displaced by this patch
						iter0->second.second += r;				//	if this vertex is displaced by other patch(smaller patches), don't move it

					//std::map<int, std::pair<int, Vec3<T>>>::iterator iter1 = vertex_patch_displace.find(vid1);
					auto iter1 = vertex_patch_displace.find(vid1);
					if (iter1 == vertex_patch_displace.end())
						vertex_patch_displace.insert(std::pair<int, std::pair<int, Vec3<T>>>(vid1, std::pair<int, Vec3<T>>(patch_idx, -r)));
					else if (iter1->second.first == patch_idx)
						iter1->second.second -= r;
				}
			}
		}

		//	displace all vertices
		//for (std::map<int, std::pair<int, Vec3<T>>>::iterator iter = vertex_patch_displace.begin(); iter != vertex_patch_displace.end(); iter++) {
		for (auto iter = vertex_patch_displace.begin(); iter != vertex_patch_displace.end(); iter++) {
			vertex[iter->first] += iter->second.second;
		}

		return small_patch_count;
	}

	/**
		Shrink small boundary segment (between joint vertices or end vertices) toward direction of making the segment shorter

		\param	boundary_segment_thre		if boundary segment is shorter than this threshould, shrink it
	*/
	int ShrinkSmallBoundarySegmentStep(T boundary_segment_thre = 1e6) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);

		//	get all short boundary segments
		std::set<Vec2i>			referenced_edge;
		std::vector<double>		segment_length;
		std::vector<Vec2i>		segment_end_vid;
		std::vector<Vec3<T>>	segment_end_displace0;
		std::vector<Vec3<T>>	segment_end_displace1;

		//	for each non-manifold, we get corresponding segment
		for (std::map<Vec2i, int>::iterator iter = non_manifold_edge_nei_face_count.begin();
			iter != non_manifold_edge_nei_face_count.end();
			iter++)
		{
			//	get a seed edge that is not referenced yet
			Vec2i edge = iter->first;
			if (referenced_edge.find(edge) != referenced_edge.end())
				continue;
			T segment_len = EdgeLength(edge.x, edge.y);
			referenced_edge.insert(edge);

			//	extend towards two directions, until meets a joint or end vertex, or detected loop
			int loop_flag = 0;
			Vec2i edge_extend[2] = { edge, Vec2i(edge.y, edge.x) };
			for (int i = 0; i < 2; i++) {
				if (loop_flag)
					break;
				while (non_manifold_boundary_vv.count(edge_extend[i].y) == 2) {		//	edge_extend[i].y is normal vertex, extend
					std::multimap<int, int>::iterator iter = non_manifold_boundary_vv.find(edge_extend[i].y);
					if (iter->second == edge_extend[i].x)
						iter++;
					edge_extend[i].x = edge_extend[i].y;
					edge_extend[i].y = iter->second;
					segment_len += EdgeLength(edge_extend[i].x, edge_extend[i].y);

					edge = edge_extend[i];
					if (edge.x > edge.y)
						mySwap(edge.x, edge.y);

					//	if edge already referenced, it forms a loop
					if (referenced_edge.find(edge) == referenced_edge.end())
						referenced_edge.insert(edge);
					else {
						loop_flag = 1;
						break;
					}
				}
			}

			if (segment_len > boundary_segment_thre || loop_flag)
				continue;

			//	we found a short boundary segment, record it
			segment_length.push_back(segment_len);
			segment_end_vid.push_back(int2(edge_extend[0].y, edge_extend[1].y));
			segment_end_displace0.push_back((vertex[edge_extend[0].x] - vertex[edge_extend[0].y]) * 0.45);
			segment_end_displace1.push_back((vertex[edge_extend[1].x] - vertex[edge_extend[1].y]) * 0.45);
		}

		//	reorder segment according to boundary length
		if (segment_length.empty())
			return 0;
		std::vector<int> index;
		utils::sortIndex(index, segment_length.begin(), segment_length.end());
		utils::reorderDataBySourceOrder(segment_length, index);
		utils::reorderDataBySourceOrder(segment_end_vid, index);
		utils::reorderDataBySourceOrder(segment_end_displace0, index);
		utils::reorderDataBySourceOrder(segment_end_displace1, index);

		//	displace vertices
		std::set<int> referenced_vertex;
		for (unsigned int i = 0; i != segment_end_vid.size(); i++) {
			if (referenced_vertex.find(segment_end_vid[i].x) == referenced_vertex.end()) {
				vertex[segment_end_vid[i].x] += segment_end_displace0[i];
				referenced_vertex.insert(segment_end_vid[i].x);
			}
			if (referenced_vertex.find(segment_end_vid[i].y) == referenced_vertex.end()) {
				vertex[segment_end_vid[i].y] += segment_end_displace1[i];
				referenced_vertex.insert(segment_end_vid[i].y);
			}
		}

		return segment_length.size();
	}

	/**
		Perform Laplacian smooth on non-manifold boundary, fixing joint vertices and end vertices

		\param	coef	linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
	*/
	void SmoothLaplacianBoundary(T coef = 1.0) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);

		std::vector<Vec3<T>> boundary_normal_vertex;
		boundary_normal_vertex.resize(non_manifold_boundary_normal_vertex.size());

		//	for each normal vertex on boundary edge, calculate average position of its neighbor, then interpolate
		for (unsigned int i = 0; i != non_manifold_boundary_normal_vertex.size(); i++) {
			int vid = non_manifold_boundary_normal_vertex[i];
			auto its = non_manifold_boundary_vv.equal_range(vid);
			Vec3<T> v_sum(0, 0, 0);
			for (auto iter = its.first; iter != its.second; iter++) {
				v_sum += vertex[iter->second];
			}
			v_sum *= 0.5;	//	for normal vertices on boundary edge, each has two neighbors
			boundary_normal_vertex[i] = (1.0 - coef) * vertex[vid] + coef * v_sum;
		}

		//	write the smoothed vertex to the mesh
		for (unsigned int i = 0; i != non_manifold_boundary_normal_vertex.size(); i++) {
			int vid = non_manifold_boundary_normal_vertex[i];
			vertex[vid] = boundary_normal_vertex[i];
		}
	}

	/**
		Perform Laplacian smooth on the mesh, exclude non-manifold edges

		\param	coef	linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
	*/
	void SmoothLaplacianFace(T coef = 1.0) {
		//	we first record all non-manifold vertices, then perform Laplacian smooth, and finally recover them
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);

		//	record vertex position on boundary
		std::vector<std::pair<int, Vec3<T>>> fixed_vertex;
		fixed_vertex.reserve(
			non_manifold_boundary_joint_vertex.size() +
			non_manifold_boundary_end_vertex.size() +
			non_manifold_boundary_normal_vertex.size());
		std::vector<int>* vid_ptr[3] = {
			&non_manifold_boundary_joint_vertex,
			&non_manifold_boundary_end_vertex,
			&non_manifold_boundary_normal_vertex
		};
		for (int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j != vid_ptr[i]->size(); j++) {
				int vid = (*vid_ptr[i])[j];
				fixed_vertex.push_back(std::pair<int, Vec3<T>>(vid, vertex[vid]));
			}
		}

		//	perform laplacian smooth on all vertices
		SmoothLaplacian(coef);

		//	recover fixed vertices
		for (unsigned int i = 0; i != fixed_vertex.size(); i++) {
			vertex[fixed_vertex[i].first] = fixed_vertex[i].second;
		}
	}

	/**
		Create all clusters (two adjacent vertices are too close) and record them in the array

		\param	dist_thre		radius of the cluster
		\return					number of detected clusters
	*/
	int CreateSmallClusters(T dist_thre) {
		std::vector<Vec3<T>>&	vertex = (*vertex_ptr);
		std::vector<int3>&		face = (*face_ptr);

		vertex_cluster_index.clear();
		vertex_cluster_index.resize(vertex.size(), -1);
		cluster_vertex.clear();
		cluster_vertex_start.clear();
		cluster_vertex_start.push_back(0);
		neighbor_cluster_pair.clear();

		T dist_thre_square = dist_thre * dist_thre;
		std::deque<int>		cluster_queue;	// queue to perform floodfill
		int cluster_index = 0;

		for (unsigned int i = 0; i != vertex.size(); i++) {
			if (vertex_cluster_index[i] != -1)	//	already in another cluster, skip
				continue;

			//	check whether the vertex is clost to its neighbor, if so, create the initial cluster
			auto its = vvf_hash.equal_range(i);
			//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
			for (auto iter = its.first; iter != its.second; iter++) {
				int nvid = iter->second.vid;
				Vec3<T> r = vertex[nvid] - vertex[i];
				//	if we found a close vertex pair
				if (r.SquareLength() < dist_thre_square) {
					//	we set the first vertex as the center of this cluster
					cluster_queue.push_back(i);
					cluster_vertex.push_back(i);
					vertex_cluster_index[i] = cluster_index;
					break;
				}
			}

			//	if this vertex is not close to its neighbor, it is not a cluster
			if (cluster_queue.empty())
				continue;

			//	now we got a seed of the cluster, extend it by floodfill
			while (!cluster_queue.empty()) {
				int vid = cluster_queue.front();
				cluster_queue.pop_front();

				//	for each vertex, check whether its neighbor can add to this cluster
				auto its = vvf_hash.equal_range(vid);
				//for (std::unordered_multimap<int, VVF>::iterator iter = its.first; iter != its.second; iter++) {
				for (auto iter = its.first; iter != its.second; iter++) {
					int nvid = iter->second.vid;
					//	if this vertex already belongs to a cluster (including itself), skip
					if (vertex_cluster_index[nvid] != -1) {
						if (vertex_cluster_index[nvid] != cluster_index)	// if two clusters are very close, record the pair
							neighbor_cluster_pair.insert(Vec2i(vertex_cluster_index[nvid], cluster_index));
						continue;
					}
					//	we have a un-labeled vertex, check whether it belongs to this cluster
					Vec3<T> r = vertex[nvid] - vertex[i];
					if (r.SquareLength() < dist_thre_square) {
						cluster_queue.push_back(nvid);
						cluster_vertex.push_back(nvid);
						vertex_cluster_index[nvid] = cluster_index;
					}
				}
			}

			//	now we have got a cluster, record it
			cluster_index++;
			cluster_vertex_start.push_back(cluster_vertex.size());
		}

		return cluster_index;
	}



protected:		//	helper functions
	/**
		Get area of the face
	*/
	inline T FaceArea(int fid) {
		Vec3<T> r1 = (*vertex_ptr)[(*face_ptr)[fid].y] - (*vertex_ptr)[(*face_ptr)[fid].x];
		Vec3<T> r2 = (*vertex_ptr)[(*face_ptr)[fid].z] - (*vertex_ptr)[(*face_ptr)[fid].x];
		return cross(r1, r2).Length() * (T)0.5;
	}

	/**
		Get edge length
	*/
	inline T EdgeLength(int vid0, int vid1) {
		Vec3<T> r = (*vertex_ptr)[vid0] - (*vertex_ptr)[vid1];
		return r.Length();
	}
};

///@}

}}}	//	namespace yz::geometry::remeshing



#endif	//	__YZ_MESH_OPTIMIZATION_H__