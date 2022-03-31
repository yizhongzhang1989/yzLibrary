/***********************************************************/
/**	\file
	\brief		Repair a Mesh
	\details	This file provides some functions to repair 
				a mesh if the mesh contain some bad features
	\author		Yizhong Zhang
	\date		9/19/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_REPAIR_H__
#define __YZ_MESH_REPAIR_H__

#include <vector>
#include <unordered_map>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_hash.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Repair Mesh Topology
*/
//	========================================

/**
	reorder vertices according to reorder_index

	Sometimes, we need to reorder vertices of a mesh. For example,
	several vertices are removed. Then we perform reorderVertices.
	The new index of each vertex is stored in reorder_index. If this
	vertex is removed, its index is negtive value. reorder_index can
	only change position of each vertex, no other usage can be made.

	\param	vertex			vertex list
	\param	face			face list
	\param	reorder_index	new index of each vertex, array size should equal to vertex
	\return					number of remaining vertices
*/
template<typename T>
int reorderVertices(std::vector<Vec3<T>>&	vertex,
					std::vector<int3>&		face,
					const std::vector<int>&	reorder_index ){
	if( vertex.size() != reorder_index.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: reorderVertices. vertex list size don't match reorder_index list size" << std::endl;
		#endif
		return vertex.size();
	}

	//	find the start position that need to be corrected
	int start = 0;
	for(int i=0; i<reorder_index.size(); i++)
		if(reorder_index[i] != i){
			start = i;
			break;
		}
	if( start == reorder_index.size() ){	//	no need to reorder
		return vertex.size();
	}

	//	update vertex
	int remove_count = 0;
	std::vector<Vec3<T>>	old_vertex;
	old_vertex.assign( vertex.begin()+start, vertex.end() );
	for(int i=start; i<reorder_index.size(); i++){
		if(reorder_index[i] >= 0)
			vertex[reorder_index[i]] = old_vertex[i-start];
		else
			remove_count ++;
	}
	vertex.resize(vertex.size() - remove_count);

	//	update face
	for(int i=0; i<face.size(); i++){
		if(face[i].x >= start)
			face[i].x = reorder_index[face[i].x];
		if(face[i].y >= start)
			face[i].y = reorder_index[face[i].y];
		if(face[i].z >= start)
			face[i].z = reorder_index[face[i].z];
	}

	return vertex.size();
}

/**
	Get index of isolated vertices from a mesh

	isolated vertex is a vertex that is not contained in any face

	\param	isolated_vertex_index	return the index of all isolated vertices in sequence
	\param	vertex_number			number of vertices
	\param	face					face list
*/
template<typename T>
int getIsolatedVerticesFromMesh(std::vector<int>&			isolated_vertex_index,
								int							vertex_number,
								const std::vector<int3>&	face){
	isolated_vertex_index.clear();

	//	check each vertex whether referenced in face
	std::vector<int>	reference;
	reference.resize(vertex_number, 0);

	for(int f=0; f<face.size(); f++){
		reference[ face[f].x ] = 1;
		reference[ face[f].y ] = 1;
		reference[ face[f].z ] = 1;
	}

	//	add isolated vertices to list
	for(int i=0; i<reference.size(); i++){
		if(!reference[i])
			isolated_vertex_index.push_back( i );
	}

	return isolated_vertex_index.size();
}

/**
	Remove Isolated Vertices from Mesh, keeping original vertex order

	Isolated vertices are vertices appear in vertex list,
	but never referenced in face. We remove those vertices
	by this function

	this function don't change the original vertex relative position
	in the array, thus space coherence is preserved

	\param	vertex		vertex list
	\param	face		face list
	\return				number of vertices removed
*/
template<typename T>
int removeIsolatedVerticesFromMesh(std::vector<Vec3<T>>&	vertex,
								   std::vector<int3>&		face){
	//	check each vertex whether referenced in face
	std::vector<int>	reference;
	reference.resize(vertex.size(), -1);

	for(int f=0; f<face.size(); f++){
		reference[ face[f].x ] = 1;
		reference[ face[f].y ] = 1;
		reference[ face[f].z ] = 1;
	}

	//	calculate reference position
	int count = 0;
	for(int i=0; i<reference.size(); i++){
		if(reference[i] == 1){
			reference[i] = count;
			count ++;
		}
	}

	if( count == vertex.size() )	//	no need to remove any vertex
		return 0;

	reorderVertices(vertex, face, reference);

	return reference.size() - count;
}

/**
	Remove Identical Duplicate Vertices from Mesh, keeping original vertex order

	Only vertices that is exactly the same (all data bits) will be merged

	\param	vertex		vertex list
	\param	face		face list
	\return				number of remaining vertices
*/
template<typename T>
int removeDuplicateVerticesFromMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face) 
{
	//	create hash table for vertices according to their exact location
	std::unordered_map<Vec3<T>, int, utils::BitwiseHasher<Vec3<T>>> vertex_hash;
	vertex_hash.reserve(vertex.size());

	//	mapping from old vertices to new vertices
	std::vector<int> mapping;
	mapping.resize(vertex.size());

	//	parse all vertices
	for (int i = 0; i<vertex.size(); i++) {
		//	map -0.0 to 0.0
		if (vertex[i].x == 0)	vertex[i].x = 0;
		if (vertex[i].y == 0)	vertex[i].y = 0;
		if (vertex[i].z == 0)	vertex[i].z = 0;

		//	check whether this vertex already exist
		//std::unordered_map<Vec3<T>, int, utils::BitwiseHasher<Vec3<T>>>::const_iterator iter = vertex_hash.find(vertex[i]);
		auto iter = vertex_hash.find(vertex[i]);

		if (iter == vertex_hash.end()) {	//	this is a new vertex
			mapping[i] = vertex_hash.size();
			vertex_hash.insert(std::pair<Vec3<T>, int>(vertex[i], mapping[i]));
			vertex[mapping[i]] = vertex[i];
		}
		else {		//	vertex already exist, merge to existing
			mapping[i] = iter->second;
		}
	}

	//	update face
	for (int i = 0; i<face.size(); i++) {
		face[i].x = mapping[face[i].x];
		face[i].y = mapping[face[i].y];
		face[i].z = mapping[face[i].z];
	}

	vertex.resize(vertex_hash.size());

	return vertex.size();
}

/**
	Remove Duplicate Vertices from Mesh, keeping original vertex order

	After calling this function, no vertices closer than threshold will exist. 
	Degenerate triangles will be removed.

	\param	vertex		vertex list
	\param	face		face list
	\param	threshold	threshold of vertice distance, smaller than which the vertices will be merged
	\return				number of remaining vertices
*/
template<typename T>
int removeDuplicateVerticesFromMesh(
	std::vector<Vec3<T>>&	vertex, 
	std::vector<int3>&		face,
	T						threshold) 
{
	if (vertex.empty())
		return 0;
	T squ_thre = threshold * threshold;

	//	calculate minimal vertex coordinate
	Vec3<T> vertex_min = vertex[0];
	for (int i = 1; i < vertex.size(); i++) {
		if (vertex[i].x < vertex_min.x)		vertex_min.x = vertex[i].x;
		if (vertex[i].y < vertex_min.y)		vertex_min.y = vertex[i].y;
		if (vertex[i].z < vertex_min.z)		vertex_min.z = vertex[i].z;
	}

	//	create hash table
	std::unordered_multimap<Vec3i, int, utils::BitwiseHasher<Vec3i>> vertex_hash;
	vertex_hash.reserve(vertex.size());
	std::vector<int> mapping;
	mapping.resize(vertex.size());

	//	neighbor voxel offset, from closest to farest
	yz::Vec3i neighbor_offset[27] = {
		{ 0,0,0 },	//	center voxel
		{ -1,0,0 },{ 1,0,0 },{ 0,-1,0 },{ 0,1,0 },{ 0,0,-1 },{ 0,0,1 },	//	face adjacent voxel
		{ 0,1,1 },{ 0,1,-1 },{ 0,-1,1 },{ 0,-1,-1 },	//	edge adjacent voxel
		{ 1,0,1 },{ 1,0,-1 },{ -1,0,1 },{ -1,0,-1 },
		{ 1,1,0 },{ 1,-1,0 },{ -1,1,0 },{ -1,-1,0 },
		{ 1,1,1 },{ -1,1,1 },{ 1,-1,1 },{ 1,1,-1 },{ -1,-1,1 },{ -1,1,-1 },{ 1,-1,-1 },{ -1,-1,-1 }	//	corner adjacent voxel
	};

	//	parse all vertices	
	for (int i = 0; i < vertex.size(); i++) {
		Vec3i idx = (vertex[i] - vertex_min) / threshold;
		int detected = 0;

		//	search vertex neighborhood to detect a close vertex
		for (int j = 0; !detected && j < 27; j++) {
			std::pair<
				std::unordered_map<Vec3i, int, utils::BitwiseHasher<Vec3i>>::iterator,
				std::unordered_map<Vec3i, int, utils::BitwiseHasher<Vec3i>>::iterator>
				range = vertex_hash.equal_range(idx + neighbor_offset[j]);
			while (!detected && range.first != range.second) {
				Vec3<T> nei_v = vertex[range.first->second];
				if ((vertex[i] - nei_v).SquareLength() < squ_thre) {
					mapping[i] = range.first->second;
					detected = 1;
				}
				range.first++;
			}
		}

		//	if no close vertices are detected, 
		if (!detected) {
			int curr_vertex_num = vertex_hash.size();
			mapping[i] = curr_vertex_num;
			vertex[curr_vertex_num] = vertex[i];
			vertex_hash.insert(std::pair<Vec3i, int>(idx, curr_vertex_num));
		}		
	}

	//	update face, remove degenerate faces
	int face_count = 0;
	for(int i=0; i<face.size(); i++){
		face[i].x = mapping[face[i].x];
		face[i].y = mapping[face[i].y];
		face[i].z = mapping[face[i].z];

		if (face[i].x != face[i].y && face[i].x != face[i].z && face[i].y != face[i].z)
			face[face_count++] = face[i];
	}

	vertex.resize(vertex_hash.size());
	face.resize(face_count);

	return vertex.size();
}


///@}


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_REPAIR_H__