/***********************************************************/
/**	\file
	\brief		Match two meshes
	\author		Yizhong Zhang
	\date		11/5/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_MATCHING_H__
#define __YZ_MESH_MATCHING_H__

#include <iostream>
#include <vector>
#include <deque>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_create_topology.h"
#include "yzLib/yz_utils/yz_reorder.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Calculate Key Points on Mesh
*/
//	========================================

/**
	calculate a unique key for each vertex on a mesh
	if a face and the corresponding vertex sequence is specified.

	The algorithm simulates the process of finding correspondence
	by hand from a single triangle.

	This function is used to match two distinct meshes. For one 
	mesh, if a face and the vertex sequence of that face is given,
	a unique number of each vertex of the mesh will be generated.
	The value depend on the topology only. It is not affected by 
	sequence of vertices and faces. So if we run this function on 
	two meshes, corresponding points can be found.

	This function is not limited to the whole mesh. If the mesh is
	seperated into several pieces, only keys on the piece containing the
	selected face will be created. The return value is the number of 
	vertices with keys. Connecting piece is defined by face connecting
	in this case.

	\param	vertex_key			return the key value of each vertex, size must be equal to vertex_number
	\param	vertex_number		number of vertices on the mesh
	\param	face				face list of the mesh
	\param	face_id				the selected face id
	\param	vertex_sequence		vertex sequence of the selected face, 
								must be the three vertices of face_id
	\param	start_key_value		the value that the key starts
	\return						number of vertices that created keys.
*/
inline int createSequentialFullKeysByGivenTriangle(std::vector<int>&		vertex_key,
												   int						vertex_number,
												   const std::vector<int3>&	face,
												   int						face_id,
												   int3						vertex_sequence,
												   int						start_key_value = 0 ){
	if( vertex_number <= 0 || face.empty() )
		return 0;
	if( vertex_key.size() != vertex_number ){
		//	this is a hard constraint, if it doesn't satisfy, we just stop the program
		std::cout << "error:createSequentialFullKeysByGivenTriangle, vertex_key.size() must equal to vertex_number " << std::endl;
		exit(0);
	}

	//	create topology
	std::vector<int2>	edge;
	std::vector<int2>	ef;
	std::vector<int3>	fe;
	int edge_number = createEdgeEFFEFromFace(edge, ef, fe, face);

	//	setup keys
	int face_number = face.size();
	std::vector<int>	face_referenced;
	std::vector<int>	vertex_referenced;
	face_referenced.resize(face_number, 0);
	vertex_referenced.resize(vertex_number, 0);
	int curr_v_key_val = start_key_value;
	vertex_key[vertex_sequence.x] = curr_v_key_val++;
	vertex_key[vertex_sequence.y] = curr_v_key_val++;
	vertex_key[vertex_sequence.z] = curr_v_key_val++;
	vertex_referenced[vertex_sequence.x] = 1;
	vertex_referenced[vertex_sequence.y] = 1;
	vertex_referenced[vertex_sequence.z] = 1;

	//	setup queue
	int count = 0;
	count ++;
	face_referenced[face_id] = 1;
	std::deque<int>	face_queue;
	face_queue.push_back(face_id);

	//	all the faces in the queue, its three corners are already keyed
	while(!face_queue.empty()){
		int f_id = face_queue.front();
		face_queue.pop_front();

		//	we start from the edge, whose opposite vertex has the minimal key value
		int index[3], vk[3] = {vertex_key[face[f_id][0]], vertex_key[face[f_id][1]], vertex_key[face[f_id][2]]};
		yz::utils::sortIndex(index, vk, vk+3);

		for(int i=0; i<3; i++){
			int e_id = fe[f_id][index[i]];	//	sequence of the edge is in index[]
			int nei_f = ef[e_id][0];
			if( f_id == nei_f )	
				nei_f = ef[e_id][1];
			if( nei_f<0 || face_referenced[nei_f] )	//	this edge is boundary, or this 
				continue;							//	face is already referenced, vertex must been keyed
			else{
				//	the face is not referenced, then push it to the queue
				int third_v = getThirdVertexOfFace(face[nei_f], edge[e_id][0], edge[e_id][1]);
				if( !vertex_referenced[third_v] ){	//	if the corresponding vertex is not keyed, key it
					vertex_key[third_v] = curr_v_key_val++;
					vertex_referenced[third_v] = 1;
				}

				count ++;
				face_referenced[nei_f] = 1;
				face_queue.push_back(nei_f);
			}
		}
	}

	return curr_v_key_val - start_key_value;
}


///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_MATCHING_H__