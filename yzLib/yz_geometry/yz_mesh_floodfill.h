/***********************************************************/
/**	\file
	\brief		Perform Mesh Floodfill
	\author		Yizhong Zhang
	\date		9/22/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_FLOODFILL_H__
#define __YZ_MESH_FLOODFILL_H__

#include <iostream>
#include <vector>
#include <deque>

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Vertex Floodfill
*/
//	========================================

/**
	perform floodfill on mesh vertex by vertex-vertex connectivity

	After calling this function, vertices are marked in vertex_mark
	array. 

	\param	vertex_mark			return floodfill label
	\param	seed_value			value of seed in vertex_mark
	\param	seed_vertex_index	index of the seed vertex
	\param	vertex_number		number of vertices
	\param	vv					vertex-vertex connectivity
	\param	vv_start			vv start position
	\return						how many vertices are marked in this floodfill
*/	
inline int floodfillMeshVertexByVV(int*			vertex_mark,
								   int			seed_value,
								   int			seed_vertex_index,
								   int			vertex_number,
								   const int*	vv,
								   const int*	vv_start ){
	if( seed_vertex_index >= vertex_number ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshVertexByVV, invalid seed vertex index" << std::endl;
		#endif
		return 0;
	}

	//	create a list implemented as double sided queue
	int count = 0;
	std::deque<int>	vertex_queue;
	vertex_mark[seed_vertex_index] = seed_value;
	vertex_queue.push_back(seed_vertex_index);		//	push the seed
	count ++;

	//	each loop, get the first vertex, and for each neighbor of the vertex,
	//	if it is not marked, add it to the queue, too. Until the queue is empty
	while(!vertex_queue.empty()){
		int v_id = vertex_queue.front();
		vertex_queue.pop_front();

		for(int j=vv_start[v_id]; j<vv_start[v_id+1]; j++){
			int nei_id = vv[j];

			if( vertex_mark[nei_id] != seed_value ){
				vertex_mark[nei_id] = seed_value;
				vertex_queue.push_back(nei_id);
				count ++;
			}
		}
	}

	//	when the queue is empty, the floodfill has completed, return marked vertex number
	return count;
}

/**
	perform floodfill on mesh vertex by vertex-vertex connectivity, with vector container

	After calling this function, vertices are marked in vertex_mark
	array. 

	\param	vertex_mark			return floodfill label, we don't resize the vector size in
								this function, so its size must match vertex_number
	\param	seed_value			value of seed in vertex_mark
	\param	seed_vertex_index	index of the seed vertex
	\param	vertex_number		number of vertices
	\param	vv					vertex-vertex connectivity
	\param	vv_start			vv start position
	\return						how many vertices are marked in this floodfill
*/	
inline int floodfillMeshVertexByVV(std::vector<int>&		vertex_mark,
								   int						seed_value,
								   int						seed_vertex_index,
								   int						vertex_number,
								   const std::vector<int>&	vv,
								   const std::vector<int>&	vv_start ){
	if( vertex_number <= 0 )
		return 0;
	if( vertex_mark.size()!=vertex_number || vv_start.size()-1!=vertex_number || vv.empty() ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshVertexByVV, invalid array size" << std::endl;
		#endif
		return 0;
	}

	return floodfillMeshVertexByVV( (int*)&vertex_mark[0], 
		seed_value, seed_vertex_index, vertex_number, (int*)&vv[0], (int*)&vv_start[0] );
}


///@}

//	========================================
///@{
/**	@name Face Floodfill
*/
//	========================================

/**
	perform floodfill on mesh face by edge-face connectivity

	After calling this function, faces are marked in face_mark
	array. 

	\param	face_mark			return floodfill label
	\param	seed_value			value of seed in face_mark
	\param	seed_face_index		index of the seed face
	\param	face_number			number of faces
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many faces are marked in this floodfill
*/	
inline int floodfillMeshFaceByEFFE(int*			face_mark,
								   int			seed_value,
								   int			seed_face_index,
								   int			face_number,
								   const int*	ef,
								   const int*	fe){
	if( seed_face_index >= face_number ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshFaceByEFFE, invalid seed face index" << std::endl;
		#endif
		return 0;
	}

	//	create a list implemented as double sided queue
	int count = 0;
	std::deque<int>	face_queue;
	face_mark[seed_face_index] = seed_value;
	face_queue.push_back(seed_face_index);		//	push the seed
	count ++;

	//	each loop, get the first face, and for each neighbor of the face,
	//	if it is not marked, add it to the queue, too. Until the queue is empty
	while(!face_queue.empty()){
		int f_id = face_queue.front();
		face_queue.pop_front();

		for(int i=0; i<3; i++){
			int nei_e = fe[f_id*3+i];
			int nei_f = ef[nei_e*2];
			if( nei_f == f_id )
				nei_f = ef[nei_e*2+1];
			if( nei_f != -1 && face_mark[nei_f] != seed_value ){	//	on the boundary, nei_f == -1
				face_mark[nei_f] = seed_value;
				face_queue.push_back(nei_f);
				count ++;
			}
		}
	}

	//	when the queue is empty, the floodfill has completed, return marked face number
	return count;
}

/**
	perform floodfill on mesh face by edge-face connectivity

	After calling this function, faces are marked in face_mark
	array. 

	\param	face_mark			return floodfill label
	\param	seed_value			value of seed in face_mark
	\param	seed_face_index		index of the seed face
	\param	face_number			number of faces
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many faces are marked in this floodfill
*/	
inline int floodfillMeshFaceByEFFE(std::vector<int>&		face_mark,
								   int						seed_value,
								   int						seed_face_index,
								   int						face_number,
								   const std::vector<int2>&	ef,
								   const std::vector<int3>&	fe ){
	if( face_number <= 0 )
		return 0;
	if( face_mark.size()!=face_number || fe.size()!=face_number || ef.empty() ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshFaceByEFFE, invalid array size" << std::endl;
		#endif
		return 0;
	}

	return floodfillMeshFaceByEFFE( (int*)&face_mark[0], 
		seed_value, seed_face_index, face_number, (int*)&ef[0], (int*)&fe[0] );

}

/**
	perform floodfill on mesh face with given edge boundary by edge-face connectivity

	After calling this function, faces are marked in face_mark array. 

	\param	face_mark			return floodfill label
	\param	seed_value			value of seed in face_mark
	\param	seed_face_index		index of the seed face
	\param	face_number			number of faces
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many faces are marked in this floodfill
*/	
inline int floodfillMeshFaceWithEdgeBoundaryByEFFE(int*			face_mark,
												   int			seed_value,
												   int			seed_face_index,
												   int			face_number,
												   const int*	boundary_edge_flag,
												   const int*	ef,
												   const int*	fe){
	if( seed_face_index >= face_number ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshFaceByEFFE, invalid seed face index" << std::endl;
		#endif
		return 0;
	}

	//	create a list implemented as double sided queue
	int count = 0;
	std::deque<int>	face_queue;
	face_mark[seed_face_index] = seed_value;
	face_queue.push_back(seed_face_index);		//	push the seed
	count ++;

	//	each loop, get the first face, and for each neighbor of the face,
	//	if it is not marked, add it to the queue, too. Until the queue is empty
	while(!face_queue.empty()){
		int f_id = face_queue.front();
		face_queue.pop_front();

		for(int i=0; i<3; i++){
			int nei_e = fe[f_id*3+i];
			if( boundary_edge_flag[nei_e] )	//	the edge is set to be boundary, go to the next edge
				continue;
			int nei_f = ef[nei_e*2];
			if( nei_f == f_id )
				nei_f = ef[nei_e*2+1];
			if( nei_f != -1 && face_mark[nei_f] != seed_value ){	//	on the boundary, nei_f == -1
				face_mark[nei_f] = seed_value;
				face_queue.push_back(nei_f);
				count ++;
			}
		}
	}

	//	when the queue is empty, the floodfill has completed, return marked face number
	return count;
}

/**
	perform floodfill on mesh face by edge-face connectivity

	After calling this function, faces are marked in face_mark
	array. 

	\param	face_mark			return floodfill label
	\param	seed_value			value of seed in face_mark
	\param	seed_face_index		index of the seed face
	\param	face_number			number of faces
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many faces are marked in this floodfill
*/	
inline int floodfillMeshFaceWithEdgeBoundaryByEFFE(std::vector<int>&		face_mark,
												   int						seed_value,
												   int						seed_face_index,
												   int						face_number,
												   const std::vector<int>&	boundary_edge_flag,
												   const std::vector<int2>&	ef,
												   const std::vector<int3>&	fe ){
	if( face_number <= 0 )
		return 0;
	if( face_mark.size()!=face_number || boundary_edge_flag.size()!=ef.size() || fe.size()!=face_number || ef.empty() ){
		#ifndef BE_QUIET
			std::cout << "error: floodfillMeshFaceByEFFE, invalid array size" << std::endl;
		#endif
		return 0;
	}

	return floodfillMeshFaceWithEdgeBoundaryByEFFE( (int*)&face_mark[0], 
		seed_value, seed_face_index, face_number, (int*)&boundary_edge_flag[0], (int*)&ef[0], (int*)&fe[0] );
}


///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_FLOODFILL_H__