/***********************************************************/
/**	\file
	\brief		Mesh Utility Functions
	\details	a gallary of mesh utility functions
	\author		Yizhong Zhang
	\date		10/18/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_UTILS_H__
#define __YZ_MESH_UTILS_H__

#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <unordered_map>
#include <set>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Dihedral Angle
*/
//	========================================

/**
	calculate angle of adjacent triangles in radius

	v1, v2 are wing vertices; v3, v4 are common vertices. 
	refer to Figure 1. in paper : 
	Simulation of Clothing with Folds and Wrinkles

		  v4			\n
		 /|\			\n
	v1	/ | \ v2		\n
		\ | /			\n
		 \|/			\n
		  v3

	\param	v1		vertex 1
	\param	v2		vertex 2
	\param	v3		vertex 3, on common edge of the two faces
	\param	v4		vertex 4, on common edge of the two faces
	\return			the inner angle of the two triangles, ranging 0 - 2PI
*/
template<typename T>
inline T calculateAngleBetweenTrianglesInRad(Vec3<T> v1, Vec3<T> v2, Vec3<T> v3, Vec3<T> v4){
	Vec3<T>	n1 = cross(v3-v1, v4-v1).Normalize();
	Vec3<T>	n2 = cross(v3-v2, v4-v2).Normalize();
	T cos_theta = dot(n1, n2);

	assert(cos_theta<=1 && cos_theta>=-1);	//	if this doesn't satisfy, the result will be NaN

	if( dot(n1, v2-v3) < 0 )	//	determine whether the angle is smaller or bigger than 180 degrees
		return acos( cos_theta );
	else
		return 2*T(YZ_PI)-acos( cos_theta );
}


/**
	calculate angle of adjacent triangles in degree

	v1, v2 are wing vertices; v3, v4 are common vertices. 
	refer to Figure 1. in paper : 
	Simulation of Clothing with Folds and Wrinkles

		  v4			\n
		 /|\			\n
	v1	/ | \ v2		\n
		\ | /			\n
		 \|/			\n
		  v3

	\param	v1		vertex 1
	\param	v2		vertex 2
	\param	v3		vertex 3, on common edge of the two faces
	\param	v4		vertex 4, on common edge of the two faces
	\return			the inner angle of the two triangles, ranging 0 - 360
*/
template<typename T>
inline T calculateAngleBetweenTrianglesInDeg(Vec3<T> v1, Vec3<T> v2, Vec3<T> v3, Vec3<T> v4){
	return calculateAngleBetweenTrianglesInRad(v1, v2, v3, v4) * 180 / T(YZ_PI);
}


///@}

//	========================================
///@{
/**	@name Extract Boundary Curves
*/
//	========================================

/**
	Get all edges that only has one adjacent face.
	An edge is represeted as int2: (x, y), x < y

	\param	boundary_edge	return the boundary edge list
	\param	face			input the face
	\return					number of boundary edges
*/
inline int getBoundaryEdges(
	std::vector<int2>&			boundary_edges,
	const std::vector<int3>&	face
) {
	boundary_edges.clear();

	//	we parse all edges two times, and record the appearance of each edge using a hash table
	std::unordered_map<Vec2i, int, utils::BitwiseHasher<Vec2i>> edge_count_hash;
	edge_count_hash.reserve(face.size() * 3);

	//	first parse, count the appearance of each edge
	for (int i = 0; i < face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			//	get the edge, index of the first vertex should be smaller
			Vec2i e(face[i][j], face[i][(j + 1) % 3]);
			if (e.x > e.y)
				mySwap(e.x, e.y);

			//	find the edge in the hash
			auto iter = edge_count_hash.find(e);
			if (iter != edge_count_hash.end())	//	if already exist, increase the count
				iter->second++;
			else								//	if not exist, insert it
				edge_count_hash.insert(std::pair<Vec2i, int>(e, 1));
		}
	}

	//	second parse, extract boundary edges
	for (int i = 0; i < face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			//	get the edge, index of the first vertex should be smaller
			Vec2i e(face[i][j], face[i][(j + 1) % 3]);
			if (e.x > e.y)
				mySwap(e.x, e.y);

			//	find the edge in the hash
			auto iter = edge_count_hash.find(e);
			if (iter->second == 1) {	//	if the edge appear only once, it is a boundary edge
				boundary_edges.push_back(iter->first);
				iter->second++;			//	increase the count, so it can only be recorded once
			}
		}
	}

	return boundary_edges.size();
}

/**
	Get all directed edges that only has one adjacent face.
	An edge is represeted as int2: (x, y). x y is the same sequence in the face.

	\param	boundary_edges_directed		return the directed boundary edge list
	\param	face						input the face
	\return								number of boundary edges
*/
inline int getBoundaryEdgesDirected(
	std::vector<int2>&			boundary_edges_directed,
	const std::vector<int3>&	face
) {
	boundary_edges_directed.clear();

	//	we parse all edges two times, and record the appearance of each edge using a hash table
	std::unordered_map<Vec2i, int, utils::BitwiseHasher<Vec2i>> edge_count_hash;
	edge_count_hash.reserve(face.size() * 3);

	//	first parse, count the appearance of each edge
	for (int i = 0; i < face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			//	get the directed edge
			Vec2i e(face[i][j], face[i][(j + 1) % 3]);

			//	find the edge in the hash
			auto iter = edge_count_hash.find(e);
			if (iter != edge_count_hash.end())	//	if already exist, increase the count
				iter->second++;
			else {								//	if not exist, we check whether the inverse direction exist
				iter = edge_count_hash.find(Vec2i(e.y, e.x));
				if (iter != edge_count_hash.end())	//	if the inverse exist, increase the count
					iter->second++;
				else								//	else, insert the edge (original direction)
					edge_count_hash.insert(std::pair<Vec2i, int>(e, 1));
			}
		}
	}

	//	second parse, extract directed boundary edges
	for (int i = 0; i < face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			//	get the edge, index of the first vertex should be smaller
			Vec2i e(face[i][j], face[i][(j + 1) % 3]);

			//	find the edge (or inverse edge) in the hash
			auto iter = edge_count_hash.find(e);
			if (iter == edge_count_hash.end())
				iter = edge_count_hash.find(Vec2i(e.y, e.x));

			if (iter->second == 1) {	//	if the edge appear only once, it is a boundary edge
				boundary_edges_directed.push_back(iter->first);
				iter->second++;			//	increase the count, so it can only be recorded once
			}
		}
	}

	return boundary_edges_directed.size();
}

/**
	get the indices of all boundary edges

	\param	boundary_edges		all the edge indices that are boundary
	\param	edge				edges of the mesh
	\param	face				faces of the mesh
	\return						number of boundary edges
*/
inline int getBoundaryEdges(std::vector<int>&			boundary_edges,
							const std::vector<int2>&	edge,
							const std::vector<int3>&	face){
	std::vector<int>	edge_mark;
	markBoundaryEdge(edge_mark, edge, face);
	
	assert(edge_mark.size() == edge.size());
	boundary_edges.clear();

	for(int i=0; i<edge_mark.size(); i++){
		if( edge_mark[i] )
			boundary_edges.push_back(i);
	}

	return boundary_edges.size();
}

/**
	get all boundary edges

	\param	boundary_edges		all the edge that are boundary
	\param	edge				edges of the mesh
	\param	face				faces of the mesh
	\return						number of boundary edges
*/
inline int getBoundaryEdges(std::vector<int2>&			boundary_edges,
							const std::vector<int2>&	edge,
							const std::vector<int3>&	face){
	std::vector<int>	boundary_edge_indices;
	getBoundaryEdges(boundary_edge_indices, edge, face);

	boundary_edges.clear();
	for(int i=0; i<boundary_edge_indices.size(); i++){
		boundary_edges.push_back( edge[boundary_edge_indices[i]] );
	}

	return boundary_edges.size();
}

/**
	Group directed edges into curves. 
	Curves are segmented by joint vertices (vertex connecting more than 2 edges), or form a loop.

	\param	curve_edgeIdx		output:	edges of each curve
	\param	curve_loop_flag		output:	whether the corresponding curve is a closed loop
	\param	directed_edge		input:	directed edges
	\return						number of detected curves
*/
inline int getCurveEdgesFromDirectedEdges(
	std::vector<std::vector<int>>&		curve_edgeIdx,
	std::vector<int>&					curve_loop_flag,
	const std::vector<int2>&			directed_edge
) {
	curve_edgeIdx.clear();
	curve_loop_flag.clear();
	if (directed_edge.empty())
		return 0;
	
	std::unordered_multimap<int, int>	edge_hash;			//	nodeID, edgeID
	std::unordered_map<int, int2>		node_degree_hash;	//	nodeID, (in, out)
	std::set<int>						joint_node;			//	nodes that (in, out) is not (1, 1)
	edge_hash.reserve(directed_edge.size());
	node_degree_hash.reserve(directed_edge.size() * 2);

	//	create edge hash and degree of each node
	for (unsigned int i = 0; i != directed_edge.size(); i++) {
		//	record this edge
		edge_hash.insert(std::pair<int, int>(directed_edge[i].x, i));

		//	record degree of each node
		auto iter0 = node_degree_hash.find(directed_edge[i].x);
		auto iter1 = node_degree_hash.find(directed_edge[i].y);

		if (iter0 == node_degree_hash.end())
			node_degree_hash.insert(std::pair<int, int2>(directed_edge[i].x, int2(0, 1))).first;
		else
			iter0->second.y++;

		if (iter1 == node_degree_hash.end())
			node_degree_hash.insert(std::pair<int, int2>(directed_edge[i].y, int2(1, 0))).first;
		else
			iter1->second.x++;
	}

	//	find all joint nodes
	for (unsigned int i = 0; i != directed_edge.size(); i++) {
		auto iter0 = node_degree_hash.find(directed_edge[i].x);
		auto iter1 = node_degree_hash.find(directed_edge[i].y);

		//	if degree of in/out is not 1, it is a joint vertex
		if (iter0->second != int2(1, 1))
			joint_node.insert(iter0->first);
		if (iter1->second != int2(1, 1))
			joint_node.insert(iter1->first);
	}

	//	parse all joint nodes and detect curves
	for (auto joint_iter = joint_node.begin(); joint_iter != joint_node.end(); joint_iter++) {
		int node_id = *joint_iter;
		int out_count = node_degree_hash.find(node_id)->second.y;

		//	get all curves starting from this node
		for (int i = 0; i < out_count; i++) {
			curve_edgeIdx.resize(curve_edgeIdx.size() + 1);
			curve_loop_flag.push_back(0);

			//	find the first edge starting from this node
			auto edge_iter = edge_hash.find(node_id);
			int next_node = directed_edge[edge_iter->second].y;
			curve_edgeIdx.back().push_back(edge_iter->second);
			edge_hash.erase(edge_iter);

			//	loop until the curve stop at the next joint
			while (joint_node.find(next_node) == joint_node.end()) {	//	not a joint node, it is a normal node with 1 out degree
				edge_iter = edge_hash.find(next_node);
				next_node = directed_edge[edge_iter->second].y;
				curve_edgeIdx.back().push_back(edge_iter->second);
				edge_hash.erase(edge_iter);
			}
		}
	}

	//	parse all remaining edges and detect curves
	for (unsigned int i = 0; i != directed_edge.size(); i++) {
		auto edge_iter = edge_hash.find(directed_edge[i].x);
		if (edge_iter == edge_hash.end())	//	this edge is already removed, skip
			continue;
		curve_edgeIdx.resize(curve_edgeIdx.size() + 1);
		curve_loop_flag.push_back(0);

		int next_node = directed_edge[edge_iter->second].y;
		curve_edgeIdx.back().push_back(edge_iter->second);
		edge_hash.erase(edge_iter);

		//	loop until the curve stop at the next joint, or forms a loop
		while (joint_node.find(next_node) == joint_node.end()) {	//	joint vertex
			edge_iter = edge_hash.find(next_node);
			if (edge_iter == edge_hash.end()) {						//	a loop
				curve_loop_flag.back() = 1;							//	we label it as loop
				break;
			}
			next_node = directed_edge[edge_iter->second].y;
			curve_edgeIdx.back().push_back(edge_iter->second);
			edge_hash.erase(edge_iter);
		}
	}

	return curve_edgeIdx.size();
}

/**
	Group directed edges into curves.
	Curves are segmented by joint vertices (vertex connecting more than 2 edges), or form a loop.

	\param	curve_edgeIdx		output:	edges of each curve
	\param	curve_loop_flag		output:	whether the corresponding curve is a closed loop
	\param	directed_edge		input:	directed edges
	\return						number of detected curves
*/
inline int getCurveEdgesFromDirectedEdges(
	std::vector<std::vector<int2>>&		curve_edge,
	std::vector<int>&					curve_loop_flag,
	const std::vector<int2>&			directed_edge
) {
	std::vector<std::vector<int>>		curve_edgeIdx;
	getCurveEdgesFromDirectedEdges(curve_edgeIdx, curve_loop_flag, directed_edge);

	curve_edge.clear();
	curve_edge.resize(curve_edgeIdx.size());

	//	parse all curves
	for (unsigned int curve_id = 0; curve_id != curve_edgeIdx.size(); curve_id++) {
		curve_edge[curve_id].reserve(curve_edgeIdx[curve_id].size());

		//	parse the curve and record vertices
		for (unsigned int i = 0; i != curve_edgeIdx[curve_id].size(); i++) {
			int edge_idx = curve_edgeIdx[curve_id][i];
			curve_edge[curve_id].push_back(directed_edge[edge_idx]);
		}
	}

	return curve_edge.size();
}

/**
	Get boundary edge curves from a mesh

	\param	curve_edge			detected curve edges
	\param	curve_loop_flag		whether each curve form a loop
	\param	face				the input mesh face
	\return						number of detected curves
*/
inline int getBoundaryCurveEdgesFromMesh(
	std::vector<std::vector<int2>>&		curve_edge,
	std::vector<int>&					curve_loop_flag,
	const std::vector<int3>&			face
) {
	std::vector<int2> directed_edge;

	getBoundaryEdgesDirected(directed_edge, face);

	return getCurveEdgesFromDirectedEdges(curve_edge, curve_loop_flag, directed_edge);
}

/**
	Get boundary vertex curves from a mesh

	\param	curve_vertex		detected curve vertices
	\param	curve_loop_flag		whether each curve form a loop
	\param	face				the input mesh face
	\return						number of detected curves
*/
inline int getBoundaryCurveVerticesFromMesh(
	std::vector<std::vector<int>>&		curve_vertex,
	std::vector<int>&					curve_loop_flag,
	const std::vector<int3>&			face
) {
	std::vector<std::vector<int2>> curve_edge;
	getBoundaryCurveEdgesFromMesh(curve_edge, curve_loop_flag, face);

	curve_vertex.resize(curve_edge.size());

	for (int i = 0; i < curve_vertex.size(); i++) {
		//	push the first vertex
		if (curve_loop_flag[i])
			curve_vertex[i].reserve(curve_edge[i].size());
		else {
			curve_vertex[i].reserve(curve_edge[i].size() + 1);
			curve_vertex[i].push_back(curve_edge[i][0].x);
		}

		//	push the remaining vertices
		for (int j = 0; j < curve_edge[i].size(); j++)
			curve_vertex[i].push_back(curve_edge[i][j].y);
	}

	return curve_vertex.size();
}

/**
	this is just a simple version, very unstable

	\todo	this function need re-implement, no duplicate vertices should be allowed in one curve
*/
inline int getBoundaryEdgeCurves(std::vector<std::vector<int>>&	edge_boundary,
								 const std::vector<int2>&		edge,
								 const std::vector<int>&		boundary_edge_flag,
								 const std::vector<int>&		vv,
								 const std::vector<int>&		ve,
								 const std::vector<int>&		vve_start ){
	std::vector<int> edge_referenced_flag;
	edge_referenced_flag.resize(edge.size(), 0);

	//	scan each edge
	for(int e_id = 0; e_id<edge.size(); e_id++){
		//	if we got a boundary edge and it is not referenced yet, 
		if( boundary_edge_flag[e_id] && !edge_referenced_flag[e_id] ){
			edge_boundary.resize( edge_boundary.size()+1 );	//	add one more element
			int v0 = edge[e_id].x;
			int v1 = edge[e_id].y;
			int e0 = e_id;
			int e1 = -1;
			do{
				//	find neighboring boundary edge
				e1 = -1;
				for(int i=vve_start[v1]; i<vve_start[v1+1]; i++){
					if( boundary_edge_flag[ve[i]] && ve[i] != e0 ){
						e1 = ve[i];
						break;
					}
				}
				assert( e1 != -1 );

				//	update v0 v1 e0 e1
				v0 = v1;
				v1 = (edge[e1].x==v0 ? edge[e1].y : edge[e1].x);
				e0 = e1;

				//	push back this edge
				edge_boundary.back().push_back(e1);
				edge_referenced_flag[e1] = 1;
			}while(e1 != e_id);
		}
	}

	return edge_boundary.size();
}


/**
	get the boundary vertex list of a mesh

	it is assumed that one vertex appear only once in a boundary curve,
	so if two boundary curve has common vertex, they are treated as 
	independent curves

	\param	vertex_boundary		return the boundary vertex list
	\param	edge				edge of a mesh
	\param	boundary_edge_flag	whether the edge is a boundary
	\param	vv					vertex-vertex connectivity
	\param	ve					vertex-edge connectivity
	\param	vve_start			index start of vv, v2
	\return						curve number
*/
inline int getBoundaryVertexCurves(std::vector<std::vector<int>>&	vertex_boundary,
								   const std::vector<int2>&			edge,
								   const std::vector<int>&			boundary_edge_flag,
								   const std::vector<int>&			vv,
								   const std::vector<int>&			ve,
								   const std::vector<int>&			vve_start ){
	std::vector<std::vector<int>>	edge_boundary;
	int curves = getBoundaryEdgeCurves(edge_boundary, edge, boundary_edge_flag, vv, ve, vve_start);
	vertex_boundary.clear();

	//	extract vertex from edge
	vertex_boundary.resize(edge_boundary.size());
	for(int j=0; j<edge_boundary.size(); j++){
		if( edge_boundary[j].size() < 3 ){	//	such boundary should never exist
			#ifndef	BE_QUIET
				std::cout << "error: getBoundaryVertexCurves, boundary length smaller than 3" << std::endl;
			#endif
			return 0;
		}
		else{
			for(int i=0; i<edge_boundary[j].size(); i++){
				int e_id0 = edge_boundary[j][i];
				int e_id1 = (i!=edge_boundary[j].size()-1 ? edge_boundary[j][i+1] : edge_boundary[j][0]);
				int v[4] = {edge[e_id0].x, edge[e_id0].y, edge[e_id1].x, edge[e_id1].y};
				//	find v that appear twice
				int v_id = ((v[0]==v[2] || v[0]==v[3]) ? v[0] : v[1]);
				vertex_boundary[j].push_back( v_id );
			}
		}
	}

	return curves;
}

/**
	get the boundary vertex list of a mesh

	\param	vertex_boundary		return the boundary vertex list
	\param	vertex_number		number of vertex of the mesh
	\param	face				face list of the mesh
	\return						curve number
*/
inline int getBoundaryVertexCurves(std::vector<std::vector<int>>&	vertex_boundary,
								   int								vertex_number,
								   const std::vector<int3>&			face){
	std::vector<int2>	edge;
	std::vector<int>	boundary_edge_flag;
	std::vector<int>	vv;
	std::vector<int>	ve;
	std::vector<int>	vve_start;
	createEdgeFromFace(edge, face);
	createVVEFromEdge(vv, ve, vve_start, vertex_number, edge);
	markBoundaryEdge(boundary_edge_flag, edge, face);

	//	we have to create edge boundary, then vertex boundary

	return getBoundaryVertexCurves(vertex_boundary, edge, boundary_edge_flag, vv, ve, vve_start);
}

/**
	get sub-curve from vertex curve list, including start and end point

	for curve containing the two points, there are two sub-curves, we want the shorter one

	\param	sub_curve		return the sub-curve from start_id to end_id
	\param	curve_list		the vertex curve list
	\param	start_id		sub curve start vertex index
	\param	end_id			sub curve end vertex index
	\return					vertex number in the sub-curve
*/
inline int getSubVertexCurveShorter(std::vector<int>&				sub_curve,
									std::vector<std::vector<int>>&	curve_list,
									int								start_id,
									int								end_id ){
	sub_curve.clear();
	int v0 = start_id, v1 = end_id;
	for(int j=0; j<curve_list.size(); j++){
		for(int i=0; i<curve_list[j].size(); i++){
			if( curve_list[j][i] == v0 ){
				//	find the start of the curve, then search for v1
				int k=i;
				int len0to1 = -1;
				int change_curve_flag = 0;
				while(curve_list[j][k] != v1){
					k++;
					len0to1 ++;
					if( len0to1 >= curve_list[j].size() ){
						//	v0 and v1 are not inside the same curve, then goto check the next curve
						change_curve_flag = 1;
						break;
					}
					if( k==curve_list[j].size() )
						k = 0;
				}
				if( change_curve_flag )
					break;	//	break for(i...), go to the next curve
				//	find the shorter path
				int len1to0 = curve_list[j].size() - len0to1;
				if( len0to1 < len1to0 ){
					while( i!=k ){
						sub_curve.push_back( curve_list[j][i] );
						i++;
						if( i==curve_list[j].size() )
							i = 0;
					}
					sub_curve.push_back(end_id);
					return sub_curve.size();
				}
				else{
					while( i!=k ){
						sub_curve.push_back( curve_list[j][i] );
						i--;
						if( i==-1 )
							i = curve_list[j].size()-1;
					}
					sub_curve.push_back(end_id);
					return sub_curve.size();
				}
			}
		}
	}

	return 0;
}


/**
	get sub-curve from vertex curve list, including start and end point

	for curve containing the two points, there are two sub-curves, we want the longer one

	\param	sub_curve		return the sub-curve from start_id to end_id
	\param	curve_list		the vertex curve list
	\param	start_id		sub curve start vertex index
	\param	end_id			sub curve end vertex index
	\return					vertex number in the sub-curve
*/
inline int getSubVertexCurveLonger(std::vector<int>&				sub_curve,
								   std::vector<std::vector<int>>&	curve_list,
								   int								start_id,
								   int								end_id ){
	sub_curve.clear();
	int v0 = start_id, v1 = end_id;
	for(int j=0; j<curve_list.size(); j++){
		for(int i=0; i<curve_list[j].size(); i++){
			if( curve_list[j][i] == v0 ){
				//	find the start of the curve, then search for v1
				int k=i;
				int len0to1 = -1;
				int change_curve_flag = 0;
				while(curve_list[j][k] != v1){
					k++;
					len0to1 ++;
					if( len0to1 >= curve_list[j].size() ){
						//	v0 and v1 are not inside the same curve, then goto check the next curve
						change_curve_flag = 1;
						break;
					}
					if( k==curve_list[j].size() )
						k = 0;
				}
				if( change_curve_flag )
					break;	//	break for(i...), go to the next curve
				//	find the shorter path
				int len1to0 = curve_list[j].size() - len0to1;
				if( len0to1 >= len1to0 ){	//	the only difference from short function is this sign
					while( i!=k ){
						sub_curve.push_back( curve_list[j][i] );
						i++;
						if( i==curve_list[j].size() )
							i = 0;
					}
					sub_curve.push_back(end_id);
					return sub_curve.size();
				}
				else{
					while( i!=k ){
						sub_curve.push_back( curve_list[j][i] );
						i--;
						if( i==-1 )
							i = curve_list[j].size()-1;
					}
					sub_curve.push_back(end_id);
					return sub_curve.size();
				}
			}
		}
	}

	return 0;
}


///@}

//	========================================
///@{
/**	@name Detect Topology Issues
*/
//	========================================

/**
	get all faces that has 2 or 3 itendital vertices

	\param	thin_face		return the boundary edge list
	\param	face			input the face
	\return					number of boundary edges
*/
inline int getThinFaces(
	std::vector<int3>&			thin_face,
	const std::vector<int3>&	face
) {
	thin_face.clear();

	for (int i = 0; i < face.size(); i++) {
		if (face[i].x == face[i].y || face[i].x == face[i].z || face[i].y == face[i].z)
			thin_face.push_back(face[i]);
	}

	return thin_face.size();
}

/**
	get all duplicate face pairs

	\param	duplicate_face_pair		return the duplicate face pairs
	\param	face					input the face
	\return							number of pairs
*/
inline int getDuplicateFacePairs(
	std::vector<int2>&			duplicate_face_pair,
	const std::vector<int3>&	face
) {
	duplicate_face_pair.clear();

	std::unordered_map<Vec3i, int, utils::BitwiseHasher<Vec3i>> face_hash;
	face_hash.reserve(face.size());

	//	parse all faces
	for (int i = 0; i < face.size(); i++) {
		//	reorder vertex indices in increasing order
		Vec3i f(face[i]);
		if (f.x > f.y)
			mySwap(f.x, f.y);
		if (f.x > f.z)
			mySwap(f.x, f.z);
		if (f.y > f.z)
			mySwap(f.y, f.z);

		//	check whether 
		std::unordered_map<Vec3i, int, utils::BitwiseHasher<Vec3i>>::iterator iter = face_hash.find(f);
		if (iter != face_hash.end())
			duplicate_face_pair.push_back(int2(iter->second, i));
		else
			face_hash.insert(std::pair<Vec3i, int>(f, i));
	}

	return duplicate_face_pair.size();
}

///@}


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_UTILS_H__