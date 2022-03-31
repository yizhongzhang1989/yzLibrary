/***********************************************************/
/**	\file
	\brief		Algorithms of Graph
	\author		Yizhong Zhang
	\date		11/4/2012
*/
/***********************************************************/
#ifndef __YZ_GRAPH_ALGORITHM_H__
#define __YZ_GRAPH_ALGORITHM_H__

#include <vector>
#include <deque>
#include <algorithm>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
///@{
/**	@name Shortest Path
*/
//	========================================

/**
	get shortest path of a undirected unweighted graph

	this is achieved by breadth first search, similar to topological sort

	\param	vertex_path			return vertices on the path, including start and end
	\param	edge_path			return edges on the path
	\param	vertex_number		number of vertices of the graph
	\param	vv					vertex - neighbor vertex
	\param	ve					vertex - neighbor edge
	\param	vve_start			start index of vv ve
	\param	start_vertex_id		id of the start vertex
	\param	end_vertex_id		id of the end vertex
	\return						the path length, equal to edge_path.size(); if unreachable, return -1
*/
inline int getGraphShortestPathFromVV(std::vector<int>&			vertex_path,
									  std::vector<int>&			edge_path,
									  int						vertex_number,
									  const std::vector<int>&	vv,
									  const std::vector<int>&	ve,
									  const std::vector<int>&	vve_start,
									  int						start_vertex_id,
									  int						end_vertex_id ){
	//	if start and end are same vertex, return directly
	if( start_vertex_id == end_vertex_id ){
		vertex_path.resize(1);
		vertex_path[0] = start_vertex_id;
		edge_path.clear();
		return 0;
	}

	//	breadth first search
	std::vector<int>	distance;
	std::vector<int>	last_vertex;
	std::vector<int>	last_edge;

	distance.resize(vertex_number, -1);
	last_vertex.resize(vertex_number, -1);
	last_edge.resize(vertex_number, -1);

	std::deque<int>	vertex_queue;
	distance[start_vertex_id] = 0;
	vertex_queue.push_back(start_vertex_id);

	//	each loop, get the first vertex, and for each neighbor of the vertex,
	//	if it is not visited, add it to the queue. Until the queue is empty
	while(!vertex_queue.empty()){
		int v_id = vertex_queue.front();
		vertex_queue.pop_front();
		int curr_dist = distance[v_id];

		//	if it is not end vertex, we scan all its neighbors
		for(int j=vve_start[v_id]; j<vve_start[v_id+1]; j++){
			int nei_id = vv[j];

			if( distance[nei_id] == -1 ){	//	the neighbor vertex is not visited
				distance[nei_id]	= curr_dist + 1;
				last_vertex[nei_id] = v_id;
				last_edge[nei_id]	= ve[j];
				vertex_queue.push_back(nei_id);
			}

			if( nei_id == end_vertex_id ){	//	we have reached the end point
				int curr_id = nei_id;
				vertex_path.clear();
				edge_path.clear();
				vertex_path.push_back(curr_id);
				while( curr_id != start_vertex_id ){
					vertex_path.push_back(last_vertex[curr_id]);
					edge_path.push_back(last_edge[curr_id]);
					curr_id = last_vertex[curr_id];
				}
				std::reverse(vertex_path.begin(), vertex_path.end());
				std::reverse(edge_path.begin(), edge_path.end());
				return edge_path.size();
			}
		}
	}

	//	going here means there is no path between start and end
	vertex_path.clear();
	edge_path.clear();
	return -1;
}




///@}

}	//	namespace yz

#endif	//	__YZ_GRAPH_ALGORITHM_H__