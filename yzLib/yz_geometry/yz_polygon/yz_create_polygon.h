/***********************************************************/
/**	\file
	\brief		Create Polygon
	\details	create polygon by various methods
	\author		Yizhong Zhang
	\date		1/17/2012
*/
/***********************************************************/
#ifndef __YZ_CREATE_POLYGON_H__
#define __YZ_CREATE_POLYGON_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"

namespace yz{  namespace geometry{  namespace polygon{

/**
	create polygon consist vertices and lines connecting them by a list of vertices

	\param	vertex			return the vertex of the created polygon
	\param	edge			return the edge of the created polygon
	\param	vertex_list		the vertex list that will become the polygon vertex sequentially
	\return					number of edges created
*/
template<typename T>
int createPolygonFromVertexList(std::vector<Vec2<T>>&		vertex, 
								std::vector<int2>&			edge,
								const std::vector<Vec2<T>>&	vertex_list){
	//	copy vertex
	vertex = vertex_list;
	
	//	create edge
	if( vertex.size() < 2 )
		edge.clear();
	else{
		edge.resize(vertex.size());
		for(int i=0; i<edge.size()-1; i++){
			edge[i].x = i;
			edge[i].y = i+1;
		}
		edge.back().x = vertex.size()-1;
		edge.back().y = 0;
	}

	return edge.size();
}

/**
	create polygon as the boundary of mesh

	\param	vertex			return vertex of the polygon
	\param	edge			return edge of the polygon
	\param	mesh_vertex		vertex of the 2d mesh
	\param	mesh_face		face of the 2d mesh
	\return					number of edges in the polygon
*/
template<typename T>
int createPolygonAs2DMeshBoundary(std::vector<Vec2<T>>&			vertex, 
								  std::vector<int2>&			edge,
								  const std::vector<Vec2<T>>&	mesh_vertex,
								  const std::vector<int3>&		mesh_face){
	//	get boundary edges of mesh
	TriMeshEdge mesh_edge;
	mesh_edge.CreateEdge(mesh_face);
	getBoundaryEdges(edge, mesh_edge.edge, mesh_face);

	//	calculate vertex mapping from mesh to polygon (we don't need extra vertices)
	std::vector<int> mapping;
	mapping.resize(mesh_vertex.size(), -1);
	for(int e=0; e<edge.size(); e++){
		mapping[edge[e].x] = 0;
		mapping[edge[e].y] = 0;
	}

	//	set vertex according to mapping
	vertex.clear();
	int count = 0;
	for(int i=0; i<mapping.size(); i++){
		if( mapping[i] != -1 ){
			mapping[i] = count;
			vertex.push_back( mesh_vertex[i] );
			count ++;
		}
	}

	//	reset edge end point index
	for(int e=0; e<edge.size(); e++){
		edge[e].x = mapping[edge[e].x];
		edge[e].y = mapping[edge[e].y];
	}

	return edge.size();
}

}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_CREATE_POLYGON_H__