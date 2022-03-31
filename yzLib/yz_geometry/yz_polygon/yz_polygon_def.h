/***********************************************************/
/**	\file
	\brief		Polygon
	\details	polygon consist of vertices and edges
	\author		Yizhong Zhang
	\date		1/17/2012
*/
/***********************************************************/
#ifndef __YZ_POLYGON_DEF_H__
#define __YZ_POLYGON_DEF_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{  namespace geometry{  namespace polygon{

/**
	Polygon consist vertices and lines connecting them
*/
template<class T>
class Polygon{
public:
	std::vector<Vec2<T>>	vertex;
	std::vector<int2>		edge;

public:
	//	constructor
	Polygon(){};

	/**
		Reset all members, clear all memory
	*/
	inline void Reset(){
		vertex.clear();
		edge.clear();
	}

	/**
		Draw the polygon
	*/
	inline int Display(int shading_mode=0){
		#ifdef YZ_gl_h
			opengl::drawMeshEdge2D(vertex, edge);
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in Polygon" << std::endl;
			return 0;
		#endif
	}

	/**
		insert a new polygon to the old polygon, forming a single polygon
	*/
	inline void AppendPolygon(Polygon<T>& polygon){
		int old_vertex_number	= vertex.size();
		int old_edge_number		= edge.size();
		vertex.insert(vertex.end(), polygon.vertex.begin(), polygon.vertex.end());
		edge.insert(edge.end(), polygon.edge.begin(), polygon.edge.end());
		for( int e=old_edge_number; e<edge.size(); e++ ){
			edge[e].x += old_vertex_number;
			edge[e].y += old_vertex_number;
		}
	}

};

}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_POLYGON_DEF_H__