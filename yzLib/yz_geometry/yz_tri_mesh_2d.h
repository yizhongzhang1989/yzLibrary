/***********************************************************/
/**	\file
	\brief		Triangle Mesh 2D
	\details	Mesh 2d is simpler than mesh 3d because mesh 2d
				don't have normal. Mesh 2D can only be displayed 
				on plane.
	\author		Yizhong Zhang
	\date		6/28/2012
*/
/***********************************************************/
#ifndef __YZ_TRI_MESH_2D__
#define __YZ_TRI_MESH_2D__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_rw.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{	namespace geometry{

/**
	base 2D triangle mesh ptr, just vertex and face pointers
*/
template<class T>
class PtrTriMesh2D{
public:
	T*		vertex_ptr;
	int*	face_ptr;

public:
	//	constructor
	PtrTriMesh2D() : vertex_ptr(NULL), face_ptr(NULL){}
	/**
		Reset all members. clear pointers, but don't touch memory
	*/
	inline void Reset(){
		vertex_ptr	= NULL;
		face_ptr	= NULL;
	}
};

/**
	2D triangle mesh using vector, just vertex and face

	2D triangle mesh is mainly used for display mesh on 2D
	canvas, so display function only show edge
*/
template<class T>
class TriMesh2D{
public:
	std::vector<Vec2<T>>	vertex;
	std::vector<int3>		face;

public:
	//	constructor
	TriMesh2D(){};

	/**
		Reset all members, clear all memory
	*/
	inline void Reset(){
		vertex.clear();
		face.clear();
	}

	/**
		Read mesh from file

		\param	file_name			file name of the 2D mesh
		\param	ignore_dimension	which dimension to ignore, 'x', 'y', 'z', default 'z'
		\return						whether read mesh succeed
	*/
	inline int ReadMeshFromFile(const char* file_name, char ignore_dimension = 'z'){
		return	readTriMeshFromFile2D(file_name, &vertex, &face, ignore_dimension );
	}

	/**
		Draw mesh edge
	*/
	inline int Display(int shading_mode=0){
		#ifdef YZ_gl_h
		opengl::drawMeshEdgeFromFace2D(vertex, face);
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in TriMesh2D" << std::endl;
			return 0;
		#endif
	}

	/**
		insert a new mesh to the old mesh, forming a single mesh
	*/
	inline void AppendTriMesh2D(TriMesh2D<T>& mesh){
		int old_vertex_number = vertex.size();
		int old_face_number = face.size();
		vertex.insert(vertex.end(), mesh.vertex.begin(), mesh.vertex.end());
		face.insert(face.end(), mesh.face.begin(), mesh.face.end());
		for( int f=old_face_number; f<face.size(); f++ ){
			face[f].x += old_vertex_number;
			face[f].y += old_vertex_number;
		}
	}

};


typedef TriMesh2D<float>	TriMesh2Df;
typedef TriMesh2D<double>	TriMesh2Dd;

}}	//	namespace yz::geometry

#endif	//	__YZ_TRI_MESH_2D__