/***********************************************************/
/**	\file
	\brief		Offset Triangle Mesh
	\author		Yizhong Zhang
	\date		10/28/2012
*/
/***********************************************************/
#ifndef __YZ_OFFSET_TRI_MESH_H__
#define __YZ_OFFSET_TRI_MESH_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_rw.h"
#include "yzLib/yz_geometry/yz_mesh_normal.h"
#include "yzLib/yz_geometry/yz_tri_mesh.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{	namespace geometry{

/**
	offset mesh of base mesh
*/
template<class T>
class OffsetTriMesh : public TriMesh<T>{
public:
	using TriMesh<T>::vertex;
	using TriMesh<T>::face;
	using TriMesh<T>::vertex_normal;
	using TriMesh<T>::face_normal;

	TriMesh<T>*	base_mesh;	///<	pointer to the base mesh
	T			offset;		///<	offset from the base mesh

public:
	OffsetTriMesh(){
		base_mesh	= NULL;
		offset		= 0;
	}

	OffsetTriMesh(TriMesh<T>& mesh, T offset_value=0){
		Setup(mesh, offset_value);
	}

	void Setup(TriMesh<T>& mesh, T offset_value=0){
		base_mesh		= &mesh;
		face			= mesh.face;
		offset			= offset_value;

		Update();
	}

	void Update(){
		if( !base_mesh )	return;

		vertex = base_mesh->vertex;
		for(int i=0; i<vertex.size(); i++){
			vertex[i] += base_mesh->vertex_normal[i] * offset;
		}
		this->CalculateNormals();
	}


};



}}	//	namespace yz::geometry

#endif	//	__YZ_OFFSET_TRI_MESH_H__