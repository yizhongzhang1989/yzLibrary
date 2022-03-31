/***********************************************************/
/**	\file
	\brief		Subdivide a Mesh
	\author		Yizhong Zhang
	\date		6/4/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_SUBD_H__
#define __YZ_MESH_SUBD_H__

#include <stdlib.h>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_create_topology.h"

namespace yz{	namespace geometry{
//	========================================
///@{
/**	@name Subdivide Triangle Mesh
*/
//	========================================
/**
	loop subd triangle mesh by vector

	This function is out of date, we used new name scheme

	Given vertex and face, subd the mesh directly.
	New vertices are placed at the center of edge

	\param	vertex		vertex vector
	\param	face		face vector
	\return				number of vertices after subd
*/
template<typename T>
inline int loopSubdTriMesh(	std::vector<Vec3<T>>&	vertex, 
							std::vector<int3>&		face){
	std::cout << "warning: loopSubdTriMesh, this function is out of date, pleace change to subdTriMeshLoop()" << std::endl;
	unsigned int old_v_num = vertex.size();
	unsigned int old_f_num = face.size();

	vertex.resize(vertex.size() + face.size()*3);
	face.resize(face.size()*4);
	int vertex_number = subdTriMeshLoop((T*)&vertex[0], old_v_num, (int*)&face[0], old_f_num);
	vertex.resize(vertex_number);

	return vertex_number;
}

/**
	loop subd triangle mesh by vector

	Given vertex and face, subd the mesh directly.
	New vertices are placed at the center of edge

	\param	vertex		vertex vector
	\param	face		face vector
	\return				number of vertices after subd
*/
template<typename T>
inline int subdTriMeshLoop(	std::vector<Vec3<T>>&	vertex, 
							std::vector<int3>&		face){
	unsigned int old_v_num = vertex.size();
	unsigned int old_f_num = face.size();

	vertex.resize(vertex.size() + face.size()*3);	//	the minimal safe size, we will resize the array later
	face.resize(face.size()*4);
	int vertex_number = subdTriMeshLoop((T*)&vertex[0], old_v_num, (int*)&face[0], old_f_num);
	vertex.resize(vertex_number);

	return vertex_number;
}

/**
	loop subd triangle mesh
	
	we perform type check inside this function, only float or double is allowed

	Given vertex and face, subd the mesh directly. 
	New vertices are placed at the center of edge

	\param	vertex			vertex list, xyz_xyz_, size unknown, 
							minimal safe size = sizeof(T)*(vertex_number*3 + face_number*9)
	\param	vertex_number	vertex number
	\param	face			face list, v0v1v2_v0v1v2_, size = sizeof(int)*(face_number*12)
	\param	face_number		face number
	\return					number of vertices after subd
*/
template<typename T>
inline int subdTriMeshLoop(T* vertex, unsigned int vertex_number, int* face, unsigned int face_number){
	//	check type of T, we only accept float or double
	if( !Is_float<T>::check_type && !Is_double<T>::check_type ){
		#ifndef	BE_QUIET
		std::cout << "error: subdTriMeshLoop, unsupported type, only accept float or double" << std::endl;
		#endif
		return 0;
	}

	//	memory check
	if( vertex==NULL || vertex_number==0 || face==NULL || face_number==0 )
		return 0;

	//	save old data
	unsigned int old_v_num = vertex_number;
	unsigned int old_f_num = face_number;

	//	create topology
	int* edge	= new int[old_f_num*9];
	int* fe		= edge	+ old_f_num*6;
	unsigned int e_num = createEdgeEFFEFromFace(edge, NULL, fe, face, old_f_num);

	//	add new vertices
	for( unsigned int i=0; i<e_num; i++ ){
		int ex3[2] = {edge[i*2]*3, edge[i*2+1]*3};
		vertex[vertex_number*3  ] = (vertex[ex3[0]  ] + vertex[ex3[1]  ]) * T(0.5);
		vertex[vertex_number*3+1] = (vertex[ex3[0]+1] + vertex[ex3[1]+1]) * T(0.5);
		vertex[vertex_number*3+2] = (vertex[ex3[0]+2] + vertex[ex3[1]+2]) * T(0.5);
		vertex_number ++;
	}

	//	add new faces, old face is not correct any more and need to be changed
	for( unsigned int i=0; i<old_f_num; i++ ){
		int v[3] = {face[i*3], face[i*3+1], face[i*3+2]};
		int e[3] = {fe[i*3], fe[i*3+1], fe[i*3+2]};

		//	sort edge sequence, so that it is on the opposite side of corresbonding v
		if(edge[e[0]*2]==v[0] || edge[e[0]*2+1]==v[0])	{int tmp=e[0]; e[0]=e[1]; e[1]=tmp;}
		if(edge[e[0]*2]==v[0] || edge[e[0]*2+1]==v[0])	{int tmp=e[0]; e[0]=e[2]; e[2]=tmp;}
		if(edge[e[1]*2]==v[1] || edge[e[1]*2+1]==v[1])	{int tmp=e[1]; e[1]=e[2]; e[2]=tmp;}

		//	insert new face
		for( int j=0; j<3; j++ ){
			int v0 = v[j];
			int v1 = old_v_num + e[(j+2)%3];
			int v2 = old_v_num + e[(j+1)%3];
			face[face_number*3  ] = v0;
			face[face_number*3+1] = v1;
			face[face_number*3+2] = v2;
			face_number ++;
		}

		//	change old face
		face[i*3  ] = old_v_num + e[0];
		face[i*3+1] = old_v_num + e[1];
		face[i*3+2] = old_v_num + e[2];
	}

	delete[] edge;

	return vertex_number;
}
///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_SUBD_H__