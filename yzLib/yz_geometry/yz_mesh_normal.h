/***********************************************************/
/**	\file
	\brief		Calculate Normal
	\author		Yizhong Zhang
	\date		6/4/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_NORMAL_H__
#define __YZ_MESH_NORMAL_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Calculate Normal on Mesh
*/
//	========================================

/**
	calculate face normal 

	\param	face_id		the index of face to calculate normal
	\param	vertex		vertex of the mesh
	\param	face		face of the mesh
	\return				the normal of the given face
*/
template<typename T>
inline Vec3<T> calculateFaceNormal(int							face_id, 
								   const std::vector<Vec3<T>>&	vertex, 
								   const std::vector<int3>&		face){
	assert(face_id>=0 && face_id<face.size());
	return calculateFaceNormal(vertex[face[face_id].x], 
		vertex[face[face_id].y], vertex[face[face_id].z]);
}

/**
	calculate face normal given three vertices

	\param	v0		position of vertex 0
	\param	v1		position of vertex 1
	\param	v2		position of vertex 2
	\return			return the normal
*/
template<typename T>
inline Vec3<T> calculateFaceNormal(Vec3<T> v0, Vec3<T> v1, Vec3<T> v2){
	return cross(v2-v1, v0-v1).Normalize();
}


/**
	Calculate face normal, with assertion that type of T is correct. 

	This function only accept T as float / double

	\param	face_normal,	face normal list ptr, format xyz_xyz_, size = face_number*3
	\param	vertex,			vertex list ptr, format xyz_xyz_
	\param	vertex_number	vertex number
	\param	face,			face list ptr, format v0v1v2_v0v1v2_
	\param	face_number		face number
*/
template<typename T>
inline void calculateFaceNormal_TypeSafe(T*			face_normal, 
										 const T*	vertex, 
										 int		vertex_number, 
										 const int* face, 
										 int		face_number){
	if( face_normal==NULL )													//	no output
		return;
	if( vertex==NULL || vertex_number==0 || face==NULL || face_number==0 )	//	no input
		return;

	for( int i=0; i<face_number; i++ ){
		int fx3[3] = {face[i*3]*3, face[i*3+1]*3, face[i*3+2]*3};
		Vec3<T> r1(	vertex[fx3[1]  ] - vertex[fx3[0]  ],	//	vector: v0 -> v1
					vertex[fx3[1]+1] - vertex[fx3[0]+1],
					vertex[fx3[1]+2] - vertex[fx3[0]+2]	);
		Vec3<T> r2(	vertex[fx3[2]  ] - vertex[fx3[0]  ],	//	vector: v0 -> v2
					vertex[fx3[2]+1] - vertex[fx3[0]+1],
					vertex[fx3[2]+2] - vertex[fx3[0]+2]	);
		Vec3<T> nor	= cross(r1, r2).Normalize();
		(*(Vec3<T>*)&face_normal[i*3]) = nor;
	}
}

/**
	calculate face normal, if T is not float/double, print error information

	\param	face_normal,	face normal list ptr, format xyz_xyz_, size = face_number*3
	\param	vertex,			vertex list ptr, format xyz_xyz_
	\param	vertex_number	vertex number
	\param	face,			face list ptr, format v0v1v2_v0v1v2_
	\param	face_number		face number
*/
template<typename T>
inline void calculateFaceNormal(T*			face_normal, 
								const T*	vertex, 
								int			vertex_number, 
								const int*	face, 
								int			face_number){
	std::cout << "invalid type, calculateFaceNormal only accept float / double" << std::endl;
}

template<> inline void calculateFaceNormal(float*		face_normal, 
										   const float* vertex, 
										   int			vertex_number, 
										   const int*	face, 
										   int			face_number){
	calculateFaceNormal_TypeSafe(face_normal, vertex, vertex_number, face, face_number);
}

template<> inline void calculateFaceNormal(double*			face_normal, 
										   const double*	vertex, 
										   int				vertex_number, 
										   const int*		face, 
										   int				face_number){
	calculateFaceNormal_TypeSafe(face_normal, vertex, vertex_number, face, face_number);
}



/**
	calculate face normal using vector

	\param	face_normal,	face normal vector
	\param	vertex,			vertex vector
	\param	face,			face vector
*/
template<typename T>
inline void calculateFaceNormal(std::vector<Vec3<T>>&		face_normal, 
								const std::vector<Vec3<T>>&	vertex, 
								const std::vector<int3>&	face){
	if( face.size() <= 0 )
		return;
	face_normal.resize( face.size() );
	calculateFaceNormal((T*)&face_normal[0], (T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}

/**
	calculate vertex normal, with assertion that type of T is correct.

	this function only accept T as float / double

	\param	vertex_normal,	vertex normal list ptr, format xyz_xyz_, size = vertex_number*3
	\param	vertex,			vertex list ptr, format xyz_xyz_
	\param	vertex_number	vertex number
	\param	face,			face list ptr, format v0v1v2_v0v1v2_
	\param	face_number		face number
*/
template<typename T>
inline void calculateVertexNormal_TypeSafe(T*			vertex_normal, 
										   const T*		vertex, 
										   int			vertex_number, 
										   const int*	face, 
										   int			face_number){
	if( vertex_normal==NULL )												//	no output
		return;
	if( vertex==NULL || vertex_number==0 || face==NULL || face_number==0 )	//	no input
		return;

	std::fill(&vertex_normal[0], &vertex_normal[vertex_number*3], 0);//memset(vertex_normal, 0, sizeof(T)*vertex_number*3);
	for( int i=0; i<face_number; i++ ){
		int fx3[3] = {face[i*3]*3, face[i*3+1]*3, face[i*3+2]*3};
		Vec3<T> r1(	vertex[fx3[1]  ] - vertex[fx3[0]  ],	//	vector: v0 -> v1
					vertex[fx3[1]+1] - vertex[fx3[0]+1],
					vertex[fx3[1]+2] - vertex[fx3[0]+2]	);
		Vec3<T> r2(	vertex[fx3[2]  ] - vertex[fx3[0]  ],	//	vector: v0 -> v2
					vertex[fx3[2]+1] - vertex[fx3[0]+1],
					vertex[fx3[2]+2] - vertex[fx3[0]+2]	);
		Vec3<T> nor	= cross(r1, r2).Normalize();

		(*(Vec3<T>*)&vertex_normal[fx3[0]]) += nor;
		(*(Vec3<T>*)&vertex_normal[fx3[1]]) += nor;
		(*(Vec3<T>*)&vertex_normal[fx3[2]]) += nor;
	}
	for( int i=0; i<vertex_number; i++ )
		(*(Vec3<T>*)&vertex_normal[i*3]).SetNormalize();	//	treate that space as vec3 and set normalize
}

/**
	calculate vertex normal, if T is not float/double, print error information

	\param	vertex_normal,	vertex normal list ptr, format xyz_xyz_, size = vertex_number*3
	\param	vertex,			vertex list ptr, format xyz_xyz_
	\param	vertex_number	vertex number
	\param	face,			face list ptr, format v0v1v2_v0v1v2_
	\param	face_number		face number
*/
template<typename T>
inline void calculateVertexNormal(T*			vertex_normal, 
								  const T*		vertex, 
								  int			vertex_number, 
								  const int*	face, 
								  int			face_number){
	std::cout << "invalid type, calculateVertexNormal only accept float / double" << std::endl;
}

template<> inline void calculateVertexNormal(float*			vertex_normal, 
											 const float*	vertex, 
											 int			vertex_number, 
											 const int*		face, 
											 int			face_number){
	calculateVertexNormal_TypeSafe(vertex_normal, vertex, vertex_number, face, face_number);
}

template<> inline void calculateVertexNormal(double*		vertex_normal, 
											 const double*	vertex, 
											 int			vertex_number, 
											 const int*		face, 
											 int			face_number){
	calculateVertexNormal_TypeSafe(vertex_normal, vertex, vertex_number, face, face_number);
}

/**
	calculate vertex normal using vector

	\param	vertex_normal,	vertex normal vector
	\param	vertex,			vertex vector
	\param	face,			face vector
*/
template<class T>
inline void calculateVertexNormal(std::vector<Vec3<T>>&			vertex_normal, 
								  const std::vector<Vec3<T>>&	vertex, 
								  const std::vector<int3>&		face){
	if (vertex.empty() || face.empty())
		return;
	vertex_normal.resize( vertex.size() );
	calculateVertexNormal((T*)&vertex_normal[0], (T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}

///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_NORMAL_H__