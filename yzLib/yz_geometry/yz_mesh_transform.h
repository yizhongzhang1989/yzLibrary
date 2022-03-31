/***********************************************************/
/**	\file
	\brief		Transform the Mesh
	\author		Yizhong Zhang
	\date		8/30/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_TRANSFORM_H__
#define __YZ_MESH_TRANSFORM_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_matrix.h"
#include "yzLib/yz_math/yz_numerical_utils.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Transform Vertices by Matrix
*/
//	========================================
/**
	Rigid transform 3D vertices

	This function is identical to combine the rotation matrix and 
	translation vector to a 4x4 matrix, then transform

	\param	vertex_ptr		pointer to vertex list, type must be basic data type float / double
	\param	vertex_number	number of vertex
	\param	rotation		rotation matrix, 3x3
	\param	translation		translation vector
*/
template<typename T1, typename T2, typename T3>
inline void transformVertices3D(T1*				vertex_ptr,
								int				vertex_number,
								Matrix3x3<T2>	rotation,
								Vec3<T3>		translation){
	for( int i=0; i<vertex_number; i++ ){
		(*(Vec3<T1>*)&vertex_ptr[i*3]) = rotation * (*(Vec3<T1>*)&vertex_ptr[i*3]) + translation;
	}
}
	
/**
	Rigid transform 3D vertices

	This function is identical to combine the rotation matrix and 
	translation vector to a 4x4 matrix, then transform

	\param	vertex		the vertex list
	\param	rotation	rotation matrix, 3x3
	\param	translation	translation vector
*/
template<typename T1, typename T2, typename T3>
inline void transformVertices(std::vector<Vec3<T1>>&	vertex,
							  Matrix3x3<T2>				rotation,
							  Vec3<T3>					translation){
	for( int i=0; i<vertex.size(); i++ )
		vertex[i] = rotation * vertex[i] + translation;
}

/**
	Rigid transform 3D vertices

	\param	vertex_ptr		the vertex list
	\param	vertex_number	number of vertices
	\param	trans_mat		transform matrix, 4x4
*/
template<typename T1, typename T2>
inline void transformVertices3D(T1*				vertex_ptr,
								int				vertex_number,
								Matrix4x4<T2>	trans_mat){
	for( int i=0; i<vertex_number; i++ ){
		(*(Vec3<T1>*)&vertex_ptr[i*3]) = trans_mat * (*(Vec3<T1>*)&vertex_ptr[i*3]);
	}
}

/**
	Rigid transform 3D vertices

	\param	vertex		the vertex list
	\param	trans_mat	transform matrix, 4x4
*/
template<typename T1, typename T2>
inline void transformVertices(std::vector<Vec3<T1>>&	vertex,
							  Matrix4x4<T2>				trans_mat){
	for( int i=0; i<vertex.size(); i++ )
		vertex[i] = trans_mat * vertex[i];
}
///@}

//	========================================
///@{
/**	@name Transform Mesh Vertices Directly
*/
//	========================================

/**
	translate each vertex by given vector

	\param	vertex_ptr		pointer to 3D vertex array
	\param	vertex_number	number of vertices
	\param	trans_x			translate vector x component
	\param	trans_y			translate vector y component
	\param	trans_z			translate vector z component
*/
template<typename T, typename E>
inline void translateVertices3D(T*		vertex_ptr,
								int		vertex_number,
								E		trans_x,
								E		trans_y,
								E		trans_z ){
	for( int i=0; i<vertex_number; i++ ){
		vertex_ptr[i*3  ] += trans_x;
		vertex_ptr[i*3+1] += trans_y;
		vertex_ptr[i*3+2] += trans_z;
	}
}

/**
	translate each vertex by given vector

	\param	vertex_ptr		pointer to 2D vertex array
	\param	vertex_number	number of vertices
	\param	trans_x			translate vector x component
	\param	trans_y			translate vector y component
*/
template<typename T, typename E>
inline void translateVertices2D(T*		vertex_ptr,
								int		vertex_number,
								E		trans_x,
								E		trans_y ){
	for( int i=0; i<vertex_number; i++ ){
		vertex_ptr[i*2  ] += trans_x;
		vertex_ptr[i*2+1] += trans_y;
	}
}

/**
	translate each vertex by given vector

	\param	vertex			vertex array
	\param	trans_vec		translation of each vertex
*/
template<typename T, typename E>
inline void translateVertices(std::vector<Vec3<T>>&	vertex,
							  Vec3<E>				trans_vec){
	if( vertex.empty() )
		return;
	translateVertices3D((T*)&vertex[0], vertex.size(), 
		trans_vec.x, trans_vec.y, trans_vec.z);
}

/**
	translate each vertex by given vector

	\param	vertex			vertex array
	\param	trans_vec		translation of each vertex
*/
template<typename T, typename E>
inline void translateVertices(std::vector<Vec2<T>>&	vertex,
							  Vec2<E>				trans_vec){
	if( vertex.empty() )
		return;
	translateVertices2D((T*)&vertex[0], vertex.size(), 
		trans_vec.x, trans_vec.y, trans_vec.z);
}

/**
	translate the mesh so the center of the mesh is original point

	\param	vertex_ptr		pointer to 3D vertex array
	\param	vertex_number	number of vertices
*/
template<typename T>
inline void translateVerticesToCenter3D(T*		vertex_ptr,
										int		vertex_number){
	Vec3<T> bb_min, bb_max;
	getAABBCoef(bb_min, bb_max, vertex_ptr, vertex_number);
	Vec3<T> center = (bb_min + bb_max) * 0.5;
	translateVertices3D(vertex_ptr, vertex_number, -center[0], -center[1], -center[2]);
}

/**
	translate the mesh so the center of the mesh is original point

	\param	vertex_ptr		pointer to 2D vertex array
	\param	vertex_number	number of vertices
*/
template<typename T>
inline void translateVerticesToCenter2D(T*		vertex_ptr,
										int		vertex_number){
	Vec2<T> bb_min, bb_max;
	getAABBCoef(bb_min, bb_max, vertex_ptr, vertex_number);
	Vec2<T> center = (bb_min + bb_max) * 0.5;
	translateVertices2D(vertex_ptr, vertex_number, -center[0], -center[1]);
}

/**
	translate the mesh so the center of the mesh is original point

	\param	vertex		vertices
*/
template<typename T>
inline void translateVerticesToCenter(std::vector<Vec3<T>>&	vertex){
	if( vertex.empty() )
		return;
	translateVerticesToCenter3D((T*)&vertex[0], vertex.size());
}

/**
	translate the mesh so the center of the mesh is original point

	\param	vertex		vertices
*/
template<typename T>
inline void translateVerticesToCenter(std::vector<Vec2<T>>&	vertex){
	if( vertex.empty() )
		return;
	translateVerticesToCenter2D((T*)&vertex[0], vertex.size());
}

/**
	Move each vertex of a mesh by given vector and scale

	\param	vertex_ptr		pointer to 3D vertex array
	\param	vertex_number	number of vertices
	\param	move_vec_ptr	moving vector of each vertex
	\param	scale			the scale of moving applied to vector
*/
template<typename T>
inline void moveVerticesByVector3D(T*		vertex_ptr,
								   int		vertex_number,
								   const T*	move_vec_ptr,
								   double	scale){
	for( int i=0; i<vertex_number*3; i++ ){
		vertex_ptr[i] += move_vec_ptr[i] * scale;
	}
}

/**
	Move each vertex of a mesh by given vector and scale

	\param	vertex_ptr		pointer to 2D vertex array
	\param	vertex_number	number of vertices
	\param	move_vec_ptr	moving vector of each vertex
	\param	scale			the scale of moving applied to vector
*/
template<typename T>
inline void moveVerticesByVector2D(T*		vertex_ptr,
								   int		vertex_number,
								   const T*	move_vec_ptr,
								   double	scale){
	for( int i=0; i<vertex_number*2; i++ ){
		vertex_ptr[i] += move_vec_ptr[i] * scale;
	}
}

/**
	Move each vertex of a mesh by given vector and scale

	\param	vertex			vertex array
	\param	move_vec		moving vector of each vertex
	\param	scale			the scale of moving applied to vector
*/
template<typename T>
inline void moveVerticesByVector(std::vector<Vec3<T>>&	vertex,
								 std::vector<Vec3<T>>&	move_vec,
								 double					scale){
	if(vertex.empty() || move_vec.size()<vertex.size())
		return;

	moveVerticesByVector3D((T*)&vertex[0], vertex.size(), (T*)&move_vec[0], scale);
}
/**
	Move each vertex of a mesh by given vector and scale

	\param	vertex			vertex array
	\param	move_vec		moving vector of each vertex
	\param	scale			the scale of moving applied to vector
*/
template<typename T>
inline void moveVerticesByVector(std::vector<Vec2<T>>&	vertex,
								 std::vector<Vec2<T>>&	move_vec,
								 double					scale){
	if(vertex.empty() || move_vec.size()<vertex.size())
		return;

	moveVerticesByVector2D((T*)&vertex[0], vertex.size(), (T*)&move_vec[0], scale);
}
/**
	Add Noise to 3D Vertices

	\param	vertex_ptr		pointer to 3D vertex array 
	\param	vertex_number	number of vertices
	\param	noise_scale		how big the noise to be
*/
template<typename T>
inline void applyNoiseToVertices3D(T*		vertex_ptr,
								   int		vertex_number,
								   double	noise_scale){
	for( int i=0; i<vertex_number*3; i++ ){
		vertex_ptr[i] += (rand0to1d() - 0.5) * noise_scale;
	}
}

/**
	Add Noise to 2D Vertices

	\param	vertex_ptr		pointer to 2D vertex array 
	\param	vertex_number	number of vertices
	\param	noise_scale		how big the noise to be
*/
template<typename T>
inline void applyNoiseToVertices2D(T*		vertex_ptr,
								   int		vertex_number,
								   double	noise_scale){
	for( int i=0; i<vertex_number*2; i++ ){
		vertex_ptr[i] += (rand0to1d() - 0.5) * noise_scale;
	}
}

/**
	Add Noise to 3D Vertices

	\param	vertex			vertex array
	\param	noise_scale		how big the noise to be
*/
template<typename T>
inline void applyNoiseToVertices(std::vector<Vec3<T>>&	vertex,
								 double					noise_scale){
	if( vertex.empty() )
		return;
	applyNoiseToVertices3D((T*)&vertex[0], vertex.size(), noise_scale);
}



/**
	Add Noise to 2D Vertices

	\param	vertex			vertex array
	\param	noise_scale		how big the noise to be
*/
template<typename T>
inline void applyNoiseToVertices(std::vector<Vec2<T>>&	vertex,
								 double					noise_scale){
	if( vertex.empty() )
		return;
	applyNoiseToVertices2D((T*)&vertex[0], vertex.size(), noise_scale);
}

///@}

//	========================================
///@{
/**	@name Mesh Deformation
*/
//	========================================

/**
	Deform the mesh by Enright Test

	\param	vertex		mesh vertex
	\param	step		step size of Enright Test
*/
template <typename T>
inline void EnrightTestStep(
	std::vector<Vec3<T>>&	vertex, 
	T						step
) {
	for (unsigned int i = 0; i != vertex.size(); i++) {
		T x = vertex[i].x;
		T y = vertex[i].y;
		T z = vertex[i].z;
		T vx = 2 * sin(YZ_PI*x)*sin(YZ_PI*x) * sin(2 * YZ_PI*y) * sin(2 * YZ_PI*z);
		T vy = -sin(2 * YZ_PI*x) * sin(YZ_PI*y)*sin(YZ_PI*y) * sin(2 * YZ_PI*z);
		T vz = -sin(2 * YZ_PI*x) * sin(2 * YZ_PI*y) * sin(YZ_PI*z)*sin(YZ_PI*z);
		vertex[i].x += vx * step;
		vertex[i].y += vy * step;
		vertex[i].z += vz * step;
	}
}

///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_TRANSFORM_H__