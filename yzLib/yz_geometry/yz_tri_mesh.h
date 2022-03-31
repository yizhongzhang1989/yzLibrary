/***********************************************************/
/**	\file
	\brief		Triangle Mesh
	\details	Data structure needed to represent a mesh for display.
				data include: vertex position, face, vertex normal, 
				face normal. Choose the right class according to your need.
	\author		Yizhong Zhang
	\date		6/1/2012
*/
/***********************************************************/
#ifndef __YZ_TRI_MESH_H__
#define __YZ_TRI_MESH_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_rw.h"
#include "yzLib/yz_geometry/yz_mesh_normal.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{	namespace geometry{

/**
	base triangle mesh ptr, just vertex and face pointers
*/
template<class T>
class PtrBaseTriMesh{
public:
	T*		vertex_ptr;
	int*	face_ptr;

public:
	//	constructor
	PtrBaseTriMesh() : vertex_ptr(NULL), face_ptr(NULL){}
	/**
		Reset all members. clear pointers, but don't touch memory
	*/
	inline void Reset(){
		vertex_ptr	= NULL;
		face_ptr	= NULL;
	}
};

/**
	flat shading triangle mesh ptr, with face normal pointer
*/
template<class T>
class PtrFlatShadingTriMesh : public PtrBaseTriMesh<T>{
public:
	T*	face_normal_ptr;

public:
	//	constructor
	PtrFlatShadingTriMesh() : PtrBaseTriMesh<T>(), face_normal_ptr(NULL){}
	/**
		Reset all members. clear pointers, but don't touch memory
	*/
	inline void Reset(){
		PtrBaseTriMesh<T>::Reset();
		face_normal_ptr	= NULL;
	}
};

/**
	smooth shading triangle mesh ptr, with vertex normal pointer
*/
template<class T>
class PtrSmoothShadingTriMesh : public PtrBaseTriMesh<T>{
public:
	T*	vertex_normal_ptr;

public:
	//	constructor
	PtrSmoothShadingTriMesh() : PtrBaseTriMesh<T>(), vertex_normal_ptr(NULL){}
	/**
		Reset all members. clear pointers, but don't touch memory
	*/
	inline void Reset(){
		PtrBaseTriMesh<T>::Reset();
		vertex_normal_ptr	= NULL;
	}
};



/**
	triangle mesh with vertex normal and face normal
*/
template<class T>
class PtrTriMesh : public PtrSmoothShadingTriMesh<T>{
public:
	T*	face_normal_ptr;

public:
	//	constructor
	PtrTriMesh() : PtrSmoothShadingTriMesh<T>(), face_normal_ptr(NULL) {}

	/**
		Reset all members. clear pointers, but don't touch memory
	*/
	inline void Reset(){
		PtrSmoothShadingTriMesh<T>::Reset();
		face_normal_ptr	= NULL;
	}
};

/**
	base triangle mesh using vector, just vertex and face
*/
template<class T>
class BaseTriMesh{
public:
	std::vector<Vec3<T>>	vertex;
	std::vector<int3>		face;

public:
	//	constructor
	BaseTriMesh(){};

	//	virtual base functions
	/**
		Reset all members, clear all memory
	*/
	virtual inline void Reset(){
		vertex.clear();
		face.clear();
	}

	/**
	*/
	virtual inline int ReadMeshFromFile(const char* file_name){	//	only vertex and triangle are loaded in BaseTriMesh
		return	readTriMeshFromFile(file_name, &vertex, &face );
	}

	virtual inline int Display(int display_mode=0) const{
		std::cout << "BaseTriMesh con't have Draw() function implemented" << std::endl;
		return 0;
	}

	virtual inline void AppendTriMesh(const BaseTriMesh<T>& mesh){
		int old_vertex_number = vertex.size();
		int old_face_number = face.size();
		vertex.insert(vertex.end(), mesh.vertex.begin(), mesh.vertex.end());
		face.insert(face.end(), mesh.face.begin(), mesh.face.end());
		for( int f=old_face_number; f<face.size(); f++ ){
			face[f].x += old_vertex_number;
			face[f].y += old_vertex_number;
			face[f].z += old_vertex_number;
		}
	}

};

/**
	flat shading triangle mesh using vector, with face normal
*/
template<class T>
class FlatShadingTriMesh : public BaseTriMesh<T>{
public:
	using BaseTriMesh<T>::vertex;
	using BaseTriMesh<T>::face;

	std::vector<Vec3<T>>	face_normal;

public:
	//	constructor
	FlatShadingTriMesh():BaseTriMesh<T>(){};

	//	overload virtual function
	virtual inline void Reset(){
		BaseTriMesh<T>::Reset();
		face_normal.clear();
	}

	virtual inline int ReadMeshFromFile(const char* file_name){
		if( BaseTriMesh<T>::ReadMeshFromFile(file_name) ){
			CalculateFaceNormal();
			return 1;
		}
		
		return 0;
	}

	virtual inline int Display(int display_mode=0) const{
		#ifdef YZ_gl_h
			opengl::drawFlatShadingTriMesh(vertex, face, face_normal);
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in FlatShadingTriMesh" << std::endl;
			return 0;
		#endif
	}

	virtual inline void AppendTriMesh(const BaseTriMesh<T>& mesh){
		BaseTriMesh<T>::AppendTriMesh( mesh );
		CalculateFaceNormal();
	}
	//	own function
	inline void CalculateFaceNormal(){
		calculateFaceNormal(face_normal, vertex, face);
	}


};

/**
	smooth shading triangle mesh using vector, with vertex normal
*/
template<class T>
class SmoothShadingTriMesh : public BaseTriMesh<T>{
public:
	using BaseTriMesh<T>::vertex;
	using BaseTriMesh<T>::face;

	std::vector<Vec3<T>>	vertex_normal;

public:
	//	constructor
	SmoothShadingTriMesh():BaseTriMesh<T>(){};

	//	overload virtual function
	virtual inline void Reset(){
		BaseTriMesh<T>::Reset();
		vertex_normal.clear();
	}

	virtual inline int ReadMeshFromFile(const char* file_name){
		if( BaseTriMesh<T>::ReadMeshFromFile(file_name) ){
			CalculateVertexNormal();
			return 1;
		}
		
		return 0;
	}

	virtual inline int Display(int display_mode=0){
		#ifdef YZ_gl_h
			opengl::drawSmoothShadingTriMesh(vertex, face, vertex_normal);
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in SmoothShadingTriMesh" << std::endl;
			return 0;
		#endif
	}

	virtual inline void AppendTriMesh(const BaseTriMesh<T>& mesh){
		BaseTriMesh<T>::AppendTriMesh( mesh );
		CalculateVertexNormal();
	}
	//	own function
	inline void CalculateVertexNormal(){
		calculateVertexNormal(vertex_normal, vertex, face);
	}
};

/**
	Triangle Mesh, contain both vertex normal and face normal

	This class inherent from SmoothShadingTriMesh, so face normal
	is added in this class.
*/
template<class T>
class TriMesh : public SmoothShadingTriMesh<T>{
public:
	using BaseTriMesh<T>::vertex;
	using BaseTriMesh<T>::face;
	using SmoothShadingTriMesh<T>::vertex_normal;

	std::vector<Vec3<T>>	face_normal;

public:
	//	constructor
	TriMesh():SmoothShadingTriMesh<T>(){}
	//	overload virtual function
	virtual inline void Reset(){
		SmoothShadingTriMesh<T>::Reset();
		face_normal.clear();
	}

	virtual inline int ReadMeshFromFile(const char* file_name){
		if( BaseTriMesh<T>::ReadMeshFromFile(file_name) ){
			CalculateNormals();
			return 1;
		}
		
		return 0;
	}

	virtual inline int Display(int display_mode= 0x1D01){	//	default value: GL_SMOOTH
		#ifdef YZ_gl_h
			if( display_mode == GL_SMOOTH )
				SmoothShadingTriMesh<T>::Display(display_mode);
			else if( display_mode == GL_FLAT )
				DisplayFlat();
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in TriMesh" << std::endl;
			return 0;
		#endif
	}

	virtual inline void AppendTriMesh(const BaseTriMesh<T>& mesh){
		BaseTriMesh<T>::AppendTriMesh( mesh );
		CalculateNormals();
	}


	//	own function
	inline void CalculateNormals(){
		CalculateFaceNormal();
		this->CalculateVertexNormal();
	}

	inline void CalculateFaceNormal(){
		calculateFaceNormal(face_normal, vertex, face);
	}

	inline int DisplayFlat() const{
		#ifdef YZ_gl_h
			opengl::drawFlatShadingTriMesh(vertex, face, face_normal);
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use DisplayFlat() in TriMesh" << std::endl;
			return 0;
		#endif
	}

};


/**
	shading triangle mesh, vertex on each face has it's own normal

	Normal arrangement is different from TriMesh. This class is used
	when the mesh contain smooth and flat part at the same time.
*/
template<class T>
class ShadingTriMesh : public BaseTriMesh<T>{
public:
	using BaseTriMesh<T>::vertex;
	using BaseTriMesh<T>::face;

	std::vector<Vec3<T>>	face_vertex_normal;		//	vertex on each face has it's own normal

public:
	//	constructor
	ShadingTriMesh():BaseTriMesh<T>(), normal_status(0){}

	//	overload virtual function
	virtual inline int ReadMeshFromFile(const char* file_name){
		std::vector<Vec3<T>> vn;
		std::vector<int3>	 vn_idx;
		if( readTriMeshFromFile(file_name, &vertex, &face, &vn, &vn_idx) ){
			if( vn.empty() )	//	the obj file don't have normal information
				CalculateFaceNormal();
			else{				//	copy read vn to face_vertex_normal
				face_vertex_normal.resize( face.size()*3 );	// each face has 3 vertices
				for( int i=0; i<face.size(); i++ ){
					face_vertex_normal[i*3  ] = vn[vn_idx[i].x];
					face_vertex_normal[i*3+1] = vn[vn_idx[i].y];
					face_vertex_normal[i*3+2] = vn[vn_idx[i].z];
				}
				normal_status = 0;
			}
			return 1;
		}

		return 0;
	}

	virtual inline void Reset(){
		BaseTriMesh<T>::Reset();
		face_vertex_normal.clear();
		normal_status = 0;
	}

	virtual inline int Display(int display_mode=0) const{
		#ifdef YZ_gl_h
			if( display_mode == GL_SMOOTH )
				CalculateVertexNormal();
			else if( display_mode == GL_FLAT )
				CalculateFaceNormal();

			glBegin( GL_TRIANGLES );
				for( unsigned int i=0; i<face.size(); i++ ){
					glNormal3f( face_vertex_normal[i*3  ].x, face_vertex_normal[i*3  ].y, face_vertex_normal[i*3  ].z );
					glVertex3f( vertex[face[i].x].x, vertex[face[i].x].y, vertex[face[i].x].z );
					glNormal3f( face_vertex_normal[i*3+1].x, face_vertex_normal[i*3+1].y, face_vertex_normal[i*3+1].z );
					glVertex3f( vertex[face[i].y].x, vertex[face[i].y].y, vertex[face[i].y].z );
					glNormal3f( face_vertex_normal[i*3+2].x, face_vertex_normal[i*3+2].y, face_vertex_normal[i*3+2].z );
					glVertex3f( vertex[face[i].z].x, vertex[face[i].z].y, vertex[face[i].z].z );
				}
			glEnd();

			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in SmoothShadingTriMesh" << std::endl;
			return 0;
		#endif
	}

	virtual inline void AppendTriMesh(const BaseTriMesh<T>& mesh){
		BaseTriMesh<T>::AppendTriMesh( mesh );
		CalculateVertexNormal();
	}
	//	own function
	inline void CalculateFaceNormal(){
		if( normal_status == 2 )
			return;

		face_vertex_normal.resize( face.size()*3 );	// each face has 3 vertices
		for( int i=0; i<face.size(); i++ ){
			Vec3<T> r1	= vertex[face[i].y] - vertex[face[i].x];
			Vec3<T> r2	= vertex[face[i].z] - vertex[face[i].x];
			Vec3<T> nor	= cross(r1, r2).Normalize();
			face_vertex_normal[i*3  ] = nor;
			face_vertex_normal[i*3+1] = nor;
			face_vertex_normal[i*3+2] = nor;
		}

		normal_status = 2;
	}

	inline void CalculateVertexNormal(){
		if( normal_status == 1 )
			return;

		//	create temp vertex normal array and calculate normal
		std::vector<Vec3<T>> vertex_normal;
		vertex_normal.resize( vertex.size() );
		memset(&vertex_normal[0], 0, sizeof(Vec3<T>)*vertex_normal.size());
		for( int i=0; i<face.size(); i++ ){
			Vec3<T> r1	= vertex[face[i].y] - vertex[face[i].x];
			Vec3<T> r2	= vertex[face[i].z] - vertex[face[i].x];
			Vec3<T> nor	= cross(r1, r2).Normalize();

			vertex_normal[face[i].x] += nor;
			vertex_normal[face[i].y] += nor;
			vertex_normal[face[i].z] += nor;
		}
		for( int i=0; i<vertex_normal.size(); i++ )
			vertex_normal[i].SetNormalize();

		//	copy vertex normal to each face vertex normal
		face_vertex_normal.resize( face.size()*3 );	// each face has 3 vertices
		for( int i=0; i<face.size(); i++ ){
			face_vertex_normal[i*3  ] = vertex_normal[face[i].x];
			face_vertex_normal[i*3+1] = vertex_normal[face[i].y];
			face_vertex_normal[i*3+2] = vertex_normal[face[i].z];
		}

		normal_status = 1;
	}

private:
	int normal_status;		//	0: same as read file;	1: vertex normal;	2: face normal
};

//	------------------------------------
//	type define
//	------------------------------------

typedef BaseTriMesh<float>				BaseTriMeshf;
typedef BaseTriMesh<double>				BaseTriMeshd;

typedef FlatShadingTriMesh<float>		FlatShadingTriMeshf;
typedef FlatShadingTriMesh<double>		FlatShadingTriMeshd;

typedef SmoothShadingTriMesh<float>		SmoothShadingTriMeshf;
typedef SmoothShadingTriMesh<double>	SmoothShadingTriMeshd;

typedef TriMesh<float>					TriMeshf;
typedef TriMesh<double>					TriMeshd;

typedef ShadingTriMesh<float>			ShadingTriMeshf;
typedef ShadingTriMesh<double>			ShadingTriMeshd;


}}	//	namespace yz::geometry

#endif	//	__YZ_TRI_MESH_H__