/***********************************************************/
/**	\file
	\brief		Tri Mesh with Texture
	\details	Tri mesh added texture and material
	\author		Yizhong Zhang
	\date		10/29/2014
*/
/***********************************************************/
#ifndef __YZ_TEXTURE_TRI_MESH_H__
#define __YZ_TEXTURE_TRI_MESH_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_tri_mesh.h"
#include "yzLib/yz_geometry/yz_mesh_texture.h"
#include "yzLib/yz_opengl/yz_texture.h"

namespace yz{	namespace geometry{

/**
	Triangle mesh with one texture and material

	Each vertex has one corresponding texture coordinate
*/
template<class T>
class SingleTextureTriMesh : 
	public TriMesh<T>, 
	public MeshTextureCoordinate<T>,
	public MeshMaterial,
	public opengl::Texture
{
public:
	using TriMesh<T>::vertex;
	using TriMesh<T>::face;
	using TriMesh<T>::vertex_normal;
	using TriMesh<T>::face_normal;
	using MeshTextureCoordinate<T>::tex_coord;

public:
	/**
		read mesh and texture from file
	*/
	virtual int ReadMeshFromFile(const char* file_name){
		std::vector<Com2<T>>	tmp_vt;
		std::vector<int3>		tmp_vt_idx;

		//	read data from the file
		int succ = yz::geometry::readTriMeshFromObjExt(
			file_name,
			&vertex,
			&face,
			(std::vector<Vec3<T>>*)NULL,
			(std::vector<int3>*)NULL,
			&tmp_vt,
			&tmp_vt_idx,
			this,
			&mtl_file_name
			);
		if( ! succ )
			return 0;

		//	create vertex coordinate
		tex_coord.resize( vertex.size() );
		if( tmp_vt_idx.size() == face.size() && !tmp_vt.empty() ){
			for(int i=0; i<tmp_vt_idx.size(); i++){
				for(int j=0; j<3; j++)
					tex_coord[face[i][j]] = tmp_vt[tmp_vt_idx[i][j]];
			}
		}

		//	create normals
		this->CalculateNormals();

		return 1;
	}

	/**
		create texture, must be called after read mesh and create window
	*/
	int CreateTexture(){
		if( map_Kd.tex_image.empty() ){
			return 0;
		}

		opengl::Texture::CreateTexture();

		SetupTexturePtr( 
			map_Kd.tex_w,
			map_Kd.tex_h,
			&map_Kd.tex_image[0],
			GL_RGB,
			GL_RGB,
			GL_UNSIGNED_BYTE);

		LoadPtrToTexture();

		return 1;
	}

	/**
		draw the mesh
	*/
	void Display(){
		GLenum	type;
		if( yz::Is_float<T>::check_type )
			type = GL_FLOAT;
		else if( yz::Is_double<T>::check_type )
			type = GL_DOUBLE;
		else{
			std::cout << "error: SingleTextureTriMesh::Display(), unsupported type" << std::endl;
			return;
		}

		if( tex_id ){
			Texture::Bind();
		}

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, type, 0, &vertex[0]);
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(type, 0, &vertex_normal[0]);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glTexCoordPointer(2, type, 0, &tex_coord[0]);

		glDrawElements(GL_TRIANGLES, face.size()*3, GL_UNSIGNED_INT, &face[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);

		if( tex_id ){
			Texture::UnBind();
		}
	}


public:
	std::string		mtl_file_name;

};

/**
	Triangle mesh with one texture and material

	Vertex on each face has a texture coordinate
*/
template<class T>
class SingleDisplayTextureTriMesh : public SingleTextureTriMesh<T>
{
public:
	using SingleTextureTriMesh<T>::vertex;
	using SingleTextureTriMesh<T>::face;
	using SingleTextureTriMesh<T>::vertex_normal;
	using SingleTextureTriMesh<T>::face_normal;
	using SingleTextureTriMesh<T>::tex_coord;
	using SingleTextureTriMesh<T>::tex_id;
	using SingleTextureTriMesh<T>::map_Kd;
	using SingleTextureTriMesh<T>::mtl_file_name;
	using SingleTextureTriMesh<T>::name;

	std::vector<yz::Vec3<T>>	display_vertex;
	std::vector<yz::Vec3<T>>	display_vertex_normal;
	std::vector<yz::int3>		display_face;
	std::vector<yz::Vec2<T>>	display_tex_coord;

public:
	/**
		read mesh and texture from file
	*/
	virtual int ReadMeshFromFile(const char* file_name){
		std::vector<Vec3<T>>	tmp_vn;
		std::vector<int3>		tmp_vn_idx;
		std::vector<Com2<T>>	tmp_vt;
		std::vector<int3>		tmp_vt_idx;

		//	read data from the file
		int succ = yz::geometry::readTriMeshFromObjExt(
			file_name,
			&vertex,
			&face,
			&tmp_vn,
			&tmp_vn_idx,
			&tmp_vt,
			&tmp_vt_idx,
			this,
			&mtl_file_name
			);
		if (!succ)
			return 0;

		//	create vertex coordinate
		if (tmp_vt_idx.size() == face.size() && !tmp_vt.empty()){
			tex_coord.resize(vertex.size());
			for (int i = 0; i < tmp_vt_idx.size(); i++){
				for (int j = 0; j < 3; j++)
					tex_coord[face[i][j]] = tmp_vt[tmp_vt_idx[i][j]];
			}
		}

		//	create normals
		this->CalculateNormals();

		//	create display vertex
		display_vertex.resize(face.size() * 3);
		for (int i = 0; i < face.size(); i++){
			for (int j = 0; j < 3; j++)
				display_vertex[i * 3 + j] = vertex[face[i][j]];
		}

		//	create display vertex normal
		display_vertex_normal.resize(display_vertex.size());
		if (tmp_vn_idx.size() == face.size() && !tmp_vn.empty()){
			for (int i = 0; i < tmp_vn_idx.size(); i++){
				for (int j = 0; j < 3; j++)
					display_vertex_normal[i * 3 + j] = tmp_vn[tmp_vn_idx[i][j]];
			}
		}
		else{
			for (int i = 0; i < face.size(); i++){
				for (int j = 0; j < 3; j++)
					display_vertex_normal[i * 3 + j] = vertex_normal[face[i][j]];
			}
		}

		//	create display face
		display_face.resize(face.size());
		for (int i = 0; i < display_face.size(); i++){
			display_face[i].x = i * 3;
			display_face[i].y = i * 3 + 1;
			display_face[i].z = i * 3 + 2;
		}

		//	create texture coordinate
		if (tmp_vt_idx.size() == face.size() && !tmp_vt.empty()){
			display_tex_coord.resize(face.size() * 3);
			for (int i = 0; i < tmp_vt_idx.size(); i++){
				for (int j = 0; j < 3; j++)
					display_tex_coord[i * 3 + j] = tmp_vt[tmp_vt_idx[i][j]];
			}
		}
		else{
			if (!tmp_vt.empty()){
				display_tex_coord.resize(face.size() * 3);
				for (int i = 0; i < tmp_vt_idx.size(); i++){
					for (int j = 0; j < 3; j++)
						display_tex_coord[face[i][j]] = tmp_vt[tmp_vt_idx[i][j]];
				}
			}
		}

		return 1;
	}

	/**
		write mesh to file

		\todo this function is not complete
	*/
	virtual int WriteMeshToFile(const char* file_name){
		std::ofstream obj(file_name);
		if (!obj.is_open()){
			return 0;
		}

		obj << "# mesh generated with SingleDisplayTextureTriMesh" << std::endl;
		obj << "mtllib " << mtl_file_name.c_str() << std::endl;

		for (int i = 0; i < vertex.size(); i++){
			obj << "v " << vertex[i].x << " " << vertex[i].y << " " << vertex[i].z << std::endl;
		}

		for (int i = 0; i < display_tex_coord.size(); i++){
			obj << "vt " << display_tex_coord[i].x << " " << display_tex_coord[i].y << std::endl;
		}

		obj << "usemtl " << name.c_str() << std::endl;
		for (int i = 0; i < face.size(); i++){
			obj << "f "
				<< face[i].x + 1 << "/" << i * 3 + 1 << " "
				<< face[i].y + 1 << "/" << i * 3 + 2 << " "
				<< face[i].z + 1 << "/" << i * 3 + 3 << std::endl;
		}

		obj.close();

		return 1;
	}

	/**
		create texture, must be called after read mesh and create window
	*/
	int CreateTexture(){
		if (map_Kd.tex_image.empty()){
			return 0;
		}

		opengl::Texture::CreateTexture();

		SetupTexturePtr(
			map_Kd.tex_w,
			map_Kd.tex_h,
			&map_Kd.tex_image[0],
			GL_RGB,
			GL_RGB,
			GL_UNSIGNED_BYTE);

		this->LoadPtrToTexture();

		return 1;
	}

	/**
		draw the mesh
	*/
	void Display(){
		GLenum	type;
		if (yz::Is_float<T>::check_type)
			type = GL_FLOAT;
		else if (yz::Is_double<T>::check_type)
			type = GL_DOUBLE;
		else{
			std::cout << "error: SingleTextureTriMesh::Display(), unsupported type" << std::endl;
			return;
		}

		if (tex_id)
			opengl::Texture::Bind();

		if (!display_vertex.empty() && !display_vertex_normal.empty()){
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, type, 0, &display_vertex[0]);
			glEnableClientState(GL_NORMAL_ARRAY);
			glNormalPointer(type, 0, &display_vertex_normal[0]);
		}

		if (!display_tex_coord.empty()){
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
			glTexCoordPointer(2, type, 0, &display_tex_coord[0]);
		}

		if (!display_face.empty()){
			glColor3f(1, 1, 1);
			glDrawElements(GL_TRIANGLES, display_face.size() * 3, GL_UNSIGNED_INT, &display_face[0]);
		}

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);

		if (tex_id)
			opengl::Texture::UnBind();
	}

};

/**
	Mesh with texture coordinate
*/
template <class T>
class TextureTriMesh :
	public TriMesh<T>,
	public MeshTextureCoordinate<T>,
	public MeshTextureFace
{
public:
	using TriMesh<T>::vertex;
	using TriMesh<T>::face;
	using MeshTextureCoordinate<T>::tex_coord;

public:
	virtual int ReadMeshFromFile(const char* file_name) {
		int succ = readTriMeshFromFile(
			file_name,
			&vertex,
			&face,
			(std::vector<Vec3<T>>*)NULL,
			(std::vector<int3>*)NULL,
			&tex_coord,
			&tex_face
		);
		if (!succ)
			return 0;

		this->CalculateNormals();
		return 1;
	}

};


}}	//	namespace yz::geometry

#endif	//	__YZ_TEXTURE_TRI_MESH_H__