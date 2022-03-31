/***********************************************************/
/**	\file
	\brief		Quad Mesh
	\details	quad only mesh
	\author		Yizhong Zhang
	\date		3/26/2019
*/
/***********************************************************/
#ifndef __YZ_QUAD_MESH__
#define __YZ_QUAD_MESH__

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
quad mesh
*/
template<class T>
class QuadMesh {
public:
	std::vector<Vec3<T>>	vertex;
	std::vector<int4>		face;
	std::vector<Vec3<T>>	vertex_normal;
	std::vector<Vec3<T>>	face_normal;

protected:
	//	variables used for display flat
	std::vector<Vec3<T>>		disp_split_vertex;
	std::vector<Vec3<T>>		disp_split_vertex_normal;
	std::vector<unsigned int>	disp_indices;

public:
	//	constructor
	QuadMesh() {};

	/**
	Reset all members, clear all memory
	*/
	virtual inline void Reset() {
		vertex.clear();
		face.clear();
		vertex_normal.clear();
		face_normal.clear();
	}

	/**
	read the mesh

	just keep quad faces
	*/
	virtual inline int ReadMeshFromFile(const char* file_name) {	//	only vertex and triangle are loaded in BaseTriMesh
		std::vector<int> f, f_start;
		readMeshFromObj(file_name, &vertex, &f, &f_start);

		face.clear();
		face.reserve(f_start.size() - 1);

		for (int i = 1; i < f_start.size(); i++) {
			int fv_count = f_start[i] - f_start[i - 1];
			if (fv_count != 4)
				continue;

			int4 curr_f(
				f[f_start[i - 1]],
				f[f_start[i - 1] + 1],
				f[f_start[i - 1] + 2],
				f[f_start[i - 1] + 3]
			);
			face.push_back(curr_f);
		}

		CalculateNormals();

		return 1;
	}

	inline void CalculateNormals() {
		CalculateFaceNormal();
		CalculateVertexNormal();
	}

	inline void CalculateFaceNormal() {
		face_normal.resize(face.size() * 3);	// each face has 3 vertices
		for (int i = 0; i<face.size(); i++) {
			Vec3<T> r1 = vertex[face[i].y] - vertex[face[i].x];
			Vec3<T> r2 = vertex[face[i].z] - vertex[face[i].x];
			Vec3<T> r3 = vertex[face[i].w] - vertex[face[i].x];
			Vec3<T> nor = (cross(r1, r2) + cross(r2, r3)).Normalize();
			face_normal[i] = nor;
		}
	}

	inline void CalculateVertexNormal() {
		//	create temp vertex normal array and calculate normal
		vertex_normal.clear();
		vertex_normal.resize(vertex.size(), Vec3<T>(0, 0, 0));
		for (int i = 0; i<face.size(); i++) {
			for (int j = 0; j < 4; j++) {
				Vec3<T> r1 = vertex[face[i][(j + 1) % 4]] - vertex[face[i][j]];
				Vec3<T> r2 = vertex[face[i][(j + 3) % 4]] - vertex[face[i][j]];
				Vec3<T> nor = cross(r1, r2).Normalize();
				vertex_normal[face[i][j]] += nor;
			}
		}
		for (int i = 0; i<vertex_normal.size(); i++)
			vertex_normal[i].SetNormalize();
	}

	virtual inline void AppendQuadMesh(const QuadMesh<T>& mesh) {
		int old_vertex_number = vertex.size();
		int old_face_number = face.size();
		vertex.insert(vertex.end(), mesh.vertex.begin(), mesh.vertex.end());
		face.insert(face.end(), mesh.face.begin(), mesh.face.end());
		for (int f = old_face_number; f<face.size(); f++) {
			face[f].x += old_vertex_number;
			face[f].y += old_vertex_number;
			face[f].z += old_vertex_number;
			face[f].w += old_vertex_number;
		}
	}

	/**
	draw the mesh with flat shading.

	We maintain a copy of vertex to display flat mesh.
	If position of vertex changes, we need to set force_refresh_flag = 1

	\param	force_refresh_flag		force update, used for dynamic mesh
	*/
	virtual int DisplayFlat(int force_refresh_flag = 0) {
		if (vertex.empty() || face.empty())
			return 1;
#ifdef YZ_gl_h
		if (disp_split_vertex.size() != face.size() * 4)
			force_refresh_flag = 1;
		if (disp_split_vertex_normal.size() != face.size() * 4)
			force_refresh_flag = 1;
		if (face_normal.size() != face.size())
			CalculateFaceNormal();

		if (force_refresh_flag) {
			disp_split_vertex.resize(face.size() * 4);
			disp_split_vertex_normal.resize(face.size() * 4);
			for (unsigned int f = 0; f != face.size(); f++) {
				for (unsigned int i = 0; i != 4; i++) {
					disp_split_vertex[f * 4 + i] = vertex[face[f][i]];
					disp_split_vertex_normal[f * 4 + i] = face_normal[f];
				}
			}
		}

		if (disp_indices.size() < face.size() * 4) {
			disp_indices.reserve(face.size() * 4);
			unsigned int count = face.size() * 4 - disp_indices.size();
			do {
				disp_indices.push_back(disp_indices.size());
			} while (--count);
		}

		opengl::DrawElements(
			GL_QUADS,
			&disp_split_vertex[0].x,
			&disp_split_vertex_normal[0].x,
			(float*)NULL,
			&disp_indices[0],
			face.size() * 4
		);

		return 1;
#else
		std::cout << "gl.h has to be included in order to use DisplayFlat() in QuadMesh" << std::endl;
		return 0;
#endif
	}

	virtual int DisplaySmooth() {
		if (vertex.empty() || face.empty())
			return 1;
#ifdef YZ_gl_h
		if (vertex_normal.size() != vertex.size())
			CalculateVertexNormal();

		if (disp_indices.size() < face.size() * 4) {
			disp_indices.reserve(face.size() * 4);
			unsigned int count = face.size() * 4 - disp_indices.size();
			do {
				disp_indices.push_back(disp_indices.size());
			} while (--count);
		}

		opengl::DrawElements(
			GL_QUADS,
			&vertex[0].x,
			&vertex_normal[0].x,
			(float*)NULL,
			&face[0].x,
			face.size() * 4
		);

		return 1;
#else
		std::cout << "gl.h has to be included in order to use DisplaySmooth() in QuadMesh" << std::endl;
		return 0;
#endif
	}

	virtual int Display(int display_mode = 0x1D01, int force_refresh_flag = 0) {	//	default value: GL_SMOOTH
		if (display_mode == 0x1D01)
			DisplaySmooth();
		else
			DisplayFlat(force_refresh_flag);

		return 0;
	}

};

typedef QuadMesh<float>		QuadMeshf;
typedef QuadMesh<double>	QuadMeshd;



}}	//	namespace yz::geometry

#endif	//	__YZ_QUAD_MESH__