/***********************************************************/
/**	\file
	\brief		Rigging to another mesh
	\details	Sometimes, we need to attach vertices to another
				mesh, then we need rigging to mesh, similar 
				to rigging to skeleton. When the base mesh deforms,
				the attached vertices deform with the mesh
	\author		Yizhong Zhang
	\date		11/12/2012
*/
/***********************************************************/
#ifndef __YZ_RIGGING_TO_MESH_H__
#define __YZ_RIGGING_TO_MESH_H__

#include <vector>
#include <algorithm>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"

namespace yz{	namespace animation{

//	========================================
///@{
/**	@name Rigging To Mesh
*/
//	========================================

/**
	Rigging Vertices to another mesh

	This class is used to calculate correspondence between two meshes,
	similar to parameterization.
	For source mesh, we calculate the nearest vertex of each vertex
	on the target mesh, and record the coordinate. So if the target mesh
	changed, we can update the source mesh according to the parameterization.

	To use this class, 
	1, call SetRigging() to calculate parameterization from source mesh to target mesh
	2, call GetVertex() to calculate new vertex position of the source mesh
*/
template<class T>
class RiggingToMesh{
public:
	struct Param{
		int			face_id;		///<	index of the rigged face
		T			coef1;			///<	coefficient point projection on edge v0-v1
		T			coef2;			///<	coefficient point projection on edge v0-v2
		yz::Vec3<T>	offset;			///<	offset of the vertex along face normal
	};

	std::vector<Param>	rigging_coef;

public:
	/**
		calculate rigging parameters of vertices to the mesh

		\param	vertex_to_rig	vertices that is to be rigged
		\param	vertex			target mesh vertex
		\param	face			target mesh face
	*/
	void SetRigging(std::vector<yz::Vec3<T>>&		vertex_to_rig,
					const std::vector<yz::Vec3<T>>&	vertex,
					const std::vector<yz::int3>&	face){
		//	calculate aabb tree
		geometry::AABBTree3D<T>	face_aabb;
		face_aabb.BuildTriangleAABBTree(vertex, face);

		//	calculate rigging parameters
		rigging_coef.resize( vertex_to_rig.size() );
		for(int i=0; i<vertex_to_rig.size(); i++){
			yz::Vec3<T>	nearest_point;
			int			face_id;
			yz::geometry::getNearestPointOnMesh(nearest_point, face_id, vertex_to_rig[i], vertex, face, face_aabb);
			rigging_coef[i].face_id = face_id;

			yz::geometry::getPointTriangleProjectionCoef(rigging_coef[i].coef1, rigging_coef[i].coef2,
				nearest_point, vertex[face[face_id].x], vertex[face[face_id].y], vertex[face[face_id].z] );

			rigging_coef[i].offset = vertex_to_rig[i] - nearest_point;
		}
	}

	/**
		get vertex position from rigging

		after rigging has been set, we can get new position of vertices

		\param	vertex_to_get	vertices that is to be rigged
		\param	vertex			target mesh vertex
		\param	face			target mesh face
	*/
	void GetVertex(	std::vector<yz::Vec3<T>>&		vertex_to_get,
					const std::vector<yz::Vec3<T>>&	vertex,
					const std::vector<yz::int3>&	face){
		//	update vertex position
		vertex_to_get.resize(rigging_coef.size());
		for(int i=0; i<vertex_to_get.size(); i++){
			int			face_id = rigging_coef[i].face_id;
			yz::Vec3<T>	v0 = vertex[face[face_id].x];
			yz::Vec3<T>	v1 = vertex[face[face_id].y];
			yz::Vec3<T>	v2 = vertex[face[face_id].z];
			yz::Vec3<T> p = v0 + (v1-v0)*rigging_coef[i].coef1 + (v2-v0)*rigging_coef[i].coef2;
			vertex_to_get[i] = p + rigging_coef[i].offset;
		}
	}

	/**
		Get any interpolated value by interpolating on the mesh

		Besides vertex position, we can get interpolated values of any parameter.
		Data type E must have E*T E+E defined

		\param	data_to_get		interpolated data
		\param	vertex_data		data on each vertex of the target mesh
		\param	face			face of the target mesh
	*/
	template<typename E>
	void GetVertexData(	std::vector<E>&					data_to_get,
						const std::vector<E>&			vertex_data,
						const std::vector<yz::int3>&	face){

		//	update vertex normal
		data_to_get.resize(rigging_coef.size());		
		for(int i=0; i<data_to_get.size(); i++){
			int			face_id = rigging_coef[i].face_id;
			yz::Vec3<T>	d0 = vertex_data[face[face_id].x];
			yz::Vec3<T>	d1 = vertex_data[face[face_id].y];
			yz::Vec3<T>	d2 = vertex_data[face[face_id].z];
			yz::Vec3<T> d = d0 + (d1-d0)*rigging_coef[i].coef1 + (d2-d0)*rigging_coef[i].coef2;
			data_to_get[i] = d;
		}
	}

	/**
		Read the rigging information from file
	*/
	int ReadRiggingFromFile(const char* file_name){
		std::ifstream file( file_name );	//	load obj failed, do not touch old data
		if(!file.is_open()){
			std::cout << "cannot open " << file_name << std::endl;
			return 0;
		}


		int size;
		file >> size;
		rigging_coef.resize(size);

		for(int i=0; i<size; i++){
			file >> rigging_coef[i].face_id >> rigging_coef[i].coef1 >> rigging_coef[i].coef2 
				>> rigging_coef[i].offset[0] >> rigging_coef[i].offset[1] >> rigging_coef[i].offset[2] ;
		}

		file.close();
	}

	/**
		Write the rigging information to file
	*/
	int WriteRiggingToFile(const char* file_name){
		std::ofstream file( file_name );
		if(!file.is_open()){
			#ifndef BE_QUIET
				std::cout << "error: RiggingToMesh::WriteRiggingToFile, cannot open " << file_name << std::endl;
			#endif
			return 0;
		}

		file << rigging_coef.size() << std::endl;

		for(int i=0; i<rigging_coef.size(); i++){
			file << rigging_coef[i].face_id << "\t"
				<< rigging_coef[i].coef1 << "\t"
				<< rigging_coef[i].coef2 << "\t"
				<< rigging_coef[i].offset[0] << "\t" 
				<< rigging_coef[i].offset[1] << "\t" 
				<< rigging_coef[i].offset[2] << std::endl;
		}

		file.close();
	}


};


///@}

}}	//	namespace yz::animation

#endif	//	__YZ_RIGGING_TO_MESH_H__