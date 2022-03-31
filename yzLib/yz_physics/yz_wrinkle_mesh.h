/***********************************************************/
/**	\file
	\brief		Position Based Dynamics
	\details	Implementation of the paper: Position based dynamics \n
				Matthias M¨¹ller, Bruno Heidelberger, Marcus Hennix, John Ratcliff \n
				3rd Workshop in Virtual Reality Interactions and Physical Simulation "VRIPHYS" (2006)
	\author		Yizhong Zhang
	\date		10/17/2012
*/
/***********************************************************/
#ifndef __YZ_WRINKLE_MESH_H__
#define __YZ_WRINKLE_MESH_H__

#include "yzLib/yz_physics/yz_position_based_dynamics.h"
#include "yzLib/yz_animation/yz_rigging_to_mesh.h"

namespace yz{	namespace physics{

/**
	Wrinkle mesh is improvement of position based dynamics

	for each vertex, we have attachment constraint to limit the position of the vertex

	see more in paper:
	Wrinkle Meshes
	Matthias Muller, Nuttapong Chentanez
*/
template<class T>
class WrinkleMesh : public PositionBasedDynamicTriMesh<T>{
public:
	using PositionBasedDynamicTriMesh<T>::vertex;
	using PositionBasedDynamicTriMesh<T>::p;

public:
	/**
		read mesh
	*/
	inline int ReadMeshFromMeshSubd(const std::vector<Vec3<T>>& vertex_attach, const std::vector<int3>& face_attach, int subd_times){
		//	create subd mesh
		std::vector<Vec3<T>> subd_vertex = vertex_attach;
		std::vector<int3> subd_face = face_attach;
		for(int i=0; i<subd_times; i++)
			geometry::subdTriMeshLoop(subd_vertex, subd_face);

		//	read mesh
		PositionBasedDynamicTriMesh<T>::ReadMeshFromMesh(subd_vertex, subd_face);

		//	create attachment
		CreateAttachment(vertex_attach, face_attach, 1);

		return 1;
	}

	/**
		create attachment
	*/
	inline int CreateAttachment(const std::vector<Vec3<T>>& vertex_attach, const std::vector<int3>& face_attach, int mode = 0){
		attach_mode = mode;
		rigging_to_mesh.SetRigging(vertex, vertex_attach, face_attach);

		int vertex_number = vertex.size();
		attach_type.resize(vertex_number);
		attach_points.resize(vertex_number);
		attach_radius.resize(vertex_number);

		for(int i=0; i<vertex_number; i++){
			attach_type[i] = attach_mode;
			attach_radius[i] = 0.01;
		}

		UpdateAttachment(vertex_attach, face_attach, mode);

		return 1;
	}

	/**
		update attach
	*/
	inline int UpdateAttachment(const std::vector<Vec3<T>>& vertex_attach, const std::vector<int3>& face_attach, int mode = 0){
		std::vector<Vec3<T>> vertex_normal_attach;
		yz::geometry::calculateVertexNormal(vertex_normal_attach, vertex_attach, face_attach);

		//	get vertex by consider the attach mesh to be a bezier surface
		attach_points.resize(rigging_to_mesh.rigging_coef.size());
		for(int i=0; i<attach_points.size(); i++){
			int face_id = rigging_to_mesh.rigging_coef[i].face_id;
			std::vector<Vec3<T>> control_points;
			geometry::remeshing::getCubicBezierTriangleControlPoints(control_points,
				vertex_attach[face_attach[face_id].x], 
				vertex_attach[face_attach[face_id].y], 
				vertex_attach[face_attach[face_id].z],
				vertex_normal_attach[face_attach[face_id].x], 
				vertex_normal_attach[face_attach[face_id].y], 
				vertex_normal_attach[face_attach[face_id].z]);

			//	read uv
			T u = rigging_to_mesh.rigging_coef[i].coef1;
			T v = rigging_to_mesh.rigging_coef[i].coef2;

			//	cubic bezier interpolate
			attach_points[i] = interpCubicBezierTriangle(&control_points[0], u, v);
		}


		//rigging_to_mesh.GetVertex(attach_points, vertex_attach, face_attach);

		if( attach_mode == 1 ){	//	one sided attach
			rigging_to_mesh.GetVertexData(attach_normal, vertex_normal_attach, face_attach);
		}

		return 1;
	}

	/**
		stretching, bending and attachment
	*/
	void ProjectConstraints(){
		//	stretching and bending, use old function
		PositionBasedDynamicTriMesh<T>::ProjectConstraints();

		//	attachment
		ProjectAttachmentConstraints();
	}

	/**
		if distance of p to attach_points is further than threshold, project it back
	*/
	void ProjectAttachmentConstraints(){
		if( p.size() != attach_points.size() || p.size() != attach_radius.size() )
			return;

		int vertex_number = p.size();
		for(int i=0; i<vertex_number; i++){
			if( attach_type[i] == 0 ){
				Vec3<T> r = p[i] - attach_points[i];
				T length = r.Length();
				if( length > attach_radius[i] ){	//	distance too far, move near
					p[i] = attach_points[i] + r * (attach_radius[i] / length);
				}
			}
			else if( attach_type[i] == 1 ){
				Vec3<T> r = p[i] - attach_points[i];
				Vec3<T> n = attach_normal[i].Normalize();
				T dot_val = dot(r, n);
				if( dot_val < 0 ){
					p[i] -= n * dot_val;
					r = p[i] - attach_points[i];
				}

				T length = r.Length();
				if( length > attach_radius[i] ){	//	distance too far, move near
					p[i] = attach_points[i] + r * (attach_radius[i] / length);
				}
			}
		}

	}

public:
	std::vector<int>		attach_type;	///<	0: normal;  1: one sided attachment
	std::vector<Vec3<T>>	attach_points;	///<	attach target, if attach type is 1, this is used
	std::vector<Vec3<T>>	attach_normal;	///<	normal of the attach point
	std::vector<T>			attach_radius;	///<	tolerance radius

	int							attach_mode;		///<	mode of attach
	animation::RiggingToMesh<T>	rigging_to_mesh;	///<	rigging to attach mesh

};



}}	//	namesapce yz::physics

#endif //	__YZ_WRINKLE_MESH_H__