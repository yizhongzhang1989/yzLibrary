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
#ifndef __YZ_POSITION_BASED_DYNAMICS_H__
#define __YZ_POSITION_BASED_DYNAMICS_H__

#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_geometry/yz_tri_mesh.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"
#include "yzLib/yz_geometry/yz_mesh_curvature.h"
#include "yzLib/yz_geometry/yz_mesh_transform.h"
#include "yzLib/yz_geometry/yz_mesh_intersection_test.h"
#include "yzLib/yz_geometry/yz_mesh_nn.h"

namespace yz{	namespace physics{


/**
	An implimentation of Position Based Dynamics

	Follow the following simulation steps:	\n
	1, [set force]						\n
	2, call SimulateStepPart1()			\n
	3, [set collisions]					\n
	4, call SimulateStepPart2()			\n
	5, [explicit solve collisions]		\n
	6, call SimulateStepPart3()			\n

	\snippet	position_based_dynamics.cpp		Position Based Dynamics Simulation
*/
template<class T>
class PositionBasedDynamicTriMesh : 
	public PhysicsBase<T>,
	public PointMass<T>,
	public PointAttachConstraint,
	public geometry::TriMesh<T>, 
	public geometry::TriMeshEdgeEFFE 
{
public:
	using PhysicsBase<T>::gravity;
	using PointMass<T>::M;
	using PointMass<T>::x;
	using PointMass<T>::v;
	using PointMass<T>::f;
	using geometry::TriMesh<T>::vertex;
	using geometry::TriMesh<T>::face;
	using geometry::TriMesh<T>::vertex_normal;
	using geometry::TriMesh<T>::face_normal;

public:
	/**
		constructor, set default parameters
	*/
	PositionBasedDynamicTriMesh(){
		SetGravity(Vec3<T>(0, -9.8, 0));

		stretch_mode		= 0;
		time_step			= 0.001;
		solver_iterations	= 5;

		density			= 1;		//	1 kg per square meter
		k_streching		= 1;		//	the cloth is rigid
		k_bending		= 1;
		k_friction		= 0;		//	default no friction
		thickness		= 0.001;	//	thickness of the mesh

		compress_threshold	= 0.1;
		k_compress_phase1	= 0.1;
		k_compress_phase2	= 1;
	}

	/**
		read the mesh to simulate
	*/
	inline int ReadMeshFromFile(const char* file_name){
		//	setup the mesh and topology
		if( !geometry::TriMesh<T>::ReadMeshFromFile(file_name) ){
			#ifndef BE_QUIET
				std::cout << "error: PositionBasedDynamicTriMesh::ReadMeshFromFile" << std::endl;
			#endif
			return 0;
		}
		geometry::applyNoiseToVertices(vertex, 1e-5);
		geometry::TriMeshEdgeEFFE::CreateTopology(vertex.size(), face);

		//	setup arrays related to physics
		SetupPhysicalValues();

		//	setup aabb tree
		face_aabb.BuildTriangleAABBTree(vertex, face, thickness);

		return 1;
	}

	/**
		read the mesh to simulate from another mesh
	*/
	inline int ReadMeshFromMesh(
		const std::vector<Vec3<T>>& vertex_data, 
		const std::vector<int3>&	face_data)
	{
		vertex	= vertex_data;
		face	= face_data;

		this->CalculateNormals();
		
		geometry::applyNoiseToVertices(vertex, 1e-5);
		geometry::TriMeshEdgeEFFE::CreateTopology(vertex.size(), face);

		//	setup arrays related to physics
		SetupPhysicalValues();

		//	setup aabb tree
		face_aabb.BuildTriangleAABBTree(vertex, face, thickness);

		return 1;
	}

	/**
		if no collision exist, this function can be called
	*/
	virtual int SimulateStep(T dt) {
		time_step = dt;

		this->SimulateStepPart1();
		this->SimulateStepPart2();
		this->SimulateStepPart3();

		return 1;
	}

	/**
		Step 1 of simulation

		Before step 1, external force should be set. 
	*/
	virtual void SimulateStepPart1(){
		for( int i=0; i<v.size(); i++ )
			v[i] = v[i] + time_step * weight[i] * f[i];

		//	damp velocity according to last collision
		//for( int i=0; i<collision_constraint.size(); i++ ){
		//	int v_id = collision_constraint[i].vertex_id;
		//	Vec3<T> n = collision_constraint[i].intersection_normal;
		//	Vec3<T>	v_n = n * dot(v[v_id], n);
		//	Vec3<T> v_t = v[v_id] - v_n;
		//	v_t *= (1 - k_friction);
		//	v[v_id] = v_t + v_n;
		//}
		ClearCollisionConstraints();	//	before each step, collision constraints are cleared
		
		for( int i=0; i<x.size(); i++ )
			p[i] = x[i] + time_step * v[i];
	}

	/**
		Step 2 of simulation

		Before step 2, collision constraint should be set.
		Original collision constraints have been removed in step 1,
		so just add all collision constraints
	*/
	virtual void SimulateStepPart2(){
		//	attach constraint
		for (int i = 0; i < x.size(); i++) {
			if (attach_flag[i]) {
				p[i] = attach_position[i];
			}
		}

		//	solver
		for(int i=0; i<solver_iterations; i++){
			ProjectConstraints();
		}
	}

	/**
		Step 3 of simulation

		Before step 3, all collisions should be resolved
	*/
	virtual void SimulateStepPart3(){
		//	update velocity and positon
		for( int i=0; i<x.size(); i++ ){
			v[i] = (p[i] - x[i]) / time_step;
			x[i] = p[i];
		}

		//	apply friction by explicitly set velocity
		for( int i=0; i<collision_constraint.size(); i++ ){
			int v_id = collision_constraint[i].vertex_id;
			Vec3<T> n = collision_constraint[i].intersection_normal;
			Vec3<T> c_v = collision_constraint[i].intersection_velocity;
			Vec3<T>	v_n = n * dot(v[v_id], n);
			Vec3<T> v_t = v[v_id] - v_n;
			c_v = c_v - n * dot(n, c_v);
			v_t = (1 - k_friction) * v_t + k_friction * c_v;
			v[v_id] = v_t + v_n;
		}

		vertex	= x;
		face_aabb.UpdateTriangleAABBTree(vertex, face, thickness);

		this->CalculateNormals();
	}

	/**
		Add a vertex to be constraint at current position. Ignore if this constraint already exist
	*/
	virtual inline void AddAttachConstraint(int vertex_id){
		if (size_t(vertex_id) >= vertex.size())			//	illegal vertex id
			return;
		PointAttachConstraint::AddAttachConstraint(vertex_id);
		attach_position[vertex_id] = x[vertex_id];
		weight[vertex_id] = 0;
	}

	/**
		Add a vertex to a given position. Set new constraint position if this vertex is already constraint
	*/
	virtual inline void AddAttachConstraint(int vertex_id, Vec3<T> attach_pos){
		if (size_t(vertex_id) >= vertex.size())			//	illegal vertex id
			return;
		PointAttachConstraint::AddAttachConstraint(vertex_id);
		attach_position[vertex_id] = attach_pos;
		weight[vertex_id] = 0;
	}

	/**
		Get the attached position of a given vertex

		\return	0:	the vertex is not attached
				-1:	the input vertex is illegal
				1:	the correct attach_position is returned
	*/
	virtual inline int GetAttachedPosition(yz::Vec3<T>& attach_pos, int vertex_id){
		if (size_t(vertex_id) >= vertex.size())			//	illegal vertex id
			return -1;
		if (!attach_flag[vertex_id])
			return -1;
		attach_pos = attach_position[vertex_id];
		return 0;
	}

	/**
		Remove the attach constraint of a given vertex
	*/
	virtual inline void RemoveAttachConstraint(int vertex_id){
		if (size_t(vertex_id) >= vertex.size())			//	illegal vertex id
			return;
		attach_flag[vertex_id] = 0;
		weight[vertex_id] = 1.0 / M[vertex_id];
	}

	/**
		Add collision constraint to a static mesh

		\param	vertex		vertex vector of the rigid mesh
		\param	face		face of the rigid mesh
		\param	face_aabb	bvh of the rigid mesh
	*/
	virtual inline void AddCollisionConstraintsWithClosedTriMesh(
		const std::vector<Vec3<T>>&		vertex,
		const std::vector<int3>&		face,
		const geometry::AABBTree3D<T>&	face_aabb)
	{
		for(int i=0; i<p.size(); i++){
			if( geometry::isPointInsideClosedMesh(p[i], vertex, face, face_aabb) ){	//	p is inside the mesh
				if( geometry::isPointInsideClosedMesh(x[i], vertex, face, face_aabb) ){	
					//	both x and p are inside the mesh
					//	find the nearesr point of p on the mesh
					yz::Vec3<T> nearest_point;
					int			face_id;
					yz::geometry::getNearestPointOnMesh(nearest_point, face_id, p[i], vertex, face, face_aabb);

					AddCollisionConstraint(i, nearest_point, geometry::calculateFaceNormal(face_id, vertex, face));
				}
				else{	//	x->p penetrate the mesh
					std::vector<yz::Vec3<T>>	intersection_points;
					std::vector<int>			face_id;
					yz::geometry::getSegmentMeshIntersectionPoints(intersection_points, face_id,
						x[i], p[i], vertex, face, face_aabb);

					if (intersection_points.size() == 1)
						AddCollisionConstraint(i, intersection_points[0], geometry::calculateFaceNormal(face_id[0], vertex, face));
				}
			}
		}
	}

	/**
		Add collision constraint to a dynamic mesh(moving or deforming)

		\param	vertex				vertex vector of the collision mesh
		\param	face				face of the collision mesh
		\param	vertex_velocity		velocity of each vertex of the collision mesh
		\param	face_aabb			bvh of the collision mesh
	*/
	virtual inline void AddCollisionConstraintsWithClosedTriMesh(
		const std::vector<Vec3<T>>&		vertex,
		const std::vector<int3>&		face,
		const std::vector<Vec3<T>>&		vertex_velocity,
		const geometry::AABBTree3D<T>&	face_aabb)
	{
		for(int i=0; i<p.size(); i++){
			if( geometry::isPointInsideClosedMesh(p[i], vertex, face, face_aabb) ){	//	p is inside the mesh

				yz::Vec3<T> nearest_point;
				int			face_id;
				yz::geometry::getNearestPointOnMesh(nearest_point, face_id, p[i], vertex, face, face_aabb);

				T coef1, coef2;
				int id0 = face[face_id].x, id1 = face[face_id].y, id2 = face[face_id].z;
				yz::geometry::getPointTriangleProjectionCoef(coef1, coef2, nearest_point, 
					vertex[id0], vertex[id1], vertex[id2]);

				yz::Vec3<T> vel = vertex_velocity[id0] + coef1 * (vertex_velocity[id1] - vertex_velocity[id0])
					+ coef2 * (vertex_velocity[id2] - vertex_velocity[id0]);

				AddCollisionConstraint(i, nearest_point, geometry::calculateFaceNormal(face_id, vertex, face), vel);


				if( geometry::isPointInsideClosedMesh(x[i], vertex, face, face_aabb) ){	
					//	both x and p are inside the mesh
					//	find the nearesr point of p on the mesh
					yz::Vec3<T> nearest_point;
					int			face_id;
					yz::geometry::getNearestPointOnMesh(nearest_point, face_id, p[i], vertex, face, face_aabb);

					T coef1, coef2;
					int id0 = face[face_id].x, id1 = face[face_id].y, id2 = face[face_id].z;
					yz::geometry::getPointTriangleProjectionCoef(coef1, coef2, nearest_point, 
						vertex[id0], vertex[id1], vertex[id2]);

					yz::Vec3<T> vel = vertex_velocity[id0] + coef1 * (vertex_velocity[id1] - vertex_velocity[id0])
						+ coef2 * (vertex_velocity[id2] - vertex_velocity[id0]);

					AddCollisionConstraint(i, nearest_point, geometry::calculateFaceNormal(face_id, vertex, face), vel);
				}
				else{	//	x->p penetrate the mesh
					std::vector<yz::Vec3<T>>	intersection_points;
					std::vector<int>			face_id;
					yz::geometry::getSegmentMeshIntersectionPoints(intersection_points, face_id,
						x[i], p[i], vertex, face, face_aabb);

					if( intersection_points.size() == 1 ){
						T coef1, coef2;
						int id0 = face[face_id[0]].x, id1 = face[face_id[0]].y, id2 = face[face_id[0]].z;
						yz::geometry::getPointTriangleProjectionCoef(coef1, coef2, intersection_points[0], 
							vertex[id0], vertex[id1], vertex[id2]);

						yz::Vec3<T> vel = vertex_velocity[id0] + coef1 * (vertex_velocity[id1] - vertex_velocity[id0])
							+ coef2 * (vertex_velocity[id2] - vertex_velocity[id0]);

						AddCollisionConstraint(i, intersection_points[0], 
							geometry::calculateFaceNormal(face_id[0], vertex, face), vel);
					}
					else if( intersection_points.size() == 0 ){
						std::cout << "intersection failed" << std::endl;
					}
				}
			}
		}
	}

	/**
		Add collision constraint to a rigid open mesh

		\param	vertex					vertex vector of the collision mesh
		\param	face					face of the collision mesh
		\param	face_aabb				bvh of the collision mesh
		\param	colliding_normal_side	which side should this collision happen. 1: same side as normal; 0: opposite side as normal
	*/
	virtual inline void AddCollisionConstraintsWithOpenTriMesh(
		const std::vector<Vec3<T>>&		vertex,
		const std::vector<int3>&		face,
		const geometry::AABBTree3D<T>&	face_aabb,
		int								colliding_normal_side = 1)
	{
		for (int i = 0; i<p.size(); i++){
			yz::Vec3<T> nearest_point_p, nearest_point_x;
			int			face_id_p, face_id_x;
			yz::geometry::getNearestPointOnMesh(nearest_point_p, face_id_p, p[i], vertex, face, face_aabb);
			yz::geometry::getNearestPointOnMesh(nearest_point_x, face_id_x, p[i], vertex, face, face_aabb);
			yz::Vec3<T> n_p = geometry::calculateFaceNormal(face_id_p, vertex, face);
			yz::Vec3<T> n_x = geometry::calculateFaceNormal(face_id_x, vertex, face);

			//	check 
			if ((p[i] - nearest_point_p).Length() < thickness * 100 && (x[i] - nearest_point_x).Length() < thickness * 100 &&
				((dot(p[i] - nearest_point_p, n_p)>0 && dot(x[i] - nearest_point_x, n_x)>0 && colliding_normal_side == 0) ||
				(dot(p[i] - nearest_point_p, n_p)<0 && dot(x[i] - nearest_point_x, n_x)<0 && colliding_normal_side == 1))){
				AddCollisionConstraint(i, nearest_point_p, n_p);
			}
			else{	//	x->p penetrate the mesh
				std::vector<yz::Vec3<T>>	intersection_points;
				std::vector<int>			face_id;
				yz::geometry::getSegmentMeshIntersectionPoints(intersection_points, face_id,
					x[i], p[i], vertex, face, face_aabb);

				if (intersection_points.size() == 1)
					AddCollisionConstraint(i, intersection_points[0],
					geometry::calculateFaceNormal(face_id[0], vertex, face));
			}
		}
	}

	/**
		Explicitly add one collision constraint
	*/
	virtual inline void AddCollisionConstraint(
		int		vertex_id, 
		Vec3<T>	intersection_point,
		Vec3<T>	intersection_normal)
	{
		CollisionConstraint coll_cons;
		coll_cons.vertex_id				= vertex_id;
		coll_cons.intersection_point	= intersection_point;
		coll_cons.intersection_normal	= intersection_normal;
		collision_constraint.push_back(coll_cons);
	}

	/**
		Explicitly add one collision constraint
	*/
	virtual inline void AddCollisionConstraint(
		int		vertex_id,
		Vec3<T>	intersection_point,
		Vec3<T>	intersection_normal,
		Vec3<T>	intersection_velocity)
	{
		CollisionConstraint coll_cons;
		coll_cons.vertex_id				= vertex_id;
		coll_cons.intersection_point	= intersection_point;
		coll_cons.intersection_normal	= intersection_normal;
		coll_cons.intersection_velocity	= intersection_velocity;
		collision_constraint.push_back(coll_cons);
	}

public:
	/**
	*/
	virtual void SetFaceAABBTreeToP(){
		face_aabb.UpdateTriangleAABBTree(p, face, thickness);
	}

	/**
		Solve streching, bending, collision constraints one time in a row

		In real simulation, this function should be called multiple times
	*/
	virtual void ProjectConstraints(){
		//	strech constraint
		T	k_streching_prime = 1.0 - pow(1.0 - k_streching, 1.0 / solver_iterations);
		T	k_compress_phase1_prime = 1.0 - pow(1.0 - k_compress_phase1, 1.0 / solver_iterations);
		T	k_compress_phase2_prime = 1.0 - pow(1.0 - k_compress_phase2, 1.0 / solver_iterations);
		for (int i = 0; i<edge.size(); i++){
			int v1 = edge[i].x;
			int v2 = edge[i].y;
			if (weight[v1] + weight[v2] < 1e-6)
				continue;
			yz::Vec3<T>	p1 = p[v1];
			yz::Vec3<T>	p2 = p[v2];
			T len = (p1 - p2).Length();
			yz::Vec3<T> dp, dp1, dp2;
			if (stretch_mode == 0)			//	same stretch and compress
				dp = k_streching_prime / (weight[v1] + weight[v2]) * (len - rest_length[i]) * (p1 - p2) / len;
			else if (stretch_mode == 1){	//	two phase compress
				if (len > rest_length[i])												//	stretched
					dp = k_streching_prime / (weight[v1] + weight[v2]) * (len - rest_length[i]) * (p1 - p2) / len;
				else if (rest_length[i] - len < rest_length[i] * compress_threshold)	//	compressed, small
					dp = k_compress_phase1_prime / (weight[v1] + weight[v2]) * (len - rest_length[i]) * (p1 - p2) / len;
				else if (rest_length[i] - len >= rest_length[i] * compress_threshold){	//	compressed, big
					T adj_length = k_compress_phase2_prime * (len - rest_length[i]);
					if (-adj_length > rest_length[i] * compress_threshold)
						adj_length = -rest_length[i] * compress_threshold;
					dp = adj_length / (weight[v1] + weight[v2]) * (p1 - p2) / len;
				}
				else
					std::cout << "this should not happen" << std::endl;
			}
			else
				std::cout << "invalud stretch mode" << std::endl;

			dp1 = -weight[v1] * dp;
			dp2 = weight[v2] * dp;

			//	apply the constraint
			p[v1] += dp1;
			p[v2] += dp2;
		}

		//	bending constraint
		T	k_bending_prime = 1.0 - pow(1.0 - k_bending, 1.0 / solver_iterations);
		for (int e = 0; e<edge.size(); e++){
			int f1 = ef[e].x;
			int f2 = ef[e].y;
			if (f1 == -1 || f2 == -1)		//	bondary edge, no constraint here
				continue;
			else{
				int v1 = edge[e].x;
				int v2 = edge[e].y;
				int v3 = yz::geometry::getThirdVertexOfFace(face[f1], v1, v2);
				int v4 = yz::geometry::getThirdVertexOfFace(face[f2], v1, v2);

				yz::Vec3<T>	p1(0, 0, 0);
				yz::Vec3<T> p2 = p[v2] - p[v1];
				yz::Vec3<T> p3 = p[v3] - p[v1];
				yz::Vec3<T> p4 = p[v4] - p[v1];

				yz::Vec3<T> n1 = yz::cross(p2, p3).Normalize();
				yz::Vec3<T> n2 = yz::cross(p2, p4).Normalize();
				T d = yz::dot(n1, n2);
				d = yz::clamp(d, -1, 1);

				yz::Vec3<T> q3 = (yz::cross(p2, n2) + yz::cross(n1, p2)*d) / yz::cross(p2, p3).Length();
				yz::Vec3<T> q4 = (yz::cross(p2, n1) + yz::cross(n2, p2)*d) / yz::cross(p2, p4).Length();
				yz::Vec3<T> q2 = -(yz::cross(p3, n2) + yz::cross(n1, p3)*d) / yz::cross(p2, p3).Length()
					- (yz::cross(p4, n1) + yz::cross(n2, p4)*d) / yz::cross(p2, p4).Length();
				yz::Vec3<T> q1 = -q2 - q3 - q4;

				T q_sum = q1.SquareLength() + q2.SquareLength() + q3.SquareLength() + q4.SquareLength();
				T w_sum = weight[v1] + weight[v2] + weight[v3] + weight[v4];
				T coef = 0;
				if (q_sum > 1e-5 && w_sum > 1e-5)	//	in case q1-4 are zero
					coef = -4 * k_bending_prime*sqrt(1 - d*d)*(acos(d) - rest_angle[e]) / (w_sum*q_sum);
				yz::Vec3<T> dp1 = weight[v1] * coef * q1;
				yz::Vec3<T> dp2 = weight[v2] * coef * q2;
				yz::Vec3<T> dp3 = weight[v3] * coef * q3;
				yz::Vec3<T> dp4 = weight[v4] * coef * q4;

				//	apply constraint
				assert(dp1.x == dp1.x);
				assert(dp2.x == dp2.x);
				assert(dp3.x == dp3.x);
				assert(dp4.x == dp4.x);
				p[v1] += dp1;
				p[v2] += dp2;
				p[v3] += dp3;
				p[v4] += dp4;
			}
		}

		//	collision constraint
		for (int i = 0; i<collision_constraint.size(); i++){
			int vertex_id = collision_constraint[i].vertex_id;
			Vec3<T> n = collision_constraint[i].intersection_normal;
			T	dot_value = dot(p[vertex_id] - collision_constraint[i].intersection_point, n);
			if (dot_value < 0){
				p[vertex_id] -= n * dot_value / n.SquareLength();
			}
		}
	}

	/**
		Setup physical parameters used in simulation
	*/
	virtual inline void SetupPhysicalValues(){
		int vertex_number = vertex.size();
		int	edge_number = edge.size();

		//	setup physical parameters
		PointMass<T>::Init(vertex_number);
		PointAttachConstraint::Init(vertex_number);
		x = vertex;
		p.resize(vertex_number);
		attach_position.resize(vertex_number);

		//	calculate mixed area
		geometry::CurvTriMesh<T>	curv_mesh;
		curv_mesh.vertex = vertex;
		curv_mesh.face = face;
		curv_mesh.CalculateNormals();
		curv_mesh.CreateTopology();
		curv_mesh.CalculateMeanCurvature();

		//	calculate mass and weight
		weight.resize(vertex_number);
		mass_sum = 0;
		for (int i = 0; i<vertex_number; i++){
			M[i] = density * curv_mesh.mixed_area[i];
			weight[i] = 1.0 / M[i];
			mass_sum += M[i];
		}

		//	calculate rest length
		rest_length.resize(edge_number);
		for (int i = 0; i<edge_number; i++){
			int v1 = edge[i].x;
			int v2 = edge[i].y;
			rest_length[i] = (vertex[v1] - vertex[v2]).Length();
		}

		//	calculate rest angle
		rest_angle.resize(edge_number);
		for (int e = 0; e<edge_number; e++){
			int f1 = ef[e].x;
			int f2 = ef[e].y;
			if (f1 == -1 || f2 == -1)		//	bondary edge
				rest_angle[e] = 0;
			else{
				yz::Vec3<T>	n1 = face_normal[f1];
				yz::Vec3<T>	n2 = -face_normal[f2];
				T dot_val = yz::dot(n1, n2);
				rest_angle[e] = acos(clamp(dot_val, T(-1), T(1)));	//	we clamp to [-1, 1] to prevent numerical error causing illegal value
			}
		}

		SetForceToBeGravity();
	}

	/**
	*/
	virtual inline void SetForceToBeGravity(){
		for (int i = 0; i<f.size(); i++){
			f[i] = M[i] * gravity;
		}
	}

	/**
		clear all attach constraints
	*/
	virtual inline void ClearAttachConstraints(){
		for (int i = 0; i < attach_flag.size(); i++) {
			if (attach_flag[i]) {
				attach_flag[i] = 0;
				weight[i] = 1.0 / M[i];
			}
		}
	}

	/**
		clear all collision constraints
	*/
	virtual inline void ClearCollisionConstraints(){
		collision_constraint.clear();
	}

	/**
		Explicitly solve collision with mesh in case constraint is not enough

		\param	vertex		vertex vector of the rigid mesh
		\param	face		face of the rigid mesh
		\param	face_aabb	bvh of the rigid mesh
	*/
	virtual inline void ExplicitSolveCollision(
		const std::vector<Vec3<T>>&		vertex,
		const std::vector<int3>&		face,
		const geometry::AABBTree3D<T>&	face_aabb)
	{
		for (int i = 0; i<p.size(); i++){
			if (geometry::isPointInsideClosedMesh(p[i], vertex, face, face_aabb)){	//	p is inside the mesh

				yz::Vec3<T> nearest_point;
				int			face_id;
				yz::geometry::getNearestPointOnMesh(nearest_point, face_id, p[i], vertex, face, face_aabb);

				yz::Vec3<T> v0 = vertex[face[face_id].x];
				yz::Vec3<T> v1 = vertex[face[face_id].y];
				yz::Vec3<T> v2 = vertex[face[face_id].z];
				yz::Vec3<T> f_nor = yz::cross(v1 - v0, v2 - v0).Normalize();
				if (yz::geometry::isPointProjectionInsideTriangle(p[i], v0, v1, v2) &&
					yz::dot(f_nor, p[i] - nearest_point) < 0){
					p[i] = nearest_point;
				}
			}
		}
	}

	/**
		Explicitly solve collision with mesh in case constraint is not enough

		\param	vertex				vertex vector of the rigid mesh
		\param	face				face of the rigid mesh
		\param	vertex_velocity		velocity of each vertex of the collision mesh
		\param	face_aabb			bvh of the rigid mesh
	*/
	virtual inline void ExplicitSolveCollision(
		const std::vector<Vec3<T>>&		vertex,
		const std::vector<int3>&		face,
		const std::vector<Vec3<T>>&		vertex_velocity,
		const geometry::AABBTree3D<T>&	face_aabb)
	{
		for (int i = 0; i<p.size(); i++){
			if (geometry::isPointInsideClosedMesh(p[i], vertex, face, face_aabb)){	//	p is inside the mesh

				yz::Vec3<T> nearest_point;
				int			face_id;
				yz::geometry::getNearestPointOnMesh(nearest_point, face_id, p[i], vertex, face, face_aabb);

				T coef1, coef2;
				int id0 = face[face_id].x, id1 = face[face_id].y, id2 = face[face_id].z;
				yz::geometry::getPointTriangleProjectionCoef(coef1, coef2, nearest_point,
					vertex[id0], vertex[id1], vertex[id2]);

				yz::Vec3<T> vel = vertex_velocity[id0] + coef1 * (vertex_velocity[id1] - vertex_velocity[id0])
					+ coef2 * (vertex_velocity[id2] - vertex_velocity[id0]);

				yz::Vec3<T> v0 = vertex[face[face_id].x];
				yz::Vec3<T> v1 = vertex[face[face_id].y];
				yz::Vec3<T> v2 = vertex[face[face_id].z];
				yz::Vec3<T> f_nor = yz::cross(v1 - v0, v2 - v0).Normalize();
				if (yz::geometry::isPointProjectionInsideTriangle(p[i], v0, v1, v2) &&
					yz::dot(f_nor, p[i] - nearest_point) < 0){
					p[i] = nearest_point;
				}
			}
		}
	}

public:
	/**
		Constraint a vertex above a given plane, which used to handle collision with another mesh
	*/
	struct CollisionConstraint{
		int		vertex_id;				///<	index of the vertex
		Vec3<T>	intersection_point;		///<	intersection point of x-p with another mesh
		Vec3<T>	intersection_normal;	///<	normal of the intersecting mesh at intersection point
		Vec3<T>	intersection_velocity;	///<	velocity of the intersected mesh at the intersecting point
	};


public:
	//	simulation parameters
	int						stretch_mode;		///<	mode of compression, 0: as stretching, 1: two phase compression
	T						time_step;			///<	time step used in simulation
	int						solver_iterations;	///<	iterations of solving constraints each step

	//	parameters of cloth
	T						density;			///<	mass per square meter
	T						k_streching;		///<	strech coefficient, 
	T						k_bending;			///<	bending coefficient
	T						k_friction;			///<	friction coefficient
	T						thickness;			///<	the thickness of cloth

	T						compress_threshold;	///<	threshold of changing from phase1 to phase 2, precentage of compression
	T						k_compress_phase1;	///<	small compression, small coefficient
	T						k_compress_phase2;	///<	big compression, big coefficient

	//	dynamic values
	std::vector<yz::Vec3<T>>	p;				///<	new location

	//	fixed values of cloth parameters
	T							mass_sum;		///<	the mass of the whole object
	std::vector<T>				weight;			///<	inverse of mass of each vertex
	std::vector<T>				rest_length;	///<	rest length of each edge
	std::vector<T>				rest_angle;		///<	rest angle of each edge, meaningless for boundary edge

	//	constraints
	std::vector<Vec3<T>>				attach_position;		///<	attach position of each vertex (if attached)
	std::vector<CollisionConstraint>	collision_constraint;	///<	collision with another mesh

	//	bvh of this mesh
	geometry::AABBTree3D<T>	face_aabb;	///<	face aabb tree of the mesh

};


}}	//	end namespace yz::physics

#endif	//	__YZ_POSITION_BASED_DYNAMICS_H__