/***********************************************************/
/**	\file
	\brief		FEM on triangle mesh
	\details	
	\author		Yizhong Zhang
	\date		1/5/2016
*/
/***********************************************************/
#ifndef __YZ_FEM_TRI_MESH_H__
#define __YZ_FEM_TRI_MESH_H__

#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_tri_mesh.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"
#include "yzLib/yz_physics/yz_physics_base.h"
#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_opengl_utils.h"
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"
#endif


namespace yz{  namespace physics{  namespace fem{

/**
	Base of FEM on tri-mesh

	Basic Guideline:	\n
	1, call ReadMeshFromFile() to read the mesh, default physical parameters will be set accordingly	\n
	2, call SimulateStep(dt) to perform simulation. \n
		If dt > 0, we perform simulation with corresponding step. dt is recommended to be no larger than 0.01	\n
		If dt < 0, we perform static simulation to the equilibrium. \n
	3, call Draw() to display the whole system.	\n

	Advanced Guideline:	\n
	All functions in this class can be overloaded with your own implementation. 
	So that when calling SimulateStep(), the simulation process can be updated directly.	\n
	1, To change physical parameters of the model, overload SetupPhysicalValues() \n
	2, To change rest shape of the model, overload PreComputation() \n
	3, To fix some vertices to certain positions, call AddAttachConstraint() and RemoveAttachConstraint() \n
	4, To add collision and other external forces, inherent AddCollisionForce() and AddExternalForce()	\n

	\TODO	currently, bending is not included in the simulator.
*/
template<class T>
class FEMTriMesh : 
	public geometry::TriMesh<T>,
	public geometry::TriMeshTopology,
	public ImplicitPhysicsBase<T>,
	public PointMass<T>,
	public PointAttachConstraint
{
public:
	using geometry::TriMesh<T>::vertex;
	using geometry::TriMesh<T>::face;
	using geometry::TriMesh<T>::vertex_normal;
	using geometry::TriMesh<T>::face_normal;
	using ImplicitPhysicsBase<T>::gravity;
	using PointMass<T>::M;
	using PointMass<T>::x;
	using PointMass<T>::v;
	using PointMass<T>::f;

public:
	/**
		Read the mesh from file

		If read succeed, setup physical values accordingly.
		The initial shape of the mesh is set to be the rest shape.

		\return		whether succeed
	*/
	int ReadMeshFromFile(const char* file_name){
		int succ = geometry::TriMesh<T>::ReadMeshFromFile(file_name);

		if (succ){	//	setup physical values and ready for simulation
			CreateTopology(vertex.size(), face);
			AllocSpace();
			SetupPhysicalValues();
			PreComputation();
		}

		return succ;
	}

	/**
		setup the mesh from existing
	*/
	int ReadMeshFromMesh(
		const std::vector<Vec3<T>>& vertex_data,
		const std::vector<int3>&	face_data)
	{
		vertex = vertex_data;
		face = face_data;
		this->CalculateNormals();

		CreateTopology(vertex.size(), face);
		AllocSpace();
		SetupPhysicalValues();
		PreComputation();

		return 1;
	}

	/**
		Simulate one time step

		\param	dt					positive time step, recommended to be no larger than 0.01.
									negative time step, quasi-static relaxiation
		\param	max_newton_steps	max Newton steps performed in simulation, default 20 steps
		\return						newton iterations actually performed
	*/
	virtual inline int SimulateStepImplicit(T dt, int max_newton_steps){
		//	initialize v and x as constant velocity motion
		//	CAUTION! the initialization of x and v in Eftychios D. Sifakis' course P34 is wrong
		prev_v = v;
		for (int i = 0; i < v.size(); i++)
			x[i] = vertex[i] + v[i] * dt;

		int newton_count = 0;
		for (int i = 0; i < max_newton_steps; i++){
			newton_count++;
			T step_size = NewtonStep(dt);
			if (step_size < 1e-8)
				break;
		}

		vertex = x;
		this->CalculateNormals();

		return newton_count;
	}

	/**
		set gravity of the system
	*/
	virtual inline void SetGravity(Vec3<T> g = Vec3<T>(0, -9.8, 0)){
		gravity = g;

		//	calculate gravity of each vertex
		geometry::CurvTriMesh<T>	curv_mesh;
		curv_mesh.vertex = vertex;
		curv_mesh.face = face;
		curv_mesh.CalculateNormals();
		curv_mesh.CreateTopology();
		curv_mesh.CalculateMeanCurvature();

		T area_sum = 0;
		for (int i = 0; i < curv_mesh.mixed_area.size(); i++)
			area_sum += curv_mesh.mixed_area[i];
		for (int i = 0; i < G.size(); i++)
			G[i] = mass * gravity * curv_mesh.mixed_area[i] / area_sum;
	}

	/**
	Add a vertex to a given position.

	\param	vertex_id			id of the vertex to be constraint
	\param	attach_position		the target position to attach the vertex
	*/
	virtual inline void AddAttachConstraint(int vertex_id) {
		AddAttachConstraint(vertex_id, vertex[vertex_id]);
	}

	/**
		Add a vertex to a given position.

		\param	vertex_id			id of the vertex to be constraint
		\param	attach_position		the target position to attach the vertex
	*/
	virtual inline void AddAttachConstraint(int vertex_id, Vec3<T> attach_position){
		if (size_t(vertex_id) >= vertex.size())			//	illegal vertex id
			return;

		attach_flag[vertex_id] = 1;
		vertex[vertex_id] = attach_position;
		x[vertex_id] = attach_position;
	}

	/**
		Detatch the attachment of a vertex

		\param	vertex_id		id of the vertex that is to remove constraint
	*/
	virtual inline void RemoveAttachConstraint(int vertex_id){
		if (size_t(vertex_id) >= vertex.size())			//	resolve illegal vertex id
			return;

		attach_flag[vertex_id] = 0;
	}

	/**
		Clear all attach constraint
	*/
	virtual inline void ClearAttachConstraint(){
		for (int i = 0; i < attach_flag.size(); i++)
			attach_flag[i] = 0;
	}

	/**
		Draw the system as a mesh
	*/
	virtual inline void Draw(){
#ifdef YZ_gl_h
		glColor3f(1, 1, 1);
		this->DisplayFlat();

		glColor3f(0, 0, 1);
		opengl::drawMeshEdgeFromFace(vertex, face, vertex_normal);

		for (int i = 0; i < vertex.size(); i++){
			if (attach_flag[i] )
				glColor3f(1, 0, 0);
			else
				glColor3f(1, 1, 0);

			opengl::drawPointAsSphere(vertex[i], 0.005);
		}
#else
		std::cout << "gl.h has to be included in order to use Draw() in FEMTriMesh" << std::endl;
		return;
#endif
	}

	virtual inline void DrawStrainPrincipleDirection(){
	}

	virtual inline void DrawStrain(){
	}

	virtual inline void DrawStress(){

	}

	virtual inline void DrawVonMisesStress(T min_color_val = 1.0, T max_color_val = 2.0){
#ifdef YZ_gl_h
		if (!CalculateVonMises())
			return;

		std::vector<uchar3>	vertex_color;
		vertex_color.resize(vertex.size());
		utils::convertToColorArraySimple(
			(unsigned char*)&vertex_color[0].x,
			&VonMises_vertex[0],
			vertex.size(),
			min_color_val,
			max_color_val
			);

		opengl::drawSmoothColorTriMesh(
			vertex,
			face,
			vertex_color);
#else
		std::cout << "gl.h has to be included in order to use DrawVonMisesStress() in FEMTriMesh" << std::endl;
		return;
#endif
	}

	virtual inline void DrawThickness(){

	}

	virtual inline void DrawGravity(){
#ifdef YZ_gl_h
		glColor3f(1, 0, 1);
		for (int i = 0; i < vertex.size(); i++){
			yz::opengl::drawArrow(vertex[i], vertex[i] + G[i]*1000, 0.001);
		}
#else
		std::cout << "gl.h has to be included in order to use DrawGravity() in FEMTriMesh" << std::endl;
		return;
#endif
	}

	virtual inline void PickingDraw(){
#ifdef YZ_gl_h
		for (int i = 0; i < vertex.size(); i++){
			opengl::setPickingIndex(i);
			opengl::drawPointAsSphere(vertex[i], 0.05);
		}
#else
		std::cout << "gl.h has to be included in order to use PickingDraw() in FEMTriMesh" << std::endl;
		return;
#endif
	}

public:
	/**
		Alloc all space required for simulation
	*/
	virtual inline void AllocSpace(){
		PointMass<T>::Init(vertex.size());
		PointAttachConstraint::Init(vertex.size());

		Dm.resize(face.size());
		Bm.resize(face.size());
		W.resize(face.size(), 0);
		Ds.resize(face.size());
		F.resize(face.size());
		E.resize(face.size());

		G.resize(vertex.size());
		prev_v.resize(vertex.size());

		//	visualization array
		VonMises_face.resize(face.size());
		VonMises_vertex.resize(vertex.size());
	}

	/**
	Setup all physical values required in simulation
	*/
	virtual inline void SetupPhysicalValues(){
		thickness = 0.001;				//	1 mm
		mass = 0.015;					//	15 g
		youngs_modulus = 1.0;			//	Young's modulus
		poisson_ratio = 0.4;			//	Poisson ratio
		gravity = Vec3<T>(0, -9.8, 0);	//	-y direction

		damping_k = 0;
		damping_mass = 0;
		bending_k = 0;

		mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
		lambda = (youngs_modulus*poisson_ratio) / ((1.0 + poisson_ratio)*(1.0 - 2.0*poisson_ratio));
		I.SetIdentity();

		//	calculate gravity of each vertex
		std::vector<T> mixed_area;
		T area_sum = geometry::calculateMixedArea(mixed_area, vertex, face);
		for (int i = 0; i < G.size(); i++){
			M[i] = mass * mixed_area[i] / area_sum;
			G[i] = M[i] * gravity;
		}
	}

	/**
	Pre-compute all variables that doesn't change
	*/
	virtual inline void PreComputation(){
		for (int i = 0; i < face.size(); i++){
			Vec3<T> r1 = vertex[face[i][1]] - vertex[face[i][0]];
			Vec3<T> r2 = vertex[face[i][2]] - vertex[face[i][0]];
			T		r1_len = r1.Length();
			Vec3<T> nr1 = r1 / r1_len;
			T		r2x = dot(r2, nr1);
			T		r2y = (r2 - nr1 * r2x).Length();

			Dm[i] = Matrix2x2<T>(r1_len, r2x, 0, r2y);
			Bm[i] = Dm[i].Inverse();
			W[i] = fabs(Dm[i].Det()) / 2.0;
		}
	}

	/**
		A Newton step

		\param	dt	time step, if dt <= 0, we perform quasi static simulation
		\return		residual: rTr
	*/
	virtual inline T NewtonStep(T dt){
		int quasi_static_flag = (dt <= 0 ? 1 : 0);
		T inv_dt = (quasi_static_flag ? 0.0 : 1.0 / dt);
		std::vector<Vec3<T>> b;
		b.resize(x.size());

		//	compute right hand side
		ComputeStrain();
		ClearForce();
		AddElasticForce();
		AddCollisionForce();
		AddExternalForce();
		AddGravity();
		if (!quasi_static_flag)
			AddDampingForce();

		for (int i = 0; i < f.size(); i++){		//	calculate right hand side
			if (attach_flag[i])		//	this vertex is constraint
				b[i] = Vec3<T>(0, 0, 0);
			else
				b[i] = f[i] + ( quasi_static_flag ? 0.0 : inv_dt * M[i] * (prev_v[i] - v[i]) );
		}

		//	compute left hand side and iterate with CG
		//	I use the same notation as wikipedia, https://en.wikipedia.org/wiki/Conjugate_gradient_method
		std::vector<Vec3<T>>	dx;
		dx.resize(x.size());

		std::vector<Vec3<T>>	r = b;
		std::vector<Vec3<T>>	p = r;
		T rTr = XTY(r, r);

		int noise_flag = 0;
		int step_count = 0;
		for (int k = 0; k < 20; k++){
			std::vector<Vec3<T>> Ap;
			std::vector<Vec3<T>> df;
			if (quasi_static_flag)  	//	in quasi-static, only elastic force exist
				ComputeElasticDifferentialForce(df, p);
			else
				ComputeElasticDampingDifferentialForce(df, p, dt);
			Ap.swap(df);
			ComputeCollisionDifferentialForce(df, p);
			for (int i = 0; i < Ap.size(); i++)
				Ap[i] += df[i];

			T inv_dt_squ = inv_dt * inv_dt;
			for (int i = 0; i < Ap.size(); i++){
				Ap[i] = inv_dt_squ * M[i] * p[i] - Ap[i];
			}

			for (int i = 0; i < attach_flag.size(); i++){	//	Add constraint
				if (attach_flag[i])
					Ap[i] = Vec3<T>(0, 0, 0);
			}

			T pTAp = XTY(p, Ap);

			T alpha = rTr / pTAp;

			if (isnan(alpha) || isinf(alpha)) {
				//	degenerate case appear, we add noise to dx
				for (int i = 0; i < dx.size(); i++) {
					if (!attach_flag[i])
						dx[i] += Vec3<T>(randNumber(0, 1e-4), randNumber(0, 1e-4), randNumber(0, 1e-4));
				}
				noise_flag = 1;
				break;
			}

			for (int i = 0; i < dx.size(); i++)
				dx[i] += alpha * p[i];
			for (int i = 0; i < r.size(); i++)
				r[i] -= alpha * Ap[i];

			step_count++;

			T old_rTr = rTr;
			rTr = XTY(r, r);

			//	exit loop if r is small enough
			if (rTr < 1e-8)
				break;

			T beta = rTr / old_rTr;

			if (beta != beta)
				break;

			for (int i = 0; i < p.size(); i++)
				p[i] = r[i] + beta * p[i];
		}	//	end of CG iteration 

			//	if the error is already very small, return directly
		T dxb = XTY(dx, b);
		if (dxb < 1e-8) {
			if (quasi_static_flag) {
				for (int i = 0; i < x.size(); i++) {
					v[i] = Vec3<T>(0, 0, 0);
					x[i] += dx[i];
				}
			}
			else {
				for (int i = 0; i < x.size(); i++) {
					v[i] += dx[i] * inv_dt;
					x[i] += dx[i];
				}
			}
			return noise_flag ? 1e3 : dxb;	//	to apply noise, the actual error is not small
		}

		//	if the error is big, perform line search to determine step
		T E_prev = TotalEnergy();
		T alpha = 1.0;

		std::vector<Vec3<T>> old_v, old_x;
		old_v.swap(v);
		old_x.swap(x);
		v.resize(old_v.size());
		x.resize(old_x.size());

		while (alpha > 1e-2) {
			for (int i = 0; i < x.size(); i++) {
				v[i] = old_v[i] + alpha * dx[i] * inv_dt;
				x[i] = old_x[i] + alpha * dx[i];
			}

			T E_curr = TotalEnergy();

			if (E_curr <= E_prev + 1e-4 * alpha * dxb)	//	we use coefficient on wikipedia
				break;

			alpha *= 0.5;
		}

		if (quasi_static_flag)
			for (int i = 0; i < x.size(); i++)
				v[i] = Vec3<T>(0, 0, 0);

		return alpha * dxb;
	}

	/**
		Compute deformation gradient and Green strain of each face
	*/
	virtual inline void ComputeStrain(){
		for (int i = 0; i < face.size(); i++){
			Ds[i] = Matrix3x2<T>(
				x[face[i][1]] - x[face[i][0]],
				x[face[i][2]] - x[face[i][0]]
				);
			F[i] = Ds[i] * Bm[i];
			E[i] = (F[i].Transpose() * F[i] - I) * 0.5;
		}
	}

	/**
		Set all force to zero
	*/
	virtual inline void ClearForce(){
		for (int i = 0; i < f.size(); i++)
			f[i] = Vec3<T>(0, 0, 0);
	}

	/**
		Compute elastic force and add to f

		\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddElasticForce(){
		for (int i = 0; i < face.size(); i++){
			Matrix3x2<T> P = Stress(i);
			Matrix3x2<T> H = P * Bm[i].Transpose() * (-W[i]);

			Vec3<T> f1(H[0][0], H[1][0], H[2][0]);
			Vec3<T> f2(H[0][1], H[1][1], H[2][1]);
			f[face[i][1]] += f1;
			f[face[i][2]] += f2;
			f[face[i][0]] += -f1 - f2;
		}

		return 1;
	}

	/**
	Compute rayleigh damping force and add to f

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddDampingForce() {
		//	add elastic damping force
		if (damping_k > 1e-5) {
			//	the equation of calculating damping force is the same as elastic differential force
			std::vector<Vec3<T>> tmp_f;
			ComputeElasticDifferentialForce(tmp_f, v);

			//	add elastic damping and mass damping
			for (int i = 0; i < f.size(); i++) {
				f[i] += damping_k * tmp_f[i];
			}
		}

		//	add mass damping force
		if (damping_mass > 1e-5) {
			for (int i = 0; i < f.size(); i++) {
				f[i] -= damping_mass * M[i] * v[i];
			}
		}

		return 1;
	}

	/**
		Compute collision force and add to f

		default no collision, but the function is here for inherent

		\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddCollisionForce(){
		return 0;
	}

	/**
		Add external force (exclude gravity) to f

		default no external force, but the function is here for inherent

		\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddExternalForce() {
		return 0;
	}

	/**
		Add gravity to f

		\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddGravity(){
		for (int i = 0; i < f.size(); i++){
			f[i] += G[i];
		}

		return 1;
	}

	/**
		add bending force to f

		equation comes from the paper: Simulation of Clothing with Folds and Wrinkles
	*/
	virtual inline int AddBendingForce(){
		for (int i = 0; i < edge.size(); i++){
			if (ef[i].x < 0 || ef[i].y < 0)
				continue;
			int vid3 = edge[i].x;
			int vid4 = edge[i].y;
			int vid1 = geometry::getThirdVertexOfFace(face[ef[i].x], vid3, vid4);
			int vid2 = geometry::getThirdVertexOfFace(face[ef[i].y], vid3, vid4);
			if (attach_flag[vid1] && attach_flag[vid2] && attach_flag[vid3] && attach_flag[vid4])
				continue;
			yz::Vec3<T> f1, f2, f3, f4;
			calculateBendingForces(
				f1, f2, f3, f4, 
				x[vid1], x[vid2], x[vid3], x[vid4],
				bending_k);
			f[vid1] += f1;
			f[vid2] += f2;
			f[vid3] += f3;
			f[vid4] += f4;
		}

		return 1;
	}

	/**
		Compute differential of elastic force

		\param	df		force difference with respect to dx
		\param	dx		input displacement of each vertex
		\return			whether df is calculated, 1: calculated. 0: not calculated
	*/
	virtual inline int ComputeElasticDifferentialForce(
		std::vector<Vec3<T>>&		df, 
		const std::vector<Vec3<T>>& dx
	){
		assert(x.size() == dx.size());
		df.clear();
		df.resize(dx.size());

		for (int i = 0; i < face.size(); i++){
			Matrix3x2<T> dP = StressDifferential(i, dx[face[i][0]], dx[face[i][1]], dx[face[i][2]]);
			Matrix3x2<T> dH = dP * Bm[i].Transpose() * (-W[i]);

			Vec3<T> df1(dH[0][0], dH[1][0], dH[2][0]);
			Vec3<T> df2(dH[0][1], dH[1][1], dH[2][1]);
			df[face[i][1]] += df1;
			df[face[i][2]] += df2;
			df[face[i][0]] += -df1 - df2;
		}

		return 1;
	}

	/**
	Compute differential of damping force

	\param	df		force difference with respect to dx
	\param	dx		input displacement of each vertex
	\param	dt		time step of simulation
	\return			whether df is calculated, 1: calculated. 0 : not calculated
	*/
	virtual inline int ComputeDampingDifferentialForce(
		std::vector<Vec3<T>>&		df,
		const std::vector<Vec3<T>>& dx,
		T							dt
	) {
		assert(x.size() == dx.size());
		df.clear();
		df.resize(dx.size());

		T inv_dt = 1.0 / dt;

		//	damping on string
		for (int i = 0; i < face.size(); i++) {
			Matrix3x2<T> dP = StressDifferential(i, dx[face[i][0]], dx[face[i][1]], dx[face[i][2]]);
			Matrix3x2<T> dH = dP * Bm[i].Transpose() * (-W[i]);

			Vec3<T> df1(dH[0][0], dH[1][0], dH[2][0]);
			Vec3<T> df2(dH[0][1], dH[1][1], dH[2][1]);
			df[face[i][1]] += df1 * damping_k * inv_dt;
			df[face[i][2]] += df2 * damping_k * inv_dt;
			df[face[i][0]] += -(df1 + df2) * damping_k * inv_dt;
		}

		//	damping on mass
		for (int i = 0; i < f.size(); i++) {
			df[i] -= damping_mass * inv_dt * M[i] * dx[i];
		}

		return 1;
	}

	/**
	Compute differential of elastic and damping force

	\param	df		force difference with respect to dx
	\param	dx		input displacement of each vertex
	\param	dt		time step of simulation
	\return			whether df is calculated, 1: calculated. 0 : not calculated
	*/
	virtual inline int ComputeElasticDampingDifferentialForce(
		std::vector<Vec3<T>>&		df,
		const std::vector<Vec3<T>>& dx,
		T							dt
	) {
		assert(x.size() == dx.size());
		df.clear();
		df.resize(dx.size());

		T inv_dt = 1.0 / dt;

		//	elastic and damping on string
		for (int i = 0; i < face.size(); i++) {
			Matrix3x2<T> dP = StressDifferential(i, dx[face[i][0]], dx[face[i][1]], dx[face[i][2]]);
			Matrix3x2<T> dH = dP * Bm[i].Transpose() * (-W[i]);

			Vec3<T> df1(dH[0][0], dH[1][0], dH[2][0]);
			Vec3<T> df2(dH[0][1], dH[1][1], dH[2][1]);
			df[face[i][1]] += df1 * (1 + damping_k * inv_dt);
			df[face[i][2]] += df2 * (1 + damping_k * inv_dt);
			df[face[i][0]] += -(df1 + df2) * (1 + damping_k * inv_dt);
		}

		//	damping on mass
		for (int i = 0; i < f.size(); i++) {
			df[i] -= damping_mass * inv_dt * M[i] * dx[i];
		}

		return 1;
	}

	/**
		Compute differential of collision force

		\param	df		force difference with respect to dx
		\param	dx		input displacement of each vertex
		\return			whether df is calculated, 1: calculated. 0: not calculated
	*/
	virtual inline int ComputeCollisionDifferentialForce(
		std::vector<Vec3<T>>&		df, 
		const std::vector<Vec3<T>>& dx
	){
		df.clear();
		df.resize(dx.size(), yz::Vec3<T>(0, 0, 0));
		return 0;
	}

	/**
		Compute differential of bending force

		\param	df		force difference with respect to dx
		\param	dx		input displacement of each vertex
		\return			whether df is calculated, 1: calculated. 0: not calculated
	*/
	virtual inline int ComputeBendingDifferentialForce(
		std::vector<Vec3<T>>&		df, 
		const std::vector<Vec3<T>>& dx
	){
		return 0;
	}

	/**
		Calculate the stress of an element. If you want to use other stress, just overload this function

		St. Venant-Kirchhoff model
	*/
	virtual inline Matrix3x2<T> Stress(const int& face_id){
		return F[face_id] * (E[face_id] * 2 * mu + I * lambda*E[face_id].Trace());
	}

	/**
		Calculate the stress differential of an element. If you want to use other stress, just overload this function

		St. Venant-Kirchhoff model

		\param	face_id		index of the face
		\param	dx0			displacement of vertex 0
		\param	dx1			displacement of vertex 1
		\param	dx2			displacement of vertex 2
	*/
	virtual inline Matrix3x2<T> StressDifferential(
		const int&		face_id,
		const Vec3<T>&	dx0,
		const Vec3<T>&	dx1,
		const Vec3<T>&	dx2)
	{
		Matrix3x2<T> dDs(dx1 - dx0, dx2 - dx0);
		Matrix3x2<T> dF = dDs * Bm[face_id];
		Matrix2x2<T> dE = (dF.Transpose() * F[face_id] + F[face_id].Transpose() * dF) * T(0.5);
		return dF * (E[face_id] * 2 * mu + I * lambda*E[face_id].Trace()) + F[face_id] * (dE * 2 * mu + I * lambda*dE.Trace());
	}

	/**
		Calculate Von Mises of F
	*/
	virtual inline int CalculateVonMises(){
		if (VonMises_face.empty() || VonMises_vertex.empty())
			return 0;
		for (int i = 0; i < face.size(); i++){
			Matrix3x2<T> P = Stress(i);
			T s1, s2;
			Vec2<T> v1, v2;
			SVD3x2(P, v1, s1, v2, s2);
			VonMises_face[i] = sqrt(s1*s1 + s2*s2 - s1*s2);
		}
		for (int i = 0; i < vertex.size(); i++){
			VonMises_vertex[i] = 0;
			for (int j = vf_start[i]; j < vf_start[i + 1]; j++){
				int nf = vf[j];
				VonMises_vertex[i] += VonMises_face[nf];
			}
			VonMises_vertex[i] /= vf_start[i + 1] - vf_start[i];
		}
		return 1;
	}

	/**
	calculate total energy of the FEM system

	to reduce numerical error in adding float point numbers, we use Kahan summation algorithm
	*/
	virtual inline T TotalEnergy() {
		ComputeStrain();

		//	spring potential energy
		T energy_sum = 0;
		T c = 0;
		for (int i = 0; i < face.size(); i++) {
			T Ei_tr = E[i].Trace();
			T energy_density = (E[i] * E[i]).Trace() * mu + Ei_tr * Ei_tr * lambda * 0.5;
			T energy = energy_density * W[i];

			T y = energy - c;
			T t = energy_sum + y;
			c = (t - energy_sum) - y;
			energy_sum = t;
		}

		//	kinetic energy and weight potential energy
		for (int i = 0; i < M.size(); i++) {
			T kinetic_energy = 0.5 * M[i] * dot(v[i], v[i]) - M[i] * dot(x[i], gravity);

			T y = kinetic_energy - c;
			T t = energy_sum + y;
			c = (t - energy_sum) - y;
			energy_sum = t;
		}

		return energy_sum;
	}

public:
	//	physical parameters
	T							thickness;				///<	thickness of the sheet			
	T							mass;					///<	mass of the whole mesh
	T							youngs_modulus;			///<	Young's Module
	T							poisson_ratio;			///<	Poisson ratio
	T							damping_k;				///<	coefficient of damping
	T							damping_mass;			///<	coefficient of damping with respect to mass
	T							bending_k;				///<	bending coefficient

protected:
	T							mu;						///<	Lame coefficient
	T							lambda;					///<	Lame coefficient
	Matrix2x2<T>				I;						///<	Identity matrix

	std::vector<Matrix2x2<T>>	Dm;						///<	vector matrix of reference mesh of each face
	std::vector<Matrix2x2<T>>	Bm;						///<	inverse of Dm of each face
	std::vector<T>				W;						///<	area of each face
	std::vector<Matrix3x2<T>>	Ds;						///<	vector matrix of current mesh of each face
	std::vector<Matrix3x2<T>>	F;						///<	deformation gradient of each face
	std::vector<Matrix2x2<T>>	E;						///<	Green strain of each face

	std::vector<Vec3<T>>		G;						///<	gravity of each vertex
	std::vector<Vec3<T>>		prev_v;					///<	velocity of the previous time step of each vertex

	std::vector<T>				VonMises_face;			///<	Von Mises stress of each face
	std::vector<T>				VonMises_vertex;		///<	Von Mises stress of each vertex
};

}}}	//	yz::physics::fem


#endif //	__YZ_FEM_TRI_MESH_H__