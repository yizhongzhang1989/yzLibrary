/***********************************************************/
/**	\file
	\brief		mass spring system
	\details	
	\author		Yizhong Zhang
	\date		3/29/2019
*/
/***********************************************************/
#ifndef __YZ_MASS_SPRING_H__
#define __YZ_MASS_SPRING_H__

#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_physics/yz_physics_base.h"

namespace yz{	namespace physics{


/**
	The mass spring system
*/
template<class T>
class MassSpring :
	public ImplicitPhysicsBase<T>,
	public PointMass<T>,
	public PointAttachConstraint
{
public:
	using ImplicitPhysicsBase<T>::gravity;
	using PointMass<T>::M;
	using PointMass<T>::x;
	using PointMass<T>::v;
	using PointMass<T>::f;

public:
	void Draw() {
#ifdef YZ_gl_h
		glColor3f(1, 1, 1);
		for (int i = 0; i < x.size(); i++)
			yz::opengl::drawPointAsCube(x[i], 0.1);

		glColor3f(1, 0, 0);
		for (int i = 0; i < spring_node.size(); i++)
			yz::opengl::drawCylinder(x[spring_node[i].x], x[spring_node[i].y], 0.02, 8);

		//	draw velocity
		glColor3f(1, 1, 0);
		for (int i = 0; i < v.size(); i++) {
			yz::opengl::drawArrow(x[i], x[i] + v[i], 0.05);
		}
#else
		std::cout << "gl.h has to be included in order to use Draw() in MassSpring" << std::endl;
		return;
#endif
	}

	void SetupMassSpring(
		const std::vector<Vec3<T>>& node,
		const std::vector<int2>&	spring)
	{
		PointMass<T>::Init(node.size());
		x = node;

		spring_node = spring;

		InitPhysicaParameters();

		PointAttachConstraint::Init(node.size());
	}

	void InitPhysicaParameters(T node_mass = 1, T spring_k = 1, T spring_damping_k = 0, T damping_M = 0) {
		//	node
		std::fill(M.begin(), M.end(), node_mass);

		spring_K.clear();
		spring_K.resize(spring_node.size());

		G.resize(M.size());
		for (int i = 0; i < M.size(); i++)
			G[i] = M[i] * gravity;

		prev_v.clear();
		prev_v.resize(v.size());

		//	spring
		rest_length.resize(spring_node.size());
		for (int i = 0; i < rest_length.size(); i++) {
			rest_length[i] = (x[spring_node[i].x] - x[spring_node[i].y]).Length();
		}

		k.clear();
		k.resize(spring_node.size(), spring_k);

		damping_k.clear();
		damping_k.resize(spring_node.size(), spring_damping_k);

		damping_mass = damping_M;
	}

	/**
	set gravity of the system
	*/
	virtual inline void SetGravity(Vec3<T> g = Vec3<T>(0, -9.8, 0)) {
		gravity = g;

		G.resize(M.size());
		for (int i = 0; i < M.size(); i++)
			G[i] = M[i] * gravity;
	}

	/**
	simulate one step using forward euler
	*/
	int SimulateStepExplicit(T dt) {
		if (prev_v.empty())
			prev_v = v;
		std::fill(f.begin(), f.end(), Vec3<T>(0, 0, 0));

		//	calculate f
		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T curr_len = r.Length();
			r /= curr_len;	//	set normalize

			//	elastic force
			Vec3<T> elas_f = k[i] * (curr_len - rest_length[i]) * r;
			f[spring_node[i].y] -= elas_f;
			f[spring_node[i].x] += elas_f;

			//	rayleigh damping force (spring)
			Vec3<T> r_v = prev_v[spring_node[i].y] - prev_v[spring_node[i].x];
			T rel_v = dot(r, r_v);
			Vec3<T> damp_f = damping_k[i] * k[i] * rel_v * r;
			f[spring_node[i].y] -= damp_f;
			f[spring_node[i].x] += damp_f;
		}

		//	rayleigh damping force (mass) and gravity
		for (int i = 0; i < f.size(); i++) {
			f[i] += G[i];
			f[i] -= damping_mass * M[i] * prev_v[i];
		}

		//	calculate v and x
		prev_v = v;
		for (int i = 0; i < x.size(); i++) {
			Vec3<T> acc = f[i] / M[i];
			if (attach_flag[i])
				acc = Vec3<T>(0, 0, 0);
			v[i] += acc * dt;
			x[i] += v[i] * dt;
		}

		return 1;
	}

	/**
	Simulate one time step

	\param	dt					positive time step, recommended to be no larger than 0.01.
								negative time step, quasi-static relaxiation
	\param	max_newton_steps	max Newton steps performed in simulation, default 20 steps
	\return						newton iterations actually performed
	*/
	int SimulateStepImplicit(T dt, int max_newton_steps) {
		//	initialize v and x as constant velocity motion
		//	CAUTION! the initialization of x and v in Eftychios D. Sifakis' course P34 is wrong
		prev_v = v;
		for (int i = 0; i < v.size(); i++)
			x[i] += v[i] * dt;

		int newton_count = 0;
		for (int i = 0; i < max_newton_steps; i++) {
			newton_count++;
			T step_size = NewtonStep(dt);
			if (step_size < 1e-8)
				break;
		}

		return newton_count;
	}

public:
	/**
	A Newton step + line search

	\param	dt	time step, if dt <= 0, we perform quasi static simulation
	\return		stopping criteria, dx * b
	*/
	virtual inline T NewtonStep(T dt) {
		int quasi_static_flag = (dt <= 0 ? 1 : 0);
		T inv_dt = (quasi_static_flag ? 0.0 : 1.0 / dt);
		std::vector<Vec3<T>> b;
		b.resize(x.size());

		CalculateSpringK();

		//	compute right hand side
		ClearForce();
		if (quasi_static_flag) {
			AddElasticForce();
			AddGravity();
			AddCollisionForce();
			AddExternalForce();
		}
		else {
			AddElasticDampingForce();
			AddGravity();
		}
		for (int i = 0; i < f.size(); i++) {		//	calculate right hand side
			if (attach_flag[i])		//	this vertex is constraint
				b[i] = Vec3<T>(0, 0, 0);
			else
				b[i] = f[i] +(quasi_static_flag ? 0.0 : inv_dt * M[i] * (prev_v[i] - v[i]));
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
		for (int k = 0; k < 20; k++) {
			//	Ap = dx * M / dtdt - df(dx)
			std::vector<Vec3<T>> Ap;
			std::vector<Vec3<T>> df;
			if (quasi_static_flag) { 	//	in quasi-static, only elastic force exist
				ComputeElasticDifferentialForce(df, p);
				Ap.swap(df);
			}
			else { 					//	in dynamic simulation, damping force exist
				ComputeElasticDampingDifferentialForce(df, p, dt);
				Ap.swap(df);
				ComputeCollisionDifferentialForce(df, p);
				for (int i = 0; i < Ap.size(); i++)
					Ap[i] += df[i];
			}
			T inv_dt_squ = inv_dt * inv_dt;
			for (int i = 0; i < Ap.size(); i++) {
				Ap[i] = inv_dt_squ * M[i] * p[i] - Ap[i];
			}

			for (int i = 0; i < attach_flag.size(); i++) {	//	Add constraint
				if (attach_flag[i])
					Ap[i] = Vec3<T>(0, 0, 0);
			}

			T pTAp = XTY(p, Ap);

			T alpha = rTr / pTAp;

			if (isnan(alpha) || isinf(alpha)) {
				//	degenerate case appear, we add noise to dx
				for (int i = 0; i < dx.size(); i++) {
					if (!attach_flag[i])
						dx[i] += Vec3<T>(randNumber(0, 1e-5), randNumber(0, 1e-5), randNumber(0, 1e-5));
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
	Set all force to zero
	*/
	virtual inline void ClearForce() {
		std::fill(f.begin(), f.end(), Vec3<T>(0, 0, 0));
	}

	/**
	Add gravity to f

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddGravity() {
		for (int i = 0; i < f.size(); i++) {
			f[i] += G[i];
		}

		return 1;
	}

	/**
	Compute elastic force and add to f

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddElasticForce() {
		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T curr_len = r.Length();
			r /= curr_len;	//	set normalize

			Vec3<T> elas_f = k[i] * r * (curr_len - rest_length[i]);
			f[spring_node[i].x] += elas_f;
			f[spring_node[i].y] -= elas_f;
		}

		return 1;
	}

	/**
	Compute rayleigh damping force and add to f

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddDampingForce() {
		Matrix3x3<T> I;
		I.SetIdentity();

		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T len = r.Length();
			r /= len;
			Matrix3x3<T> rTr(
				r[0] * r[0], r[0] * r[1], r[0] * r[2],
				r[1] * r[0], r[1] * r[1], r[1] * r[2],
				r[2] * r[0], r[2] * r[1], r[2] * r[2]);

			Vec3<T> r_v = v[spring_node[i].y] - v[spring_node[i].x];
			Vec3<T> damp_f = -damping_k[i] * k[i] * ((1 - rest_length[i] / len)*(I - rTr) + rTr) * r_v;
			f[spring_node[i].x] -= damp_f;
			f[spring_node[i].y] += damp_f;
		}

		//	rayleigh damping force (mass)
		for (int i = 0; i < f.size(); i++) {
			f[i] -= damping_mass * M[i] * v[i];
		}

		return 1;
	}

	/**
	Compute forces and add to f, combine elastic and damping to reduce computational cost

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddElasticDampingForce() {		
		Matrix3x3<T> I;
		I.SetIdentity();

		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T curr_len = r.Length();
			r /= curr_len;	//	set normalize

			//	elastic force
			Vec3<T> elas_f = k[i] * (curr_len - rest_length[i]) * r;
			f[spring_node[i].x] += elas_f;
			f[spring_node[i].y] -= elas_f;

			//	rayleigh damping force (spring)
			Vec3<T> r_v = v[spring_node[i].y] - v[spring_node[i].x];
			Vec3<T> damp_f = -damping_k[i] * spring_K[i] * r_v;
			f[spring_node[i].x] -= damp_f;
			f[spring_node[i].y] += damp_f;
		}

		//	rayleigh damping force on mass
		for (int i = 0; i < f.size(); i++) {
			f[i] -= damping_mass * M[i] * v[i];
		}

		return 1;
	}

	/**
	Compute collision force and add to f

	default no collision, but the function is here for inherent

	\return		whether force added to f, 1: added. 0: not added
	*/
	virtual inline int AddCollisionForce() {
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
	Compute differential of elastic force

	equations of force Jacobian can be found here:
	http://blog.mmacklin.com/2012/05/04/implicitsprings/

	\param	df		force difference with respect to dx
	\param	dx		input displacement of each vertex
	\return			whether df is calculated, 1: calculated. 0: not calculated
	*/
	virtual inline int ComputeElasticDifferentialForce(
		std::vector<Vec3<T>>&		df, 
		const std::vector<Vec3<T>>& dx
	) {
		df.clear();
		df.resize(dx.size());

		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = dx[spring_node[i].y] - dx[spring_node[i].x];
			df[spring_node[i].x] += spring_K[i] * r;
			df[spring_node[i].y] -= spring_K[i] * r;
		}

		return 1;
	}

	/**
	Compute differential of damping force
	*/
	virtual inline int ComputeDampingDifferentialForce(
		std::vector<Vec3<T>>&		df, 
		const std::vector<Vec3<T>>& dx, 
		T							dt
	) {
		df.clear();
		df.resize(dx.size());

		T inv_dt = 1.0 / dt;

		//	damping on string
		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = dx[spring_node[i].y] - dx[spring_node[i].x];
			df[spring_node[i].x] += spring_K[i] * r * damping_k[i] * inv_dt;
			df[spring_node[i].y] -= spring_K[i] * r * damping_k[i] * inv_dt;
		}

		//	damping on mass
		for (int i = 0; i < f.size(); i++) {
			df[i] -= damping_mass * inv_dt * M[i] * dx[i];
		}

		return 1;
	}

	/**
	Compute differential of all forces
	*/
	virtual inline int ComputeElasticDampingDifferentialForce(
		std::vector<Vec3<T>>&		df,
		const std::vector<Vec3<T>>& dx,
		T							dt
	) {
		df.clear();
		df.resize(dx.size());

		T inv_dt = 1.0 / dt;

		//	elastic and damping on string
		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = dx[spring_node[i].y] - dx[spring_node[i].x];
			df[spring_node[i].x] += spring_K[i] * r * (damping_k[i] * inv_dt + 1);
			df[spring_node[i].y] -= spring_K[i] * r * (damping_k[i] * inv_dt + 1);
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
	virtual inline int ComputeCollisionDifferentialForce(std::vector<Vec3<T>>& df, const std::vector<Vec3<T>>& dx) {
		df.clear();
		df.resize(dx.size(), yz::Vec3<T>(0, 0, 0));
		return 0;
	}

	/**
		Calculate K (force jacobian) of each spring

		To improve rubustness, we enforce each 3x3 matrix to be semi-positive definite
	*/
	void CalculateSpringK() {
		Matrix3x3<T> I;
		I.SetIdentity();

		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T curr_len = r.Length();
			r /= curr_len;

			Matrix3x3<T> rTr(
				r[0] * r[0], r[0] * r[1], r[0] * r[2],
				r[1] * r[0], r[1] * r[1], r[1] * r[2],
				r[2] * r[0], r[2] * r[1], r[2] * r[2]);
			Matrix3x3<T> D = k[i] * ((1 - rest_length[i] / curr_len)*(I - rTr) + rTr);
			spring_K[i] = D;
			EnforceSemiPositiveDefinite(spring_K[i]);
		}
	}

	/**
		calculate total energy of the mass spring system

		to reduce numerical error in adding float point numbers, we use Kahan summation algorithm
	*/
	inline T TotalEnergy() {
		//	spring potential energy
		T energy_sum = 0;
		T c = 0;
		for (int i = 0; i < spring_node.size(); i++) {
			Vec3<T> r = x[spring_node[i].y] - x[spring_node[i].x];
			T len_diff = r.Length() - rest_length[i];
			T spring_energy = 0.5 * k[i] * len_diff * len_diff;

			T y = spring_energy - c;
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
	//	spring parameter
	std::vector<int2>			spring_node;			///<	nodes connected by spring
	std::vector<T>				rest_length;			///<	rest length of each spring
	std::vector<T>				k;						///<	k of each spring
	std::vector<T>				damping_k;				///<	Rayleigh damping coef of each spring, associate with K
	T							damping_mass;			///<	Rayleigh damping coef of each mass, associate with M

protected:
	std::vector<Matrix3x3<T>>	spring_K;				///<	K matrix of each spring (jacobian matrix of force)
	std::vector<Vec3<T>>		G;						///<	gravity of each node
	std::vector<Vec3<T>>		prev_v;					///<	velocity of the previous time step of each vertex
};


}}	//	end namespace yz::physics

#endif	//	__YZ_MASS_SPRING_H__