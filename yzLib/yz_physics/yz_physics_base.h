/***********************************************************/
/**	\file
	\brief		Base of physical simulation
	\details	including virtual base for simulation, and basic physical parameters
	\author		Yizhong Zhang
	\date		3/29/2019
*/
/***********************************************************/
#ifndef __YZ_PHYSICS_BASE_H__
#define __YZ_PHYSICS_BASE_H__

#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace physics{


/**
	The virtual base of physical simulation
*/
template<class T>
class PhysicsBase {
public:
	/**
	simulate a step
	*/
	virtual int SimulateStep(T dt) = 0;

	/**
	set gravity, default value is gravity on earth
	*/
	virtual inline void SetGravity(Vec3<T> g = Vec3<T>(0, -9.8, 0)) {
		gravity = g;
	}

	/**
	set the environment to be weighless
	*/
	virtual inline void SetWeightlessness() {
		gravity = Vec3<T>(0, 0, 0);
	}

public:
	Vec3<T>	gravity;		///<	gravity of the system
};


/**
	The virtual base of implicit physical simulation

	helper functions are provided in this class
*/
template<class T>
class ImplicitPhysicsBase : public PhysicsBase<T> {
public:
	/**
	simulate a step using implicit method, with default 20 iterations
	*/
	virtual int SimulateStep(T dt) {
		return SimulateStepImplicit(dt, 20);
	}

	/**
	simulate a step, with iterations
	*/
	virtual int SimulateStepImplicit(T dt, int max_newton_steps) = 0;

public:
	/**
	Compute vector dot product

	\param	X		left vector
	\param	Y		right vector
	\return			the dot product of X*Y
	*/
	virtual inline T XTY(const std::vector<T>& X, const std::vector<T>& Y) {
		assert(X.size() == Y.size());

		T sum = 0;
		for (int i = 0; i < X.size(); i++) {
			sum += X[i] * Y[i];
		}

		return sum;
	}

	/**
	Compute vector dot product

	\param	X		left vector
	\param	Y		right vector
	\return			the dot product of X*Y
	*/
	virtual inline T XTY(const std::vector<Vec3<T>>& X, const std::vector<Vec3<T>>& Y) {
		assert(X.size() == Y.size());

		T sum = 0;
		for (int i = 0; i < X.size(); i++) {
			sum += X[i].x*Y[i].x + X[i].y*Y[i].y + X[i].z*Y[i].z;
		}

		return sum;
	}

	/**
	SVD decomposition on 3x2 matrix, get the planar part

	\param	mat			the input 3x2 matrix
	\param	eigen_val	2x2 matrix SIGMA, diagonal matrix of eigen values
	\param	eigen_vec	2x2 matrix V, corresponding matrix of eigen vectors
	*/
	virtual inline void SVD3x2(Matrix3x2<T> mat, Com2<T>& eigen_vec1, T& eigen_val1, Com2<T>& eigen_vec2, T& eigen_val2) {
		Matrix2x2<T> ATA = mat.Transpose() * mat;
		ATA.Eigen(eigen_vec1, eigen_val1, eigen_vec2, eigen_val2);
		eigen_val1 = sqrt(eigen_val1);
		eigen_val2 = sqrt(eigen_val2);
	}

	/**
	SVD decomposition on 3x2 matrix, get the planar part, in matrix format

	\param	mat			the input 3x2 matrix
	\param	SIGMA		3x2 matrix SIGMA, the last row is zero
	\param	VT			2x2 matrix VT, transpose of V
	*/
	virtual inline void SVD3x2(Matrix3x2<T> mat, Matrix3x2<T>& SIGMA, Matrix2x2<T>& VT) {
		Com2<T> eigen_vec1, eigen_vec2;
		T eigen_val1, eigen_val2;
		Matrix2x2<T> ATA = mat.Transpose() * mat;
		ATA.Eigen(eigen_vec1, eigen_val1, eigen_vec2, eigen_val2);
		SIGMA[0][0] = sqrt(eigen_val1);
		SIGMA[1][1] = sqrt(eigen_val2);
		SIGMA[0][1] = SIGMA[1][0] = SIGMA[2][0] = SIGMA[2][1] = 0;
		VT[0][0] = eigen_vec1[0];
		VT[0][1] = eigen_vec1[1];
		VT[1][0] = eigen_vec2[0];
		VT[1][1] = eigen_vec2[1];
	}

	/**
	SVD decomposition on 3x2 matrix, get the planar part, in matrix format

	\param	mat			the input 3x2 matrix
	\param	SIGMA		2x2 matrix SIGMA, remove the last empty row
	\param	VT			2x2 matrix VT, transpose of V
	*/
	virtual inline void SVD3x2(Matrix3x2<T> mat, Matrix2x2<T>& SIGMA, Matrix2x2<T>& VT) {
		Com2<T> eigen_vec1, eigen_vec2;
		T eigen_val1, eigen_val2;
		Matrix2x2<T> ATA = mat.Transpose() * mat;
		ATA.Eigen(eigen_vec1, eigen_val1, eigen_vec2, eigen_val2);
		SIGMA[0][0] = sqrt(eigen_val1);
		SIGMA[1][1] = sqrt(eigen_val2);
		SIGMA[0][1] = SIGMA[1][0] = 0;
		VT[0][0] = eigen_vec1[0];
		VT[0][1] = eigen_vec1[1];
		VT[1][0] = eigen_vec2[0];
		VT[1][1] = eigen_vec2[1];
	}

	/**
	Enforce a 3x3 symmetric matrix to be semi-positive definite

	The input matrix must be symmetric, or the calculation will have problem

	\param	mat		the matrix to enforce
	\return			eigen values of the input matrix
	*/
	virtual inline Vec3<T> EnforceSemiPositiveDefinite(Matrix3x3<T>& mat) const {
		//	calculate eigen values
		Vec3<T> S;
		int is_symm = mat.EigenSymm(S[0], S[1], S[2]);
		if (!is_symm)
			return Vec3<T>(-1, -1, -1);
		if (S[2] >= 0)
			return S;

		Matrix3x3d UM, SM;
		mat.EigenSymm(UM, S);
		SM.SetAsDiagonal(S);
		for (int i = 0; i < 3; i++)
			if (SM[i][i] < 0)
				SM[i][i] = 0;
		mat = UM * SM * UM.Transpose();		//	UM is rotational matrix, transpose = inverse
		return S;
	}

};


/**
	physical parameters of point mass system
*/
template<class T>
class PointMass {
public:
	std::vector<T>				M;		///<	mass of each node

	std::vector<Vec3<T>>		x;		///<	position of each node
	std::vector<Vec3<T>>		v;		///<	velocity of each node
	std::vector<Vec3<T>>		f;		///<	force on each node

public:
	virtual void Init(int node_size) {
		M.clear();
		x.clear();
		v.clear();
		f.clear();
		M.resize(node_size);
		x.resize(node_size);
		v.resize(node_size);
		f.resize(node_size);
	}
};


/**
	constraint status of points
*/
class PointAttachConstraint {
public:
	std::vector<int>		attach_flag;	///<	the attach status of each node

public:
	virtual void Init(int node_size) {
		attach_flag.clear();
		attach_flag.resize(node_size, 0);
	}

	virtual inline void AddAttachConstraint(unsigned int node_id, int flag = 1) {
		if (node_id >= attach_flag.size())
			return;

		attach_flag[node_id] = flag;
	}
};


}}	//	end namespace yz::physics

#endif	//	__YZ_PHYSICS_BASE_H__
