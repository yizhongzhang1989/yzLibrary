/***********************************************************/
/**	\file
	\brief		point point icp
	\details	iterative closest point
	\author		Yizhong Zhang
	\date		12/7/2018
*/
/***********************************************************/
#ifndef __YZ_POINT_POINT_H__
#define __YZ_POINT_POINT_H__

#include "yzLib_config.h"
#ifdef yzLib_ENABLE_Eigen
#	include <Eigen/Dense>
#endif

namespace yz {	namespace geometry {	namespace icp {

/**
calculate a rigid transform that transform source point cloud to match the target

Kabsch algorithm is used for least square

\param	rot_M		rotation part of the transform matrix
\param	trans_V		translation part of the transform matrix
\param	source		the source point cloud
\param	target		the target point cloud
\return				whether the calculating succeed
*/
template <typename T>
int calculateTransformPointToPoint(
	Matrix3x3<T>& rot_M,
	Vec3<T>& trans_V,
	const std::vector<Vec3<T>>& source,
	const std::vector<Vec3<T>>& target
) {
	rot_M.SetIdentity();
	trans_V = Vec3<T>(0, 0, 0);

	if (source.size() < 3)	//	at least 3 points are required
		return 0;
	if (source.size() != target.size())
		return 0;

	//	calculate center
	Vec3<T> source_center(0, 0, 0), target_center(0, 0, 0);
	for (int i = 0; i < source.size(); i++) {
		source_center += source[i];
		target_center += target[i];
	}
	source_center /= source.size();
	target_center /= source.size();

	Matrix3x3<T> Rot;
	for (int i = 0; i < source.size(); i++) {
		Vec3<T> vl = source[i] - source_center;
		Vec3<T> vr = target[i] - target_center;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				Rot[j][k] += vr[j] * vl[k];
			}
		}
	}

	//	perform SVD to remove scale and screw 
#ifdef yzLib_ENABLE_Eigen
	Eigen::Matrix<double, 3, 3, Eigen::RowMajor> A;
	for (int i = 0; i < 9; i++)
		A.data()[i] = Rot[0][i];
	Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	A = svd.matrixU() * svd.matrixV().transpose();
	Eigen::DiagonalMatrix<double, 3> diag(1, 1, A.determinant());
	A = svd.matrixU() * diag * svd.matrixV().transpose();
	for (int i = 0; i < 9; i++)
		Rot[0][i] = A.data()[i];
#else
	yz::Matrix3x3<T> RRT = Rot * Rot.Transpose();
	yz::Matrix3x3<T> U, SIGMA, VT;
	yz::Vec3<T> S;
	RRT.EigenSymm(U, S);
	SIGMA.SetAsDiagonal(sqrt(S[0]), sqrt(S[1]), sqrt(S[2]));
	VT = (U * SIGMA).Inverse() * Rot;
	Rot = U * VT;
#endif

	//	calculate transform
	rot_M = Rot;
	trans_V = target_center - Rot * source_center;

	return 1;
}

/**
calculate a rigid transform that transform source point cloud to match the target

Kabsch algorithm is used for least square

\param	rot_M		rotation part of the transform matrix
\param	trans_V		translation part of the transform matrix
\param	source		the source point cloud
\param	target		the target point cloud
\param	weight		weight of each point
\return				whether the calculating succeed
*/
template <typename T>
int calculateTransformPointToPoint(
	Matrix3x3<T>& rot_M,
	Vec3<T>& trans_V,
	const std::vector<Vec3<T>>& source,
	const std::vector<Vec3<T>>& target,
	const std::vector<T>& weight
) {
	rot_M.SetIdentity();
	trans_V = Vec3<T>(0, 0, 0);

	if (source.size() < 3)	//	at least 3 points are required
		return 0;
	if (source.size() != target.size() ||
		source.size() != weight.size())
		return 0;

	//	calculate center
	Vec3<T> source_center(0, 0, 0), target_center(0, 0, 0);
	T weight_sum = 0;
	for (int i = 0; i < source.size(); i++) {
		source_center += source[i] * weight[i];
		target_center += target[i] * weight[i];
		weight_sum += weight[i];
	}
	source_center /= weight_sum;
	target_center /= weight_sum;

	//	create R
	Matrix3x3<T> Rot;
	for (int i = 0; i < source.size(); i++) {
		Vec3<T> vl = source[i] - source_center;
		Vec3<T> vr = target[i] - target_center;
		T w = weight[i];
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				Rot[j][k] += vr[j] * vl[k] * w;
			}
		}
	}

	//	perform SVD to remove scale and screw 
#ifdef yzLib_ENABLE_Eigen
	Eigen::Matrix<double, 3, 3, Eigen::RowMajor> A;
	for (int i = 0; i < 9; i++)
		A.data()[i] = Rot[0][i];
	Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	A = svd.matrixU() * svd.matrixV().transpose();
	Eigen::DiagonalMatrix<double, 3> diag(1, 1, A.determinant());
	A = svd.matrixU() * diag * svd.matrixV().transpose();
	for (int i = 0; i < 9; i++)
		Rot[0][i] = A.data()[i];
#else
	yz::Matrix3x3<T> RRT = Rot * Rot.Transpose();
	yz::Matrix3x3<T> U, SIGMA, VT;
	yz::Vec3<T> S;
	RRT.EigenSymm(U, S);
	SIGMA.SetAsDiagonal(sqrt(S[0]), sqrt(S[1]), sqrt(S[2]));
	VT = (U * SIGMA).Inverse() * Rot;
	Rot = U * VT;
#endif

	//	calculate transform
	rot_M = Rot;
	trans_V = target_center - Rot * source_center;

	return 1;
}

/**
calculate a rigid transform that transform source point cloud to match the target

\param	trans		the transform matrix
\param	source		the source point cloud
\param	target		the target point cloud
\return				whether the calculating succeed
*/
template <typename T>
int calculateTransformPointToPoint(
	Matrix4x4<T>& trans,
	const std::vector<Vec3<T>>& source,
	const std::vector<Vec3<T>>& target
) {
	Matrix3x3<T>	rot_M;
	Vec3<T>			trans_V;

	int ret = calculateTransformPointToPoint(rot_M, trans_V, source, target);
	if (ret)
		trans = Matrix4x4<T>(rot_M, trans_V);

	return ret;
}

/**
calculate a rigid transform that transform source point cloud to match the target

\param	trans		the transform matrix
\param	source		the source point cloud
\param	target		the target point cloud
\param	weight		weight of each point
\return				whether the calculating succeed
*/
template <typename T>
int calculateTransformPointToPoint(
	Matrix4x4<T>& trans,
	const std::vector<Vec3<T>>& source,
	const std::vector<Vec3<T>>& target,
	const std::vector<T>& weight
) {
	Matrix3x3<T>	rot_M;
	Vec3<T>			trans_V;

	int ret = calculateTransformPointToPoint(rot_M, trans_V, source, target, weight);
	if (ret)
		trans = Matrix4x4<T>(rot_M, trans_V);

	return ret;
}

}	}	}

#endif	//	__YZ_POINT_POINT_H__