/***********************************************************/
/**	\file
	\brief		Fitting
	\author		Yizhong Zhang
	\date		11/7/2019
*/
/***********************************************************/
#ifndef __YZ_FITTING_H__
#define __YZ_FITTING_H__

#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Fit plane
*/
//	========================================

/**
fit plane using point cloud, local coordinate is represented using right hand coordinate

\param	plane_org		local coordinate origin of the plane
\param	plane_x			x axis of the plane
\param	plane_y			y axis of the plane
\param	plane_z			z axis of the plane
\param	eigen_xyz		eigen value on each axis
\param	point			input the point cloud
\return					1 on succeed, 0 on fail
*/
template<typename T>
int fitPlaneFromPoints(
	Vec3<T>&					plane_org,
	Vec3<T>&					plane_x,
	Vec3<T>&					plane_y,
	Vec3<T>&					plane_z,
	Vec3<T>&					eigen_xyz,
	const std::vector<Vec3<T>>& point
) {
	if (point.size() < 3)	//	at least 3 points required
		return 0;

	//	calculate mean of the points
	Vec3<T> mean(0, 0, 0);
	for (int i = 0; i < point.size(); i++)
		mean += point[i];
	mean /= point.size();

	//	setup the matrix
	Matrix3x3<T> M(0, 0, 0, 0, 0, 0, 0, 0, 0);
	for (int i = 0; i < point.size(); i++) {
		Vec3<T> r = point[i] - mean;
		M[0][0] += r[0] * r[0];
		M[0][1] += r[0] * r[1];
		M[0][2] += r[0] * r[2];
		M[1][1] += r[1] * r[1];
		M[1][2] += r[1] * r[2];
		M[2][2] += r[2] * r[2];
	}
	M[1][0] = M[0][1];
	M[2][0] = M[0][2];
	M[2][1] = M[1][2];

	//	SVD
	Matrix3x3<T> U;
	Vec3<T> S;
	int ret = M.EigenSymm(U, S);

	//	write the plane
	plane_org = mean;
	plane_x = Vec3<T>(U[0][0], U[1][0], U[2][0]);
	plane_y = Vec3<T>(U[0][1], U[1][1], U[2][1]);
	plane_z = Vec3<T>(U[0][2], U[1][2], U[2][2]);
	eigen_xyz = Vec3<T>(sqrt(S[0]), sqrt(S[1]), sqrt(S[2]));

	//	check right hand coordinate
	if (dot(cross(plane_x, plane_y), plane_z) < 0)
		plane_z = -plane_z;

	return 1;
}

/**
fit plane using point cloud

\param	plane_cen		center of the plane
\param	plane_nor		normal of the plane
\param	point			input the point cloud
\return					1 on succeed, 0 on fail
*/
template<typename T>
int fitPlaneFromPoints(
	Vec3<T>&					plane_cen,
	Vec3<T>&					plane_nor,
	const std::vector<Vec3<T>>& point
) {
	Vec3<T>	plane_org, plane_x, plane_y, plane_z, eigen_xyz;

	int succ = fitPlaneFromPoints(
		plane_org, plane_x, plane_y, plane_z, eigen_xyz, point);

	if (succ) {
		plane_cen = plane_org;
		plane_nor = plane_z;
	}

	return succ;
}

/**
RANSAC fit plane using point cloud, local coordinate is represented using right hand coordinate

\param	plane_org		local coordinate origin of the plane
\param	plane_x			x axis of the plane
\param	plane_y			y axis of the plane
\param	plane_z			z axis of the plane
\param	eigen_xyz		eigen value on each axis
\param	inlier_mask		whether the point is inlier when fitting plane
\param	point			input the point cloud
\param	dist_thre		distance threshold to the plane
\param	max_RANSAC_iterations	max RANSAC iterations
\return					number of inliers
*/
template<typename T>
int fitPlaneFromPointsRANSAC(
	Vec3<T>&					plane_org,
	Vec3<T>&					plane_x,
	Vec3<T>&					plane_y,
	Vec3<T>&					plane_z,
	Vec3<T>&					eigen_xyz,
	std::vector<int>&			inlier_mask,
	const std::vector<Vec3<T>>& point,
	T							dist_thre,
	int							max_RANSAC_iterations = 20
) {
	if (point.size() <= 3) {
		int succ = fitPlaneFromPoints(
			plane_org, plane_x, plane_y, plane_z, eigen_xyz, point);
		if (succ) {
			inlier_mask.clear();
			inlier_mask.resize(point.size(), 1);
		}
		return succ;
	}

	//	each iteration, randomly select 3 points, then fit a plane
	std::srand(point.size());
	int curr_best_fit_num = 0;
	std::vector<Vec3<T>> test_point;
	std::vector<int> test_mask;

	for (int ransac_idx = 0; ransac_idx < max_RANSAC_iterations; ransac_idx++) {
		test_point.clear();
		test_mask.clear();
		test_mask.resize(point.size(), 0);

		//	select seed points
		while (test_point.size() < 3) {
			int index = std::rand() % point.size();
			while (test_mask[index] == 1)
				index = std::rand() % point.size();
			test_mask[index] = 1;
			test_point.push_back(point[index]);
		}

		//	initial fit
		Vec3<T> test_plane_org, test_plane_x, test_plane_y, test_plane_z, test_eigen_xyz;
		fitPlaneFromPoints(
			test_plane_org,
			test_plane_x,
			test_plane_y,
			test_plane_z,
			test_eigen_xyz,
			test_point);

		//	loop until cannot find more inlier
		for (int inc_loop = 0; inc_loop < 20; inc_loop++) {
			int more_inlier_flag = 0;
			for (int i = 0; i < point.size(); i++) {
				if (!test_mask[i] && fabs(dot(point[i] - test_plane_org, test_plane_z)) < dist_thre) {
					test_mask[i] = 1;
					test_point.push_back(point[i]);
					more_inlier_flag = 1;
				}
			}

			if (!more_inlier_flag)
				break;

			fitPlaneFromPoints(
				test_plane_org,
				test_plane_x,
				test_plane_y,
				test_plane_z,
				test_eigen_xyz,
				test_point);
		}

		//	exclude outlier
		test_point.clear();
		for (int i = 0; i < point.size(); i++) {
			if (fabs(dot(point[i] - test_plane_org, test_plane_z)) < dist_thre) {
				test_point.push_back(point[i]);
				test_mask[i] = 1;
			}
			else
				test_mask[i] = 0;
		}

		//	check whether improve
		if (test_point.size() > curr_best_fit_num) {
			fitPlaneFromPoints(
				test_plane_org,
				test_plane_x,
				test_plane_y,
				test_plane_z,
				test_eigen_xyz,
				test_point);

			curr_best_fit_num = test_point.size();
			plane_org = test_plane_org;
			plane_x = test_plane_x;
			plane_y = test_plane_y;
			plane_z = test_plane_z;
			eigen_xyz = test_eigen_xyz;
			inlier_mask.swap(test_mask);
		}

		if (curr_best_fit_num == point.size())
			break;
	}

	return curr_best_fit_num;
}

///@}


}}	//	namespace yz::geometry

#endif	//	__YZ_SPHERE_H__