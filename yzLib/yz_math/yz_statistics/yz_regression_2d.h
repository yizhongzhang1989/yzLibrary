/***********************************************************/
/**	\file
	\brief		regression of 2d data
	\author		Yizhong Zhang
	\date		8/27/2013
*/
/***********************************************************/
#ifndef __YZ_REGRESSION_2D_H__
#define __YZ_REGRESSION_2D_H__

#include <vector>
#include "yzLib/yz_math/yz_matrix.h"

namespace yz{	namespace statistics{

//	========================================
///@{
/**	@name	Fitting Straight Line 2D
*/
//	========================================

/**
	line fitting from 2D points, the line is y = k x + b

	the error is calculated in y direction, so k cannot be too big

	\param	k			return the slope of the line
	\param	b			return the y-intercept of the line
	\param	xy_ptr		pointer to the array of points, data arranged in xy_xy_xy_ format
	\param	point_num	number of points in the array
	\return				whether regression succeed
*/
template<typename T>
inline int fitLineLeastSquare2D(T& k, T& b, const T* xy_ptr, int point_num){
	if( point_num < 2 )
		return 0;

	Matrix2x2<PROMOTE_T_TO_FLOAT>	mat(0, 0, 0, 0);
	Vec2<PROMOTE_T_TO_FLOAT>		rhs(0, 0);

	//	fill the matrix
	for(int i=0; i<point_num; i++){
		T x = xy_ptr[i*2];
		T y = xy_ptr[i*2+1];
		mat.data[0][0] += x * x;
		mat.data[0][1] += x;
		rhs.x += x * y;
		rhs.y += y;
	}
	mat.data[1][0] = mat.data[0][1];
	mat.data[1][1] = point_num;

	Vec2<PROMOTE_T_TO_FLOAT>	res = mat.Inverse() * rhs;
	k = res[0];
	b = res[1];
	return 1;
}

/**
	line fitting from 2D points, the line is y = k x + b

	the error is calculated in y direction, so k cannot be too big

	\param	k		return the slope of the line
	\param	b		return the y-intercept of the line
	\param	xy		list of the points
	\return			whether regression succeed
*/
template<typename T>
inline int fitLineLeastSquare(T& k, T& b, const std::vector<Vec2<T>>& xy){
	if( xy.empty() )
		return 0;

	return fitLineLeastSquare2D(k, b, (T*)&xy[0], xy.size());
}

/**
	line fitting from 2D points, the line is x = a y + b

	the error is calculated in x direction, so parameter a cannot be too big

	\param	a			return the inverse of slope of the line
	\param	b			return the x-intercept of the line
	\param	xy_ptr		pointer to the array of points, data arranged in xy_xy_xy_ format
	\param	point_num	number of points in the array
	\return				whether regression succeed
*/
template<typename T>
inline int fitLineLeastSquareX2D(T& a, T& b, const T* xy_ptr, int point_num){
	if( point_num < 2 )
		return 0;

	Matrix2x2<PROMOTE_T_TO_FLOAT>	mat(0, 0, 0, 0);
	Vec2<PROMOTE_T_TO_FLOAT>		rhs(0, 0);

	//	fill the matrix
	for(int i=0; i<point_num; i++){
		T x = xy_ptr[i*2];
		T y = xy_ptr[i*2+1];
		mat.data[0][0] += y * y;
		mat.data[0][1] += y;
		rhs.x += x * y;
		rhs.y += x;
	}
	mat.data[1][0] = mat.data[0][1];
	mat.data[1][1] = point_num;

	Vec2<PROMOTE_T_TO_FLOAT>	res = mat.Inverse() * rhs;
	a = res[0];
	b = res[1];
	return 1;
}

/**
	line fitting from 2D points, the line is x = a y + b

	the error is calculated in x direction, so parameter a cannot be too big

	\param	a		return the inverse of slope of the line
	\param	b		return the x-intercept of the line
	\param	xy		list of the points
	\return			whether regression succeed
*/
template<typename T>
inline int fitLineLeastSquareX(T& a, T& b, const std::vector<Vec2<T>>& xy){
	if( xy.empty() )
		return 0;

	return fitLineLeastSquareX2D(a, b, (T*)&xy[0], xy.size());
}

/**
	line fitting from 2D points with PCA

	for PCA, the error the perpendicular distance from each point to the line

	\param	mean_x		return center of the points, x coordinate
	\param	mean_y		return center of the points, y coordinate
	\param	pa_x		return principle axis, x
	\param	pa_y		return principle axis, y
	\param	xy_ptr		pointer to the array of points, data arranged in xy_xy_xy_ format
	\param	point_num	number of points in the array
	\return				whether regression succeed
*/
template<typename T>
inline int fitLinePCA2D(T& mean_x, T& mean_y, T& pa_x, T& pa_y, const T* xy_ptr, int point_num){
	if( point_num < 2 )
		return 0;

	//	calculate mean
	PROMOTE_T_TO_FLOAT cx = 0, cy = 0;
	for(int i=0; i<point_num; i++){
		cx += xy_ptr[i*2];
		cy += xy_ptr[i*2+1];
	}
	cx /= point_num;
	cy /= point_num;

	//	calculate covariance matrix
	Matrix2x2<PROMOTE_T_TO_FLOAT>	cov(0, 0, 0, 0);
	for(int i=0; i<point_num; i++){
		PROMOTE_T_TO_FLOAT x_dif = xy_ptr[i*2] - cx;
		PROMOTE_T_TO_FLOAT y_dif = xy_ptr[i*2+1] - cy;
		cov.data[0][0] += x_dif * x_dif;
		cov.data[0][1] += x_dif * y_dif;
		cov.data[1][1] += y_dif * y_dif;
	}
	cov.data[1][0] = cov.data[0][1];
	cov /= point_num - 1;

	//	calculate eigen values
	Com2<PROMOTE_T_TO_FLOAT>	evec1, evec2;
	PROMOTE_T_TO_FLOAT			eval1, eval2;
	int eigens = cov.Eigen(evec1, eval1, evec2, eval2);
	if( !eigens )
		return 0;

	//	write result
	mean_x	= cx;
	mean_y	= cy;
	pa_x	= evec1[0];
	pa_y	= evec1[1];

	return 1;
}

/**
	line fitting from 2D points with PCA

	for PCA, the error the perpendicular distance from each point to the line

	\param	mean			return center of the points
	\param	principle_axis	return principle axis
	\param	xy				list of the points
	\return					whether regression succeed
*/
template<typename T>
inline int fitLinePCA(Vec2<T>& mean, Vec2<T>& principle_axis, const std::vector<Vec2<T>>& xy){
	if( xy.empty() )
		return 0;

	return fitLinePCA2D(mean.x, mean.y, principle_axis.x, principle_axis.y, (T*)&xy[0], xy.size());
}




/**
	line fitting from weighted 2D points with PCA

	for PCA, the error the perpendicular distance from each point to the line

	\param	mean_x		return center of the points, x coordinate
	\param	mean_y		return center of the points, y coordinate
	\param	pa_x		return principle axis, x
	\param	pa_y		return principle axis, y
	\param	xy_ptr		pointer to the array of points, data arranged in xy_xy_xy_ format
	\param	w_ptr		pointer to the array of weight
	\param	point_num	number of points in the array
	\return				whether regression succeed
*/
template<typename T>
inline int fitLineWeightedPCA2D(T& mean_x, T& mean_y, T& pa_x, T& pa_y, const T* xy_ptr, const T* w_ptr, int point_num){
	if( point_num < 2 )
		return 0;

	//	calculate mean
	PROMOTE_T_TO_FLOAT cx = 0, cy = 0, w_sum = 0;
	for(int i=0; i<point_num; i++){
		cx += xy_ptr[i*2] * w_ptr[i];
		cy += xy_ptr[i*2+1] * w_ptr[i];
		w_sum += w_ptr[i];
	}
	cx /= w_sum;
	cy /= w_sum;

	//	calculate covariance matrix
	Matrix2x2<PROMOTE_T_TO_FLOAT>	cov(0, 0, 0, 0);
	for(int i=0; i<point_num; i++){
		PROMOTE_T_TO_FLOAT x_dif = xy_ptr[i*2] - cx;
		PROMOTE_T_TO_FLOAT y_dif = xy_ptr[i*2+1] - cy;
		cov.data[0][0] += x_dif * x_dif * w_ptr[i];
		cov.data[0][1] += x_dif * y_dif * w_ptr[i];
		cov.data[1][1] += y_dif * y_dif * w_ptr[i];
	}
	cov.data[1][0] = cov.data[0][1];
	cov /= point_num - 1;

	//	calculate eigen values
	Com2<PROMOTE_T_TO_FLOAT>	evec1, evec2;
	PROMOTE_T_TO_FLOAT			eval1, eval2;
	int eigens = cov.Eigen(evec1, eval1, evec2, eval2);
	if( !eigens )
		return 0;

	//	write result
	mean_x	= cx;
	mean_y	= cy;
	pa_x	= evec1[0];
	pa_y	= evec1[1];

	return 1;
}

/**
	line fitting from weighted 2D points with PCA

	for PCA, the error the perpendicular distance from each point to the line

	\param	mean			return center of the points
	\param	principle_axis	return principle axis
	\param	xy				list of the points
	\param	w				weight of each point
	\return					whether regression succeed
*/
template<typename T>
inline int fitLineWeightedPCA(Vec2<T>& mean, Vec2<T>& principle_axis, const std::vector<Vec2<T>>& xy, const std::vector<T>& w){
	if( xy.empty() || xy.size()!=w.size() )
		return 0;

	return fitLineWeightedPCA2D(mean.x, mean.y, principle_axis.x, principle_axis.y, (T*)&xy[0], (T*)&w[0], xy.size());
}

///@}


}}	//	namespace yz::statistics

#endif	//	__YZ_STATISTICS_H__