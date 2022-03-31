/***********************************************************/
/**	\file
	\brief		Matrix Type for Graphics
	\details	Matrix 2x2, 3x3, 4x4. 
				They are implemented as transform matrix defined
				in OpenGL coordinate
	\author		Yizhong Zhang
	\date		5/19/2012
*/
/***********************************************************/
#ifndef __YZ_MATRIX_H__
#define __YZ_MATRIX_H__

#include <assert.h>
#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_numerical_utils.h"

namespace yz{

template<class T> class Quaternion;
template<class T> class Matrix3x2;
template<class T> class Matrix2x3;

//	========================================
//	Matrix Type
//	========================================
/**
	2x2 Matrix
*/
template< class T >
class Matrix2x2{
public:
	T	data[2][2];

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Matrix2x2(const T* val_ptr){
		data[0][0] = val_ptr[0];
		data[0][1] = val_ptr[1];
		data[1][0] = val_ptr[2];
		data[1][1] = val_ptr[3];
	}
	Matrix2x2(T val00=0, T val01=0, T val10=0, T val11=0){
		data[0][0] = val00;
		data[0][1] = val01;
		data[1][0] = val10;
		data[1][1] = val11;
	}
	template<typename TYPE> Matrix2x2(const Matrix2x2<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
	}
	Matrix2x2(const Com2<T> vec1, const Com2<T> vec2){
		data[0][0] = vec1.x;
		data[0][1] = vec2.x;
		data[1][0] = vec1.y;
		data[1][1] = vec2.y;
	}

	//	operators
	inline T* operator[] (const int id){
		assert(id==0 || id==1);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else
			return data[1];
	}

	inline const T* operator[] (const int id) const{
		assert(id==0 || id==1);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else
			return data[1];
	}

	inline T& operator()(const int x, const int y){
		return data[x][y];
	}

	inline const T& operator()(const int x, const int y) const{
		return data[x][y];
	}

	inline Matrix2x2 operator - () const{
		Matrix2x2 result(*this);
		result.data[0][0] = - data[0][0];
		result.data[0][1] = - data[0][1];
		result.data[1][0] = - data[1][0];
		result.data[1][1] = - data[1][1];
		return result;
	}

	template<typename TYPE> inline Matrix2x2& operator = (const Matrix2x2<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		return *this;
	}

	template<typename TYPE> inline Matrix2x2& operator += (const Matrix2x2<TYPE> mat){
		data[0][0] += mat.data[0][0];
		data[0][1] += mat.data[0][1];
		data[1][0] += mat.data[1][0];
		data[1][1] += mat.data[1][1];
		return *this;
	}

	template<typename TYPE> inline Matrix2x2& operator -= (const Matrix2x2<TYPE> mat){
		data[0][0] -= mat.data[0][0];
		data[0][1] -= mat.data[0][1];
		data[1][0] -= mat.data[1][0];
		data[1][1] -= mat.data[1][1];
		return *this;
	}

	template<typename TYPE> inline Matrix2x2& operator *= (const Matrix2x2<TYPE> mat){
		T new_data[2][2];
		new_data[0][0] = data[0][0]*mat.data[0][0] + data[0][1]*mat.data[1][0];
		new_data[0][1] = data[0][0]*mat.data[0][1] + data[0][1]*mat.data[1][1];
		new_data[1][0] = data[1][0]*mat.data[0][0] + data[1][1]*mat.data[1][0];
		new_data[1][1] = data[1][0]*mat.data[0][1] + data[1][1]*mat.data[1][1];
		data[0][0] = new_data[0][0];
		data[0][1] = new_data[0][1];
		data[1][0] = new_data[1][0];
		data[1][1] = new_data[1][1];
		return *this;
	}

	template<typename TYPE> inline Matrix2x2& operator *= (const TYPE val){
		data[0][0] *= val;
		data[0][1] *= val;
		data[1][0] *= val;
		data[1][1] *= val;
		return *this;
	}

	template<typename TYPE> inline Matrix2x2& operator /= (const TYPE val){
		data[0][0] /= val;
		data[0][1] /= val;
		data[1][0] /= val;
		data[1][1] /= val;
		return *this;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator + (const Matrix2x2<TYPE> mat) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator - (const Matrix2x2<TYPE> mat) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix2x2<TYPE> mat) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix2x3<TYPE> mat) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result;
		result[0][0] = data[0][0] * mat[0][0] + data[0][1] * mat[1][0];
		result[0][1] = data[0][0] * mat[0][1] + data[0][1] * mat[1][1];
		result[1][0] = data[0][0] * mat[1][0] + data[0][1] * mat[1][2];
		result[1][0] = data[1][0] * mat[0][0] + data[1][1] * mat[1][0];
		result[1][1] = data[1][0] * mat[0][1] + data[1][1] * mat[1][1];
		result[1][2] = data[1][0] * mat[1][0] + data[1][1] * mat[1][2];
		return result;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator * (const Vec2<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0]*vec.x + data[0][1]*vec.y;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0]*vec.x + data[1][1]*vec.y;
		return Vec2<TYPE_PROMOTE(T, TYPE)>(new_x, new_y);
	}

	/**
		Return an 2x2 identity matrix
	*/
	inline Matrix2x2 Identity() const{
		Matrix2x2 result(1, 0, 0, 1);
		return result;
	}

	/**
		Return transposed matrix, but don't change old data
	*/
	inline Matrix2x2 Transpose() const{
		Matrix2x2 result(*this);
		result.SetTranspose();
		return result;
	}

	/**
		Return inversed matrix, but don't change old data
	*/
	inline Matrix2x2 Inverse() const{
		Matrix2x2 result(*this);
		result.SetInverse();
		return result;
	}

	//	set matrix
	inline Matrix2x2& SetIdentity(){
		data[0][0] = data[1][1] = 1;
		data[0][1] = data[1][0] = 0;
		return *this;
	}

	inline Matrix2x2& SetTranspose(){
		mySwap(data[0][1], data[1][0]);
		return *this;
	}

	inline Matrix2x2& SetInverse(){
		T det = Det();
		if( det != 0 ){
			mySwap(data[0][0], data[1][1]);
			data[0][0] /= det;
			data[0][1] /= -det;
			data[1][0] /= -det;
			data[1][1] /= det;
		}
		else{
			std::cout << "not invertable" << std::endl;
		}
		return *this;
	}

	/**
		Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix2x2& SetAsDiagonal(T diag1, T diag2) {
		data[0][0] = diag1;
		data[1][1] = diag2;
		data[0][1] = data[1][0] = 0;
		return *this;
	}

	/**
		Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix2x2& SetAsDiagonal(Com2<T> diag) {
		return SetAsDiagonal(diag[0], diag[1]);
	}

	/**
		Set the matrix to be a rotation matrix, angle_rad(rad) counter-clockwise
	*/
	inline Matrix2x2& SetRotationRad(YZ_REAL angle_rad){
		YZ_REAL cos_a = cos(angle_rad);
		YZ_REAL sin_a = sin(angle_rad);
		data[0][0] = data[1][1] = cos_a;
		data[0][1] = - sin_a;
		data[1][0] = sin_a;
		return *this;
	}

	/**
		Set the matrix to be a rotation matrix, angle_deg(degrees) counter-clockwise
	*/
	inline Matrix2x2& SetRotationDeg(YZ_REAL angle_deg){
		return SetRotationRad( angle_deg*YZ_PI/180 );
	}

	//	matrix properties
	inline void RetriveColumn(Com2<T>& vec1, Com2<T>& vec2) const {
		vec1.x = data[0][0];
		vec1.y = data[1][0];
		vec2.x = data[0][1];
		vec2.y = data[1][1];
	}

	inline int IsSymmetric() const {
		if (fabs(data[0][1] - data[1][0]) > 1e-6)
			return 0;
		return 1;
	}

	inline T Det() const{
		return data[0][0]*data[1][1] - data[0][1]*data[1][0];
	}

	inline T Trace() const{
		return data[0][0] + data[1][1];
	}

	/**
		get eigen vectors and eigen values of the matrix

		\param	eigen_vec1	return the eigen vector with bigger eigen value
		\param	eigen_val1	return the bigger eigen value
		\param	eigen_vec2	return the eigen vector with smaller eigen value
		\param	eigen_val2	return the smaller eigen value
		\return				the number of eigen values
	*/
	template<typename TYPE> 
	inline int Eigen(Com2<TYPE>& eigen_vec1, TYPE& eigen_val1, Com2<TYPE>& eigen_vec2, TYPE& eigen_val2) const{
		double	tra	= Trace();
		double	det	= Det();
		double	root = tra*tra*0.25 - det;
		if( root < -1e-6 )			//	no eigen values
			return 0;
		else if (root <= 0){		//	due to numerical error, we assume two equal eigen values in this case
			eigen_val1 = tra * 0.5;	//	in exact mathematics, this case will lead to two complex eigen values
			eigen_val2 = tra * 0.5;
			eigen_vec1[0] = 1;
			eigen_vec1[1] = 0;
			eigen_vec2[0] = 0;
			eigen_vec2[1] = 1;
			return 1;
		}

		double	eval1	= tra * 0.5 + sqrt(root);
		double	eval2	= tra - eval1;
		Vec2<double>	evec1, evec2;
		if( data[1][0] < -1e-5 || data[1][0] > 1e-5 ){
			evec1[0] = eval1 - data[1][1];
			evec1[1] = data[1][0];
			evec2[0] = eval2 - data[1][1];
			evec2[1] = data[1][0];
			evec1.SetNormalize();
			evec2.SetNormalize();
		}
		else if( data[0][1] < -1e-5 || data[0][1] > 1e-5 ){
			evec1[0] = data[0][1];
			evec1[1] = eval1 - data[0][0];
			evec2[0] = data[0][1];
			evec2[1] = eval2 - data[0][0];
			evec1.SetNormalize();
			evec2.SetNormalize();
		}
		else if( fabs(data[0][0]) > fabs(data[1][1]) ){
			evec1[0] = 1;
			evec1[1] = 0;
			evec2[0] = 0;
			evec2[1] = 1;
		}
		else{
			evec1[0] = 0;
			evec1[1] = 1;
			evec2[0] = 1;
			evec2[1] = 0;
		}

		eigen_vec1 = evec1;
		eigen_vec2 = evec2;
		eigen_val1 = eval1;
		eigen_val2 = eval2;

		if( root < 1e-6 )	//	two eigen values are almost the same
			return 1;
		else
			return 2;
	}
};

/**
	3x3 Matrix
*/
template< class T >
class Matrix3x3{
public:
	T	data[3][3];

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Matrix3x3(const T* val_ptr){
		data[0][0] = val_ptr[0];
		data[0][1] = val_ptr[1];
		data[0][2] = val_ptr[2];
		data[1][0] = val_ptr[3];
		data[1][1] = val_ptr[4];
		data[1][2] = val_ptr[5];
		data[2][0] = val_ptr[6];
		data[2][1] = val_ptr[7];
		data[2][2] = val_ptr[8];
	}
	Matrix3x3(T val00=0, T val01=0, T val02=0, T val10=0, T val11=0, T val12=0, T val20=0, T val21=0, T val22=0){
		data[0][0] = val00;
		data[0][1] = val01;
		data[0][2] = val02;
		data[1][0] = val10;
		data[1][1] = val11;
		data[1][2] = val12;
		data[2][0] = val20;
		data[2][1] = val21;
		data[2][2] = val22;
	}
	template<typename TYPE> Matrix3x3(const Matrix3x3<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
		data[2][2] = mat.data[2][2];
	}
	template<typename TYPE> Matrix3x3(const Quaternion<TYPE> quaternion){
		this->SetRotation(quaternion);
	}
	Matrix3x3(const Com3<T> vec1, const Com3<T> vec2, const Com3<T> vec3){
		data[0][0] = vec1.x;
		data[0][1] = vec2.x;
		data[0][2] = vec3.x;
		data[1][0] = vec1.y;
		data[1][1] = vec2.y;
		data[1][2] = vec3.y;
		data[2][0] = vec1.z;
		data[2][1] = vec2.z;
		data[2][2] = vec3.z;
	}

	//	operators
	inline T* operator[] (const int id){
		assert(id==0 || id==1 || id==2);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else if( id==1 )
			return data[1];
		else
			return data[2];
	}

	inline const T* operator[] (const int id) const{
		assert(id==0 || id==1 || id==2);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else if( id==1 )
			return data[1];
		else
			return data[2];
	}

	inline T& operator()(const int x, const int y){
		return data[x][y];
	}

	inline const T& operator()(const int x, const int y) const{
		return data[x][y];
	}

	inline Matrix3x3 operator - () const{
		Matrix3x3 result(*this);
		result.data[0][0] = - data[0][0];
		result.data[0][1] = - data[0][1];
		result.data[0][2] = - data[0][2];
		result.data[1][0] = - data[1][0];
		result.data[1][1] = - data[1][1];
		result.data[1][2] = - data[1][2];
		result.data[2][0] = - data[2][0];
		result.data[2][1] = - data[2][1];
		result.data[2][2] = - data[2][2];
		return result;
	}

	template<typename TYPE> inline Matrix3x3& operator = (const Matrix3x3<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
		data[2][2] = mat.data[2][2];
		return *this;
	}

	template<typename TYPE> inline Matrix3x3& operator = (const Quaternion<TYPE> quaternion){
		return SetRotation(quaternion);
	}

	template<typename TYPE> inline Matrix3x3& operator += (const Matrix3x3<TYPE> mat){
		data[0][0] += mat.data[0][0];
		data[0][1] += mat.data[0][1];
		data[0][2] += mat.data[0][2];
		data[1][0] += mat.data[1][0];
		data[1][1] += mat.data[1][1];
		data[1][2] += mat.data[1][2];
		data[2][0] += mat.data[2][0];
		data[2][1] += mat.data[2][1];
		data[2][2] += mat.data[2][2];
		return *this;
	}

	template<typename TYPE> inline Matrix3x3& operator -= (const Matrix3x3<TYPE> mat){
		data[0][0] -= mat.data[0][0];
		data[0][1] -= mat.data[0][1];
		data[0][2] -= mat.data[0][2];
		data[1][0] -= mat.data[1][0];
		data[1][1] -= mat.data[1][1];
		data[1][2] -= mat.data[1][2];
		data[2][0] -= mat.data[2][0];
		data[2][1] -= mat.data[2][1];
		data[2][2] -= mat.data[2][2];
		return *this;
	}

	template<typename TYPE> inline Matrix3x3& operator *= (const Matrix3x3<TYPE> mat){
		T new_data[3][3];
		new_data[0][0] = data[0][0]*mat.data[0][0] + data[0][1]*mat.data[1][0] + data[0][2]*mat.data[2][0];
		new_data[0][1] = data[0][0]*mat.data[0][1] + data[0][1]*mat.data[1][1] + data[0][2]*mat.data[2][1];
		new_data[0][2] = data[0][0]*mat.data[0][2] + data[0][1]*mat.data[1][2] + data[0][2]*mat.data[2][2];
		new_data[1][0] = data[1][0]*mat.data[0][0] + data[1][1]*mat.data[1][0] + data[1][2]*mat.data[2][0];
		new_data[1][1] = data[1][0]*mat.data[0][1] + data[1][1]*mat.data[1][1] + data[1][2]*mat.data[2][1];
		new_data[1][2] = data[1][0]*mat.data[0][2] + data[1][1]*mat.data[1][2] + data[1][2]*mat.data[2][2];
		new_data[2][0] = data[2][0]*mat.data[0][0] + data[2][1]*mat.data[1][0] + data[2][2]*mat.data[2][0];
		new_data[2][1] = data[2][0]*mat.data[0][1] + data[2][1]*mat.data[1][1] + data[2][2]*mat.data[2][1];
		new_data[2][2] = data[2][0]*mat.data[0][2] + data[2][1]*mat.data[1][2] + data[2][2]*mat.data[2][2];
		data[0][0] = new_data[0][0];
		data[0][1] = new_data[0][1];
		data[0][2] = new_data[0][2];
		data[1][0] = new_data[1][0];
		data[1][1] = new_data[1][1];
		data[1][2] = new_data[1][2];
		data[2][0] = new_data[2][0];
		data[2][1] = new_data[2][1];
		data[2][2] = new_data[2][2];
		return *this;
	}

	template<typename TYPE> inline Matrix3x3& operator *= (const TYPE val){
		data[0][0] *= val;
		data[0][1] *= val;
		data[0][2] *= val;
		data[1][0] *= val;
		data[1][1] *= val;
		data[1][2] *= val;
		data[2][0] *= val;
		data[2][1] *= val;
		data[2][2] *= val;
		return *this;
	}

	template<typename TYPE> inline Matrix3x3& operator /= (const TYPE val){
		data[0][0] /= val;
		data[0][1] /= val;
		data[0][2] /= val;
		data[1][0] /= val;
		data[1][1] /= val;
		data[1][2] /= val;
		data[2][0] /= val;
		data[2][1] /= val;
		data[2][2] /= val;
		return *this;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator + (const Matrix3x3<TYPE> mat) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator - (const Matrix3x3<TYPE> mat) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix3x3<TYPE> mat) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix3x2<TYPE> mat) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result;
		result[0][0] = data[0][0] * mat[0][0] + data[0][1] * mat[1][0] + data[0][2] * mat[2][0];
		result[0][1] = data[0][0] * mat[0][1] + data[0][1] * mat[1][1] + data[0][2] * mat[2][1];
		result[1][0] = data[1][0] * mat[0][0] + data[1][1] * mat[1][0] + data[1][2] * mat[2][0];
		result[1][1] = data[1][0] * mat[0][1] + data[1][1] * mat[1][1] + data[1][2] * mat[2][1];
		result[2][0] = data[2][0] * mat[0][0] + data[2][1] * mat[1][0] + data[2][2] * mat[2][0];
		result[2][1] = data[2][0] * mat[0][1] + data[2][1] * mat[1][1] + data[2][2] * mat[2][1];
		return result;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator * (const Vec3<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0]*vec.x + data[0][1]*vec.y + data[0][2]*vec.z;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0]*vec.x + data[1][1]*vec.y + data[1][2]*vec.z;
		TYPE_PROMOTE(T, TYPE) new_z = data[2][0]*vec.x + data[2][1]*vec.y + data[2][2]*vec.z;
		return Vec3<TYPE_PROMOTE(T, TYPE)>(new_x, new_y, new_z);
	}

	/**
		Return an 3x3 identity matrix
	*/	
	inline Matrix3x3 Identity() const{
		Matrix3x3 result(1, 0, 0, 0, 1, 0, 0, 0, 1);
		return result;
	}

	/**
		Return transposed matrix, but don't change old data
	*/
	inline Matrix3x3 Transpose() const{
		Matrix3x3 result(*this);
		result.SetTranspose();
		return result;
	}

	/**
		Return inversed matrix, but don't change old data
	*/
	inline Matrix3x3 Inverse() const{
		Matrix3x3 result(*this);
		result.SetInverse();
		return result;
	}

	//	set matrix
	inline Matrix3x3& SetIdentity(){
		data[0][0] = data[1][1] = data[2][2] = 1;
		data[0][1] = data[0][2] = data[1][0] = data[1][2] = data[2][0] = data[2][1] = 0;
		return *this;
	}

	inline Matrix3x3& SetTranspose(){
		mySwap(data[0][1], data[1][0]);
		mySwap(data[0][2], data[2][0]);
		mySwap(data[1][2], data[2][1]);
		return *this;
	}

	inline Matrix3x3& SetInverse(){
		T det = Det();
		if( det != 0 ){
			T inv[3][3];
			inv[0][0] = data[1][1]*data[2][2]-data[2][1]*data[1][2];
			inv[0][1] = data[0][2]*data[2][1]-data[2][2]*data[0][1];
			inv[0][2] = data[0][1]*data[1][2]-data[1][1]*data[0][2];
			inv[1][0] = data[1][2]*data[2][0]-data[2][2]*data[1][0];
			inv[1][1] = data[0][0]*data[2][2]-data[2][0]*data[0][2];
			inv[1][2] = data[0][2]*data[1][0]-data[1][2]*data[0][0];
			inv[2][0] = data[1][0]*data[2][1]-data[2][0]*data[1][1];
			inv[2][1] = data[0][1]*data[2][0]-data[2][1]*data[0][0];
			inv[2][2] = data[0][0]*data[1][1]-data[1][0]*data[0][1];
			data[0][0] = inv[0][0] / det;
			data[0][1] = inv[0][1] / det;
			data[0][2] = inv[0][2] / det;
			// data[0][3] = inv[0][3] / det;
			// data[0][4] = inv[0][4] / det;
			// data[0][5] = inv[0][5] / det;
			// data[0][6] = inv[0][6] / det;
			// data[0][7] = inv[0][7] / det;
			// data[0][8] = inv[0][8] / det;
			data[1][0] = inv[0][3] / det;
			data[1][1] = inv[0][4] / det;
			data[1][2] = inv[0][5] / det;
			data[2][0] = inv[0][6] / det;
			data[2][1] = inv[0][7] / det;
			data[2][2] = inv[0][8] / det;
		}
		else{
			std::cout << "not invertable" << std::endl;
		}
		return *this;
	}

	/**
		Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix3x3& SetAsDiagonal(T diag1, T diag2, T diag3) {
		data[0][0] = diag1;
		data[1][1] = diag2;
		data[2][2] = diag3;
		data[0][1] = data[0][2] = data[1][0] = data[1][2] = data[2][0] = data[2][1] = 0;
		return *this;
	}

	/**
		Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix3x3& SetAsDiagonal(Com3<T> diag) {
		return SetAsDiagonal(diag[0], diag[1], diag[2]);
	}

	/**
		Set the matrix according to quaternion

		\author	Chen Cao
	*/
	template<typename TYPE> 
	inline Matrix3x3& SetRotation(Quaternion<TYPE> quaternion){
		quaternion.SetNormalize();

		// Cal rotation matrix
		T w = quaternion.w;
		T x = quaternion.x;
		T y = quaternion.y;
		T z = quaternion.z;

		T w2 = w * w;
		T x2 = x * x;
		T y2 = y * y;
		T z2 = z * z;

		data[0][0] = w2 + x2 - y2 - z2;
		data[1][1] = w2 - x2 + y2 - z2;
		data[2][2] = w2 - x2 - y2 + z2;

		T tmp1 = x * y;
		T tmp2 = z * w;
		data[0][1] = 2 * (tmp1 - tmp2);
		data[1][0] = 2 * (tmp1 + tmp2);

		tmp1 = x * z;
		tmp2 = y * w;
		data[2][0] = 2 * (tmp1 - tmp2);
		data[0][2] = 2 * (tmp1 + tmp2);

		tmp1 = y * z;
		tmp2 = x * w;
		data[1][2] = 2 * (tmp1 - tmp2);
		data[2][1] = 2 * (tmp1 + tmp2);

		return *this;
	}

	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	cos_angle	cos(rotation_angle) value
		\param	sin_angle	sin(rotation_angle) value
		\return				*this
	*/
	inline Matrix3x3& SetRotation(Com3<YZ_REAL> axis, YZ_REAL cos_angle, YZ_REAL sin_angle){
		axis = Vec3<YZ_REAL>(axis).Normalize();

		data[0][0] = cos_angle + axis.x*axis.x*(1-cos_angle);
		data[0][1] = axis.x*axis.y*(1-cos_angle) - axis.z*sin_angle;
		data[0][2] = axis.x*axis.z*(1-cos_angle) + axis.y*sin_angle;
		data[1][0] = axis.y*axis.x*(1-cos_angle) + axis.z*sin_angle;
		data[1][1] = cos_angle + axis.y*axis.y*(1-cos_angle);
		data[1][2] = axis.y*axis.z*(1-cos_angle) - axis.x*sin_angle;
		data[2][0] = axis.x*axis.z*(1-cos_angle) - axis.y*sin_angle;
		data[2][1] = axis.y*axis.z*(1-cos_angle) + axis.x*sin_angle;
		data[2][2] = cos_angle + axis.z*axis.z*(1-cos_angle);

		return *this;
	}

	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	angle_rad	rotation angle in rad
		\return				*this
	*/
	inline Matrix3x3& SetRotationRad(Com3<YZ_REAL> axis, YZ_REAL angle_rad){
		YZ_REAL cos_angle = cos(angle_rad);
		YZ_REAL sin_angle = sin(angle_rad);

		return SetRotation(axis, cos_angle, sin_angle);
	}

	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	angle_deg	rotation angle in degree
		\return				*this
	*/
	inline Matrix3x3& SetRotationDeg(Com3<YZ_REAL> axis, YZ_REAL angle_deg){
		return SetRotationRad( axis, angle_deg*YZ_PI/180 );
	}

	/**
		Set the matrix to be a scale matrix

		\param	x		scale value in x direction
		\param	y		scale value in y direction
		\param	z		scale value in z direction
		\return			*this
	*/
	inline Matrix3x3& SetScale(T x, T y, T z){
		return SetAsDiagonal(x, y, z);
	}

	/**
		Set the matrix to be a scale matrix

		\param	scale	scale value in all directions
		\return			*this
	*/
	inline Matrix3x3& SetScale(T scale){
		return SetScale(scale, scale, scale);
	}

	/**
		Set the matrix to be a scale matrix

		\param	scale	scale value in xyz directions
		\return			*this
	*/
	inline Matrix3x3& SetScale(Com3<T> scale){
		return SetScale(scale.x, scale.y, scale.z);
	}

	//	matrix properties
	inline void RetriveColumn(Com3<T>& vec1, Com3<T>& vec2, Com3<T>& vec3) const {
		vec1.x = data[0][0];
		vec1.y = data[1][0];
		vec1.z = data[2][0];
		vec2.x = data[0][1];
		vec2.y = data[1][1];
		vec2.z = data[2][1];
		vec3.x = data[0][2];
		vec3.y = data[1][2];
		vec3.z = data[2][2];
	}

	inline int IsSymmetric() const {
		if (fabs(data[0][1] - data[1][0]) > 1e-6)
			return 0;
		if (fabs(data[0][2] - data[2][0]) > 1e-6)
			return 0;
		if (fabs(data[1][2] - data[2][1]) > 1e-6)
			return 0;
		return 1;
	}

	inline T Det() const{
		return		data[0][0]*data[1][1]*data[2][2]
				-	data[0][0]*data[1][2]*data[2][1]
				-	data[0][1]*data[1][0]*data[2][2]
				+	data[0][1]*data[1][2]*data[2][0]
				+	data[0][2]*data[1][0]*data[2][1]
				-	data[0][2]*data[1][1]*data[2][0];
	}

	inline T Trace() const{
		return data[0][0] + data[1][1] + data[2][2];
	}

	/**
	Calculate eigen values of 3x3 symmetric matrix

	The algorithm comes from wikipedia, https://en.wikipedia.org/wiki/Eigenvalue_algorithm
	eigen_val1 >= eigen_val2 >> eigen_val3

	\return		number of eigen values (should be 3. for non-symmetric matrix, return 0)
	*/
	template<typename TYPE>
	inline int EigenSymm(TYPE& eigen_val1, TYPE& eigen_val2, TYPE& eigen_val3) const {
		if (!IsSymmetric())	//	this algorithm only deal with symmetric matrix
			return 0;

		double p1 = data[0][1] * data[0][1] + data[0][2] * data[0][2] + data[1][2] * data[1][2];
		if (p1 < 1e-6) {
			eigen_val1 = data[0][0];
			eigen_val2 = data[1][1];
			eigen_val3 = data[2][2];
			if (eigen_val1 < eigen_val2)
				mySwap(eigen_val1, eigen_val2);
			if (eigen_val1 < eigen_val3)
				mySwap(eigen_val1, eigen_val3);
			if (eigen_val2 < eigen_val3)
				mySwap(eigen_val2, eigen_val3);
		}
		else {
			double q = Trace() / 3;
			double p2 = (data[0][0] - q) * (data[0][0] - q) +
				(data[1][1] - q) * (data[1][1] - q) +
				(data[2][2] - q) * (data[2][2] - q) + 2 * p1;
			double p = sqrt(p2 / 6);
			Matrix3x3<double> B = *this;
			B.data[0][0] -= q;
			B.data[1][1] -= q;
			B.data[2][2] -= q;
			B /= p;
			double r = B.Det() * 0.5;

			double phi = 0;
			if (r <= -1.0)
				phi = YZ_PI / 3;
			else if (r < 1.0)
				phi = acos(r) / 3;

			eigen_val1 = q + 2 * p * cos(phi);
			eigen_val3 = q + 2 * p * cos(phi + YZ_PI * 2 / 3);
			eigen_val2 = 3 * q - eigen_val1 - eigen_val3;
		}

		return 3;
	}

	/**
	Eigen decomposition of 3x3 symmetric matrix
	(*this) = U * S.AsDiagonal() * U.Inverse()

	\return		number of eigen values (should be 3. return 0 if matrix is not symmetric)
	*/
	template<typename TYPE>
	inline int EigenSymm(Matrix3x3<TYPE>& U, Vec3<TYPE>& S) const {
		int res = EigenSymm(S[0], S[1], S[2]);
		if (res == 0)	//	non-symmetric matrix, calculate failed
			return 0;

		//	three identical eigen values, U could be any rotational matrix, simply set it to I
		if (S[0] - S[1] < 1e-6 && S[1] - S[2] < 1e-6) {	//	S is already sorted from biggest to smallest
			U.SetIdentity();
			return res;
		}

		//	two identical eigen values, first calculate the unique vector, the other two just need to be any orthogonal pairs
		int unique_id = -1;
		if (S[0] - S[1] < 1e-6)
			unique_id = 2;
		else if (S[1] - S[2] < 1e-6)
			unique_id = 0;
		if (unique_id >= 0) {
			//	calculate eigen vector associated with the unique eigen value
			Matrix3x3<TYPE> C = (*this - Identity() * S[(unique_id + 1) % 3]) * (*this - Identity() * S[(unique_id + 2) % 3]);
			Vec3<TYPE> V[3];
			C.RetriveColumn(V[0], V[1], V[2]);
			int idx = 0;
			TYPE max_squ_len = V[0].SquareLength();
			for (int j = 1; j < 3; j++) {
				TYPE squ_len = V[j].SquareLength();
				if (squ_len > max_squ_len) {
					max_squ_len = squ_len;
					idx = j;
				}
			}
			Vec3<TYPE> unique_V = V[idx].Normalize();

			//	calculate orthogonal vector pairs, that also orthogonal to unique_V
			int max_dir = 0;
			max_squ_len = unique_V[max_dir] * unique_V[max_dir];
			for (int i = 1; i < 3; i++) {
				TYPE squ_len = unique_V[i] * unique_V[i];
				if (squ_len > max_squ_len) {
					max_squ_len = squ_len;
					max_dir = i;
				}
			}
			Vec3<TYPE> ortho_V1(0, 0, 0);
			ortho_V1[(max_dir + 1) % 3] = 1;
			Vec3<TYPE> ortho_V2 = cross(unique_V, ortho_V1);
			ortho_V1 = cross(unique_V, ortho_V2);

			if (unique_id == 0)
				U = Matrix3x3<TYPE>(unique_V, ortho_V1, ortho_V2);
			else
				U = Matrix3x3<TYPE>(ortho_V1, ortho_V2, unique_V);

			return res;
		}

		//	three distinct eigen values
		Matrix3x3<TYPE> cross_term[3];
		for (int i = 0; i < 3; i++)
			cross_term[i] = (*this - Identity() * S[i]);

		for (int i = 0; i < 3; i++) {
			Matrix3x3<TYPE> C = cross_term[(i + 1) % 3] * cross_term[(i + 2) % 3];
			Vec3<TYPE> V[3];
			C.RetriveColumn(V[0], V[1], V[2]);
			int idx = 0;
			TYPE max_squ_len = V[0].SquareLength();
			for (int j = 1; j < 3; j++) {
				TYPE squ_len = V[j].SquareLength();
				if (squ_len > max_squ_len) {
					max_squ_len = squ_len;
					idx = j;
				}
			}
			V[idx].SetNormalize();
			U[0][i] = V[idx][0];
			U[1][i] = V[idx][1];
			U[2][i] = V[idx][2];
		}

		return res;
	}

};

/**
	4x4 Matrix
*/
template< class T >
class Matrix4x4{
public:
	T	data[4][4];

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Matrix4x4(){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] = 0;
	}
	Matrix4x4(const T* val_ptr){
		for( int i=0; i<16; i++ )
			data[0][i] = val_ptr[i];
	}
	template<typename TYPE> Matrix4x4(const Matrix4x4<TYPE> mat){
		for( int i=0; i<16; i++ )
			data[0][i] = mat.data[0][i];
	}
	template<typename TYPE> Matrix4x4(const Matrix3x3<TYPE> mat){
		//	treate them as transform matrix
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[0][3] = 0;
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
		data[1][3] = 0;
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
		data[2][2] = mat.data[2][2];
		data[2][3] = 0;
		data[3][0] = data[3][1] = data[3][2] = 0;
		data[3][3] = 1;
	}
	template<typename TYPE> Matrix4x4(const Quaternion<TYPE> quaternion){
		SetRotation(quaternion);
	}
	template<typename T1, typename T2> Matrix4x4(const Matrix3x3<T1> rot_mat, const Vec3<T2> trans_vec){
		data[0][0] = rot_mat.data[0][0];
		data[0][1] = rot_mat.data[0][1];
		data[0][2] = rot_mat.data[0][2];
		data[0][3] = trans_vec.x;
		data[1][0] = rot_mat.data[1][0];
		data[1][1] = rot_mat.data[1][1];
		data[1][2] = rot_mat.data[1][2];
		data[1][3] = trans_vec.y;
		data[2][0] = rot_mat.data[2][0];
		data[2][1] = rot_mat.data[2][1];
		data[2][2] = rot_mat.data[2][2];
		data[2][3] = trans_vec.z;
		data[3][0] = data[3][1] = data[3][2] = 0;
		data[3][3] = 1;
	}

	Matrix4x4(const Com4<T> vec1, const Com4<T> vec2, const Com4<T> vec3, const Com4<T> vec4){
		data[0][0] = vec1.x;
		data[0][1] = vec2.x;
		data[0][2] = vec3.x;
		data[0][3] = vec4.x;
		data[1][0] = vec1.y;
		data[1][1] = vec2.y;
		data[1][2] = vec3.y;
		data[1][3] = vec4.y;
		data[2][0] = vec1.z;
		data[2][1] = vec2.z;
		data[2][2] = vec3.z;
		data[2][3] = vec4.z;
		data[3][0] = vec1.w;
		data[3][1] = vec2.w;
		data[3][2] = vec3.w;
		data[3][3] = vec4.w;
	}

	//	operators
	inline T* operator[] (const int id){
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else if( id==1 )
			return data[1];
		else if( id==2 )
			return data[2];
		else
			return data[3];
	}

	inline const T* operator[] (const int id) const{
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return data[0];
		else if( id==1 )
			return data[1];
		else if( id==2 )
			return data[2];
		else
			return data[3];
	}

	inline T& operator()(const int x, const int y){
		return data[x][y];
	}

	inline const T& operator()(const int x, const int y) const{
		return data[x][y];
	}

	inline Matrix4x4 operator - () const{
		Matrix4x4 result(*this);
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				result.data[j][i] = -data[j][i];
		return result;
	}

	template<typename TYPE> inline Matrix4x4& operator = (const Matrix4x4<TYPE> mat){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] = mat.data[j][i];
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator = (const Matrix3x3<TYPE> mat){
		//	treate them as transform matrix
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[0][3] = 0;
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
		data[1][3] = 0;
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
		data[2][2] = mat.data[2][2];
		data[2][3] = 0;
		data[3][0] = 0;
		data[3][1] = 0;
		data[3][2] = 0;
		data[3][3] = 1;
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator = (const Quaternion<TYPE> quaternion){
		return SetRotation(quaternion);
	}

	template<typename TYPE> inline Matrix4x4& operator += (const Matrix4x4<TYPE> mat){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] += mat.data[j][i];
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator -= (const Matrix4x4<TYPE> mat){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] -= mat.data[j][i];
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator *= (const Matrix4x4<TYPE> mat){
		T new_data[4][4];
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				new_data[j][i] = 0;
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				for (int k = 0; k < 4; k++)
					new_data[j][i] += data[j][k] * mat.data[k][i];
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] = new_data[j][i];
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator *= (const TYPE val){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] *= val;
		return *this;
	}

	template<typename TYPE> inline Matrix4x4& operator /= (const TYPE val){
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[j][i] /= val;
		return *this;
	}

	template<typename TYPE> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator + (const Matrix4x4<TYPE> mat) const{
		Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += mat;
		return result;
	}

	template<typename TYPE> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator - (const Matrix4x4<TYPE> mat) const{
		Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= mat;
		return result;
	}

	template<typename TYPE> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix4x4<TYPE> mat) const{
		Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= mat;
		return result;
	}

	template<typename TYPE> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	template<typename TYPE> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator * (const Vec4<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0]*vec.x + data[0][1]*vec.y + data[0][2]*vec.z + data[0][3]*vec.w;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0]*vec.x + data[1][1]*vec.y + data[1][2]*vec.z + data[1][3]*vec.w;
		TYPE_PROMOTE(T, TYPE) new_z = data[2][0]*vec.x + data[2][1]*vec.y + data[2][2]*vec.z + data[2][3]*vec.w;
		TYPE_PROMOTE(T, TYPE) new_w = data[3][0]*vec.x + data[3][1]*vec.y + data[3][2]*vec.z + data[3][3]*vec.w;
		return Vec4<TYPE_PROMOTE(T, TYPE)>(new_x, new_y, new_z, new_w);
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator * (const Vec3<TYPE> vec) const{
		return (*this) * Vec4<TYPE>(vec, 1);
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> TransformVector(const Vec3<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0]*vec.x + data[0][1]*vec.y + data[0][2]*vec.z + data[0][3];
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0]*vec.x + data[1][1]*vec.y + data[1][2]*vec.z + data[1][3];
		TYPE_PROMOTE(T, TYPE) new_z = data[2][0]*vec.x + data[2][1]*vec.y + data[2][2]*vec.z + data[2][3];
		return Vec3<TYPE_PROMOTE(T, TYPE)>(new_x, new_y, new_z);
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> RotateVector(const Vec3<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0]*vec.x + data[0][1]*vec.y + data[0][2]*vec.z;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0]*vec.x + data[1][1]*vec.y + data[1][2]*vec.z;
		TYPE_PROMOTE(T, TYPE) new_z = data[2][0]*vec.x + data[2][1]*vec.y + data[2][2]*vec.z;
		return Vec3<TYPE_PROMOTE(T, TYPE)>(new_x, new_y, new_z);
	}

	/**
		Return an 4x4 identity matrix
	*/	
	inline Matrix4x4 Identity() const{
		Matrix4x4 result(*this);
		result.SetIdentity();
		return result;
	}

	/**
		Return transposed matrix, but don't change old data
	*/
	inline Matrix4x4 Transpose() const{
		Matrix4x4 result(*this);
		result.SetTranspose();
		return result;
	}

	/**
		Return inversed matrix, but don't change old data
	*/
	inline Matrix4x4 Inverse() const{
		Matrix4x4 result(*this);
		result.SetInverse();
		return result;
	}

	//	set matrix
	inline Matrix4x4& SetIdentity(){
		for( int j=0; j<4; j++ )
			for( int i=0; i<4; i++ )
				data[j][i] = (j==i);
		return *this;
	}

	inline Matrix4x4& SetTranspose(){
		for( int j=0; j<4; j++ )
			for( int i=j+1; i<4; i++ )
				mySwap(data[j][i], data[i][j]);
		return *this;
	}

	inline Matrix4x4& SetInverse(){
		T inv[16], det;

		inv[0] = data[1][1]  * data[2][2] * data[3][3] -               
			data[1][1]  * data[2][3] * data[3][2] -               
			data[2][1]  * data[1][2]  * data[3][3] +               
			data[2][1]  * data[1][3]  * data[3][2] +              
			data[3][1] * data[1][2]  * data[2][3] -               
			data[3][1] * data[1][3]  * data[2][2];      

		inv[4] = -data[1][0]  * data[2][2] * data[3][3] +                
			data[1][0]  * data[2][3] * data[3][2] +                
			data[2][0]  * data[1][2]  * data[3][3] -                
			data[2][0]  * data[1][3]  * data[3][2] -                
			data[3][0] * data[1][2]  * data[2][3] +                
			data[3][0] * data[1][3]  * data[2][2];      

		inv[8] = data[1][0]  * data[2][1] * data[3][3] -               
			data[1][0]  * data[2][3] * data[3][1] -               
			data[2][0]  * data[1][1] * data[3][3] +               
			data[2][0]  * data[1][3] * data[3][1] +               
			data[3][0] * data[1][1] * data[2][3] -               
			data[3][0] * data[1][3] * data[2][1];      

		inv[12] = -data[1][0]  * data[2][1] * data[3][2] +                 
			data[1][0]  * data[2][2] * data[3][1] +                
			data[2][0]  * data[1][1] * data[3][2] -                 
			data[2][0]  * data[1][2] * data[3][1] -                 
			data[3][0] * data[1][1] * data[2][2] +                 
			data[3][0] * data[1][2] * data[2][1];      

		inv[1] = -data[0][1]  * data[2][2] * data[3][3] +                
			data[0][1]  * data[2][3] * data[3][2] +                
			data[2][1]  * data[0][2] * data[3][3] -                
			data[2][1]  * data[0][3] * data[3][2] -                
			data[3][1] * data[0][2] * data[2][3] +                
			data[3][1] * data[0][3] * data[2][2];      

		inv[5] = data[0][0]  * data[2][2] * data[3][3] -               
			data[0][0]  * data[2][3] * data[3][2] -               
			data[2][0]  * data[0][2] * data[3][3] +               
			data[2][0]  * data[0][3] * data[3][2] +               
			data[3][0] * data[0][2] * data[2][3] -               
			data[3][0] * data[0][3] * data[2][2];      

		inv[9] = -data[0][0]  * data[2][1] * data[3][3] +                
			data[0][0]  * data[2][3] * data[3][1] +                
			data[2][0]  * data[0][1] * data[3][3] -                
			data[2][0]  * data[0][3] * data[3][1] -                
			data[3][0] * data[0][1] * data[2][3] +                
			data[3][0] * data[0][3] * data[2][1];      

		inv[13] = data[0][0]  * data[2][1] * data[3][2] -                
			data[0][0]  * data[2][2] * data[3][1] -                
			data[2][0]  * data[0][1] * data[3][2] +                
			data[2][0]  * data[0][2] * data[3][1] +                
			data[3][0] * data[0][1] * data[2][2] -                
			data[3][0] * data[0][2] * data[2][1];      

		inv[2] = data[0][1]  * data[1][2] * data[3][3] -               
			data[0][1]  * data[1][3] * data[3][2] -               
			data[1][1]  * data[0][2] * data[3][3] +               
			data[1][1]  * data[0][3] * data[3][2] +               
			data[3][1] * data[0][2] * data[1][3] -               
			data[3][1] * data[0][3] * data[1][2];      

		inv[6] = -data[0][0]  * data[1][2] * data[3][3] +                
			data[0][0]  * data[1][3] * data[3][2] +                
			data[1][0]  * data[0][2] * data[3][3] -                
			data[1][0]  * data[0][3] * data[3][2] -                
			data[3][0] * data[0][2] * data[1][3] +                
			data[3][0] * data[0][3] * data[1][2];      

		inv[10] = data[0][0]  * data[1][1] * data[3][3] -                
			data[0][0]  * data[1][3] * data[3][1] -                
			data[1][0]  * data[0][1] * data[3][3] +                
			data[1][0]  * data[0][3] * data[3][1] +                
			data[3][0] * data[0][1] * data[1][3] -                
			data[3][0] * data[0][3] * data[1][1];      

		inv[14] = -data[0][0]  * data[1][1] * data[3][2] +                 
			data[0][0]  * data[1][2] * data[3][1] +                 
			data[1][0]  * data[0][1] * data[3][2] -                 
			data[1][0]  * data[0][2] * data[3][1] -                 
			data[3][0] * data[0][1] * data[1][2] +                 
			data[3][0] * data[0][2] * data[1][1];      

		inv[3] = -data[0][1] * data[1][2] * data[2][3] +                
			data[0][1] * data[1][3] * data[2][2] +                
			data[1][1] * data[0][2] * data[2][3] -                
			data[1][1] * data[0][3] * data[2][2] -                
			data[2][1] * data[0][2] * data[1][3] +                
			data[2][1] * data[0][3] * data[1][2];      

		inv[7] = data[0][0] * data[1][2] * data[2][3] -               
			data[0][0] * data[1][3] * data[2][2] -               
			data[1][0] * data[0][2] * data[2][3] +               
			data[1][0] * data[0][3] * data[2][2] +               
			data[2][0] * data[0][2] * data[1][3] -               
			data[2][0] * data[0][3] * data[1][2];      

		inv[11] = -data[0][0] * data[1][1] * data[2][3] +                 
			data[0][0] * data[1][3] * data[2][1] +                 
			data[1][0] * data[0][1] * data[2][3] -                 
			data[1][0] * data[0][3] * data[2][1] -                 
			data[2][0] * data[0][1] * data[1][3] +                 
			data[2][0] * data[0][3] * data[1][1];      

		inv[15] = data[0][0] * data[1][1] * data[2][2] -                
			data[0][0] * data[1][2] * data[2][1] -                
			data[1][0] * data[0][1] * data[2][2] +                
			data[1][0] * data[0][2] * data[2][1] +                
			data[2][0] * data[0][1] * data[1][2] -                
			data[2][0] * data[0][2] * data[1][1];      

		det = data[0][0] * inv[0] + data[0][1] * inv[4] + data[0][2] * inv[8] + data[0][3] * inv[12];

		if( det != 0 ){
			for (int i = 0; i < 16; i++)         
				data[0][i] = inv[i] / det;
		}
		else{
			std::cout << "not invertable" << std::endl;
		}

		return *this;
	}

	/**
	Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix4x4& SetAsDiagonal(T diag1, T diag2, T diag3, T diag4) {
		data[0][0] = diag1;
		data[0][1] = 0;
		data[0][2] = 0;
		data[0][3] = 0;
		data[1][0] = 0;
		data[1][1] = diag2;
		data[1][2] = 0;
		data[1][3] = 0;
		data[2][0] = 0;
		data[2][1] = 0;
		data[2][2] = diag3;
		data[2][3] = 0;
		data[3][0] = 0;
		data[3][1] = 0;
		data[3][2] = 0;
		data[3][3] = diag4;
		return *this;
	}
	/**
	Set the matrix to be diagonal matrix, with given diagonal values
	*/
	inline Matrix4x4& SetAsDiagonal(Com4<T> diag) {
		return SetAsDiagonal(diag[0], diag[1], diag[2], diag[3]);
	}
	/**
		set rotation matrix by quaternion
	*/
	template<typename TYPE>
	inline Matrix4x4& SetRotation(const Quaternion<TYPE> quaternion){
		Matrix3x3<T> mat3(quaternion);
		(*this) = mat3;
		return *this;
	}
	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	cos_angle	cos(rotation_angle) value
		\param	sin_angle	sin(rotation_angle) value
		\return				*this
	*/
	inline Matrix4x4& SetRotation(Com3<YZ_REAL> axis, YZ_REAL cos_angle, YZ_REAL sin_angle){
		*this = Matrix3x3<T>().SetRotation(axis, cos_angle, sin_angle);
		return *this;
	}
	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	angle_rad	rotation angle in rad
		\return				*this
	*/
	inline Matrix4x4& SetRotationRad(Com3<YZ_REAL> axis, YZ_REAL angle_rad){
		*this = Matrix3x3<T>().SetRotationRad(axis, angle_rad);
		return *this;
	}
	/**
		Set the matrix to be a rotation matrix along axis, counter-clockwise

		\param	axis		rotation axis
		\param	angle_deg	rotation angle in degree
		\return				*this
	*/
	inline Matrix4x4& SetRotationDeg(Com3<YZ_REAL> axis, YZ_REAL angle_deg){
		return SetRotationRad( axis, angle_deg*YZ_PI/180 );
	}
	/**
		Set the matrix to be a scale matrix

		\param	x		scale value in x direction
		\param	y		scale value in y direction
		\param	z		scale value in z direction
		\return			*this
	*/
	inline Matrix4x4& SetScale(T x, T y, T z){
		SetAsDiagonal(x, y, z, 1);
		return *this;
	}
	/**
		Set the matrix to be a scale matrix

		\param	scale	scale value in all directions
		\return			*this
	*/
	inline Matrix4x4& SetScale(T scale){
		return SetScale(scale, scale, scale);
	}
	/**
		Set the matrix to be a scale matrix

		\param	scale	scale value in xyz directions
		\return			*this
	*/
	inline Matrix4x4& SetScale(Com3<T> scale){
		return SetScale(scale.x, scale.y, scale.z);
	}
	/**
		Set the matrix to be a translation matrix

		\param	vec		translate value in xyz directions
		\return			*this
	*/
	inline Matrix4x4& SetTranslation(Com3<T> vec){
		SetIdentity();
		data[0][3] = vec.x;
		data[1][3] = vec.y;
		data[2][3] = vec.z;
		return *this;
	}
	/**
		Set the matrix to be a translation matrix

		\param	x		translate value in x direction
		\param	y		translate value in y direction
		\param	z		translate value in z direction
		\return			*this
	*/
	inline Matrix4x4& SetTranslation(T x, T y, T z){
		return SetTranslation( Vec3<T>(x, y, z) );
	}
	/**
		Set the matrix to be a transform matrix as OpenGL LookAt

		\param	eye		eye position
		\param	center	the center of the eye looking at
		\param	up		up direction of eye
		\return			*this
	*/
	inline Matrix4x4& SetLookAt(Com3<T> eye, Com3<T> center, Com3<T> up){
		//	set the matrix to be the same as calling gluLookAt, without transpose
		typedef YZ_REAL real;

		Vec3<real> f(center.x-eye.x, center.y-eye.y, center.z-eye.z);
		Vec3<real> s = cross(f, Vec3<real>(up));
		Vec3<real> u = cross(s, f);

		f.SetNormalize();
		s.SetNormalize();
		u.SetNormalize();

		real rot_mat[16] = {
			s.x,	s.y,	s.z,	0,
			u.x,	u.y,	u.z,	0,
			-f.x,	-f.y,	-f.z,	0,
			0,		0,		0,		1	};

		real trans_mat[16]	= {
			1,	0,	0,	-eye.x,
			0,	1,	0,	-eye.y,
			0,	0,	1,	-eye.z,
			0,	0,	0,	1	};

		Matrix4x4<real> rotM(rot_mat);
		Matrix4x4<real> transM(trans_mat);

		*this = rotM * transM;
		return *this;
	}
	/**
		set the matrix to be frustum projection
	*/
	inline Matrix4x4& SetFrustum(T left, T right, T bottom, T top, T znear, T zfar){
		T temp, temp2, temp3, temp4;
		temp = 2.0 * znear;
		temp2 = right - left;
		temp3 = top - bottom;
		temp4 = zfar - znear;
		data[0][0] = temp / temp2;
		data[0][1] = 0.0;
		data[0][2] = 0.0;
		data[0][3] = 0.0;
		data[1][0] = 0.0;
		data[1][1] = temp / temp3;
		data[1][2] = 0.0;
		data[1][3] = 0.0;
		data[2][0] = (right + left) / temp2;
		data[2][1] = (top + bottom) / temp3;
		data[2][2] = (-zfar - znear) / temp4;
		data[2][3] = -1.0;
		data[3][0] = 0.0;
		data[3][1] = 0.0;
		data[3][2] = (-temp * zfar) / temp4;
		data[3][3] = 0.0;
		return *this;
	}
	/**
		set the matrix to be perspective projection
	*/
	inline Matrix4x4& SetPerspective(T fovy_deg, T aspect, T znear, T zfar){
		T ymax, xmax;
		T temp, temp2, temp3, temp4;
		ymax = znear * tanf(fovy_deg * YZ_PI / 360.0);
		//ymin = -ymax;
		//xmin = -ymax * aspectRatio;
		xmax = ymax * aspect;

		SetFrustum(-xmax, xmax, -ymax, ymax, znear, zfar);

		return *this;
	}

	//	matrix properties
	inline void RetriveColumn(Com4<T>& vec1, Com4<T>& vec2, Com4<T>& vec3, Com4<T>& vec4) const {
		vec1.x = data[0][0];
		vec1.y = data[1][0];
		vec1.z = data[2][0];
		vec1.w = data[3][0];
		vec2.x = data[0][1];
		vec2.y = data[1][1];
		vec2.z = data[2][1];
		vec2.w = data[3][1];
		vec3.x = data[0][2];
		vec3.y = data[1][2];
		vec3.z = data[2][2];
		vec3.w = data[3][2];
		vec4.x = data[0][3];
		vec4.y = data[1][3];
		vec4.z = data[2][3];
		vec4.w = data[3][3];
	}

	inline int IsSymmetric() const {
		if (fabs(data[0][1] - data[1][0]) > 1e-6)
			return 0;
		if (fabs(data[0][2] - data[2][0]) > 1e-6)
			return 0;
		if (fabs(data[0][3] - data[3][0]) > 1e-6)
			return 0;
		if (fabs(data[1][2] - data[2][1]) > 1e-6)
			return 0;
		if (fabs(data[1][3] - data[3][1]) > 1e-6)
			return 0;
		if (fabs(data[2][3] - data[3][2]) > 1e-6)
			return 0;
		return 1;
	}

	inline T Det() const{
		T tmp[4];

		tmp[0] = data[1][1]  * data[2][2] * data[3][3] -               
			data[1][1]  * data[2][3] * data[3][2] -               
			data[2][1]  * data[1][2]  * data[3][3] +               
			data[2][1]  * data[1][3]  * data[3][2] +              
			data[3][1] * data[1][2]  * data[2][3] -               
			data[3][1] * data[1][3]  * data[2][2];      

		tmp[1] = -data[1][0]  * data[2][2] * data[3][3] +                
			data[1][0]  * data[2][3] * data[3][2] +                
			data[2][0]  * data[1][2]  * data[3][3] -                
			data[2][0]  * data[1][3]  * data[3][2] -                
			data[3][0] * data[1][2]  * data[2][3] +                
			data[3][0] * data[1][3]  * data[2][2];      

		tmp[2] = data[1][0]  * data[2][1] * data[3][3] -               
			data[1][0]  * data[2][3] * data[3][1] -               
			data[2][0]  * data[1][1] * data[3][3] +               
			data[2][0]  * data[1][3] * data[3][1] +               
			data[3][0] * data[1][1] * data[2][3] -               
			data[3][0] * data[1][3] * data[2][1];      

		tmp[3] = -data[1][0]  * data[2][1] * data[3][2] +                 
			data[1][0]  * data[2][2] * data[3][1] +                
			data[2][0]  * data[1][1] * data[3][2] -                 
			data[2][0]  * data[1][2] * data[3][1] -                 
			data[3][0] * data[1][1] * data[2][2] +                 
			data[3][0] * data[1][2] * data[2][1];  

		return data[0][0] * tmp[0] + data[0][1] * tmp[1] + data[0][2] * tmp[2] + data[0][3] * tmp[3];
	}

	inline T Trace() const{
		return data[0][0] + data[1][1] + data[2][2] + data[3][3];
	}
	/**
		Get Rotation Matrix from the full transform matrix
	*/
	inline Matrix3x3<YZ_REAL> RotationMatrix(){
		Matrix3x3<YZ_REAL> rot_mat(
			data[0][0], data[0][1], data[0][2], 
			data[1][0], data[1][1], data[1][2], 
			data[2][0], data[2][1], data[2][2]);
		return rot_mat;
	}

	/**
		Get Translation Vector from the full transform matrix
	*/
	inline Vec3<YZ_REAL> TranslationVector(){
		return Vec3<YZ_REAL>(data[0][3], data[1][3], data[2][3]);
	}
};

/**
	3x2 Matrix

	Non-square matrix is often used in 2D-3D mapping
*/
template< class T >
class Matrix3x2{
public:
	T	data[3][2];

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Matrix3x2(const T* val_ptr){
		data[0][0] = val_ptr[0];
		data[0][1] = val_ptr[1];
		// data[0][2] = val_ptr[2];
		// data[0][3] = val_ptr[3];
		// data[0][4] = val_ptr[4];
		// data[0][5] = val_ptr[5];
		data[1][0] = val_ptr[2];
		data[1][1] = val_ptr[3];
		data[2][0] = val_ptr[4];
		data[2][1] = val_ptr[5];
	}
	Matrix3x2(T val00 = 0, T val01 = 0, T val10 = 0, T val11 = 0, T val20 = 0, T val21 = 0){
		data[0][0] = val00;
		data[0][1] = val01;
		data[1][0] = val10;
		data[1][1] = val11;
		data[2][0] = val20;
		data[2][1] = val21;
	}
	template<typename TYPE> Matrix3x2(const Matrix3x2<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
	}
	Matrix3x2(const Com3<T> vec1, const Com3<T> vec2){
		data[0][0] = vec1.x;
		data[0][1] = vec2.x;
		data[1][0] = vec1.y;
		data[1][1] = vec2.y;
		data[2][0] = vec1.z;
		data[2][1] = vec2.z;
	}

	//	operators
	inline T* operator[] (const int id){
		assert(id == 0 || id == 1 || id == 2);	//	id must be valid, or the last value is returned
		if (!id)
			return data[0];
		else if (id == 1)
			return data[1];
		else
			return data[2];
	}

	inline const T* operator[] (const int id) const{
		assert(id == 0 || id == 1 || id == 2);	//	id must be valid, or the last value is returned
		if (!id)
			return data[0];
		else if (id == 1)
			return data[1];
		else
			return data[2];
	}

	inline T& operator()(const int x, const int y){
		return data[x][y];
	}

	inline const T& operator()(const int x, const int y) const{
		return data[x][y];
	}

	inline Matrix3x2 operator - () const{
		Matrix3x2 result(*this);
		result.data[0][0] = -data[0][0];
		result.data[0][1] = -data[0][1];
		result.data[1][0] = -data[1][0];
		result.data[1][1] = -data[1][1];
		result.data[2][0] = -data[2][0];
		result.data[2][1] = -data[2][1];
		return result;
	}

	template<typename TYPE> inline Matrix3x2& operator = (const Matrix3x2<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[2][0] = mat.data[2][0];
		data[2][1] = mat.data[2][1];
		return *this;
	}

	template<typename TYPE> inline Matrix3x2& operator += (const Matrix3x2<TYPE> mat){
		data[0][0] += mat.data[0][0];
		data[0][1] += mat.data[0][1];
		data[1][0] += mat.data[1][0];
		data[1][1] += mat.data[1][1];
		data[2][0] += mat.data[2][0];
		data[2][1] += mat.data[2][1];
		return *this;
	}

	template<typename TYPE> inline Matrix3x2& operator -= (const Matrix3x2<TYPE> mat){
		data[0][0] -= mat.data[0][0];
		data[0][1] -= mat.data[0][1];
		data[1][0] -= mat.data[1][0];
		data[1][1] -= mat.data[1][1];
		data[2][0] -= mat.data[2][0];
		data[2][1] -= mat.data[2][1];
		return *this;
	}

	template<typename TYPE> inline Matrix3x2& operator *= (const Matrix2x2<TYPE> mat){
		T new_data[3][2];
		new_data[0][0] = data[0][0] * mat.data[0][0] + data[0][1] * mat.data[1][0];
		new_data[0][1] = data[0][0] * mat.data[0][1] + data[0][1] * mat.data[1][1];
		new_data[1][0] = data[1][0] * mat.data[0][0] + data[1][1] * mat.data[1][0];
		new_data[1][1] = data[1][0] * mat.data[0][1] + data[1][1] * mat.data[1][1];
		new_data[2][0] = data[2][0] * mat.data[0][0] + data[2][1] * mat.data[1][0];
		new_data[2][1] = data[2][0] * mat.data[0][1] + data[2][1] * mat.data[1][1];
		data[0][0] = new_data[0][0];
		data[0][1] = new_data[0][1];
		data[1][0] = new_data[1][0];
		data[1][1] = new_data[1][1];
		data[2][0] = new_data[2][0];
		data[2][1] = new_data[2][1];
		return *this;
	}

	template<typename TYPE> inline Matrix3x2& operator *= (const TYPE val){
		data[0][0] *= val;
		data[0][1] *= val;
		data[1][0] *= val;
		data[1][1] *= val;
		data[2][0] *= val;
		data[2][1] *= val;
		return *this;
	}

	template<typename TYPE> inline Matrix3x2& operator /= (const TYPE val){
		data[0][0] /= val;
		data[0][1] /= val;
		data[1][0] /= val;
		data[1][1] /= val;
		data[2][0] /= val;
		data[2][1] /= val;
		return *this;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator + (const Matrix3x2<TYPE> mat) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator - (const Matrix3x2<TYPE> mat) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix2x3<TYPE> mat) const{
		Matrix3x3<TYPE_PROMOTE(T, TYPE)> result;
		result.data[0][0] = data[0][0] * mat.data[0][0] + data[0][1] * mat.data[1][0];
		result.data[0][1] = data[0][0] * mat.data[0][1] + data[0][1] * mat.data[1][1];
		result.data[0][2] = data[0][0] * mat.data[0][2] + data[0][1] * mat.data[1][2];
		result.data[1][0] = data[1][0] * mat.data[0][0] + data[1][1] * mat.data[1][0];
		result.data[1][1] = data[1][0] * mat.data[0][1] + data[1][1] * mat.data[1][1];
		result.data[1][2] = data[1][0] * mat.data[0][2] + data[1][1] * mat.data[1][2];
		result.data[2][0] = data[2][0] * mat.data[0][0] + data[2][1] * mat.data[1][0];
		result.data[2][1] = data[2][0] * mat.data[0][1] + data[2][1] * mat.data[1][1];
		result.data[2][2] = data[2][0] * mat.data[0][2] + data[2][1] * mat.data[1][2];
		return result;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix2x2<TYPE> mat) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= mat;
		return result;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator * (const Vec2<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0] * vec.x + data[0][1] * vec.y;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0] * vec.x + data[1][1] * vec.y;
		TYPE_PROMOTE(T, TYPE) new_z = data[2][0] * vec.x + data[2][1] * vec.y;
		return Vec3<TYPE_PROMOTE(T, TYPE)>(new_x, new_y, new_z);
	}

	/**
		Return transposed matrix, become 2x3 matrix
	*/
	inline Matrix2x3<T> Transpose() const{
		Matrix2x3<T> result;
		result.data[0][0] = data[0][0];
		result.data[1][0] = data[0][1];
		result.data[0][1] = data[1][0];
		result.data[1][1] = data[1][1];
		result.data[0][2] = data[2][0];
		result.data[1][2] = data[2][1];
		return result;
	}

};

/**
	2x3 Matrix

	Non-square matrix is often used in 2D-3D mapping
*/
template< class T >
class Matrix2x3{
public:
	T	data[2][3];

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Matrix2x3(const T* val_ptr){
		data[0][0] = val_ptr[0];
		data[0][1] = val_ptr[1];
		data[0][2] = val_ptr[2];
		// data[0][3] = val_ptr[3];
		// data[0][4] = val_ptr[4];
		// data[0][5] = val_ptr[5];
		data[1][0] = val_ptr[3];
		data[1][1] = val_ptr[4];
		data[1][2] = val_ptr[5];
	}
	Matrix2x3(T val00 = 0, T val01 = 0, T val02 = 0, T val10 = 0, T val11 = 0, T val12 = 0){
		data[0][0] = val00;
		data[0][1] = val01;
		data[0][2] = val02;
		data[1][0] = val10;
		data[1][1] = val11;
		data[1][2] = val12;
	}
	template<typename TYPE> Matrix2x3(const Matrix2x3<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
	}
	Matrix2x3(const Com2<T> vec1, const Com2<T> vec2, const Com2<T> vec3){
		data[0][0] = vec1.x;
		data[0][1] = vec2.x;
		data[0][2] = vec3.x;
		data[1][0] = vec1.y;
		data[1][1] = vec2.y;
		data[1][2] = vec3.y;
	}

	//	operators
	inline T* operator[] (const int id){
		assert(id == 0 || id == 1);	//	id must be valid, or the last value is returned
		if (!id)
			return data[0];
		else
			return data[1];
	}

	inline const T* operator[] (const int id) const{
		assert(id == 0 || id == 1);	//	id must be valid, or the last value is returned
		if (!id)
			return data[0];
		else
			return data[1];
	}

	inline T& operator()(const int x, const int y){
		return data[x][y];
	}

	inline const T& operator()(const int x, const int y) const{
		return data[x][y];
	}

	inline Matrix2x3 operator - () const{
		Matrix2x3 result(*this);
		result.data[0][0] = -data[0][0];
		result.data[0][1] = -data[0][1];
		result.data[0][2] = -data[0][2];
		result.data[1][0] = -data[1][0];
		result.data[1][1] = -data[1][1];
		result.data[1][2] = -data[1][2];
		return result;
	}

	template<typename TYPE> inline Matrix2x3& operator = (const Matrix2x3<TYPE> mat){
		data[0][0] = mat.data[0][0];
		data[0][1] = mat.data[0][1];
		data[0][2] = mat.data[0][2];
		data[1][0] = mat.data[1][0];
		data[1][1] = mat.data[1][1];
		data[1][2] = mat.data[1][2];
		return *this;
	}

	template<typename TYPE> inline Matrix2x3& operator += (const Matrix2x3<TYPE> mat){
		data[0][0] += mat.data[0][0];
		data[0][1] += mat.data[0][1];
		data[0][2] += mat.data[0][2];
		data[1][0] += mat.data[1][0];
		data[1][1] += mat.data[1][1];
		data[1][2] += mat.data[1][2];
		return *this;
	}

	template<typename TYPE> inline Matrix2x3& operator -= (const Matrix2x3<TYPE> mat){
		data[0][0] -= mat.data[0][0];
		data[0][1] -= mat.data[0][1];
		data[0][2] -= mat.data[0][2];
		data[1][0] -= mat.data[1][0];
		data[1][1] -= mat.data[1][1];
		data[1][2] -= mat.data[1][2];
		return *this;
	}

	template<typename TYPE> inline Matrix2x3& operator *= (const Matrix3x3<TYPE> mat){
		T new_data[2][3];
		new_data[0][0] = data[0][0] * mat.data[0][0] + data[0][1] * mat.data[1][0] + data[0][2] * mat.data[2][0];
		new_data[0][1] = data[0][0] * mat.data[0][1] + data[0][1] * mat.data[1][1] + data[0][2] * mat.data[2][1];
		new_data[0][2] = data[0][0] * mat.data[0][2] + data[0][1] * mat.data[1][2] + data[0][2] * mat.data[2][2];
		new_data[1][0] = data[1][0] * mat.data[0][0] + data[1][1] * mat.data[1][0] + data[1][2] * mat.data[2][0];
		new_data[1][1] = data[1][0] * mat.data[0][1] + data[1][1] * mat.data[1][1] + data[1][2] * mat.data[2][1];
		new_data[1][2] = data[1][0] * mat.data[0][2] + data[1][1] * mat.data[1][2] + data[1][2] * mat.data[2][2];
		data[0][0] = new_data[0][0];
		data[0][1] = new_data[0][1];
		data[0][2] = new_data[0][2];
		data[1][0] = new_data[1][0];
		data[1][1] = new_data[1][1];
		data[1][2] = new_data[1][2];
		return *this;
	}

	template<typename TYPE> inline Matrix2x3& operator *= (const TYPE val){
		data[0][0] *= val;
		data[0][1] *= val;
		data[0][2] *= val;
		data[1][0] *= val;
		data[1][1] *= val;
		data[1][2] *= val;
		return *this;
	}

	template<typename TYPE> inline Matrix2x3& operator /= (const TYPE val){
		data[0][0] /= val;
		data[0][1] /= val;
		data[0][2] /= val;
		data[1][0] /= val;
		data[1][1] /= val;
		data[1][2] /= val;
		return *this;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator + (const Matrix2x3<TYPE> mat) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator - (const Matrix2x3<TYPE> mat) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix3x2<TYPE> mat) const{
		Matrix2x2<TYPE_PROMOTE(T, TYPE)> result;
		result.data[0][0] = data[0][0] * mat.data[0][0] + data[0][1] * mat.data[1][0] + data[0][2] * mat.data[2][0];
		result.data[0][1] = data[0][0] * mat.data[0][1] + data[0][1] * mat.data[1][1] + data[0][2] * mat.data[2][1];
		result.data[1][0] = data[1][0] * mat.data[0][0] + data[1][1] * mat.data[1][0] + data[1][2] * mat.data[2][0];
		result.data[1][1] = data[1][0] * mat.data[0][1] + data[1][1] * mat.data[1][1] + data[1][2] * mat.data[2][1];
		return result;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator * (const Matrix3x3<TYPE> mat) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= mat;
		return result;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator * (const Vec3<TYPE> vec) const{
		TYPE_PROMOTE(T, TYPE) new_x = data[0][0] * vec.x + data[0][1] * vec.y + data[0][2] * vec.z;
		TYPE_PROMOTE(T, TYPE) new_y = data[1][0] * vec.x + data[1][1] * vec.y + data[1][2] * vec.z;
		return Vec2<TYPE_PROMOTE(T, TYPE)>(new_x, new_y);
	}

	/**
		Return transposed matrix, become 3x2 matrix
	*/
	inline Matrix3x2<T> Transpose() const{
		Matrix3x2<T> result;
		result.data[0][0] = data[0][0];
		result.data[1][0] = data[0][1];
		result.data[2][0] = data[0][2];
		result.data[0][1] = data[1][0];
		result.data[1][1] = data[1][1];
		result.data[2][1] = data[1][2];
		return result;
	}
};


//	========================================
///@{
/**	@name Non-member functions of Matrix
*/
//	========================================
template<typename TYPE, class T> inline Matrix2x2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Matrix2x2<T> mat){
	Matrix2x2<TYPE_PROMOTE(T, TYPE)> result(mat);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Matrix3x3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Matrix3x3<T> mat){
	Matrix3x3<TYPE_PROMOTE(T, TYPE)> result(mat);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Matrix4x4<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Matrix4x4<T> mat){
	Matrix4x4<TYPE_PROMOTE(T, TYPE)> result(mat);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Matrix3x2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Matrix3x2<T> mat){
	Matrix3x2<TYPE_PROMOTE(T, TYPE)> result(mat);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Matrix2x3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Matrix2x3<T> mat){
	Matrix2x3<TYPE_PROMOTE(T, TYPE)> result(mat);
	result *= val;
	return result;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Matrix2x2<T> mat){
	for( int j=0; j<2; j++ ){
		for( int i=0; i<2; i++ )
			stream << mat.data[j][i] << '\t';
		stream << std::endl;
	}
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Matrix3x3<T> mat){
	for( int j=0; j<3; j++ ){
		for( int i=0; i<3; i++ )
			stream << mat.data[j][i] << '\t';
		stream << std::endl;
	}
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Matrix4x4<T> mat){
	for( int j=0; j<4; j++ ){
		for( int i=0; i<4; i++ )
			stream << mat.data[j][i] << '\t';
		stream << std::endl;
	}
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Matrix3x2<T> mat){
	for (int j = 0; j<3; j++){
		for (int i = 0; i<2; i++)
			stream << mat.data[j][i] << '\t';
		stream << std::endl;
	}
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Matrix2x3<T> mat){
	for (int j = 0; j<2; j++){
		for (int i = 0; i<3; i++)
			stream << mat.data[j][i] << '\t';
		stream << std::endl;
	}
	return stream;
}


///@}

//	========================================
//	type rename
//	========================================
typedef Matrix2x2<char>					Matrix2x2c;
typedef Matrix2x2<unsigned char>		Matrix2x2uc;
typedef Matrix2x2<short>				Matrix2x2s;
typedef Matrix2x2<unsigned short>		Matrix2x2us;
typedef Matrix2x2<int>					Matrix2x2i;
typedef Matrix2x2<unsigned int>			Matrix2x2ui;
typedef Matrix2x2<long>					Matrix2x2l;
typedef Matrix2x2<unsigned long>		Matrix2x2ul;
typedef Matrix2x2<long long>			Matrix2x2ll;
typedef Matrix2x2<unsigned long long>	Matrix2x2ull;
typedef Matrix2x2<float>				Matrix2x2f;
typedef Matrix2x2<double>				Matrix2x2d;
typedef Matrix2x2<long double>			Matrix2x2ld;

typedef Matrix3x3<char>					Matrix3x3c;
typedef Matrix3x3<unsigned char>		Matrix3x3uc;
typedef Matrix3x3<short>				Matrix3x3s;
typedef Matrix3x3<unsigned short>		Matrix3x3us;
typedef Matrix3x3<int>					Matrix3x3i;
typedef Matrix3x3<unsigned int>			Matrix3x3ui;
typedef Matrix3x3<long>					Matrix3x3l;
typedef Matrix3x3<unsigned long>		Matrix3x3ul;
typedef Matrix3x3<long long>			Matrix3x3ll;
typedef Matrix3x3<unsigned long long>	Matrix3x3ull;
typedef Matrix3x3<float>				Matrix3x3f;
typedef Matrix3x3<double>				Matrix3x3d;
typedef Matrix3x3<long double>			Matrix3x3ld;

typedef Matrix4x4<char>					Matrix4x4c;
typedef Matrix4x4<unsigned char>		Matrix4x4uc;
typedef Matrix4x4<short>				Matrix4x4s;
typedef Matrix4x4<unsigned short>		Matrix4x4us;
typedef Matrix4x4<int>					Matrix4x4i;
typedef Matrix4x4<unsigned int>			Matrix4x4ui;
typedef Matrix4x4<long>					Matrix4x4l;
typedef Matrix4x4<unsigned long>		Matrix4x4ul;
typedef Matrix4x4<long long>			Matrix4x4ll;
typedef Matrix4x4<unsigned long long>	Matrix4x4ull;
typedef Matrix4x4<float>				Matrix4x4f;
typedef Matrix4x4<double>				Matrix4x4d;
typedef Matrix4x4<long double>			Matrix4x4ld;

typedef Matrix3x2<char>					Matrix3x2c;
typedef Matrix3x2<unsigned char>		Matrix3x2uc;
typedef Matrix3x2<short>				Matrix3x2s;
typedef Matrix3x2<unsigned short>		Matrix3x2us;
typedef Matrix3x2<int>					Matrix3x2i;
typedef Matrix3x2<unsigned int>			Matrix3x2ui;
typedef Matrix3x2<long>					Matrix3x2l;
typedef Matrix3x2<unsigned long>		Matrix3x2ul;
typedef Matrix3x2<long long>			Matrix3x2ll;
typedef Matrix3x2<unsigned long long>	Matrix3x2ull;
typedef Matrix3x2<float>				Matrix3x2f;
typedef Matrix3x2<double>				Matrix3x2d;
typedef Matrix3x2<long double>			Matrix3x2ld;

typedef Matrix2x3<char>					Matrix2x3c;
typedef Matrix2x3<unsigned char>		Matrix2x3uc;
typedef Matrix2x3<short>				Matrix2x3s;
typedef Matrix2x3<unsigned short>		Matrix2x3us;
typedef Matrix2x3<int>					Matrix2x3i;
typedef Matrix2x3<unsigned int>			Matrix2x3ui;
typedef Matrix2x3<long>					Matrix2x3l;
typedef Matrix2x3<unsigned long>		Matrix2x3ul;
typedef Matrix2x3<long long>			Matrix2x3ll;
typedef Matrix2x3<unsigned long long>	Matrix2x3ull;
typedef Matrix2x3<float>				Matrix2x3f;
typedef Matrix2x3<double>				Matrix2x3d;
typedef Matrix2x3<long double>			Matrix2x3ld;

}	//	end namespace yz

#endif	//	__YZ_MATRIX_H__