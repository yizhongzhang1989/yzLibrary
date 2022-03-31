/***********************************************************/
/**	\file
	\brief		Quaternion
	\details	This file contain definition of quaternion 
				and corresponding utility functions, such as
				interpolation, resampling and so on
	\author		Yizhong Zhang
	\date		10/11/2012
*/
/***********************************************************/
#ifndef __YZ_QUATERNION_H__
#define __YZ_QUATERNION_H__

#include <assert.h>
#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_matrix.h"

namespace yz{

/**
	Quaternion
*/
template<class T>
class Quaternion{
public:
	T	w, x, y, z;

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Quaternion(const T w_val=1, const T x_val=0, const T y_val=0, const T z_val=0): w(w_val), x(x_val), y(y_val), z(z_val){}
	Quaternion(const T *val_ptr) : w(val_ptr[0]), x(val_ptr[1]), y(val_ptr[2]), z(val_ptr[3]){}
	template<typename TYPE> Quaternion(const Quaternion<TYPE> rhs){
		w = rhs.w;
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}
	template<typename TYPE> Quaternion(const Com4<TYPE> vec4)
		: w(vec4.x), x(vec4.y), y(vec4.z), z(vec4.w){
	}

	template<typename TYPE> Quaternion(const Matrix3x3<TYPE> mat){
		this->Set(mat);
	}

	//	operators
	inline T& operator[] (const int id){
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return w;
		else if( id==1 )
			return x;
		else if( id==2 )
			return y;
		else
			return z;
	}

	inline const T& operator[] (const int id) const{
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return w;
		else if( id==1 )
			return x;
		else if( id==2 )
			return y;
		else
			return z;
	}

	inline Quaternion operator - () const{
		Quaternion result(-x, -y, -z, -w);
		return result;
	}

	template<typename TYPE> inline Quaternion& operator = (const Quaternion<TYPE> rhs){
		w = rhs.w;
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator = (const Com4<TYPE> rhs){
		//	be cautious about the sequence
		w = rhs.x;
		x = rhs.y;
		y = rhs.z;
		z = rhs.w;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator += (const Quaternion<TYPE> rhs){
		w += rhs.w;
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator -= (const Quaternion<TYPE> rhs){
		w -= rhs.w;
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator *= (const Quaternion<TYPE> rhs){
		T new_w = w*rhs.w - x*rhs.x - y*rhs.y - z*rhs.z;
		T new_x = w*rhs.x + x*rhs.w + y*rhs.z - z*rhs.y;
		T new_y = w*rhs.y - x*rhs.z + y*rhs.w + z*rhs.x;
		T new_z = w*rhs.z + x*rhs.y - y*rhs.x + z*rhs.w;
		w = new_w;
		x = new_x;
		y = new_y;
		z = new_z;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator *= (const TYPE scale){
		w *= scale;
		x *= scale;
		y *= scale;
		z *= scale;
		return *this;
	}

	template<typename TYPE> inline Quaternion& operator /= (const TYPE scale){
		w /= scale;
		x /= scale;
		y /= scale;
		z /= scale;
		return *this;
	}

	template<typename TYPE> inline Quaternion<TYPE_PROMOTE(T, TYPE)> operator + (const Quaternion<TYPE> rhs) const{
		Quaternion<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += rhs;
		return result;
	}
	template<typename TYPE> inline Quaternion<TYPE_PROMOTE(T, TYPE)> operator - (const Quaternion<TYPE> rhs) const{
		Quaternion<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= rhs;
		return result;
	}
	template<typename TYPE> inline Quaternion<TYPE_PROMOTE(T, TYPE)> operator * (const Quaternion<TYPE> rhs) const{
		Quaternion<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= rhs;
		return result;
	}

	template<typename TYPE> inline Quaternion<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE scale) const{
		Quaternion<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= scale;
		return result;
	}
	template<typename TYPE> inline Quaternion<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE scale) const{
		Quaternion<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= scale;
		return result;
	}
	/**
		get the negative quaternion

		negative and original quaternion represent the same rotation
	*/
	inline Quaternion Negative() const{
		Quaternion neg(-w, -x, -y, -z);
		return neg;
	}

	/**
		get the quaternion to negative

		negative and original quaternion represent the same rotation
	*/
	inline Quaternion& SetNegative(){
		w = -w;
		x = -x;
		y = -y;
		z = -z;
		return *this;
	}

	/**
		get the conjugate quaternion

		conjugate quaternion represents the inverse rotation
	*/
	inline Quaternion Conjugate() const{
		Quaternion conj(w, -x, -y, -z);
		return conj;
	}

	/**
		Set the quaternion to conjugate

		conjugate quaternion represents the inverse rotation
	*/
	inline Quaternion& SetConjugate(){
		x = -x;
		y = -y;
		z = -z;
		return *this;
	}

	/**
		get normalized quaternion
	*/
	inline Quaternion<PROMOTE_T_TO_FLOAT> Normalize() const{
		Quaternion<TYPE_PROMOTE(T, float)> result(*this);
		result.SetNormalize();
		return result;
	}

	/**
		Normalize the quaternion
	*/
	inline Quaternion& SetNormalize(){
		TYPE_PROMOTE(T, float) inv_length = 1.0 / sqrt(w*w + x*x + y*y + z*z);
		w *= inv_length;
		x *= inv_length;
		y *= inv_length;
		z *= inv_length;
		return *this;
	}

	/**
		Get the square norm of quaternion
	*/
	inline YZ_REAL SquareNorm() const{
		return w*w + x*x + y*y + z*z;
	}

	/**
		Get the norm of quaternion
	*/
	inline YZ_REAL Norm() const{
		return sqrt( SquareNorm() );
	}

	/**
		negative quaternion and original quaternion represent
		the same rotation matrix, we set the quaternion to 
		the form that is numerically closer to the target

		\param	q	the target quaternion
	*/
	template<typename TYPE> 
	inline Quaternion& SetCloser(Quaternion<TYPE> q){
		YZ_REAL cur_dist = (q-(*this)).SquareNorm();
		YZ_REAL neg_dist = (q+(*this)).SquareNorm();
		if( neg_dist < cur_dist )
			SetNegative();
		return (*this);
	}

	/**
		Set quaternion from rotation matrix

		In debug mode, we check whether the matrix is a rotation 
		matrix by checking determinant equal to 1. In release mode,
		we don't do this check, you must make sure that the matrix
		is legal, or the result will be incorrect

		\author	Chen Cao
	*/
	template<typename TYPE> 
	inline Quaternion& Set(const Matrix3x3<TYPE> mat){
		//	determinant of rotation matrix must be 1
		assert( fabs(mat.Det() -1 ) < 1e-3 );

		TYPE trace = mat(0, 0) + mat(1, 1) + mat(2, 2);
		if( trace > 0 ){
			T s = 0.5 / sqrt( trace + 1 );
			w = 0.25 / s;
			x = ( mat(2,1) - mat(1,2) ) * s;
			y = ( mat(0,2) - mat(2,0) ) * s;
			z = ( mat(1,0) - mat(0,1) ) * s;
		}
		else{
			if ( mat(0,0) > mat(1,1) && mat(0,0) > mat(2,2) ) {
				T s = 2 * sqrt( 1 + mat(0,0) - mat(1,1) - mat(2,2) );
				w = ( mat(2,1) - mat(1,2) ) / s;
				x = 0.25 * s;
				y = ( mat(0,1) + mat(1,0) ) / s;
				z = ( mat(0,2) + mat(2,0) ) / s;
			}
			else if ( mat(1,1) > mat(2,2)) {
				T s = 2 * sqrt( 1 + mat(1,1) - mat(0,0) - mat(2,2) );
				w = ( mat(0,2) - mat(2,0) ) / s;
				x = ( mat(0,1) + mat(1,0) ) / s;
				y = 0.25 * s;
				z = ( mat(1,2) + mat(2,1) ) / s;
			} 
			else {
				T s = 2 * sqrtf( 1 + mat(2,2) - mat(0,0) - mat(1,1) );
				w = ( mat(1,0) - mat(0,1) ) / s;
				x = ( mat(0,2) + mat(2,0) ) / s;
				y = ( mat(1,2) + mat(2,1) ) / s;
				z = 0.25 * s;
			}
		}

		SetNormalize();

		return *this;
	}

	/**
		Set the quaternion by axis and angle, counter-clockwise

		\param	axis		rotation axis
		\param	angle_rad	rotation angle in rad
		\return				*this
	*/
	inline Quaternion& SetRad(Vec3<YZ_REAL> axis, TYPE_PROMOTE(T, float) angle_rad){
		TYPE_PROMOTE(T, float) cos_half = cos(angle_rad*0.5);
		TYPE_PROMOTE(T, float) sin_half = sin(angle_rad*0.5);

		axis.SetNormalize();
		w = cos_half;
		x = axis.x * sin_half;
		y = axis.y * sin_half;
		z = axis.z * sin_half;

		//	the quaternion is normalized
		return *this;
	}

	/**
		Set the quaternion by axis and angle, counter-clockwise

		\param	axis		rotation axis
		\param	angle_deg	rotation angle in degree
		\return				*this
	*/
	inline Quaternion& SetDeg(Vec3<YZ_REAL> axis, TYPE_PROMOTE(T, float) angle_deg){
		return SetRad( axis, angle_deg*YZ_PI/180 );
	}

	/**
		return the rotation axis
	*/
	inline Vec3<YZ_REAL> Axis() const{
		Vec3<YZ_REAL> axis(x, y, z);
		return axis.Normalize();
	}

	/**
		return the rotation angle in radius
	*/
	inline YZ_REAL AngleRad() const{
		return acos(w) * 2;
	}

	/**
		return the rotation angle in degree
	*/
	inline YZ_REAL AngleDeg() const{
		return AngleRad() * 180 / YZ_REAL(YZ_PI);
	}

	/**
	*/

};

//	========================================
///@{
/**	@name Print Data Container
*/
//	========================================
template<class T> inline std::ostream& operator << (std::ostream& stream, const Quaternion<T> rhs){
	stream << "("	<< TYPE_PROMOTE(T, int)(rhs.w) << ", " 
					<< TYPE_PROMOTE(T, int)(rhs.x) << ", " 
					<< TYPE_PROMOTE(T, int)(rhs.y) << ", " 
					<< TYPE_PROMOTE(T, int)(rhs.z) << ")";
	return stream;
}
///@}

//	========================================
///@{
/**	@name Non-menber Quaternion Functions
*/
//	========================================
template<class T1, class T2> inline PROMOTE_T1_T2 dot(const Quaternion<T1> a, const Quaternion<T2> b){
	return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

///@}

//	========================================
///@{
/**	@name 1D Interpolation
*/
//	========================================

/**
	Linear Euler Interpolation of quatations

	\param	q0		quaternion 0
	\param	q1		quaternion 1
	\param	t		interpolation position, should be 0 ~ 1
	\return			interplated quaternion (normalized)
*/
template <typename T1, typename T2>
inline Quaternion<PROMOTE_T1_T2_TO_FLOAT> interpLerp(const Quaternion<T1> q0, const Quaternion<T2> q1, double t){
	Quaternion<PROMOTE_T1_T2_TO_FLOAT> q = q0 * (1.0 - t) + q1 * t;
	return q;
}

/**
	Slerp Interpolation of quatations

	this function come from: \n
	http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/

	\param	qa		quaternion 0
	\param	qb		quaternion 1
	\param	t		interpolation position, should be 0 ~ 1
	\return			interplated quaternion (normalized)
*/
template <typename T1, typename T2>
inline Quaternion<PROMOTE_T1_T2_TO_FLOAT> interpSlerp(const Quaternion<T1> qa, const Quaternion<T2> qb, double t){
	Quaternion<PROMOTE_T1_T2_TO_FLOAT> qm;
	// Calculate angle between them.
	PROMOTE_T1_T2_TO_FLOAT cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
	// if qa=qb or qa=-qb then theta = 0 and we can return qa
	if (abs(cosHalfTheta) >= 1.0 - 1e-6){
		qm.w = qa.w;qm.x = qa.x;qm.y = qa.y;qm.z = qa.z;
		return qm;
	}

	// Calculate temporary values.
	PROMOTE_T1_T2_TO_FLOAT halfTheta = acos(cosHalfTheta);
	PROMOTE_T1_T2_TO_FLOAT sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);

	// if theta = 180 degrees then result is not fully defined
	// we could rotate around any axis normal to qa or qb
	if (fabs(sinHalfTheta) < 1e-5){
		qm.w = (qa.w * 0.5 + qb.w * 0.5);
		qm.x = (qa.x * 0.5 + qb.x * 0.5);
		qm.y = (qa.y * 0.5 + qb.y * 0.5);
		qm.z = (qa.z * 0.5 + qb.z * 0.5);
		return qm;
	}

	PROMOTE_T1_T2_TO_FLOAT ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
	PROMOTE_T1_T2_TO_FLOAT ratioB = sin(t * halfTheta) / sinHalfTheta; 
	//calculate Quaternion.
	qm.w = (qa.w * ratioA + qb.w * ratioB);
	qm.x = (qa.x * ratioA + qb.x * ratioB);
	qm.y = (qa.y * ratioA + qb.y * ratioB);
	qm.z = (qa.z * ratioA + qb.z * ratioB);

	return qm;
}

///@}

//	========================================
///@{
/**	@name 1D Resample
*/
//	========================================

/**
	Resampling quaternions using Linear Euler Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleLerp(Quaternion<T1>* des, int des_size, const Quaternion<T2>* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		des[i] = interpLerp(src[i_coef], src[i_coef+1], coef);
	}
	des[des_size-1] = src[src_size-1];
}



/**
	Resampling quaternions using Spherical Linear Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleSlerp(Quaternion<T1>* des, int des_size, const Quaternion<T2>* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		des[i] = interpSlerp(src[i_coef], src[i_coef+1], coef);
	}
	des[des_size-1] = src[src_size-1];
}



///@}

typedef Quaternion<float>	Quaternionf;
typedef Quaternion<double>	Quaterniond;

}	//	end namespace yz

#endif	//	__YZ_QUATERNION_H__