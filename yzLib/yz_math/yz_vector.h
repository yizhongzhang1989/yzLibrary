/***********************************************************/
/**	\file
	\brief		Vector Type for Graphics
	\details	just implemented Vec2, Vec3, Vec4, functions of 
				them are slightly different from vector in mathematics
	\author		Yizhong Zhang
	\date		5/14/2012
*/
/***********************************************************/
#ifndef __YZ_VECTOR_H__
#define __YZ_VECTOR_H__

#include <assert.h>
#include <iostream>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

//	========================================
//	Combination of data: Com2, Com3, Com4
//	========================================
/**
	Two Data Container
*/
template<class T>
class Com2{
public:
	T	x, y;

	static const int dims;
	typedef T Type;
	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Com2(const T x_val = 0, const T y_val = 0) :x(x_val), y(y_val) {}
	Com2(const T *val_ptr) :x(*val_ptr), y(*(val_ptr + 1)) {}
	template<typename TYPE> Com2(const Com2<TYPE> c):x(c.x), y(c.y) {}
	inline T& operator[] (const int id){
		assert(id==0 || id==1);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else
			return y;
	}

	inline const T& operator[] (const int id) const{
		assert(id==0 || id==1);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else
			return y;
	}

	inline bool operator == (const Com2<T> vec) const{
		return vec.x == x && vec.y == y;
	}

	inline bool operator != (const Com2<T> vec) const {
		return vec.x != x || vec.y != y;
	}

	inline bool operator > (const Com2<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? y > vec.y : false);
	}

	inline bool operator < (const Com2<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? y < vec.y : false);
	}

	inline bool operator >= (const Com2<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? y >= vec.y : false);
	}

	inline bool operator <= (const Com2<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? y <= vec.y : false);
	}
};

/**
	Three Data Container
*/
template<class T>
class Com3{
public:
	T	x, y, z;

	static const int dims;
	typedef T Type;
	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Com3(const T x_val = 0, const T y_val = 0, const T z_val = 0) :x(x_val), y(y_val), z(z_val) {}
	Com3(const T *val_ptr) :x(*val_ptr), y(*(val_ptr + 1)), z(*(val_ptr + 2)) {}
	template<typename TYPE> Com3(const Com3<TYPE> c) : x(c.x), y(c.y), z(c.z) {}
	template<typename TYPE> Com3(const Com2<TYPE> c, const T z_val = 0) : x(c.x), y(c.y), z(z_val) {}
	inline T& operator[] (const int id){
		assert(id==0 || id==1 || id==2);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else if( id==1 )
			return y;
		else
			return z;
	}

	inline const T& operator[] (const int id) const{
		assert(id==0 || id==1 || id==2);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else if( id==1 )
			return y;
		else
			return z;
	}

	inline bool operator == (const Com3<T> vec) const{
		return vec.x == x && vec.y == y && vec.z == z;
	}

	inline bool operator != (const Com3<T> vec) const {
		return  vec.x != x || vec.y != y || vec.z != z;
	}

	inline bool operator > (const Com3<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? Com2<T>(y, z) > Com2<T>(vec.y, vec.z) : false);
	}

	inline bool operator < (const Com3<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? Com2<T>(y, z) < Com2<T>(vec.y, vec.z) : false);
	}

	inline bool operator >= (const Com3<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? Com2<T>(y, z) >= Com2<T>(vec.y, vec.z) : false);
	}

	inline bool operator <= (const Com3<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? Com2<T>(y, z) <= Com2<T>(vec.y, vec.z) : false);
	}
};

/**
	Four Data Container
*/
template<class T>
class Com4{
public:
	T	x, y, z, w;

	static const int dims;
	typedef T Type;
	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	Com4(const T x_val = 0, const T y_val = 0, const T z_val = 0, const T w_val = 0) :x(x_val), y(y_val), z(z_val), w(w_val) {}
	Com4(const T *val_ptr) :x(*val_ptr), y(*(val_ptr + 1)), z(*(val_ptr + 2)), w(*(val_ptr + 3)) {}
	template<typename TYPE> Com4(const Com4<TYPE> c) : x(c.x), y(c.y), z(c.z), w(c.w) {}
	template<typename TYPE> Com4(const Com3<TYPE> c, const T w_val = 0) : x(c.x), y(c.y), z(c.z), w(w_val) {}
	inline T& operator[] (const int id){
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else if( id==1 )
			return y;
		else if( id==2 )
			return z;
		else
			return w;
	}

	inline const T& operator[] (const int id) const{
		assert(id==0 || id==1 || id==2 || id==3);	//	id must be valid, or the last value is returned
		if( !id )
			return x;
		else if( id==1 )
			return y;
		else if( id==2 )
			return z;
		else
			return w;
	}

	inline bool operator == (const Com4<T> vec) const{
		return vec.x == x && vec.y == y && vec.z == z && vec.w == w;
	}

	inline bool operator != (const Com4<T> vec) const {
		return  vec.x != x || vec.y != y || vec.z != z || vec.w != w;
	}

	inline bool operator > (const Com4<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? Com3<T>(y, z, w) > Com3<T>(vec.y, vec.z, vec.w) : false);
	}

	inline bool operator < (const Com4<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? Com3<T>(y, z, w) < Com3<T>(vec.y, vec.z, vec.w) : false);
	}

	inline bool operator >= (const Com4<T> vec) const {
		return x > vec.x ? true : (x == vec.x ? Com3<T>(y, z, w) >= Com3<T>(vec.y, vec.z, vec.w) : false);
	}

	inline bool operator <= (const Com4<T> vec) const {
		return x < vec.x ? true : (x == vec.x ? Com3<T>(y, z, w) <= Com3<T>(vec.y, vec.z, vec.w) : false);
	}
};

//	========================================
///@{
/**	@name Print Data Container
*/
//	========================================
template<class T> inline std::ostream& operator << (std::ostream& stream, const Com2<T> com){
	stream << "("	<< TYPE_PROMOTE(T, int)(com.x) << ", " 
					<< TYPE_PROMOTE(T, int)(com.y) << ")";
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Com3<T> com){
	stream << "("	<< TYPE_PROMOTE(T, int)(com.x) << ", " 
					<< TYPE_PROMOTE(T, int)(com.y) << ", " 
					<< TYPE_PROMOTE(T, int)(com.z) << ")";
	return stream;
}

template<class T> inline std::ostream& operator << (std::ostream& stream, const Com4<T> com){
	stream << "("	<< TYPE_PROMOTE(T, int)(com.x) << ", " 
					<< TYPE_PROMOTE(T, int)(com.y) << ", " 
					<< TYPE_PROMOTE(T, int)(com.z) << ", " 
					<< TYPE_PROMOTE(T, int)(com.w) << ")";
	return stream;
}
///@}


//	========================================
//	Vector: Vec2, Vec3, Vec4
//	========================================
/**
	Vector 2D
*/
template<class T>
class Vec2 : public Com2<T>{
public:
	using Com2<T>::x;
	using Com2<T>::y;
	using Com2<T>::dims;
	using typename Com2<T>::Type;
	using typename Com2<T>::YZ_REAL;
	
	//	constructors
	Vec2(const T x_val=0, const T y_val=0):Com2<T>(x_val, y_val){}
	Vec2(const T *val_ptr):Com2<T>(val_ptr){}
	template<typename TYPE> Vec2(const Com2<TYPE> v):Com2<T>(v){}

	//	operators
	inline Vec2 operator - () const{
		Vec2 result(-x, -y);
		return result;
	}

	template<typename TYPE> inline Vec2& operator = (const Vec2<TYPE> vec){
		x = vec.x;
		y = vec.y;
		return *this;
	}

	template<typename TYPE> inline Vec2& operator += (const Vec2<TYPE> vec){
		x += vec.x;
		y += vec.y;
		return *this;
	}

	template<typename TYPE> inline Vec2& operator -= (const Vec2<TYPE> vec){
		x -= vec.x;
		y -= vec.y;
		return *this;
	}

	template<typename TYPE> inline Vec2& operator *= (const TYPE val){
		x *= val;
		y *= val;
		return *this;
	}

	template<typename TYPE> inline Vec2& operator /= (const TYPE val){
		x /= val;
		y /= val;
		return *this;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator + (const Vec2<TYPE> vec) const{
		Vec2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += vec;
		return result;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator - (const Vec2<TYPE> vec) const{
		Vec2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= vec;
		return result;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Vec2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Vec2<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	/**
		Return normalized vector, but don't change old data
	*/
	inline Vec2<YZ_REAL> Normalize() const{
		Vec2<YZ_REAL> result(*this);
		result.SetNormalize();
		return result;
	}

	/**
		Return rotated vector angle_rad(rad) counter-clockwise to original vector, but don't change old data
	*/
	inline Vec2<YZ_REAL> RotateRad(const YZ_REAL angle_rad) const{
		Vec2<YZ_REAL> result(*this);
		result.SetRotateRad(angle_rad);
		return result;
	}

	/**
		Return rotated vector angle_deg(degrees) counter-clockwise to original vector, but don't change old data
	*/
	inline Vec2<YZ_REAL> RotateDeg(const YZ_REAL angle_deg) const{
		Vec2<YZ_REAL> result(*this);
		result.SetRotateDeg(angle_deg);
		return result;
	}

	//	set vectors
	/**
		Normalize the vector
	*/
	inline Vec2& SetNormalize(){
		YZ_REAL length = Length();
		(*this) /= length;
		return *this;
	}

	/**
		Rotate the vector angle_rad(rad) counter-clockwise
	*/
	inline Vec2& SetRotateRad(const YZ_REAL angle_rad){
		YZ_REAL cos_a = cos(angle_rad);
		YZ_REAL sin_a = sin(angle_rad);
		YZ_REAL x_new = cos_a * x - sin_a * y;
		YZ_REAL y_new = sin_a * x + cos_a * y;
		x = x_new;
		y = y_new;
		return *this;
	}

	/**
		Rotate the vector angle_deg(degrees) counter-clockwise
	*/
	inline Vec2& SetRotateDeg(const YZ_REAL angle_deg){
		return SetRotateRad( angle_deg*YZ_PI/180 );
	}

	//	vector properties
	inline YZ_REAL SquareLength() const{
		return x*x + y*y;
	}

	inline YZ_REAL Length() const{
		return sqrt(SquareLength());
	}

};

/**
	Vector 3D
*/
template<class T>
class Vec3 : public Com3<T>{
public:
	using Com3<T>::x;
	using Com3<T>::y;
	using Com3<T>::z;
	using Com3<T>::dims;
	using typename Com3<T>::Type;
	using typename Com3<T>::YZ_REAL;

	//	constructors
	Vec3(const T x_val=0, const T y_val=0, const T z_val=0):Com3<T>(x_val, y_val, z_val){}
	Vec3(const T *val_ptr):Com3<T>(val_ptr){}

	template<typename TYPE> Vec3(const Com4<TYPE> v):Com3<T>(v.x, v.y, v.z){}
	template<typename TYPE> Vec3(const Com3<TYPE> v):Com3<T>(v){}
	template<typename TYPE> Vec3(const Com2<TYPE> v, const T z_val=0):Com3<T>(v, z_val){}

	//	operators
	inline Vec3 operator - () const{
		Vec3 result(-x, -y, -z);
		return result;
	}

	template<typename TYPE> inline Vec3& operator = (const Vec3<TYPE> vec){
		x = vec.x;
		y = vec.y;
		z = vec.z;
		return *this;
	}

	template<typename TYPE> inline Vec3& operator += (const Vec3<TYPE> vec){
		x += vec.x;
		y += vec.y;
		z += vec.z;
		return *this;
	}

	template<typename TYPE> inline Vec3& operator -= (const Vec3<TYPE> vec){
		x -= vec.x;
		y -= vec.y;
		z -= vec.z;
		return *this;
	}

	template<typename TYPE> inline Vec3& operator *= (const TYPE val){
		x *= val;
		y *= val;
		z *= val;
		return *this;
	}

	template<typename TYPE> inline Vec3& operator /= (const TYPE val){
		x /= val;
		y /= val;
		z /= val;
		return *this;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator + (const Vec3<TYPE> vec) const{
		Vec3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += vec;
		return result;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator - (const Vec3<TYPE> vec) const{
		Vec3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= vec;
		return result;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Vec3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Vec3<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	/**
		Return normalized vector, but don't change old data
	*/
	inline Vec3<YZ_REAL> Normalize() const{
		Vec3<YZ_REAL> result(*this);
		result.SetNormalize();
		return result;
	}

	/**
		Return rotated vector angle_rad(rad) around axis
		counter-clockwise to original vector, but don't change old data
	*/
	inline Vec3<YZ_REAL> RotateRad(const Vec3<YZ_REAL> axis, const YZ_REAL angle_rad) const{
		Vec3<YZ_REAL> result(*this);
		result.SetRotateRad(axis, angle_rad);
		return result;
	}

	/**
		Return rotated vector angle_deg(degrees) around axis
		counter-clockwise to original vector, but don't change old data
	*/
	inline Vec3<YZ_REAL> RotateDeg(const Vec3<YZ_REAL> axis, const YZ_REAL angle_deg) const{
		Vec3<YZ_REAL> result(*this);
		result.SetRotateDeg(axis, angle_deg);
		return result;
	}

	/**
		rotate the vector by a quaternion

		p' = q * p * q'
	*/
	inline Vec3<YZ_REAL> RotateQuaternion(const YZ_REAL qw, const YZ_REAL qx, const YZ_REAL qy, const YZ_REAL qz){
		Vec3<YZ_REAL> result(*this);
		result.SetRotateQuaternion(qw, qx, qy, qz);
		return result;
	}
	//	set vectors
	/**
		Normalize the vector
	*/
	inline Vec3& SetNormalize(){
		YZ_REAL length = Length();
		(*this) /= length;
		return *this;
	}

	/**
		Return rotated vector angle_rad(rad) around axis counter-clockwise to original vector
	*/
	inline Vec3& SetRotateRad(Vec3<YZ_REAL> axis, const YZ_REAL angle_rad){
		axis.SetNormalize();

		YZ_REAL cos_a = cos(angle_rad);
		YZ_REAL sin_a = sin(angle_rad);

		YZ_REAL a00 = cos_a + axis.x*axis.x*(1-cos_a);
		YZ_REAL a01 = axis.x*axis.y*(1-cos_a) - axis.z*sin_a;
		YZ_REAL a02 = axis.x*axis.z*(1-cos_a) + axis.y*sin_a;
		YZ_REAL a10 = axis.y*axis.x*(1-cos_a) + axis.z*sin_a;
		YZ_REAL a11 = cos_a + axis.y*axis.y*(1-cos_a);
		YZ_REAL a12 = axis.y*axis.z*(1-cos_a) - axis.x*sin_a;
		YZ_REAL a20 = axis.x*axis.z*(1-cos_a) - axis.y*sin_a;
		YZ_REAL a21 = axis.y*axis.z*(1-cos_a) + axis.x*sin_a;
		YZ_REAL a22 = cos_a + axis.z*axis.z*(1-cos_a);

		YZ_REAL x_new = a00 * x + a01 * y + a02 * z;
		YZ_REAL y_new = a10 * x + a11 * y + a12 * z;
		YZ_REAL z_new = a20 * x + a21 * y + a22 * z;

		x = x_new;
		y = y_new;
		z = z_new;

		return *this;
	}

	/**
		Return rotated vector angle_deg(degrees) around axis counter-clockwise to original vector
	*/	
	inline Vec3& SetRotateDeg(const Vec3<YZ_REAL> axis, const YZ_REAL angle_deg){
		return SetRotateRad(axis, angle_deg*YZ_PI/180);
	}

	/**
		set the vector to be rotated by a quaternion

		p' = q * p * q'
	*/
	inline Vec3& SetRotateQuaternion(const YZ_REAL qw, const YZ_REAL qx, const YZ_REAL qy, const YZ_REAL qz){
		YZ_REAL new_w = -qx*x - qy*y - qz*z;
		YZ_REAL new_x = qw*x + qy*z - qz*y;
		YZ_REAL new_y = qw*y - qx*z + qz*x;
		YZ_REAL new_z = qw*z + qx*y - qy*x;

		x = -new_w*qx + new_x*qw - new_y*qz + new_z*qy;
		y = -new_w*qy + new_x*qz + new_y*qw - new_z*qx;
		z = -new_w*qz - new_x*qy + new_y*qx + new_z*qw;

		return *this;
	}
	//	vector properties
	inline YZ_REAL SquareLength() const{
		return x*x + y*y + z*z;
	}

	inline YZ_REAL Length() const{
		return sqrt(SquareLength());
	}

};

/**
	Vector 4D
*/
template<class T>
class Vec4 : public Com4<T>{
public:
	using Com4<T>::x;
	using Com4<T>::y;
	using Com4<T>::z;
	using Com4<T>::w;
	using Com4<T>::dims;
	using typename Com4<T>::Type;
	using typename Com4<T>::YZ_REAL;

	//	constructors
	Vec4(const T x_val=0, const T y_val=0, const T z_val=0, const T w_val=0):Com4<T>(x_val, y_val, z_val, w_val){}
	Vec4(const T *val_ptr):Com4<T>(val_ptr){}
	template<typename TYPE> Vec4(const Com4<TYPE> v):Com4<TYPE>(v){}
	template<typename TYPE> Vec4(const Com3<TYPE> v, const T w_val = 0)  { this->x = v.x, this->y = v.y, this->z = v.z, this->w = w_val;}

	//	operators
	inline Vec4 operator - () const{
		Vec4 result(-x, -y, -z, -w);
		return result;
	}

	template<typename TYPE> inline Vec4& operator = (const Vec4<TYPE> vec){
		x = vec.x;
		y = vec.y;
		z = vec.z;
		w = vec.w;
		return *this;
	}

	template<typename TYPE> inline Vec4& operator += (const Vec4<TYPE> vec){
		x += vec.x;
		y += vec.y;
		z += vec.z;
		w += vec.w;
		return *this;
	}

	template<typename TYPE> inline Vec4& operator -= (const Vec4<TYPE> vec){
		x -= vec.x;
		y -= vec.y;
		z -= vec.z;
		w -= vec.w;
		return *this;
	}

	template<typename TYPE> inline Vec4& operator *= (const TYPE val){
		x *= val;
		y *= val;
		z *= val;
		w *= val;
		return *this;
	}

	template<typename TYPE> inline Vec4& operator /= (const TYPE val){
		x /= val;
		y /= val;
		z /= val;
		w /= val;
		return *this;
	}

	template<typename TYPE> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator + (const Vec4<TYPE> vec) const{
		Vec4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += vec;
		return result;
	}

	template<typename TYPE> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator - (const Vec4<TYPE> vec) const{
		Vec4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= vec;
		return result;
	}

	template<typename TYPE> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val) const{
		Vec4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator / (const TYPE val) const{
		Vec4<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	/**
		Return normalized vector, but don't change old data
	*/
	inline Vec4<YZ_REAL> Normalize() const{
		Vec4<YZ_REAL> result(*this);
		result.SetNormalize();
		return result;
	}


	//	set vectors
	/**
		Normalize the vector
	*/
	inline Vec4& SetNormalize(){
		YZ_REAL length = Length();
		(*this) /= length;
		return *this;
	}


	//	vector properties
	inline YZ_REAL SquareLength() const{
		return x*x + y*y + z*z + w*w;
	}

	inline YZ_REAL Length() const{
		return sqrt(SquareLength());
	}

};

//	========================================
///@{
/**	@name Non-menber Vector Functions
*/
//	========================================
template<typename TYPE, class T> inline Vec2<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Vec2<T> vec){
	Vec2<TYPE_PROMOTE(T, TYPE)> result(vec);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Vec3<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Vec3<T> vec){
	Vec3<TYPE_PROMOTE(T, TYPE)> result(vec);
	result *= val;
	return result;
}

template<typename TYPE, class T> inline Vec4<TYPE_PROMOTE(T, TYPE)> operator * (const TYPE val, const Vec4<T> vec){
	Vec4<TYPE_PROMOTE(T, TYPE)> result(vec);
	result *= val;
	return result;
}

template<class T1, class T2> inline PROMOTE_T1_T2 dot(const Vec2<T1> a, const Vec2<T2> b){
	return a.x * b.x + a.y * b.y;
}

template<class T1, class T2> inline PROMOTE_T1_T2 dot(const Vec3<T1> a, const Vec3<T2> b){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<class T1, class T2> inline PROMOTE_T1_T2 dot(const Vec4<T1> a, const Vec4<T2> b){
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}



template<class T> inline PROMOTE_T_TO_FLOAT cross(const Vec2<T> a, const Vec2<T> b){
	return a.x * b.y - a.y * b.x;
}

template<class T1, class T2> inline Vec3<TYPE_PROMOTE(T1, T2)> cross(const Vec3<T1> a, const Vec3<T2> b){
	Vec3<TYPE_PROMOTE(T1, T2)> result;
	result.x = a.y * b.z - a.z * b.y;
	result.y = a.z * b.x - a.x * b.z;
	result.z = a.x * b.y - a.y * b.x;
	return result;
}

///@}

//	========================================
//	Static member
//	========================================
template <class T> const int Com2<T>::dims = 2;
template <class T> const int Com3<T>::dims = 3;
template <class T> const int Com4<T>::dims = 4;


//	========================================
//	Type Rename
//	========================================
typedef Com2<char>					char2;
typedef Com2<unsigned char>			uchar2;
typedef	Com2<short>					short2;
typedef	Com2<unsigned short>		ushort2;
typedef Com2<int>					int2;
typedef Com2<unsigned int>			uint2;
typedef Com2<long>					long2;
typedef Com2<unsigned long>			ulong2;
typedef Com2<long long>				longlong2;
typedef Com2<unsigned long long>	ulonglong2;
typedef Com2<float>					float2;
typedef Com2<double>				double2;
typedef Com2<long double>			longdouble2;

typedef Com3<char>					char3;
typedef Com3<unsigned char>			uchar3;
typedef	Com3<short>					short3;
typedef	Com3<unsigned short>		ushort3;
typedef Com3<int>					int3;
typedef Com3<unsigned int>			uint3;
typedef Com3<long>					long3;
typedef Com3<unsigned long>			ulong3;
typedef Com3<long long>				longlong3;
typedef Com3<unsigned long long>	ulonglong3;
typedef Com3<float>					float3;
typedef Com3<double>				double3;
typedef Com3<long double>			longdouble3;

typedef Com4<char>					char4;
typedef Com4<unsigned char>			uchar4;
typedef	Com4<short>					short4;
typedef	Com4<unsigned short>		ushort4;
typedef Com4<int>					int4;
typedef Com4<unsigned int>			uint4;
typedef Com4<long>					long4;
typedef Com4<unsigned long>			ulong4;
typedef Com4<long long>				longlong4;
typedef Com4<unsigned long long>	ulonglong4;
typedef Com4<float>					float4;
typedef Com4<double>				double4;
typedef Com4<long double>			longdouble4;

typedef Vec2<char>					Vec2c;
typedef Vec2<unsigned char>			Vec2uc;
typedef Vec2<short>					Vec2s;
typedef Vec2<unsigned short>		Vec2us;
typedef Vec2<int>					Vec2i;
typedef Vec2<unsigned int>			Vec2ui;
typedef Vec2<long>					Vec2l;
typedef Vec2<unsigned long>			Vec2ul;
typedef Vec2<long long>				Vec2ll;
typedef Vec2<unsigned long long>	Vec2ull;
typedef Vec2<float>					Vec2f;
typedef Vec2<double>				Vec2d;
typedef Vec2<long double>			Vec2ld;

typedef Vec3<char>					Vec3c;
typedef Vec3<unsigned char>			Vec3uc;
typedef Vec3<short>					Vec3s;
typedef Vec3<unsigned short>		Vec3us;
typedef Vec3<int>					Vec3i;
typedef Vec3<unsigned int>			Vec3ui;
typedef Vec3<long>					Vec3l;
typedef Vec3<unsigned long>			Vec3ul;
typedef Vec3<long long>				Vec3ll;
typedef Vec3<unsigned long long>	Vec3ull;
typedef Vec3<float>					Vec3f;
typedef Vec3<double>				Vec3d;
typedef Vec3<long double>			Vec3ld;

typedef Vec4<char>					Vec4c;
typedef Vec4<unsigned char>			Vec4uc;
typedef Vec4<short>					Vec4s;
typedef Vec4<unsigned short>		Vec4us;
typedef Vec4<int>					Vec4i;
typedef Vec4<unsigned int>			Vec4ui;
typedef Vec4<long>					Vec4l;
typedef Vec4<unsigned long>			Vec4ul;
typedef Vec4<long long>				Vec4ll;
typedef Vec4<unsigned long long>	Vec4ull;
typedef Vec4<float>					Vec4f;
typedef Vec4<double>				Vec4d;
typedef Vec4<long double>			Vec4ld;

}	//	end namespace yz

#endif	//	__YZ_VECTOR_H__