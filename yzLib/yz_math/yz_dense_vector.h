/***********************************************************/
/**	\file
	\brief		Dense Vector
	\details	dense vector is used to represent vector of any
				dimension. It can operate with dense matrix or 
				sparse matrix. 
	\author		Yizhong Zhang
	\date		9/10/2012
*/
/***********************************************************/
#ifndef __YZ_DENSE_VECTOR_H__
#define __YZ_DENSE_VECTOR_H__

#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_math_setting.h"

namespace yz{

/**
	Dense vector

	dense vector can represent vector of any dimension using 
	std::vector storing all the data. 
*/
template<class T>
class DenseVector{
public:
	std::vector<T>	value;

	typedef PROMOTE_T_TO_FLOAT YZ_REAL;

public:
	//	constructors
	DenseVector(){}

	DenseVector(int dim_num){
		value.resize(dim_num, 0);
	}

	template<typename TYPE> DenseVector(const TYPE* data, int dim_num){
		value.resize(dim_num);
		value.assign(data, data+dim_num);
	}

	DenseVector(const DenseVector<T>& vec){
		if( &vec != &(*this) )
			value.assign(vec.value.begin(), vec.value.end());
	}

	template<typename TYPE> DenseVector(const DenseVector<TYPE>& vec){
		value.resize(vec.Dim());
		value.assign(vec.value.begin(), vec.value.end());
	}

	/**
		Reset the vector, clear everything
	*/
	inline void Reset(){
		value.clear();
	}

	/**
		Set dimension of the vector
	*/
	inline void SetDimension(int dim_num){
		assert(value.size() == 0);
		assert(dim_num > 0);
		value.resize(dim_num, 0);
	}

	/**
		Set Zero
	*/
	inline void SetZero(){
		//memset(&value[0], 0, sizeof(T)*value.size());
		value.assign(value.size(), 0);
	}

	/**
		Set all elements to the same value
	*/
	template<typename TYPE> inline void SetValue(TYPE val){
		for( int i=0; i<Dim(); i++ )
			value[i] = val;
	}

	/**
		Copy the value from an array to the vector
	*/
	template<typename TYPE> inline void SetValue(TYPE* data){
		value.assign(data, data + Dim());
	}

	//	operators
	inline T& operator[] (int id){
		assert(id>=0 && id<value.size());
		return value[id];
	}

	inline const T& operator[] (int id) const{
		assert(id>=0 && id<value.size());
		return value[id];
	}

	inline T& operator() (int id){
		assert(id>=0 && id<value.size());
		return value[id];
	}

	inline const T& operator() (int id) const{
		assert(id>=0 && id<value.size());
		return value[id];
	}

	inline DenseVector<T> operator - () const{
		DenseVector<T> result(*this);
		for( int i=0; i<result.Dim(); i++ )
			result.value[i] = -result.value[i];
		return result;
	}

	inline DenseVector& operator = (const DenseVector<T>& vec){
		if( &vec != &(*this) )
			value.assign(vec.value.begin(), vec.value.end());
		return *this;
	}

	template<typename TYPE> inline DenseVector& operator = (const DenseVector<TYPE>& vec){
		value.resize(vec.Dim());
		value.assign(vec.value.begin(), vec.value.end());
		return *this;
	}

	template<typename TYPE> inline DenseVector& operator += (const DenseVector<TYPE>& vec){
		AssureDimLegal(vec);
		for( int i=0; i<Dim(); i++ )
			value[i] += vec.value[i];
		return *this;
	}

	template<typename TYPE> inline DenseVector& operator -= (const DenseVector<TYPE>& vec){
		AssureDimLegal(vec);
		for( int i=0; i<Dim(); i++ )
			value[i] -= vec.value[i];
		return *this;
	}

	template<typename TYPE> inline DenseVector& operator *= (TYPE val){
		for( int i=0; i<Dim(); i++ )
			value[i] *= val;
		return *this;
	}

	template<typename TYPE> inline DenseVector& operator /= (TYPE val){
		for( int i=0; i<Dim(); i++ )
			value[i] /= val;
		return *this;
	}

	template<typename TYPE> inline DenseVector<TYPE_PROMOTE(T, TYPE)> operator + (const DenseVector<TYPE>& vec) const{
		DenseVector<TYPE_PROMOTE(T, TYPE)> result(*this);
		result += vec;
		return result;
	}

	template<typename TYPE> inline DenseVector<TYPE_PROMOTE(T, TYPE)> operator - (const DenseVector<TYPE>& vec) const{
		DenseVector<TYPE_PROMOTE(T, TYPE)> result(*this);
		result -= vec;
		return result;
	}

	template<typename TYPE> inline DenseVector<TYPE_PROMOTE(T, TYPE)> operator * (TYPE val) const{
		DenseVector<TYPE_PROMOTE(T, TYPE)> result(*this);
		result *= val;
		return result;
	}

	template<typename TYPE> inline DenseVector<TYPE_PROMOTE(T, TYPE)> operator / (TYPE val) const{
		DenseVector<TYPE_PROMOTE(T, TYPE)> result(*this);
		result /= val;
		return result;
	}

	/**
		Return normalized vector, but don't change old data
	*/
	inline DenseVector<YZ_REAL> Normalize() const{
		DenseVector<YZ_REAL> result(*this);
		result.SetNormalize();
		return result;
	}

	/**
		normalize this vector
	*/
	inline DenseVector& SetNormalize(){
		YZ_REAL length = Length();
		(*this) /= length;
		return *this;
	}

	//	vector properties
	/**
		return the dimension of the vector
	*/
	inline int Dim() const{
		return value.size();
	}

	/**
		return the sum of all elements
	*/
	inline T Sum() const{
		T sum = 0;
		for( int i=0; i<Dim(); i++ )
			sum += value[i];
		return sum;
	}

	inline YZ_REAL SquareLength() const{
		T sqr_sum = 0;
		for( int i=0; i<Dim(); i++ )
			sqr_sum += value[i] * value[i];
		return sqr_sum;
	}

	inline YZ_REAL Length() const{
		return sqrt(SquareLength());
	}
protected:
	/**
		assure dimension of the vector to calculate is the same with the current vector

		\param	vec		the vector to check
	*/
	template<typename TYPE> inline void AssureDimLegal(const DenseVector<TYPE>& vec){
		if( Dim() != vec.Dim() ){
			std::cout << "Dimension of vector illegal, vector dimension don't match" << std::endl;
			exit(0);
		}
	}

};

//	========================================
///@{
/**	@name Non-menber Vector Functions
*/
//	========================================
template<class T1, class T2> 
inline PROMOTE_T1_T2 dot(const DenseVector<T1>& den_vec1, const DenseVector<T2>& den_vec2){
	if( den_vec1.Dim() != den_vec2.Dim() ){
		std::cout << "error:dot, dimension of two dense vectors doesn't match" << std::endl;
		exit(0);
	}
	PROMOTE_T1_T2	sum = 0;
	for(int i=0; i<den_vec1.Dim(); i++)
		sum += den_vec1[i] * den_vec2[i];

	return sum;
}
///@}

}	//	end namespace yz

#endif	//	__YZ_DENSE_VECTOR_H__