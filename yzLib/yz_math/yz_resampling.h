/***********************************************************/
/**	\file
	\brief		Resampling
	\author		Yizhong Zhang
	\date		7/5/2012
*/
/***********************************************************/
#ifndef __YZ_RESAMPLING_H__
#define __YZ_RESAMPLING_H__

#include <assert.h>
#include "yzLib/yz_math/yz_interpolation.h"

namespace yz{

//	========================================
///@{
/**	@name 1D Resample
*/
//	========================================
/**
	Resampling using Nearest Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleNearest(T1* des, int des_size, const T2* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		des[i] = interpNearest(src[i_coef], src[i_coef+1], coef);
	}
	des[des_size-1] = src[src_size-1];
}

/**
	Resampling using Linear Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleLinear(T1* des, int des_size, const T2* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		des[i] = interpLinear(src[i_coef], src[i_coef+1], coef);
	}
	des[des_size-1] = src[src_size-1];
}

/**
	Resampling using Cosine Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleCosine(T1* des, int des_size, const T2* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		des[i] = interpCosine(src[i_coef], src[i_coef+1], coef);
	}
	des[des_size-1] = src[src_size-1];
}

/**
	Resampling using Cubic Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleCubic(T1* des, int des_size, const T2* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		T2	src_l = (i_coef==0 ? src[0] : src[i_coef-1]);
		T2	src_r = (i_coef==src_size-2 ? src[i_coef+1] : src[i_coef+2]);
		des[i] = interpCubic(src_l, src[i_coef], src[i_coef+1], src_r, coef);
	}
	des[des_size-1] = src[src_size-1];
}

/**
	Resampling using Cutmull-Rom Spline Interpolation, from src to des

	\param	des				destination resampling array
	\param	des_size		destination resampling array size
	\param	src				source value array
	\param	src_size		source value array size, minimal value: 2, or the program will cause error
*/
template<typename T1, typename T2>
inline void resampleCatmullRom(T1* des, int des_size, const T2* src, int src_size){
	assert(src_size>1 && des_size>0);

	for( int i=0; i<des_size-1; i++ ){
		double	coef = double(i) / (des_size-1) * (src_size-1);
		int		i_coef = coef;
		coef -= i_coef;
		T2	src_l = (i_coef==0 ? src[0] : src[i_coef-1]);
		T2	src_r = (i_coef==src_size-2 ? src[i_coef+1] : src[i_coef+2]);
		des[i] = interpCatmullRom(src_l, src[i_coef], src[i_coef+1], src_r, coef);
	}
	des[des_size-1] = src[src_size-1];
}

///@}

//	========================================
///@{
/**	@name 2D Resample
*/
//	========================================
/**
	Resampling between arrays using Nearest Interpolation 2D, from src to des

	\param	des				destination 2D array
	\param	des_w			width of destination array
	\param	des_h			height of destination array
	\param	src				source value array
	\param	src_w			width of source array
	\param	src_h			height of source array
*/
template<typename T1, typename T2>
inline void resampleNearest(T1* des, int des_w, int des_h, const T2* src, int src_w, int src_h){
	assert(src_w>1 && src_h>1 && des_w>0 && des_h>0);

	for( int j=0; j<des_h; j++ ){
		double	ty = double(j)/(des_h-1) * (src_h-1);
		int i_ty = ty;
		ty -= i_ty;
		if( j==des_h-1 )	ty = 1, i_ty = src_h-2;
		for( int i=0; i<des_w-1; i++ ){
			double	tx = double(i)/(des_w-1) * (src_w-1);
			int		i_tx = tx;
			tx -= i_tx;
			des[j*des_w+i] = interpNearest(
				src[i_ty*src_w+i_tx], src[i_ty*src_w+i_tx+1], 
				src[(i_ty+1)*src_w+i_tx], src[(i_ty+1)*src_w+i_tx+1], 
				tx, ty);
		}
		des[j*des_w+des_w-1] = interpNearest(
			src[(i_ty+1)*src_w-2], src[(i_ty+1)*src_w-1], 
			src[(i_ty+2)*src_w-2], src[(i_ty+2)*src_w-1], 
			1, ty);
	}

}

/**
	Resampling between arrays using Bilinear Interpolation, from src to des

	\param	des				destination 2D array
	\param	des_w			width of destination array
	\param	des_h			height of destination array
	\param	src				source value array
	\param	src_w			width of source array
	\param	src_h			height of source array
*/
template<typename T1, typename T2>
inline void resampleBilinear(T1* des, int des_w, int des_h, const T2* src, int src_w, int src_h){
	assert(src_w>1 && src_h>1 && des_w>0 && des_h>0);

	for( int j=0; j<des_h; j++ ){
		double	ty = double(j)/(des_h-1) * (src_h-1);
		int i_ty = ty;
		ty -= i_ty;
		if( j==des_h-1 )	ty = 1, i_ty = src_h-2;
		for( int i=0; i<des_w-1; i++ ){
			double	tx = double(i)/(des_w-1) * (src_w-1);
			int		i_tx = tx;
			tx -= i_tx;
			des[j*des_w+i] = interpBilinear(
				src[i_ty*src_w+i_tx], src[i_ty*src_w+i_tx+1], 
				src[(i_ty+1)*src_w+i_tx], src[(i_ty+1)*src_w+i_tx+1], 
				tx, ty);
		}
		des[j*des_w+des_w-1] = interpBilinear(
			src[(i_ty+1)*src_w-2], src[(i_ty+1)*src_w-1], 
			src[(i_ty+2)*src_w-2], src[(i_ty+2)*src_w-1], 
			1, ty);
	}
}

/**
	Resampling between arrays using Bicubic Interpolation, from src to des

	\param	des				destination 2D array
	\param	des_w			width of destination array
	\param	des_h			height of destination array
	\param	src				source value array
	\param	src_w			width of source array
	\param	src_h			height of source array
*/
template<typename T1, typename T2>
inline void resampleBicubic(T1* des, int des_w, int des_h, const T2* src, int src_w, int src_h){
	assert(src_w>1 && src_h>1 && des_w>0 && des_h>0);

	for( int j=0; j<des_h; j++ ){
		double	ty = double(j)/(des_h-1) * (src_h-1);
		int i_ty = ty;
		ty -= i_ty;
		if( j==des_h-1 )	ty = 1, i_ty = src_h-2;
		int i_u = (i_ty==0 ? 0 : i_ty-1);
		int i_d = (i_ty==src_h-2 ? src_h-1 : i_ty+2);

		for( int i=0; i<des_w; i++ ){
			double tx = double(i)/(des_w-1) * (src_w-1);
			int i_tx = tx;
			tx -= i_tx;
			if(i == des_w-1)	tx = 1, i_tx = src_w-2;
			int i_l = (i_tx==0 ? 0 : i_tx-1);
			int i_r = (i_tx==src_w-2 ? src_w-1 : i_tx+2);
			des[j*des_w+i] = interpBicubic(
				src[i_u*src_w+i_l],			src[i_u*src_w+i_tx],		src[i_u*src_w+i_tx+1],		src[i_u*src_w+i_r],
				src[i_ty*src_w+i_l],		src[i_ty*src_w+i_tx],		src[i_ty*src_w+i_tx+1],		src[i_ty*src_w+i_r],
				src[(i_ty+1)*src_w+i_l],	src[(i_ty+1)*src_w+i_tx],	src[(i_ty+1)*src_w+i_tx+1],	src[(i_ty+1)*src_w+i_r],
				src[i_d*src_w+i_l],			src[i_d*src_w+i_tx],		src[i_d*src_w+i_tx+1],		src[i_d*src_w+i_r],
				tx, ty);
		}
	}
}


///@}

}	//	namespace yz

#endif	//	__YZ_RESAMPLING_H__