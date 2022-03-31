/***********************************************************/
/**	\file
	\brief		Clipping
	\author		Yizhong Zhang
	\date		9/18/2012
*/
/***********************************************************/
#ifndef __YZ_CLIPPING_H__
#define __YZ_CLIPPING_H__

#include <iostream>
#include "yzLib/yz_math/yz_math_setting.h"
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

/**
	Clipper using Cohen-Sutherland Algorithm

	To use the clipper \n
	1, create an instance of this class	with type T and dimension N. 
	Pass clipping box to the class by constructor or set explicitly. \n
	2, call Clip(p0, p1) to perform clipping. Clipped points are stored
	in the original array.
*/
template<class T, int N>
class ClipperCohenSutherland{
public:
	//	constructor
	ClipperCohenSutherland(T bbmin[N], T bbmax[N]){
		SetClippingBox(bbmin, bbmax);
	}

	/**
		set the clipping box explicitly, with security check

		bbmin is expected to be smaller than bbmax in any dimension

		\param	bbmin	bonding box min
		\param	bbmax	bonding box max
	*/
	inline void SetClippingBox(T bbmin[N], T bbmax[N]){
		//	check 
		for(int i=0; i<N; i++){
			if( bbmin[i] > bbmax[i] ){
				#ifndef BE_QUIET
				std::cout << "error: set clipping box, min > max" << std::endl;
				#endif
			}
			bmin[i] = bbmin[i];
			bmax[i] = bbmax[i];
		}	
	}

	/**
		Clip the given point to bonding box

		write cliped point position to the original array directly

		\param	p0		point 0 position
		\param	p1		point 1 position
		\return			bit 0:	whether line segment intersect the box \n
						bit 1:	whether p0 moved \n
						bit 2:	whether p1 moved \n
						if the line segment don't intersect the box, the value should be 0, \n
						if neigher p0 nor p1 is moved, this indicates that the line segment is inside the box, return 1
	*/
	inline int Clip(T p0[N], T p1[N]){
		// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
		int p0_move_flag = 0;
		int p1_move_flag = 0;

		int outcode0 = ComputeOutCode( p0 );
		int outcode1 = ComputeOutCode( p1 );

		bool accept = false;

		while (true) {
			if (!(outcode0 | outcode1)) {	// Bitwise OR is 0. Trivially accept and get out of loop
				accept = true;
				break;
			} 
			else if (outcode0 & outcode1) { // Bitwise AND is not 0. Trivially reject and get out of loop
				break;
			} 
			else {
				// failed both tests, so calculate the line segment to clip
				// from an outside point to an intersection with clip edge
				T p[N];

				// At least one endpoint is outside the clip rectangle; pick it.
				int outcodeOut = outcode0 ? outcode0 : outcode1;

				// Now find the intersection point;
				for( int i=0; i<N; i++ ){
					if(outcodeOut & GetCode(i, 0)){			//	below the clip plane of dimension i
						p[i] = bmin[i];
						PROMOTE_T_TO_FLOAT	t = PROMOTE_T_TO_FLOAT(p[i] - p0[i]) / (p1[i] - p0[i]);
						for(int j=0; j<N; j++)
							if(i != j)	p[j] = p0[j] + t * (p1[j] - p0[j]);
						break;
					}
					else if(outcodeOut & GetCode(i, 1)){	//	above the clip plane of dimension i
						p[i] = bmax[i];
						PROMOTE_T_TO_FLOAT	t = PROMOTE_T_TO_FLOAT(p[i] - p0[i]) / (p1[i] - p0[i]);
						for(int j=0; j<N; j++)
							if(i != j)	p[j] = p0[j] + t * (p1[j] - p0[j]);
						break;
					}
				}

				// Now we move outside point to intersection point to clip
				// and get ready for next pass.
				if (outcodeOut == outcode0) {
					for( int i=0; i<N; i++ )
						p0[i] = p[i];
					outcode0 = ComputeOutCode(p0);
					p0_move_flag = 0x02;
				} 
				else {
					for( int i=0; i<N; i++ )
						p1[i] = p[i];
					outcode1 = ComputeOutCode(p1);
					p1_move_flag = 0x04;
				}
			}
		}
		if (accept)		//	the line cross the screen
			return p0_move_flag | p1_move_flag | 0x01;
		else
			return 0;
	}
public:
	//	the bonding box
	T bmin[N], bmax[N];

protected:
	/**
		get code of point, point are stored in the array in sequential dimension
	*/
	inline int ComputeOutCode(T p[N]){
		int code = 0;

		for( int i=0; i<N; i++ ){
			if( p[i] < bmin[i] )
				code |= GetCode(i, 0);
			else if( p[i] > bmax[i] )
				code |= GetCode(i, 1);
		}

		return code;
	}

	/**
		Get Cohen-Sutherland code of this dimension and side

		\param	dim		which dimension of the window
		\param	side	0: min side, 1: max side
		\return			return Cohen-Sutherland space code
	*/
	inline int GetCode(int dim, int side){
		return (0x01 << (dim*2+side));
	}
};

//	========================================
///@{
/**	@name Cohen-Sutherland Clipping
*/
//	========================================

/**
	clipping by Cohen-Sutherland Clipping 2D

	\param	seg_v0	line segment end point 0
	\param	seg_v1	line segment end point 1
	\param	bb_min	clipping box min
	\param	bb_max	clipping box max
	\return			bit 0:	whether line segment intersect the box \n
					bit 1:	whether end point 0 moved \n
					bit 2:	whether end point 1 moved \n
					if the line segment don't intersect the box, the value should be 0, \n
					if neigher end points are moved, this indicates that the line segment is inside the box, return 1
*/
template<typename T>
inline int clipCohenSutherland(Com2<T>& seg_v0, 
							   Com2<T>& seg_v1, 
							   Com2<T>	bb_min, 
							   Com2<T>	bb_max){
	ClipperCohenSutherland<T, 2> clipper((T*)&bb_min, (T*)&bb_max);
	return clipper.Clip((T*)&seg_v0, (T*)&seg_v1);
}

/**
	clipping by Cohen-Sutherland Clipping

	\param	seg_v0	line segment end point 0
	\param	seg_v1	line segment end point 1
	\param	bb_min	clipping box min
	\param	bb_max	clipping box max
	\return			bit 0:	whether line segment intersect the box \n
					bit 1:	whether end point 0 moved \n
					bit 2:	whether end point 1 moved \n
					if the line segment don't intersect the box, the value should be 0, \n
					if neigher end points are moved, this indicates that the line segment is inside the box, return 1
*/
template<typename T>
inline int clipCohenSutherland(Com3<T>& seg_v0, 
							   Com3<T>& seg_v1, 
							   Com3<T>	bb_min, 
							   Com3<T>	bb_max){
	ClipperCohenSutherland<T, 3> clipper((T*)&bb_min, (T*)&bb_max);
	return clipper.Clip((T*)&seg_v0, (T*)&seg_v1);
}

///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_CLIPPING_H__