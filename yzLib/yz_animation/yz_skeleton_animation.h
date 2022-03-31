/***********************************************************/
/**	\file
	\brief		Skeleton and Animation
	\author		Yizhong Zhang
	\date		6/28/2012
*/
/***********************************************************/
#ifndef __YZ_SKELETON_ANIMATION_H__
#define __YZ_SKELETON_ANIMATION_H__

#include <iostream>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_animation/yz_skeleton.h"
#include "yzLib/yz_animation/yz_skeleton_utils.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_animation_opengl_utils.h"	//	to display the skeleton
#endif

namespace yz{	namespace animation{

/**
	Class of skeleton animation

	this class contain one skeleton and corresponding transformer frames.
	currently, the data are read from .bvh file.
*/
template<class T>
class SkeletonAnimation{
public:
	Skeleton<T>						skeleton;			///<	skeleton, deform as you like
	SkeletonTransformerData<T>		transformer_data;	///<	deformation frames data

public:
	/**
		read skeleton and frames data from .bvh file

		\param	bvh_file_name	file name of the .bvh file
		\param	meter_per_unit	actual length in meters of one unit in .bvh file, also scale value
		\param	angle_unit		unit of the angle in .bvh file, 'd' degree,  'r' rad
		\return					whether read file succeed
	*/
	inline int ReadSkeletonAndFramesFromBVH(const char* bvh_file_name, double meter_per_unit=1, char angle_unit='d'){
		int ret = readMocapFromBVH(bvh_file_name, &rest_skeleton, &motion_data, meter_per_unit, angle_unit);
		if( ret ){
			skeleton = rest_skeleton;
			transformer_data.SetTransformerFromMotion(motion_data);
		}
		return ret;
	}

	/**
		Reset frame interval used in transformer data

		if we want to use different frame interval, call this 
		function to interpolate the transformer data
	*/
	int ResetFrameInterval(T target_interval){
		SkeletonMotionQuaternionData<T> quaternion_data;
		quaternion_data.Set(motion_data);

		SkeletonMotionQuaternionData<T> new_quaternion_data;
		resampleSkeletonMotion(new_quaternion_data, target_interval, quaternion_data);

		transformer_data.Set(new_quaternion_data);
		return 1;
	}

	/**
		get transformed skeleton of a given frame

		\param	frame_id	id of the frame, if invalid id passed, return 0
		\return				whether transform succeed
	*/
	inline int TransformSkeletonToFrame(int frame_id){
		if( frame_id<0 )								return 0;
		if( frame_id>=transformer_data.FrameNumber() )	return 0;

		skeleton = rest_skeleton;
		return transformSkeleton(skeleton, transformer_data.frame[frame_id]);
	}

	/**
		set the skeleton to rest pose
	*/
	inline void TransformSkeletonToRest(){
		skeleton = rest_skeleton;
	}

	/**
		Draw Skeleton

		if OpenGL is not enabled, print error information

		\param	display_mode	0: stick,  1: solid bones,  2: wire bones
		\param	stick_radius	radius of stick if display_mode is 0(stick)
		\return					whether display succeed
	*/
	inline int Display(int display_mode=0, float stick_radius=0.01f) const{
		return skeleton.Display(display_mode, stick_radius);
	}

public:
	Skeleton<T>				rest_skeleton;	///<	rest pose of the skeleton, never try to change it 
	SkeletonMotionData<T>	motion_data;	///<	original motion data, never try to change it
};
/**
	\example	character_animation.cpp

	This file shows who to create character animation
*/



}}	//	namespace yz::animation

#endif	//	__YZ_SKELETON_ANIMATION_H__