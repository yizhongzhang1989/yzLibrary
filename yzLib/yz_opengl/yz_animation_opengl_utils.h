/***********************************************************/
/**	\file
	\brief		OpenGL Utilities of Animation
	\author		Yizhong Zhang
	\date		7/1/2012
*/
/***********************************************************/
#ifndef __YZ_ANIMATION_OPENGL_UTILS_H__
#define __YZ_ANIMATION_OPENGL_UTILS_H__

#include "yzLib/yz_setting.h"

#if !(	defined(YZ_glut_h) || defined(YZ_freeglut_h) )
#	error yz_animation_opengl_utils.h must be included after glut.h or freeglut.h
#endif

#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_animation/yz_skeleton.h"
#include "yzLib/yz_animation/yz_skeleton_utils.h"
#include "yzLib/yz_opengl/yz_opengl_utils.h"
#include "yzLib/yz_opengl/yz_vector_opengl_utils.h"

namespace yz{	namespace opengl{

template<class T>	class Skeleton;

//	========================================
///@{
/**	@name Draw Skeleton
*/
//	========================================

/**
	Draw the skeleton as sticks

	If the skeleton is not continuous, then it is possible that
	the shape displayed doesn't not fully match the skeleton shape
	at discontinuous places
*/
template<typename T>
inline void drawSkeletonAsSticks(const animation::Skeleton<T>& skeleton, 
								 float stick_radius	= 0.01f, 
								 float sphere_red	= 1,
								 float sphere_green	= 0,
								 float sphere_blue	= 0,
								 float stick_red	= 1,
								 float stick_green	= 1,
								 float stick_blue	= 0 ){
	if( skeleton.bone.empty() )
		return;

	std::vector<Vec3<T>> joint;
	animation::getSkeletonJointPosition(joint, skeleton);

	glColor3f(sphere_red, sphere_green, sphere_blue);
	for(int i=0; i<joint.size(); i++)
		drawPointAsSphere(joint[i], stick_radius*2);

	glColor3f(stick_red, stick_green, stick_blue);
	for(int i=0; i<skeleton.bone.size(); i++)
		drawCylinder(joint[i+1], joint[i+1] - skeleton.bone[i].bone_vector, stick_radius);
}

/**
	Draw the skeleton as flat solid bones
*/
template<typename T>
inline void drawSkeletonAsFlatSolidBones(const animation::Skeleton<T>& skeleton){
	if( skeleton.bone.empty() )
		return;

	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++)
		drawFlatSolidBone(node[i], skeleton.bone[i].bone_vector, skeleton.bone[i].bone_front);
}

/**
	Draw the skeleton as color cube bones
*/
template<typename T>
inline void drawSkeletonAsColorCubeBones(const animation::Skeleton<T>& skeleton){
	if( skeleton.bone.empty() )
		return;

	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++){
		drawColorCubeBone(node[i], skeleton.bone[i].bone_vector, skeleton.bone[i].bone_front);
	}
}
/**
	Draw the skeleton as solid bones
*/
template<typename T>
inline void drawSkeletonAsSolidBones(const animation::Skeleton<T>& skeleton){
	if( skeleton.bone.empty() )
		return;

	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++)
		drawSolidBone(node[i], skeleton.bone[i].bone_vector, skeleton.bone[i].bone_front);
}

/**
	Draw the skeleton as solid bones, color of each bone is specified explicitly

	\param	skeleton		the skeleton to display
	\param	bone_color		color of each bone
*/
template<typename T>
inline void drawSkeletonAsSolidBones(const animation::Skeleton<T>&	skeleton,
									 const std::vector<float3>&		bone_color){
	if( skeleton.bone.size() > bone_color.size() ){
		#ifndef BE_QUIET
			std::cout << "error: drawSkeletonAsSolidBones, bone_color size don't match skeleton" << std::endl;
		#endif
		return;
	}
	if( skeleton.bone.empty() )
		return;

	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++){
		glColor3f(bone_color[i].x, bone_color[i].y, bone_color[i].z);
		drawSolidBone(node[i], skeleton.bone[i].bone_vector, skeleton.bone[i].bone_front);
	}
}

/**
	Draw the skeleton as wire bones
*/
template<typename T>
inline void drawSkeletonAsWireBones(const animation::Skeleton<T>& skeleton){
	if( skeleton.bone.empty() )
		return;
	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++)
		drawWireBone(node[i], skeleton.bone[i].bone_vector, skeleton.bone[i].bone_front);
}

/**
	draw name of each bone of the skeleton
*/
template<typename T>
inline void drawSkeletonBoneName(const animation::Skeleton<T>& skeleton){
	if( skeleton.bone.empty() )
		return;
	std::vector<Vec3<T>> node;
	animation::getSkeletonNodePosition(node, skeleton);

	for(int i=0; i<skeleton.bone.size(); i++){
		Vec3<T> center = node[i] + skeleton.bone[i].bone_vector*0.5;
		drawString(skeleton.bone[i].name, center[0], center[1], center[2]);
	}
}


///@}

}}	//	end namespace yz::opengl


#endif	//	__YZ_ANIMATION_OPENGL_UTILS_H__