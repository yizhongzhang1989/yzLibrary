/***********************************************************/
/**	\file
	\brief		util functions of skeleton
	\author		Yizhong Zhang
	\date		6/29/2012
*/
/***********************************************************/
#ifndef __YZ_SKELETON_UTILS_H__
#define __YZ_SKELETON_UTILS_H__

#include <vector>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_quaternion.h"
#include "yzLib/yz_animation/yz_skeleton.h"
#include "yzLib/yz_animation/yz_skin.h"

namespace yz{	namespace animation{

template<class T>	class Bone;
template<class T>	class Skeleton;
template<class T>	class BoneMotion;
template<class T>	class BoneMotionQuaternion;
template<class T>	class SkeletonMotion;
template<class T>	class SkeletonMotionQuaternion;
template<class T>	class SkeletonMotionData;
template<class T>	class SkeletonMotionQuaternionData;
template<class T>	class BoneTransformer;
template<class T>	class SkeletonTransformer;
template<class T>	class SkeletonTransformerData;

//	========================================
///@{
/**	@name Calculate transformer for skeleton
*/
//	========================================
/**
	Transform relative local transformer to absolute global transformer

	by default, the transformer apply to each bone in local coordinate
	system, which means to transform one bone, you must apply same 
	transform of its parents, then this local transform. The purpose 
	of this function is to calculate global transformer of each bone.

	\param	transformer		the local transformer, it will be changed to 
							global after calling this function
	\param	skeleton		skeleton of the transformer
	\return					whether the data is correst, if not, don't touch old data
*/
template<typename T>
int setTransformerToGlobal(SkeletonTransformer<T>& transformer, const Skeleton<T>& skeleton){
	if( skeleton.bone.size() != transformer.transformer.size() )
		return 0;

	if( transformer.global_flag )		//	the transformer is already global
		return 1;

	for( int i=0; i<skeleton.bone.size(); i++ ){
		int f_id = skeleton.bone[i].parent;
		if( f_id != -1 ){	//	multiply local rotation with its parents'
			transformer.transformer[i].rotation = 
				transformer.transformer[f_id].rotation * transformer.transformer[i].rotation;
		}
	}

	transformer.global_flag = 1;

	return 1;
}

/**
	calculate scale of each bone from source to target skeleton

	source * scale = target

	\param	scale				return scale of each bone
	\param	source_skeleton		the source skeleton
	\param	target_skeleton		the target skeleton
	\return						1: success, 0: failed
*/
template<typename T>
int getScaleOfEachBone(	std::vector<T>&			scale,
						const Skeleton<T>&		source_skeleton, 
						const Skeleton<T>&		target_skeleton){
	if( source_skeleton.bone.size() != target_skeleton.bone.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: getScaleOfEachBone, skeleton size don't match" << std::endl;
		#endif
		return 0;
	}

	scale.resize(source_skeleton.bone.size());
	for(int i=0; i<target_skeleton.bone.size(); i++){
		if(target_skeleton.bone[i].bone_vector.Length()<1e-6 && source_skeleton.bone[i].bone_vector.Length()<1e-6 )
			scale[i] = 1;
		else if(source_skeleton.bone[i].bone_vector.Length()<1e-6)
			scale[i] = 0;
		else
			scale[i] = target_skeleton.bone[i].bone_vector.Length() / source_skeleton.bone[i].bone_vector.Length();
	}	

	return 1;
}

/**
	Get scale transformer for the skeleton

	If the input transformer is global, we return global transformer, and vice versa.

	\param	transformer		return the transformer, we don't touch global status
	\param	skeleton		the skeleton to transform
	\param	scale			scale of each bone of the skeleton, 1 if don't scale
	\return					0: failed, 1: succeed
*/
template<typename T>
int getScaleTransformer(SkeletonTransformer<T>&	transformer, 
						const Skeleton<T>&		skeleton, 
						const std::vector<T>&	scale){
	if( skeleton.bone.size() != scale.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: getScaleTransformer, skeleton scale size don't match" << std::endl;
		#endif
		return 0;
	}

	//	create motion of scale
	int is_global = transformer.global_flag;
	SkeletonMotion<T>	skeleton_motion;
	skeleton_motion.motion.resize(skeleton.bone.size());
	for(int i=0; i<skeleton.bone.size(); i++){
		if( is_global )	//	the transformer is global, we calculate scale of each bone independently
			skeleton_motion.motion[i].scale = scale[i];
		else{							//	the transformer is local, we have to remove scale from its parent
			int parent_id = skeleton.bone[i].parent;

			if( parent_id != -1 )
				skeleton_motion.motion[i].scale = 1.0 / skeleton_motion.motion[parent_id].scale;
			else
				skeleton_motion.motion[i].scale = 1;

			skeleton_motion.motion[i].scale *= scale[i];
		}
	}

	transformer.SetTransformerFromMotion(skeleton_motion);
	transformer.global_flag = is_global;

	return 1;
}

/**
	get rotation transformer from source to target skeleton, 
	no translation of the root

	\param	rotation			return the rotation transformer
	\param	source_skeleton		the source skeleton
	\param	target_skeleton		the target skeleton
	\return						1: success, 0: failed
*/
template<typename T>
int getRotationTransformer(SkeletonTransformer<T>&	rotation, 
						   const Skeleton<T>&		source_skeleton, 
						   const Skeleton<T>&		target_skeleton){
	if( source_skeleton.bone.size() != target_skeleton.bone.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: getRotationTransformer, skeleton size don't match" << std::endl;
		#endif
		return 0;
	}

	rotation.global_flag = 1;

	int bones = target_skeleton.bone.size();
	rotation.transformer.resize(bones);

	for(int i=0; i<bones; i++){
		yz::Vec3<T>	r0 = source_skeleton.bone[i].bone_vector;
		yz::Vec3<T>	r1 = target_skeleton.bone[i].bone_vector;
		if( r0.Length()<1e-5 || r1.Length()<1e-5 )	//	if the bone is really short(a virtual bone), there is no need to rotate in global transformation
			rotation.transformer[i].rotation.SetIdentity();
		else{	//	calculate transformation regardless of scale
			assert(source_skeleton.bone[i].bone_front.SquareLength() > 0.99 &&	//	bone front is supposed to be with length 1
				target_skeleton.bone[i].bone_front.SquareLength() > 0.99 );
			r0.SetNormalize();
			r1.SetNormalize();
			yz::Vec3<T> bf0 = source_skeleton.bone[i].bone_front.Normalize();
			yz::Vec3<T> bf1 = target_skeleton.bone[i].bone_front.Normalize();

			//	bone vector should be perpendicular to bone front, if not, we just rotate the bone directly
			if( fabs(dot(r0, bf0)) < 1e-5 && fabs(dot(r1, bf1)) < 1e-5 ){
				Matrix3x3<T> source_mat(r0, bf0, cross(r0, bf0));
				Matrix3x3<T> target_mat(r1, bf1, cross(r1, bf1));
				rotation.transformer[i].rotation = target_mat * source_mat.Inverse();
			}
			else{	//	bone vector is not correct
				std::cout << "warning: getRotationTransformer, bone front not perpendicular to bone vector" << std::endl;
				yz::Vec3<T>	axis = yz::cross(r0, r1);
				T			angle = yz::angleRadBetweenVectors(r0, r1);
				rotation.transformer[i].rotation.SetRotationRad(axis, angle);
			}
		}
	}

	return 1;
}

/**
	get global transformer from source to target skeleton with rotation and translation,
	but scale is not included in the transformer

	\param	transformer			return the global transformer
	\param	source_skeleton		the source skeleton
	\param	target_skeleton		the target skeleton
	\return						1: success, 0: failed
*/
template<typename T>
int getTransformerWithoutScale(SkeletonTransformer<T>&	transformer, 
							   const Skeleton<T>&		source_skeleton, 
							   const Skeleton<T>&		target_skeleton){
	if( source_skeleton.bone.size() != target_skeleton.bone.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: getTransformer, skeleton size don't match" << std::endl;
		#endif
		return 0;
	}

	//	calculate rotation matrix
	transformer.global_flag = 1;
	getRotationTransformer(transformer, source_skeleton, target_skeleton);

	//	set translation
	transformer.transformer[0].translation = target_skeleton.bone[0].offset - source_skeleton.bone[0].offset;
	
	return 1;
}

/**
	get global transformer from source to target skeleton,

	\param	transformer			return the global transformer
	\param	source_skeleton		the source skeleton
	\param	target_skeleton		the target skeleton
	\return						1: success, 0: failed
*/
template<typename T>
int getTransformer(SkeletonTransformer<T>&	transformer, 
				   const Skeleton<T>&		source_skeleton, 
				   const Skeleton<T>&		target_skeleton){
	if( source_skeleton.bone.size() != target_skeleton.bone.size() ){
		#ifndef	BE_QUIET
			std::cout << "error: getTransformer, skeleton size don't match" << std::endl;
		#endif
		return 0;
	}

	//	calculate scale and rotation matrix
	std::vector<T>	scale;
	getScaleOfEachBone(scale, source_skeleton, target_skeleton);

	SkeletonTransformer<T>	scale_trans;
	scale_trans.global_flag = 1;
	getScaleTransformer(scale_trans, source_skeleton, scale);

	transformer.global_flag = 1;
	getRotationTransformer(transformer, source_skeleton, target_skeleton);

	//	make into transformer
	assert(transformer.transformer.size() == scale_trans.transformer.size());
	for(int i=0; i<transformer.transformer.size(); i++){
		transformer.transformer[i].rotation *= scale_trans.transformer[i].rotation;
		
		if( i==0 )	//	add translation for root
			transformer.transformer[i].translation = target_skeleton.bone[0].offset - source_skeleton.bone[0].offset;
	}
	
	return 1;
}

///@}

//	========================================
///@{
/**	@name Skeleton Utilities
*/
//	========================================

/**
	get number of root of the skeleton

	A skeleton may contain several roots, this means
	the structure may contain several independent hierarchy.
	In this function, we count the number of roots.

	\param	skeleton		the skeleton to count root
	\return					root number
*/
template<typename T>
int getSkeletonRootNumber(const Skeleton<T>&	skeleton){
	int root_num = 0;
	for( int i=0; i<skeleton.bone.size(); i++ ){
		if( skeleton.bone[i].parent == -1 )		//	the root node
			root_num ++;
	}
	return root_num;
}

/**
	get the nearest point from a point to skeleton

	\param	bone_id			id of the bone that the nearest point lies
	\param	point_on_bone	nearest point on the bone
	\param	distance		distance from point to nearest point on the skeleton
	\param	skeleton		the skeleton
	\param	point			the point
	\return					whether succeed, 0 if the skeleton don't have bones
*/
template<typename T>
int getNearestPointOnSkeleton(int&					bone_id,
							  Vec3<T>&				point_on_bone,
							  T&					distance,
							  const Skeleton<T>&	skeleton,
							  Vec3<T>				point ){
	bone_id	= -1;
	distance = 1e6;

	std::vector<Vec3<T>> node;
	getSkeletonNodePosition(node, skeleton);

	return getNearestPointOnSkeleton(bone_id, point_on_bone, distance, skeleton, node, point);
}

/**
	get the nearest point from a point to skeleton

	this function assumes that start point of each bone have been given

	\param	bone_id			id of the bone that the nearest point lies
	\param	point_on_bone	nearest point on the bone
	\param	distance		distance from point to nearest point on the skeleton
	\param	skeleton		the skeleton
	\param	node			start point of each bone
	\param	point			the point
	\return					whether succeed, 0 if the skeleton don't have bones
*/
template<typename T>
int getNearestPointOnSkeleton(int&							bone_id,
							  Vec3<T>&						point_on_bone,
							  T&							distance,
							  const Skeleton<T>&			skeleton,
							  const std::vector<Vec3<T>>&	node,
							  Vec3<T>						point ){
	assert(skeleton.bone.size() == node.size());

	bone_id	= -1;
	distance = 1e6;

	//	scan each bone
	for( int j=0; j<skeleton.bone.size(); j++ ){
		Vec3<T> bone_v0 = node[j];
		Vec3<T> bone_v1 = bone_v0 + skeleton.bone[j].bone_vector;

		//	calculate distance to this bone and get nearest point on the bone8
		Vec3<T>	curr_nearest_v;
		geometry::getNearestPointOnSegment(curr_nearest_v, point, bone_v0, bone_v1);
		T dist_to_this_bone = (point - curr_nearest_v).Length();
		if( dist_to_this_bone < distance ){	//	if more near, update
			bone_id			= j;
			point_on_bone	= curr_nearest_v;
			distance		= dist_to_this_bone;
		}
	}

	if( bone_id == -1 )
		return 0;
	else
		return 1;
}

/**
	get start position of each bone on the skeleton

	This function is simply a rename of getSkeletonNodePosition().
	We do this because this Bone Start is more clear than Node,
	but that function comes ealier, so we keep both

	\param	bone_start	return the start position of each bone
	\param	skeleton	the skeleton, must satisfy hierarchy order
*/
template<typename T>
void getSkeletonBoneStartPosition(std::vector<Vec3<T>>& bone_start, const Skeleton<T>& skeleton){
	getSkeletonNodePosition(bone_start, skeleton);
}

/**
	get end position of each bone on the skeleton

	\param	bone_end	return the end position of each bone
	\param	skeleton	the skeleton, must satisfy hierarchy order
*/
template<typename T>
void getSkeletonBoneEndPosition(std::vector<Vec3<T>>& bone_end, const Skeleton<T>& skeleton){
	getSkeletonNodePosition(bone_end, skeleton);

	for( int i=0; i<skeleton.bone.size(); i++ ){
		bone_end[i] += skeleton.bone[i].bone_vector;
	}	
}

/**
	get the position of each node of the skeleton

	the position of each node is the start position of each bone

	\param	node		return the start position of each bone
	\param	skeleton	the skeleton, must satisfy hierarchy order
*/
template<typename T>
void getSkeletonNodePosition(std::vector<Vec3<T>>& node, const Skeleton<T>& skeleton){
	node.resize(skeleton.bone.size());
	for( int i=0; i<skeleton.bone.size(); i++ ){
		if( skeleton.bone[i].parent == -1 )		//	the root node
			node[i] = skeleton.bone[i].offset;
		else							//	parent's end plus bone vector
			node[i] = node[skeleton.bone[i].parent] + skeleton.bone[i].offset;
	}
}

/**
	Get position of each joint on the skeleton

	The skeleton has several limits: \n
	1,	the skeleton contains just one root	\n
	2,	the skeleton must be connected, no space between bones. 
		This means the start of bone is the end of its parent

	the first joint is the root, the following are end position of each bone

	\param	joint		return the position of each joint
	\param	skeleton	the skeleton to calculate joint position
*/
template<typename T>
void getSkeletonJointPosition(std::vector<Vec3<T>>& joint, const Skeleton<T>& skeleton){
	assert( getSkeletonRootNumber(skeleton) == 1 );		//	only one root is allowed
	assert( skeleton.bone[0].parent == -1 );			//	the first bone should be the root

	joint.resize(skeleton.bone.size()+1);
	joint[0] = skeleton.bone[0].offset;
	for( int i=0; i<skeleton.bone.size(); i++ ){	//	parent of root is -1, and already stored in the first joint
		joint[i+1] = joint[skeleton.bone[i].parent+1] + skeleton.bone[i].bone_vector;
	}
}

/**
	set skeleton according to joint position

	The original skeleton only hold the connectivity information. 
	Length of the each bone will be changed.

	The skeleton has several limits: \n
	1,	the skeleton contains just one root	\n
	2,	the skeleton must be connected, no space between bones. 
		This means the start of bone is the end of its parent

	the first joint is the root, the following are end position of each bone

	\param	skeleton	return the skeleton
	\param	joint		the position of each joint
*/
template<typename T>
void setSkeletonByJointPosition(Skeleton<T>& skeleton, const std::vector<Vec3<T>>& joint){
	assert( skeleton.bone.size()+1 == joint.size() );	//	skeleton must match joint
	assert( getSkeletonRootNumber(skeleton) == 1 );		//	only one root is allowed
	assert( skeleton.bone[0].parent == -1 );			//	the first bone should be the root

	//	set bone vector
	for(int i=0; i<skeleton.bone.size(); i++){
		int pid = skeleton.bone[i].parent;	//	if pid == -1, then joint[pid+1] is the root offset
		skeleton.bone[i].bone_vector = joint[i+1] - joint[pid+1];
	}

	//	set bone offset
	for(int i=0; i<skeleton.bone.size(); i++){
		int pid = skeleton.bone[i].parent;
		if( pid == -1 )
			skeleton.bone[i].offset = joint[i];
		else
			skeleton.bone[i].offset = skeleton.bone[pid].bone_vector;
	}

	//	calculate bone front
	setSkeletonBoneFront(skeleton);
}

/**
	set bone front for bone_vector that is already calculated
*/
template<typename T>
void setSkeletonBoneFront(Skeleton<T>& skeleton){
	for( int i=0; i<skeleton.bone.size(); i++ ){
		//	calculate bone_front, since bvh don't have front information, 
		//	we set front to be on the plane of bone_vector and +z axis.
		//	if bone_vector is close to z axis, then change to +y axis
		Vec3<T> axis(0, 0, 1);
		if( skeleton.bone[i].bone_vector.Length() < 1e-6 )	//	the bone is just a temp bone, don't have length
			skeleton.bone[i].bone_front = axis;
		else{
			skeleton.bone[i].bone_front = skeleton.bone[i].bone_vector.Normalize();
			Vec3<T> nor = cross(skeleton.bone[i].bone_front, axis);
			if( nor.Length() < 1e-3 ){
				axis = Vec3<T>(0, 1, 0);
				nor = cross(skeleton.bone[i].bone_front, axis);
			}
			skeleton.bone[i].bone_front.SetRotateDeg(nor, 90);
		}
	}
}
/**
	translate joint position together with its children
*/
template<typename T>
int translateJointWithChildren(std::vector<Vec3<T>>&	joint, 
							   const Skeleton<T>&		skeleton, 
							   int						joint_id,
							   Vec3<T>					trans_vec){
	if( joint_id == 0 ){				//	translate root
		for(int i=0; i<joint.size(); i++)
			joint[i] += trans_vec;
		return joint.size();
	}

	if( joint_id >= joint.size() )		//	translate nothing
		return 0;

	//	only part of joints are moved
	std::vector<int> flag;
	flag.resize(skeleton.bone.size(), 0);
	joint[joint_id] += trans_vec;
	flag[joint_id-1] = 1;
	for(int bid=joint_id-1; bid<skeleton.bone.size(); bid++){
		int pid = skeleton.bone[bid].parent;
		if( flag[pid] ){
			joint[bid+1] += trans_vec;
			flag[bid] = 1;
		}
	}
}
/**
	Transform the skeleton according to transformer

	parameter transformer is a copy of the old transformer, 
	because its value will be modified during calculation

	\param	skeleton		the skeleton to transform
	\param	transformer		the transformer class, select local and global automatically
	\return					whether the transform data is correst, if not, don't touch old data
*/
template<typename T>
int transformSkeleton(Skeleton<T>& skeleton, SkeletonTransformer<T> transformer){
	if( skeleton.bone.size() != transformer.transformer.size() )
		return 0;

	if( !transformer.global_flag ){	//	convert transformer to global
		setTransformerToGlobal(transformer, skeleton);
	}

	for( int i=0; i<skeleton.bone.size(); i++ ){
		int f_id = skeleton.bone[i].parent;

		//	rotation
		skeleton.bone[i].bone_vector = transformer.transformer[i].rotation * skeleton.bone[i].bone_vector;
		skeleton.bone[i].bone_front = transformer.transformer[i].rotation * skeleton.bone[i].bone_front;

		//	offset rotation as parent
		if( f_id != -1 )
			skeleton.bone[i].offset = transformer.transformer[f_id].rotation * skeleton.bone[i].offset;

		//	translation
		skeleton.bone[i].offset += transformer.transformer[i].translation;
	}

	return 1;
}

/**
	Transform the skeleton according to transformer

	transformer is a copy of the old transformer, because its value will be modified during calculation

	\param	skeleton		the skeleton at rest pose
	\param	transformer		the transformer class
	\param	vertex			mesh vertex rigged to the rest pose skeleton
	\param	rigging			rigging information of each vertex to the bone
	\return					whether the transform data is correst, if not, don't touch old data
*/
template<typename T>
int transformSkeletonAndMeshLBS(Skeleton<T>&				skeleton, 
								std::vector<Vec3<T>>&		vertex,
								SkeletonTransformer<T>		transformer,
								const RiggingLBS<T>&		rigging){
	if( skeleton.bone.size() != transformer.transformer.size() )	//	assure skeleton match transformer
		return 0;
	if( vertex.size() != rigging.start.size()-1 )					//	assure mesh match rigging weight
		return 0;

	//	calculate global transformer
	if( !transformer.global_flag )
		setTransformerToGlobal(transformer, skeleton);

	//	transform skeleton
	transformSkeleton(skeleton, transformer);

	//	calculate current bone start position
	std::vector<Vec3<T>>	bone_start;
	getSkeletonBoneStartPosition(bone_start, skeleton);

	//	transform each vertex
#ifdef ENABLE_OPENMP
#	pragma omp parallel for
#endif
	for( int i=0; i<vertex.size(); i++ ){
		Vec3<T> new_pos(0, 0, 0);
		for(int j=rigging.start[i]; j<rigging.start[i+1]; j++){
			int bone_id = rigging.bone_weight[j].bone_id;
			Vec3<T>	r	= rigging.bone_weight[j].offset;
			T	weight	= rigging.bone_weight[j].weight;

			r = transformer.transformer[bone_id].rotation * r;
			new_pos += (bone_start[bone_id] + r) * weight;	//	global position * weight
		}

		vertex[i] = new_pos;
	}

	return 1;
}

/**
	resample skeleton motion data represented using quaternion

	Frame rate of skeleton motion data are fixed. If we want to achieve
	different frame rate, we have to interpolate the data. The best method 
	to interpolate rotation is to use quaternion. Given the target interval
	of each frame, we will give you the resampled data

	\param	target_motion_data		return resampled data in quaternion format
	\param	target_interval			interval of target frames
	\param	source_motion_data		the original motion data in quaternion format
	\return							how many frames generated
*/
template<typename T>
int resampleSkeletonMotion(SkeletonMotionQuaternionData<T>&		target_motion_data,
						   T									target_interval,
						   SkeletonMotionQuaternionData<T>&		source_motion_data){
	//	the source data is not able to resample
	if( source_motion_data.frame.size() < 2 ){
		target_motion_data = source_motion_data;
		return target_motion_data.frame.size();
	}

	int source_size = source_motion_data.frame.size();
	T	total_time	= source_motion_data.frame_interval * (source_size-1);
	int target_size = (total_time / target_interval) + 1;
	target_motion_data.frame_interval = target_interval;
	target_motion_data.frame.clear();
	target_motion_data.frame.resize(target_size);

	SkeletonMotionQuaternion<T>& frame0 = source_motion_data.frame[0];
	int bone_number = frame0.motion.size();
	for(int i=0; i<target_motion_data.frame.size(); i++){
		target_motion_data.frame[i].motion.resize(bone_number);
	}
	Quaternion<T>* source_quaternion = new Quaternion<T>[source_size];
	Quaternion<T>* target_quaternion = new Quaternion<T>[target_size];
	T*	source_data = new T[source_size];
	T*	target_data = new T[target_size];

	//	for each bone
	for(int bone_id=0; bone_id<bone_number; bone_id++){
		//	for each rotation channel
		for(int i=0; i<source_size; i++)
			source_quaternion[i] = source_motion_data.frame[i].motion[bone_id].quaternion;
		resampleSlerp(target_quaternion, target_size, source_quaternion, source_size);
		for(int i=0; i<target_size; i++)
			target_motion_data.frame[i].motion[bone_id].quaternion = target_quaternion[i];

		//	for each translation channel
		for(int j=0; j<3; j++){
			for(int i=0; i<source_size; i++)
				source_data[i] = source_motion_data.frame[i].motion[bone_id].translation[j];
			resampleCatmullRom(target_data, target_size, source_data, source_size);
			for(int i=0; i<target_size; i++)
				target_motion_data.frame[i].motion[bone_id].translation[j] = target_data[i];
		}

		//	scale channel
		for(int i=0; i<source_size; i++)
			source_data[i] = source_motion_data.frame[i].motion[bone_id].scale;
		resampleCatmullRom(target_data, target_size, source_data, source_size);
		for(int i=0; i<target_size; i++)
			target_motion_data.frame[i].motion[bone_id].scale = target_data[i];
	}

	delete[] source_quaternion;
	delete[] target_quaternion;
	delete[] source_data;
	delete[] target_data;

	return target_size;
}


///@}

//	========================================
///@{
/**	@name Morph Utilities
*/
//	========================================

/**
*/



///@}

}}	//	namespace yz::animation

#endif	//	__YZ_SKELETON_UTILS_H__