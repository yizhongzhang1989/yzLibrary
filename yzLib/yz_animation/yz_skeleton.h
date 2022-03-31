/***********************************************************/
/**	\file
	\brief		skeleton and related data structure
	\author		Yizhong Zhang
	\date		6/29/2012
*/
/***********************************************************/
#ifndef __YZ_SKELETON_H__
#define __YZ_SKELETON_H__

#include <vector>
#include <string>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_matrix.h"
#include "yzLib/yz_math/yz_quaternion.h"
#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_animation_opengl_utils.h"
#endif

namespace yz{	namespace animation{

/**
	A single bone

	A single bone has 6 DOFs. Bones in a well defined skeleton 
	should connect each other without interval. However, interval
	is also acceptable. So offset is useful if the skeleton is not
	connected.

	bone_vector is the direction of the bone, while bone_front is
	the orientation of the bone. When transform the bone, both bone_vector
	and bone_front are transformed(rotate and scale)
*/
template<class T>
class Bone{
public:
	std::string	name;			///<	name of the bone, expected to be unique for each bone
	int			parent;			///<	parent of this bone, -1 for the root
	Vec3<T>		offset;			///<	offset of the bone's origin compared to its parent's origin
	Vec3<T>		bone_vector;	///<	vector from bone origin to bone end, indicating the bone length and direction
	Vec3<T>		bone_front;		///<	bone orientation, should be perpendicular to bone_vector
};

/**
	Skeleton, a hierarchy of bones

	skeleton is simply an array of bones, connected together by
	parent relationship. bone_vector of root (parent==-1) is the
	vector of original point to the real-root of the skeleton
*/
template<class T>
class Skeleton{
public:
	std::vector<Bone<T>>	bone;

public:
	/**
		read skeleton from .bvh file

		\param	bvh_file_name	file name of the .bvh file
		\param	meter_per_unit	actual length in meters of one unit in .bvh file, also scale value
		\return					whether read file succeed
	*/
	inline int ReadSkeletonFromBVH(const char* bvh_file_name, double meter_per_unit=1){
		//return readMocapFromBVH(bvh_file_name, this, (SkeletonMotionData<T>*)NULL, meter_per_unit, 'd');
		return readMocapFromBVH(bvh_file_name, this, NULL, meter_per_unit, 'd');
	}

	/**
		Draw Skeleton

		if OpenGL is not enabled, print error information

		\param	display_mode	0: stick,  1: solid bones,  2: wire bones
		\param	stick_radius	radius of stick if display_mode is 0(stick)
		\return					whether display succeed
	*/
	inline int Display(int display_mode=0, float stick_radius=0.01f) const{
		#ifdef YZ_gl_h
			if( display_mode == 0 ){
				opengl::drawSkeletonAsSticks(*this, stick_radius);
			}
			else if( display_mode == 1 ){
				//glColor3f(0, 1, 1);
				opengl::drawSkeletonAsSolidBones(*this);
			}
			else if( display_mode == 2 ){
				//glColor3f(0, 1, 1);
				opengl::drawSkeletonAsWireBones(*this);
			}
			else if( display_mode == 3 ){
				//glColor3f(0, 1, 1);
				opengl::drawSkeletonAsFlatSolidBones(*this);
			}
			else if( display_mode == 4 ){
				opengl::drawSkeletonAsColorCubeBones(*this);
			}
			return 1;
		#else
			std::cout << "gl.h has to be included in order to use Display() in Skeleton" << std::endl;
			return 0;
		#endif
	}
};

/**
	direct motion data

	rotation angle and sequence, translation, scale

	Rotation saves the angle of rotation around x,y,z axis,
	while the sequence is saved in member: sequence. \n
	sequence has 8 bits: \n
	7 6 | 5 4 | 3 2 | 1 0 \n
	(7 6) indicate the degree of freedoms\n
	(5 4) indicate the order of z rotation in the sequence, 
	for example, if (5 4)=1, Rot_Z appear left most in the matrix. \n
	Similar to (3 2) and (1 0), which is Rot_Y and Rot_X \n
	Example: (1 1 | 1 0 | 1 1 | 0 1) indicate the transform matrix
	is cross of three rotation matrixes, from left to right, they are
	rot_x, rot_z, rot_y, so the matrix is (Rot_X x Rot_Z x Rot_Y)

	CAUTION!!! sequence order is the order appear in the matrix from left, 
	not order of rotation take place. (should be inverse)

	member function are provided to deal with sequence
*/
template<class T>
class BoneMotion{
public:
	Vec3<T>			rotation;	///<	Euler angles of each axis
	unsigned char	sequence;	///<	record the sequence of the rotation
	Vec3<T>			translation;///<	translation vector
	T				scale;		///<	scaler of the bone vector

public:
	BoneMotion() : sequence(0), scale(1){}

	/**
		\param	axis_id			0: x,  1: y,  2: z
		\param	angle_in_rad	angle of rotation
		\return					whether add succeed
	*/
	inline int AddRotation(int axis_id, T angle_in_rad){
		if( axis_id!=0 && axis_id!=1 && axis_id!=2 )
			return 0;
		int rot_count = GetSequenceBlockValue(3);
		if(rot_count == 3)	return 0;
		int guard = GetSequenceBlockValue(axis_id);
		if(guard)	return 0;
		rot_count ++;
		SetSequenceBlockValue(axis_id, rot_count);
		SetSequenceBlockValue(3, rot_count);
		if(axis_id == 0)		rotation.x = angle_in_rad;
		else if(axis_id == 1)	rotation.y = angle_in_rad;
		else					rotation.z = angle_in_rad;
		return 1;
	}

	/**
		get rotation matrix, exclude scale
	*/
	inline Matrix3x3<T> GetRotationMatrix() const{
		Matrix3x3<T>	mat, rot;
		mat.SetIdentity();
		int rot_count = GetSequenceBlockValue(3);
		for( int i=1; i<=rot_count; i++ ){
			if(		i == GetSequenceBlockValue(0) )
				rot.SetRotationRad(Vec3i(1, 0, 0), rotation.x);
			else if(i == GetSequenceBlockValue(1) )
				rot.SetRotationRad(Vec3i(0, 1, 0), rotation.y);
			else
				rot.SetRotationRad(Vec3i(0, 0, 1), rotation.z);
			mat *= rot;
		}

		return mat;
	}

	/**
		This function return rotation matrix that treate each eular angle as local rotation

		Local rotation means when one rotation happen, the rotation axis of other rotations rotate, too.
		In this case, we have to calculate rotation from right to left and calculate current rotation
		axis every time a rotation happens.
	*/
	inline Matrix3x3<T> GetRotationMatrixInLocalCoord() const{
		Matrix3x3<T>	mat, rot;
		mat.SetIdentity();
		int rot_count = GetSequenceBlockValue(3);
		Vec3<T> rot_axis[3] = { int3(1,0,0), int3(0,1,0), int3(0,0,1) };
		for( int i=rot_count; i>=1; i-- ){
			if(		i == GetSequenceBlockValue(0) )
				rot.SetRotationRad(rot_axis[0], rotation.x);
			else if(i == GetSequenceBlockValue(1) )
				rot.SetRotationRad(rot_axis[1], rotation.y);
			else
				rot.SetRotationRad(rot_axis[2], rotation.z);

			//	rotate axis
			rot_axis[0] = rot * rot_axis[0];
			rot_axis[1] = rot * rot_axis[1];
			rot_axis[2] = rot * rot_axis[2];

			mat = rot * mat;
		}

		return mat;
	}


	/**
		transform the rotation sequence to transform matrix

		transform matrix include rotation and scale
	*/
	inline Matrix3x3<T> GetTransformMatrix() const{
		Matrix3x3<T>	mat = GetRotationMatrix();

		mat.data[0][0] *= scale;
		mat.data[1][1] *= scale;
		mat.data[2][2] *= scale;

		return mat;
	}

	/**
		get quaternion from eular angles
	*/
	inline Quaternion<T> GetQuaternion() const{
		Matrix3x3<T>	mat = GetRotationMatrix();
		Quaternion<T>	q(mat);
		return q;
	}

private:
	/**
		get value from sequence by block id

		block id should be 0 ~ 3, block sequence is 
		(3 | 2 | 1 | 0) with two bits each
	*/
	inline int GetSequenceBlockValue(int block_id) const{
		unsigned char mask = (0x03 << (block_id*2) );
		return (sequence & mask) >> (block_id*2);
	}

	/**
		set value to sequence by block id

		block_id and value should be 0 ~ 3, block sequence is 
		(3 | 2 | 1 | 0) with two bits each
	*/
	inline void SetSequenceBlockValue(int block_id, unsigned char value){
		unsigned char mask = (0x03 << (block_id*2) );
		value = (value << (block_id*2)) & mask;
		sequence &= ~mask;
		sequence |= value;
	}
};

/**
	bone motion data represented by quaternion
*/
template<class T>
class BoneMotionQuaternion{
public:
	Quaternion<T>	quaternion;	///<	quaternion representing rotation
	Vec3<T>			translation;///<	translation vector
	T				scale;		///<	scaler of the bone vector

public:
	BoneMotionQuaternion() : scale(1){}

	/**
		get transform matrix from quaternion and scale

		transform matrix include rotation and scale
	*/
	inline Matrix3x3<T> GetTransformMatrix() const{
		Matrix3x3<T>	mat(quaternion);

		mat.data[0][0] *= scale;
		mat.data[1][1] *= scale;
		mat.data[2][2] *= scale;

		return mat;
	}
	/**
		set by bone motion
	*/
	template<typename TYPE>
	void Set(const BoneMotion<TYPE>& bone_motion){
		quaternion	= bone_motion.GetQuaternion();
		translation	= bone_motion.translation;
		scale		= bone_motion.scale;
	}
};

/**
	motion of an skeleton
*/
template<class T>
class SkeletonMotion{
public:
	std::vector<BoneMotion<T>> motion;	///<	motion array of skeleton
};

/**
	motion of an skeleton by quaternion
*/
template<class T>
class SkeletonMotionQuaternion{
public:
	std::vector<BoneMotionQuaternion<T>> motion;	///<	motion array of skeleton

public:
	/**
		set quaternion motion by direct motion
	*/
	template<typename TYPE>
	void Set(const SkeletonMotion<TYPE>& skeleton_motion){
		motion.resize( skeleton_motion.motion.size() );
		for(int i=0; i<motion.size(); i++)
			motion[i].Set( skeleton_motion.motion[i] );
	}
};

/**
	skeleton motion of several frames
*/
template<class T>
class SkeletonMotionData{
public:
	std::vector<SkeletonMotion<T>>	frame;
	T								frame_interval;
};

/**
	skeleton motion of several frames by quaternion
*/
template<class T>
class SkeletonMotionQuaternionData{
public:
	std::vector<SkeletonMotionQuaternion<T>>	frame;
	T											frame_interval;

public:
	/**
		Set quaternin frame data by direct motion data
	*/
	template<typename TYPE>
	void Set(const SkeletonMotionData<TYPE>& data){
		frame_interval = data.frame_interval;
		frame.resize( data.frame.size() );
		for(int i=0; i<frame.size(); i++)
			frame[i].Set( data.frame[i] );

		//	we check each quaternion sequence to eleminate jumping
		if( frame.empty() )
			return;
		int bone_number = frame[0].motion.size();
		for(int b_id = 0; b_id < bone_number; b_id++){	//	for each bone
			for(int i=1; i<frame.size(); i++){			//	for each frame after the first frame
				//	choose quaternion that is closer than previous frame
				frame[i].motion[b_id].quaternion.SetCloser( frame[i-1].motion[b_id].quaternion );
			}
		}
	}
};

/**
	The transform applied on a bone

	For each bone, there is a corresponding bone transformation, 
	represented by a rotation matrix (with scale), a translation 
	vector . To apply the transform to the bone, first rotate 
	the bone_vector by rotation, then add the translation vector 
	to the bone_vector. bone_front is only rotated by rotation.

	since bones in skeleton are constraint by each other, translation
	should only apply on root.
*/
template<class T>
class BoneTransformer{
public:
	Matrix3x3<T>	rotation;		///<	rotation matrix with scale
	Vec3<T>			translation;	///<	translation vector
public:
	/**
		constructor, set rotation matrix to be identity matrix and scale to 1
	*/
	BoneTransformer(){
		rotation.SetIdentity();
	}

	/**
		set the transformer
	*/
	void SetTransformerFromMotion(const BoneMotion<T>& bone_motion){
		rotation	= bone_motion.GetTransformMatrix();
		translation	= bone_motion.translation;
	}

	/**
		set the transformer from bone motion in local coordinate

		local bone motion means eular angle are used for local coordinate 
	*/
	void SetTransformerFromMotionInLocalCoord(const BoneMotion<T>& bone_motion_local){
		rotation	= bone_motion_local.GetRotationMatrixInLocalCoord();
		translation	= bone_motion_local.translation;
	}

	/**
		set the transformer
	*/
	void Set(const BoneMotionQuaternion<T>& bone_motion_quaternion){
		rotation	= bone_motion_quaternion.GetTransformMatrix();
		translation	= bone_motion_quaternion.translation;
	}
};

/**
	Transformers of an array of bones

	for a given bone array, the transformer array contain
	an array of transformer with each transformer apply on
	corresponding bone
*/
template<class T>
class SkeletonTransformer{
public:
	std::vector<BoneTransformer<T>>	transformer;	///<	transformer array of skeleton
	int	global_flag;								///<	whether this transformer is global transformer;

public:
	//	constructor
	SkeletonTransformer() : global_flag(0){}

	/**
		set the transformer

		skeleton motion data are local, so transformer is local
	*/
	void SetTransformerFromMotion(const SkeletonMotion<T>& skeleton_motion){
		transformer.resize( skeleton_motion.motion.size() );
		for( int i=0; i<transformer.size(); i++ )
			transformer[i].SetTransformerFromMotion( skeleton_motion.motion[i] );
		global_flag = 0;
	}

	/**
		set the transformer using data that motion are represented in local coordinate

		skeleton motion data are local, so transformer is local. This local is not the same
		as local coordinate
	*/
	void SetTransformerFromMotionInLocalCoord(const SkeletonMotion<T>& skeleton_motion){
		transformer.resize( skeleton_motion.motion.size() );
		for( int i=0; i<transformer.size(); i++ )
			transformer[i].SetTransformerFromMotionInLocalCoord( skeleton_motion.motion[i] );
		global_flag = 0;
	}

	/**
		set the transformer
	*/
	void Set(const SkeletonMotionQuaternion<T>& skeleton_motion_quaternion){
		transformer.resize( skeleton_motion_quaternion.motion.size() );
		for( int i=0; i<transformer.size(); i++ )
			transformer[i].Set( skeleton_motion_quaternion.motion[i] );
		global_flag = 0;
	}
};

/**
	bone array transformer of a lot of frames
*/
template<class T>
class SkeletonTransformerData{
public:
	std::vector<SkeletonTransformer<T>>	frame;
	T									frame_interval;

public:
	/**
		how many frames are contained
	*/
	int FrameNumber(){
		return frame.size();
	}

	/**
		set the transformer
	*/
	void SetTransformerFromMotion(const SkeletonMotionData<T>& motion_data){
		frame_interval = motion_data.frame_interval;
		frame.resize( motion_data.frame.size() );
		for( int i=0; i<frame.size(); i++ )
			frame[i].SetTransformerFromMotion( motion_data.frame[i] );
	}

	/**
		set the transformer data
	*/
	void Set(const SkeletonMotionQuaternionData<T>& motion_data){
		frame_interval = motion_data.frame_interval;
		frame.resize( motion_data.frame.size() );
		for( int i=0; i<frame.size(); i++ )
			frame[i].Set( motion_data.frame[i] );
	}
};


}}	//	namespace yz::animation

#endif	//	__YZ_SKELETON_H__