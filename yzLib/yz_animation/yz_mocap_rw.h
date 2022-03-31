/***********************************************************/
/**	\file
	\brief		motion capture data r & w
	\author		Yizhong Zhang
	\date		6/28/2012
*/
/***********************************************************/
#ifndef __YZ_MOCAP_RW_H__
#define __YZ_MOCAP_RW_H__

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_matrix.h"
#include "yzLib/yz_utils/yz_string_utils.h"
#include "yzLib/yz_animation/yz_skeleton.h"

namespace yz{	namespace animation{

/**
	read mocap data from .bvh file

	Error check of the file is weak, so make sure the format
	
of input file is correct.

	bvh file contain hierarchy and motion data, we can read
	part or both

	since the function don't know the unit used in the file, 
	it is left for user to specify the unit used in the file. 
	the created skeleton use unit of meter

	\param	bvh_file_name		file name of the .bvh file
	\param	skeleton			the read skeleton, if you don't need, set to NULL
	\param	motion_data			frames of skeleton motion, if you don't need, set to NULL
	\param	meter_per_unit		specify the length unit used in the file, how many meters per unit
	\param	angle_unit			unit used for angle used in the file, 'd'/'D' Degree ; 'r'/'R' Rad
	\return						whether read succeed
*/
template<typename T>
int readMocapFromBVH(const char*			bvh_file_name, 
					 Skeleton<T>*			skeleton		= NULL,
					 SkeletonMotionData<T>*	motion_data		= NULL,
					 double					meter_per_unit	= 1,
					 char					angle_unit		= 'd'){
	std::ifstream bvh( bvh_file_name );
	if(!bvh.is_open()){
		#ifndef BE_QUIET
		std::cout << "error: readMocapFromBVH, cannot open " << bvh_file_name << std::endl;
		#endif
		return 0;
	}

	readMocapFromBVH(bvh, skeleton, motion_data, meter_per_unit, angle_unit);

	bvh.close();
	return 1;
}

/**
	read mocap data from stream

	Error check of the file is weak, so make sure the format
	of input file is correct.

	bvh file contain hierarchy and motion data, we can read
	part or both

	since the function don't know the unit used in the file, 
	it is left for user to specify the unit used in the file. 
	the created skeleton use unit of meter

	\param	bvh					stream of the .bvh, can be either kind of stream
	\param	skeleton			the read skeleton, if you don't need, set to NULL
	\param	motion_data			frames of skeleton motion, if you don't need, set to NULL
	\param	meter_per_unit		specify the length unit used in the file, how many meters per unit
	\param	angle_unit			unit used for angle used in the file, 'd'/'D' Degree ; 'r'/'R' Rad
	\return						whether read succeed
*/
template<typename T>
int readMocapFromBVH(std::istream&			bvh,
					 Skeleton<T>*			skeleton		= NULL,
					 SkeletonMotionData<T>*	motion_data		= NULL,
					 double					meter_per_unit	= 1,
					 char					angle_unit		= 'd'){
	if( skeleton==NULL && motion_data==NULL )	//	no output
		return 0;

	//	setup tmp array
	Skeleton<T> tmp_skeleton;
	std::vector<int2> motion_list;
	enum{
		Xrotation,
		Yrotation,
		Zrotation,
		Xposition,
		Yposition,
		Zposition
	};

	//	parse the file
	enum{
		NOT_SET,
		HIERARCHY,
		MOTION
	};	
	int status = 0;
	std::string cmd, skip;
	std::stack<int> parent_stack;
	while( !bvh.eof() ){
		bvh >> cmd;

		if( status == NOT_SET ){	//	no explicit status, accept everything
			if( cmd == "HIERARCHY" ){
				status = HIERARCHY;
			}
			else if( cmd == "MOTION" ){
				status = MOTION;
				break;
			}
		}
		else if( status == HIERARCHY ){
			if( cmd=="ROOT" || cmd=="JOINT" ){	//	a new bone
				Bone<T> bone;
				bone.bone_vector = Vec3<T>(0, 0, 0);
				std::getline(bvh, bone.name);	//	get root name, space in name is accepted
				utils::removePrefixSuffixBlankCharacters(bone.name);
				std::getline(bvh, skip);		//	skip "{"
				bvh >> skip;					//	skip "OFFSET"
				T x, y, z;
				bvh >> x >> y >> z;				//	read offset
				if( cmd=="ROOT" ){
					bone.offset = Vec3<T>(x, y, z);
					bone.parent	= -1;
				}
				else{
					if( parent_stack.empty() ){
						#ifndef BE_QUIET
							std::cout << "hierarchy error" << std::endl;
						#endif
						return 0;
					}
					bone.offset = Vec3<T>(x, y, z);
					bone.parent = parent_stack.top();;
				}
				parent_stack.push( tmp_skeleton.bone.size() );	//	push new added bone to parent stack
				tmp_skeleton.bone.push_back(bone);				//	so the top of the stack is the parent of new added bone
			}
			else if( cmd=="End" ){
				std::getline(bvh, skip);	//	skip "Site"
				std::getline(bvh, skip);	//	skip "{"
				bvh >> skip;				//	skip "OFFSET"
				T x, y, z;
				bvh >> x >> y >> z;			//	read offset
				if( parent_stack.empty() ){
					#ifndef BE_QUIET
						std::cout << " hierarchy error" << std::endl;
					#endif
					return 0;
				}
				int i = parent_stack.top();
				tmp_skeleton.bone[i].bone_vector = Vec3<T>(x, y, z);
				bvh >> skip;				//	skip "}"
			}
			else if( cmd=="CHANNELS" ){
				int num;
				bvh >> num;		//	get channels
				int2 motion;
				motion.x = tmp_skeleton.bone.size() - 1;
				for( int i=0; i<num; i++ ){
					bvh >> cmd;
					if( cmd == "Xposition" )
						motion.y = Xposition;
					else if( cmd == "Yposition" )
						motion.y = Yposition;
					else if( cmd == "Zposition" )
						motion.y = Zposition;
					else if( cmd == "Xrotation" )
						motion.y = Xrotation;
					else if( cmd == "Yrotation" )
						motion.y = Yrotation;
					else if( cmd == "Zrotation" )
						motion.y = Zrotation;
					else{
						#ifndef BE_QUIET
							std::cout << "invalid CHANNEL key word: " << cmd << std::endl;
						#endif
						return 0;
					}
					motion_list.push_back(motion);
				}
			}
			else if( cmd=="}" ){
				if( parent_stack.empty() ){
					#ifndef BE_QUIET
						std::cout << "hierarchy error" << std::endl;
					#endif
					return 0;
				}
				parent_stack.pop();
				if( parent_stack.empty() )	//	we have parsed the whole ROOT
					status = 0;
			}
			else{
				#ifndef BE_QUIET
					std::cout << "unexpedted command : " << cmd << std::endl;
				#endif
				return 0;
			}
		}
	}

	//	calculate bone_vector according to offset
	std::vector<int> child_count;
	child_count.resize(tmp_skeleton.bone.size(), 0);
	for( int i=0; i<tmp_skeleton.bone.size(); i++ ){
		int p_id = tmp_skeleton.bone[i].parent;
		if( p_id == -1 )	//	root
			continue;
		tmp_skeleton.bone[p_id].bone_vector += tmp_skeleton.bone[i].offset;
		child_count[p_id] ++;
	}
	for( int i=0; i<tmp_skeleton.bone.size(); i++ ){
		if( child_count[i] ){
			tmp_skeleton.bone[i].bone_vector *= (1.0/child_count[i]);
		}
	}

	//	if output frame data
	T scale = meter_per_unit;
	if( motion_data ){
		motion_data->frame.clear();

		T num;
		int frames;
		bvh >> skip;		//	skip "Frames:"
		bvh >> frames;
		bvh >> skip;		//	skip "Frame"
		bvh >> skip;		//	skip "Time:"
		bvh >> motion_data->frame_interval;

		motion_data->frame.resize(frames);
		for( int i=0; i<frames; i++ ){		//	read frame data
			motion_data->frame[i].motion.resize(tmp_skeleton.bone.size());
			for( int j=0; j<motion_list.size(); j++ ){
				bvh >> num;
				int bone_id = motion_list[j].x;
				switch(motion_list[j].y){
					case Xrotation:
					case Yrotation:
					case Zrotation:
						if( angle_unit=='d' || angle_unit=='D' )	num = deg2rad(num);
						if( !motion_data->frame[i].motion[bone_id].AddRotation(motion_list[j].y, num) ){
							#ifndef BE_QUIET
								std::cout << "duplicate rotation of bone " << bone_id << std::endl;
							#endif
							return 0;
						}
						break;
					case Xposition:
						motion_data->frame[i].motion[bone_id].translation.x = num * scale;
						break;
					case Yposition:
						motion_data->frame[i].motion[bone_id].translation.y = num * scale;
						break;
					case Zposition:
						motion_data->frame[i].motion[bone_id].translation.z = num * scale;
						break;
				}
			}
		}
	}

	//	if output bone array
	if( skeleton ){
		for( int i=0; i<tmp_skeleton.bone.size(); i++ ){
			tmp_skeleton.bone[i].bone_vector *= scale;
			tmp_skeleton.bone[i].offset *= scale;
			//	calculate bone_front, since bvh don't have front information, 
			//	we set front to be on the plane of bone_vector and +z axis.
			//	if bone_vector is close to z axis, then change to +y axis
			Vec3<T> axis(0, 0, 1);
			if( tmp_skeleton.bone[i].bone_vector.Length() < 1e-6 )	//	the bone is just a temp bone, don't have length
				tmp_skeleton.bone[i].bone_front = axis;
			else{
				tmp_skeleton.bone[i].bone_front = tmp_skeleton.bone[i].bone_vector.Normalize();
				Vec3<T> nor = cross(tmp_skeleton.bone[i].bone_front, axis);
				if( nor.Length() < 1e-3 ){
					axis = Vec3<T>(0, 1, 0);
					nor = cross(tmp_skeleton.bone[i].bone_front, axis);
				}
				tmp_skeleton.bone[i].bone_front.SetRotateDeg(nor, 90);
			}
		}
		*skeleton = tmp_skeleton;
	}

	return 1;
}

/**
	\todo	implement this function
*/
template<typename T>
int writeMocapToBVH(std::ostream&			bvh,
					Skeleton<T>*			skeleton		= NULL,
					SkeletonMotionData<T>*	motion_data		= NULL,
					double					meter_per_unit	= 1,
					char					angle_unit		= 'd'){

}
/**
	write skeleton to .bvh file

	.bvh file contail both hierarchy and motion, this function 
	exports hierarchy only. So we use dafault motion sequence

	\param	bvh_file_name		file name of the .bvh file
	\param	skeleton			the read skeleton, if you don't need, set to NULL
	\param	meter_per_unit		specify the length unit used in the file, how many meters per unit
	\return						whether read succeed
*/
template<typename T>
int writeSkeletonToBVH(	const char*		bvh_file_name, 
						Skeleton<T>&	skeleton,
						double			meter_per_unit	= 1){
	//	open file
	std::ofstream bvh( bvh_file_name );
	if(!bvh.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: writeSkeletonToBVH, cannot open " << bvh_file_name << std::endl;
		#endif
		return 0;
	}

	//	write the skeleton hierarchy to the file
	writeSkeletonToBVH(bvh, skeleton, meter_per_unit);

	//	write motion to the file, since this file has only skeleton, so there are ZERO frames
	bvh << "MOTION" << std::endl;
	bvh << "Frames: 0" << std::endl;
	bvh << "Frame Time: 0.033333" << std::endl;

	bvh.close();
}

/**
	write skeleton to bvh stream

	.bvh file contail both hierarchy and motion, this function 
	exports hierarchy only. So we use dafault motion sequence, RotZ*RotY*RotX

	\param	bvh					file stream 
	\param	skeleton			the read skeleton, if you don't need, set to NULL
	\param	meter_per_unit		specify the length unit used in the file, how many meters per unit
	\return						whether read succeed
*/
template<typename T>
int writeSkeletonToBVH(	std::ostream&	bvh, 
						Skeleton<T>&	skeleton,
						double			meter_per_unit	= 1){
	bvh << "HIERARCHY" << std::endl;

	double scale = 1.0 / meter_per_unit;

	//	parse the skeleton
	std::stack<int> parent_stack;
	for(int bone_id=0; bone_id<skeleton.bone.size(); bone_id++){
		int parent_id = skeleton.bone[bone_id].parent;
		if( parent_id == -1 ){	//	root node
			bvh << "ROOT\t" << skeleton.bone[bone_id].name << std::endl;
			bvh << "{" << std::endl;
			parent_stack.push(bone_id);
			for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
			bvh << "OFFSET\t" << skeleton.bone[bone_id].offset[0] * scale << '\t'
				<< skeleton.bone[bone_id].offset[1] * scale << '\t'
				<< skeleton.bone[bone_id].offset[2] * scale << std::endl;
			for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
			bvh << "CHANNELS 6 Xposition Yposition Zposition Zrotation Yrotation Xrotation" << std::endl;
		}
		else{
			assert(!parent_stack.empty());
			if(parent_id != parent_stack.top()){	//	reached a end of the hierarchy
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "End Site" << std::endl;
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "{" << std::endl;
				for(int i=0; i<parent_stack.size()+1; i++)	bvh << "\t";
				bvh << "OFFSET\t" << skeleton.bone[parent_stack.top()].bone_vector[0] * scale << '\t'
					<< skeleton.bone[parent_stack.top()].bone_vector[1] * scale << '\t'
					<< skeleton.bone[parent_stack.top()].bone_vector[2] * scale << std::endl;
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "}" << std::endl;

				while(!parent_stack.empty() && parent_stack.top()!=parent_id){
					parent_stack.pop();
					for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
					bvh << "}" << std::endl;
				}
			}

			if(!parent_stack.empty() && parent_stack.top()==parent_id){	//	move hierarchy further
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "JOINT\t" << skeleton.bone[bone_id].name << std::endl;
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "{" << std::endl;
				parent_stack.push(bone_id);
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "OFFSET\t" << skeleton.bone[bone_id].offset[0] * scale <<	 '\t'
					<< skeleton.bone[bone_id].offset[1] * scale << '\t'
					<< skeleton.bone[bone_id].offset[2] * scale << std::endl;
				for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
				bvh << "CHANNELS 3 Zrotation Yrotation Xrotation" << std::endl;
			}
		}
	}

	//	sweep the end
	for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
	bvh << "End Site" << std::endl;
	for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
	bvh << "{" << std::endl;
	for(int i=0; i<parent_stack.size()+1; i++)	bvh << "\t";
	bvh << "OFFSET\t" << skeleton.bone[parent_stack.top()].bone_vector[0] * scale << '\t'
		<< skeleton.bone[parent_stack.top()].bone_vector[1] * scale << '\t'
		<< skeleton.bone[parent_stack.top()].bone_vector[2] * scale << std::endl;
	for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
	bvh << "}" << std::endl;

	while(!parent_stack.empty()){
		parent_stack.pop();
		for(int i=0; i<parent_stack.size(); i++)	bvh << "\t";
		bvh << "}" << std::endl;
	}

	return 1;
}


}}	//	namespace yz::animation

#endif	//	__YZ_MOCAP_RW_H__
