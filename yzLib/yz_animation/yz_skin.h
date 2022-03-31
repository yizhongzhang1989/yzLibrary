/***********************************************************/
/**	\file
	\brief		Skin of Character
	\author		Yizhong Zhang
	\date		9/21/2012
*/
/***********************************************************/
#ifndef __YZ_SKIN_H__
#define __YZ_SKIN_H__

#include <vector>
#include <fstream>
#include <string>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_reorder.h"

namespace yz{	namespace animation{

/**
	rigging information of a single vertex
*/
template<class T>
class BoneWeightLBS{
public:
	int		bone_id;	///<	id of the binded bone 
	Vec3<T>	offset;		///<	offset to the start position of the bone
	T		weight;		///<	weight to the bone

	BoneWeightLBS(int bone=0, T weight_of_bone=0) : bone_id(bone), weight(weight_of_bone) {}
};

/**
	Rigging information of a mesh

	bone id, offset vector and weight of each vertex are stored in 
	bone_weight	vector, vertex by vertex. It is possible that 
	a single vertex is binded to several bones, then 
	start store the start position and end position of 
	this vertex in bone_weight array

	if we want bone weight of vertex i, lookup bone_weight with
	index from start[i] to start[i+1]
*/
template<class T>
class RiggingLBS{
public:
	std::vector<BoneWeightLBS<T>>	bone_weight;	///<	bone, offset and corresponding weight of each vertex
	std::vector<int>				start;			///<	bone weight pair start position in the list of each vertex

public:
	/**
		clear data
	*/
	void Reset(){
		bone_weight.clear();
		start.clear();
	}

	/**
		Read rigging data from file
	*/
	int ReadRiggingFromFile(const char* file_name){
		std::ifstream file( file_name );	//	load obj failed, do not touch old data
		if(!file.is_open()){
			std::cout << "cannot open " << file_name << std::endl;
			return 0;
		}

		//	read start
		int size;
		file >> size;
		start.resize(size);
		for(int i=0; i<size; i++)
			file >> start[i];

		//	read bone weight
		file >> size;
		bone_weight.resize(size);
		for(int i=0; i<size; i++){
			file >> bone_weight[i].bone_id >> bone_weight[i].weight 
				>> bone_weight[i].offset[0] >> bone_weight[i].offset[1] >> bone_weight[i].offset[2];
		}

		file.close();
		return 1;
	}

	/**
		Write rigging data to file
	*/
	int WriteRiggingToFile(const char* file_name){
		std::ofstream file( file_name );
		if(!file.is_open()){
			#ifndef BE_QUIET
				std::cout << "error: RiggingLBS::WriteRiggingToFile, cannot open " << file_name << std::endl;
			#endif
			return 0;
		}

		//	write start
		file << start.size() << std::endl;
		for(int i=0; i<start.size(); i++)
			file << start[i] << std::endl;
		file << std::endl;

		//	write bone weight
		file << bone_weight.size() << std::endl;
		for(int i=0; i<bone_weight.size(); i++){
			file << bone_weight[i].bone_id << "\t" << bone_weight[i].weight << "\t" 
				<< bone_weight[i].offset[0] << "\t" << bone_weight[i].offset[1] << "\t" 
				<< bone_weight[i].offset[2] << std::endl;
		}

		file.close();
		return 1;
	}

	/**
		limit the bone number that a vertex can rig to

		it is possible for a vertex to rig to a lot of bones, 
		by calling this function, we can reduce the number
		of bones each vertex is rigged to

		\param	max_bone_number		the target number of bones
		\return						the remaining bone_weight number
	*/
	int SetMaxBoneNumberOfEachVertex(int max_bone_number){
		std::vector<BoneWeightLBS<T>>	tmp_bone_weight;
		std::vector<int>				tmp_start;

		int vertex_number = start.size() - 1;
	
		tmp_start.push_back(0);

		//	scan each vertex
		for(int i=0; i<vertex_number; i++){
			int bone_number = start[i+1] - start[i];
			if( bone_number <= max_bone_number ){	//	bone number is no more than max
				for(int j=start[i]; j<start[i+1]; j++)
					tmp_bone_weight.push_back(bone_weight[j]);
				tmp_start.push_back(tmp_bone_weight.size());
			}
			else{	//	sort bone weight and save bigger weights only
				std::vector<T>		weight;
				std::vector<int>	index;
				weight.resize(bone_number);
				for(int j=start[i], k=0; j<start[i+1]; j++, k++){
					weight[k] = bone_weight[j].weight;
				}
				utils::sortIndex(index, weight.begin(), weight.end(), std::greater<T>());
				//	normalize the biggest several weights
				T weight_sum = 0;
				for(int j=0; j<max_bone_number; j++)
					weight_sum += weight[index[j]];
				//	add biggest weights
				int bone_start_id = start[i];
				for(int j=0; j<max_bone_number; j++){
					tmp_bone_weight.push_back(bone_weight[index[j]+bone_start_id]);
					tmp_bone_weight.back().weight /= weight_sum;
				}
				tmp_start.push_back(tmp_bone_weight.size());
			}
		}

		//	set data
		bone_weight = tmp_bone_weight;
		start		= tmp_start;

		return bone_weight.size();
	}
	/**
		reduce number of bones by set threshold

		It is possible for a vertex to rig to a lot of bones, 
		by calling this function, we can reduce the number
		of bones each vertex is rigged to. 
		If the threshold is 0.97, we only save bones that 
		make up 97% of the total weight

		It is recommended to call this function only once after
		the calculation of bone heat, because the second time 
		to call this function, the base number of weight has changed

		\param	threshold		threshold value that make up of the total weight, legal value is between 0 and 1
		\return					the remaining bone_weight number
	*/
	int ReduceBoneNumberByThreshold(T threshold = 0.97){
		if( threshold >= 1 || threshold <= 0 )	//	illegal value, do nothing on the data
			return bone_weight.size();

		std::vector<BoneWeightLBS<T>>	tmp_bone_weight;
		std::vector<int>				tmp_start;

		int vertex_number = start.size() - 1;
	
		tmp_start.push_back(0);

		//	scan each vertex
		for(int i=0; i<vertex_number; i++){
			int bone_number = start[i+1] - start[i];

			//	sort bone weight and save bigger weights only
			std::vector<T>		weight;
			std::vector<int>	index;
			weight.resize(bone_number);
			for(int j=start[i], k=0; j<start[i+1]; j++, k++){
				weight[k] = bone_weight[j].weight;
			}
			utils::sortIndex(index, weight.begin(), weight.end(), std::greater<T>());
			//	normalize the biggest several weights
			T weight_sum = 0;
			int save_bone_number = 0;
			do{
				weight_sum += weight[index[save_bone_number]];
				save_bone_number ++;
			}while(weight_sum < threshold && save_bone_number<bone_number);
			//	add biggest weights
			int bone_start_id = start[i];
			for(int j=0; j<save_bone_number; j++){
				tmp_bone_weight.push_back( bone_weight[index[j]+bone_start_id] );
				tmp_bone_weight.back().weight /= weight_sum;
			}
			tmp_start.push_back(tmp_bone_weight.size());
		}

		//	set data
		bone_weight = tmp_bone_weight;
		start		= tmp_start;

		return bone_weight.size();
	}

	/**
		retarget weight of bones

		Some time we create temp bones for rigging, but no longer need them
		after weight are calculated. This function is designed for remove bones
		that is no longer needed. We can add the weight to other bones.

		It is highly recommended to create temp bones at the end of the skeleton,
		so that we don't need to change the bone id of existing bones.

		\param	bone_retarget		the target position of each bone, set to 
									-1 it we don't want to change that bone
		\param	bone_start_pos		start point of each bone at rest pose
		\return						the remaining bone_weight number
	*/
	int RetargetBoneWeight(	const std::vector<int>&			bone_retarget,
							const std::vector<Vec3<T>>&		bone_start_pos ){
		if( bone_retarget.size() != bone_start_pos.size() ){
			#ifndef	BE_QUIET
				std::cout << "error: RiggingLBS::RetargetBoneWeight, bone_retarget size doesn't match bone_start_pos size" << std::endl;
			#endif
			return -1;
		}

		std::vector<BoneWeightLBS<T>>	tmp_bone_weight;
		std::vector<int>				tmp_start;

		int vertex_number = start.size() - 1;	
		tmp_start.push_back(0);

		for(int i=0; i<vertex_number; i++){
			//	first parse, retarget bones in the original array
			for(int j=start[i]; j<start[i+1]; j++){
				int bone_id = bone_weight[j].bone_id;
				if( bone_retarget[bone_id]!=-1 && bone_retarget[bone_id]!=bone_id ){	//	this bone need to be retargeted
					int target_bone_id = bone_retarget[bone_id];
					int k = start[i];
					for(; k<start[i+1]; k++){	//	search for target bone
						if( bone_weight[k].bone_id == target_bone_id ){	//	the vertex has been rigged to the target bone,
							bone_weight[k].weight += bone_weight[j].weight;	//	we simply add the weight to the target bone
							bone_weight[j].weight = 0;
							break;
						}
					}
					if( k == start[i+1] ){		//	vertex is not rigged to target bone, we change bone_id to the target
						bone_weight[j].bone_id = target_bone_id;
						bone_weight[j].offset += bone_start_pos[bone_id] - bone_start_pos[target_bone_id];
					}
				}
			}
			//	second parse, copy to tmp_bone_weight
			for(int j=start[i]; j<start[i+1]; j++){
				if( bone_weight[j].weight > 1e-6 )
					tmp_bone_weight.push_back(bone_weight[j]);
			}

			tmp_start.push_back(tmp_bone_weight.size());
		}

		bone_weight = tmp_bone_weight;
		start		= tmp_start;
		return bone_weight.size();
	}
};


/**
	A single vector for morphing, contain vertex index and a vector
*/
template<class T>
class MorphVector : public Vec3<T>{
public:
	int			vertex_index;	///<	index of the vertex to morph

public:
	MorphVector(int idx=-1, T x=0, T y=0, T z=0):vertex_index(idx), Vec3<T>(x, y, z){}
	MorphVector(int idx, Vec3<T> v):vertex_index(idx), Vec3<T>(v){}
};

/**
	morph vector of a mesh
*/
template<class T>
class MorphMesh{
public:
	std::string					name;			///<	name of this morph
	std::vector<MorphVector<T>>	morph_vector;	///<	morphing vector of each vertex

public:
	/**
		set morph vector from displacement vector

		For example, given two meshes with same mesh topology,
		the vectors between each corresponding vertices are displacement vector,
		this can be treated as morphing vector

		\param	name		name of the morph
		\param	vec			the displacement vector 
		\param	threshold	the threshold of displacement, length smaller than
							this threshold, we don't add to morph vector
		\return				1 on success
	*/
	int SetMorph(const char* name, std::vector<Vec3<T>>& vec, T threshold){
		this->name = name;
		morph_vector.clear();

		T squ_thre = threshold * threshold;
		for(int i=0; i<vec.size(); i++){
			if( vec[i].SquareLength() > squ_thre ){
				morph_vector.push_back(MorphVector<T>(i, vec[i]));
			}
		}
		return 1;
	}

	/**
		translate vertices by morph vectors

		\param	vertex		the vertex list
		\param	scale		scale of morph vector
	*/
	void TranslateVerticesByMorph(std::vector<Vec3<T>>& vertex, T scale){
		assert(vertex.size());

		for(int i=0; i<morph_vector.size(); i++){
			int vid = morph_vector[i].vertex_index;
			vertex[vid] += morph_vector[i] * scale;
		}
	}
};

/**
	morph of a whole skin

	morph of a skin is linear combination of several morph vectors
*/
template<class T>
class MorphSkin{
public:
	std::vector<T>				morph_coef;	///<	coefficient of each morph parameter
	std::vector<MorphMesh<T>>	morph_data;	///<	offset vector of each vertex at rest pose

public:
	/**
		set the number of morph parameters
	*/
	void SetMorphParameterNumber(int param_number){
		morph_coef.resize(param_number);
		morph_data.resize(param_number);
	}

	/**
		set a single morph vector by given morph vector index

		\param	id			index of the morph vector
		\param	name		name of the morph
		\param	vec			the displacement vector
		\param	threshold	threshold over which the displacement vector will take into account
		\return				1 on succeed
	*/
	int SetMorphData(int id, const char* name, std::vector<Vec3<T>>& vec, T threshold){
		assert(morph_coef.size() == morph_data.size());
		assert(morph_coef.size() > id);

		morph_data[id].SetMorph(name, vec, threshold);

		return 1;
	}

	/**
		Translate vertex list by the morph data

		\param	vertex		the vertices to be translated
	*/
	void TranslateVertices(std::vector<Vec3<T>>& vertex){
		assert(morph_coef.size() == morph_data.size());

		for(int i=0; i<morph_data.size(); i++){
			if( fabs(morph_coef[i]) > 1e-5 )
				morph_data[i].TranslateVerticesByMorph(vertex, morph_coef[i]);
		}
	}
	/**
		Calculate the overall morph vector

		Morph of a mesh is a linear combination of several morph vectors,
		with a given weight of each vector. This function is used to
		calculate the weighted sum of the morph vector

		\param	vec		return the overall morph vector
	*/
	void GetMorphVector(std::vector<Vec3<T>>& vec){
		assert(morph_coef.size() == morph_data.size());

		vec.clear();
		vec.resize( morph_data.size(), Vec3<T>(0, 0, 0) );
		for(int i=0; i<morph_data.size(); i++){
			if( fabs(morph_coef[i]) > 1e-5 )
				morph_data[i].TranslateVerticesByMorph(vec, morph_coef[i]);
		}
	}
	/**
		read morph data from file
	*/
	int ReadMorphFromFile(const char* file_name){
		std::ifstream file( file_name );	//	load obj failed, do not touch old data
		if(!file.is_open()){
			std::cout << "cannot open " << file_name << std::endl;
			return 0;
		}

		//	write number of morph data
		morph_coef.clear();
		morph_data.clear();
		int num;
		file >> num;
		morph_coef.resize(num);
		morph_data.resize(num);

		//	write each data
		for(int i=0; i<num; i++){
			int data_size;
			file >> morph_data[i].name >> data_size >> morph_coef[i];
			morph_data[i].morph_vector.resize(data_size);

			for(int j=0; j<data_size; j++){
				file >> morph_data[i].morph_vector[j].vertex_index
					>> morph_data[i].morph_vector[j].x
					>> morph_data[i].morph_vector[j].y
					>> morph_data[i].morph_vector[j].z;
			}
		}

		file.close();
		return 1;
	}

	/**
		write morph data to file
	*/
	int WriteMorphToFile(const char* file_name){
		std::ofstream file( file_name );
		if(!file.is_open()){
			#ifndef BE_QUIET
				std::cout << "error: MorphSkin::WriteMorphToFile, cannot open " << file_name << std::endl;
			#endif
			return 0;
		}

		//	write number of morph data
		int num = morph_coef.size();
		assert(num == morph_data.size());
		file << num << std::endl;

		//	write each data
		for(int i=0; i<num; i++){
			int data_size = morph_data[i].morph_vector.size();
			file << morph_data[i].name << '\t' << data_size << '\t' << morph_coef[i] << std::endl;

			for(int j=0; j<data_size; j++){
				file << morph_data[i].morph_vector[j].vertex_index << '\t' 
					<< morph_data[i].morph_vector[j].x << '\t' 
					<< morph_data[i].morph_vector[j].y << '\t' 
					<< morph_data[i].morph_vector[j].z << std::endl;
			}
		}

		file.close();

		return 1;
	}

};



}}	//	namespace yz::animation

#endif	//	__YZ_SKIN_H__