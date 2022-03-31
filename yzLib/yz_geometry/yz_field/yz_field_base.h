/***********************************************************/
/**	\file
	\brief		Field Base
	\author		Yizhong Zhang
	\date		12/13/2017
*/
/***********************************************************/
#ifndef __YZ_FIELD_BASE_H__
#define __YZ_FIELD_BASE_H__

#include <iostream>
#include <fstream>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_nn.h"

namespace yz{	namespace geometry{		namespace field{

/**
	Base class of field
*/
template<class T>
class FieldBase{
public:
	FieldBase(){
		dim.x = dim.y = dim.z = 0;
	}

	/**
		Get the voxel index in the field array
	*/
	inline int GetVoxelID(int x, int y, int z){
		return (z * dim.y + y) * dim.x + x;
	}

	/**
		Get the sdf of a voxel
	*/
	inline T GetData(int x, int y, int z){
		return data[(z * dim.y + y) * dim.x + x];
	}

public:
	Vec3ui			dim;	///<	dimension of the field
	std::vector<T>	data;	///<	data of the field
};


}}}	//	namespace yz::geometry::field

#endif	//	__YZ_FIELD_BASE_H__