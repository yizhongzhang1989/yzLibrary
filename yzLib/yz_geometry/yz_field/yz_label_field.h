/***********************************************************/
/**	\file
	\brief		Label Field
	\author		Yizhong Zhang
	\date		12/18/2017
*/
/***********************************************************/
#ifndef __YZ_LABEL_FIELD_H__
#define __YZ_LABEL_FIELD_H__

#include <iostream>
#include <fstream>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_field/yz_field_base.h"

namespace yz{	namespace geometry{		namespace field{

/**
	A field of labels, represented using integer
*/
class LabelField : public FieldBase<int> {
public:
	/**
	setup data of the volume

	\param	dim_x		dimension in x direction
	\param	dim_y		dimension in y direction
	\param	dim_z		dimension in z direction
	*/
	void SetupVolume(unsigned int dim_x, unsigned int dim_y, unsigned int dim_z) {
		this->dim = uint3(dim_x, dim_y, dim_z);
		this->data.resize(this->dim.x * this->dim.y * this->dim.z);
	}

	/**
	read the field from file in binary format

	the file is arranged in the following format:	\n
	dim				3 int							\n
	data			dim.x*dim.y*dim.z double		\n
	*/
	int ReadFromBinaryFile(const char* file_name) {
		std::ifstream file(file_name, std::ifstream::binary);
		if (!file.is_open()) {
			std::cout << "error: LabelField::ReadFromBinaryFile, open file failed" << std::endl;
			return 0;
		}

		//	read parameters
		file.read((char*)&this->dim[0], sizeof(unsigned int) * 3);

		//	read data
		this->data.resize(this->dim.x* this->dim.y* this->dim.z);
		if (!this->data.empty()) {
			file.read((char*)&this->data[0], sizeof(int)* this->data.size());
		}

		file.close();

		return 1;
	}

	/**
	write the field to file in binary format

	file format is the same as ReadFromBinaryFile()
	*/
	int WriteToBinaryFile(const char* file_name) {
		std::ofstream file(file_name, std::ofstream::binary);
		if (!file.is_open()) {
			std::cout << "error: LabelField::WriteToBinaryFile, open file failed" << std::endl;
			return 0;
		}

		//	write parameters
		file.write((char*)&dim[0], sizeof(unsigned int) * 3);

		//	write data
		if (!data.empty()) {
			file.write((char*)&data[0], sizeof(int)*data.size());
		}

		return 1;
	}

};


}}}	//	namespace yz::geometry::field

#endif	//	__YZ_LABEL_FIELD_H__