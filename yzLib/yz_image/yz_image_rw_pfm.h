/***********************************************************/
/**	\file
	\brief		Read & Write PFM Image
	\details	Read .pfm files
	\author		Yizhong Zhang
	\date		10/22/2019
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RW_PFM_H__
#define __YZ_IMAGE_RW_PFM_H__

#ifdef _MSC_VER
#	include <io.h>
#else
#	include <sys/io.h>
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_utils/yz_string_utils.h"

namespace yz{ namespace image{

/**
read a pfm image from file

\param	pfm_file_name	image file name
\param	image_ptr		pointer to the image. memory will be allocated inside this function
\param	width			the width of the image
\param	height			the height of the image
\param	channels		number of channels of the image, must be 1 or 3
\param	scale			scale value (positive)
\return					1 succeed, 0 read failed
*/
inline int readPfmFromFile(
	const char*		pfm_file_name,
	float* &		image_ptr,
	int&			width,
	int&			height,
	int&			channels,
	float&			scale
) {
	std::ifstream pfm(pfm_file_name, std::ios::binary);
	if (!pfm.is_open()) {
		std::cout << "cannot open " << pfm_file_name << ", read pfm failed" << std::endl;
		return 0;
	}

	// init variables 
	std::string identifier;	// "Pf" = grayscale (1 band), "PF" = color (3 band)
	int w, h;				// width and height of the image
	float scalef;			// scale factor

	if (!(pfm >> identifier >> w >> h >> scalef)) {
		std::cout << "read " << pfm_file_name << " failed" << std::endl;
		return 0;
	}

	int bands = 0;
	if (identifier == "Pf")
		bands = 1;
	else if (identifier == "PF")
		bands = 3;
	else {
		std::cout << "illegal identifier: " << identifier << ", read " << pfm_file_name << " failed" << std::endl;
		return 0;
	}

	channels = bands;
	scale = fabs(scalef);

	// determine endianness 
	int intval = 1;
	char *uval = (char *)&intval;
	int littleEndianMachine = (uval[0] == 1);

	int littleEndianFile = (scalef < 0);
	int needSwap = (littleEndianFile != littleEndianMachine);

	// skip SINGLE newline character after reading third arg
	char c = pfm.get();
	if (c == '\r')      // <cr> in some files before newline
		c = pfm.get();
	if (c != '\n') {
		if (c == ' ' || c == '\t' || c == '\r') {
			std::cout << "error: readPfmFromFile, newline expected";
			return 0;
		}
		else {
			std::cout << "error: readPfmFromFile, whitespace expected";
			return 0;
		}
	}

	//	create image space
	width = w;
	height = h;
	int image_size = w * h * bands;
	image_ptr = new float[image_size];

	//	read image
	pfm.read((char*)image_ptr, image_size * 4);
	if (needSwap) {
		for (int i = 0; i < image_size; i++) {
			char* ptr = (char*)&image_ptr[i];
			char tmp;
			tmp = ptr[0]; ptr[0] = ptr[3]; ptr[3] = tmp;
			tmp = ptr[1]; ptr[1] = ptr[2]; ptr[2] = tmp;
		}
	}

	return 1;
}

}}	//	end namespace yz::image

#endif	//	__YZ_IMAGE_RW_PFM_H__