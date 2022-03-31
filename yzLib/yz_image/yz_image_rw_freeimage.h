/***********************************************************/
/**	\file
	\brief		Read & Write Image via FreeImage library
	\details	FreeImage can help read a lot kinds of files
				user must make sure FreeImage.h is included before
				including this file, or report error \n
	\author		Yizhong Zhang
	\date		6/20/2016
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RW_FREEIMAGE_H__
#define __YZ_IMAGE_RW_FREEIMAGE_H__

#ifndef YZ_FreeImage_h
#	error yz_image_rw_freeimage.h must be included after FreeImage.h
#endif

#include <io.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_utils/yz_string_utils.h"

namespace yz{ namespace image{

//	========================================
///@{
/**	@name Read / Write Image by FreeImage
*/
//	========================================

/**
	Read image from hard disk

	If FreeImage is included, then we can read all image format
	that FreeImage support. If not, only bmp image can be read.

	If read image failed, don't touch old data.

	The minimal bits per pixel is 8, if bpp of image data is less
	than 8, convert the image to grayscale. (FreeImage required)

	\param	image_file_name	file name of the image
	\param	image_ptr		pointer of the image. new space will be allocated. \n
							Pixels are stored in RGB sequence. \n
							Storage from top-left to bottom-right without padding
	\param	width			width of image in pixels
	\param	height			height of image in pixels
	\param	target_bpp		target_depth of read image raw data. valid value: 0/8/24/32	\n
							0:	default value, keep original depth of the file \n
							8:	set the return raw data to gray scale \n
							24:	set the return raw data to RGB \n
							32:	set the return raw data to RGBA
	\return					bits per pixel, 0: read image failed \n
							It is possible thet the return value doesn't match target_bpp, 
							the return value is assumed to be the correct bits per pixel value
*/
inline int readImageFromFile(const char* image_file_name, unsigned char* &image_ptr, int& width, int& height, int target_bpp=0){
	//	read using API provided by FreeImage
	//	If user included FreeImage.h not modified by Yizhong, then _WINDOWS_ macro
	//	is introduced without including <windows.h> which will cause mistake.
	//	This problem is not solved, just print error information
#	ifdef YZ_windows_h
#		ifndef _INC_WINDOWS
#			error	_WINDOWS_ introduced by FreeImage.h without include <windows.h>, \
					fix this by include <windows.h> before FreeImage.h. \
					If you are sure this doesn't introduce any problem, \
					remove this code segment from yzLib
#		endif
#	endif

	//	check file existance
	if( _access( image_file_name, 0 ) ){
		std::cout << image_file_name << " doesn't exist, read file failed" << std::endl;
		return 0;
	}

	FREE_IMAGE_FORMAT fif = FreeImage_GetFileType(image_file_name, 0);
	if(fif == FIF_UNKNOWN)
		fif = FreeImage_GetFIFFromFilename(image_file_name);
	FIBITMAP *img;
	if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)){
		img = FreeImage_Load(fif, image_file_name);
		if( !img ){
			#ifndef BE_QUIET
				std::cout << "read " << image_file_name << " failed" << std::endl;
			#endif
			return 0;
		}
	}
	else{
		#ifndef BE_QUIET
			std::cout << "yz::readImageFromFile() doesn't support reading ." 
				<< utils::getFileExtensionFromString(image_file_name) << " image, read file failed" << std::endl;
		#endif
		return 0;
	}

	//	if target bits per pixel is set, then set the image to match target
	if( target_bpp != 0 && target_bpp != FreeImage_GetBPP(img) ){
		FIBITMAP* tmp_img = NULL;
		if( target_bpp == 8 )
			tmp_img = FreeImage_ConvertToGreyscale(img);
		else if( target_bpp == 24 )
			tmp_img = FreeImage_ConvertTo24Bits(img);
		else if( target_bpp == 32 )
			tmp_img = FreeImage_ConvertTo32Bits(img);

		if( tmp_img ){
			FreeImage_Unload(img);
			img = tmp_img;
		}	//	if this failed, convertion is just failed, we keep original format
	}

	//	if bits per pixel is smaller than 8, promote to 8
	if( FreeImage_GetBPP(img) < 8 ){
		FIBITMAP* tmp_img = FreeImage_ConvertToGreyscale(img);
		if( tmp_img ){
			FreeImage_Unload(img);
			img = tmp_img;
		}
		else{
			#ifndef BE_QUIET
			std::cout << image_file_name << " is " << FreeImage_GetBPP(img) << " bits per pixel, "
				<< "fail to promote to 8 bits per pixel, read file failed" << std::endl;
			#endif
			FreeImage_Unload(img);
			return 0;
		}
	}

	width	= FreeImage_GetWidth(img);
	height	= FreeImage_GetHeight(img);
	int bpp		= FreeImage_GetBPP(img);
	int pitch	= FreeImage_GetPitch(img);
	int bytes_pp = bpp / 8;

	image_ptr = new unsigned char[width*height*bytes_pp];
	FreeImage_ConvertToRawBits(image_ptr, img, width*bytes_pp, bpp, 0, 0, 0, TRUE);

	//	swap R and B
	if( bytes_pp == 3 || bytes_pp == 4 ){
		for( int j=0; j<height; j++ ){
			unsigned char* ptr = image_ptr + j * width * bytes_pp;
			for( int i=0; i<width; i++ ){
				unsigned char tmp = ptr[i*bytes_pp  ];
				ptr[i*bytes_pp  ] = ptr[i*bytes_pp+2];
				ptr[i*bytes_pp+2] = tmp;
			}
		}
	}


	FreeImage_Unload(img);
	return bpp;
}

/**
	read the image to a normalized float array

	\param	image_file_name	file name of the image
	\param	image_ptr		pointer of the image. new space will be allocated. \n
							Pixels are stored in RGB sequence. \n
							Storage from top-left to bottom-right without padding
	\param	width			width of image in pixels
	\param	height			height of image in pixels
	\param	target_bpp		target_depth of read image raw data. valid value: 0/8/24/32	\n
							0:	default value, keep original depth of the file \n
							8:	set the return raw data to gray scale \n
							24:	set the return raw data to RGB \n
							32:	set the return raw data to RGBA
	\return					bits per pixel, 0: read image failed \n
							It is possible thet the return value doesn't match target_bpp, 
							the return value is assumed to be the correct bits per pixel value

*/
inline int readImageFromFile(const char* image_file_name, float* &image_ptr, int& width, int& height, int target_bpp=0){
	unsigned char* tmp_image_ptr;
	int bpp = readImageFromFile(image_file_name, tmp_image_ptr, width, height, target_bpp);
	if( bpp ){
		int size = width * height * bpp/8;
		image_ptr = new float[size];
		for( int i=0; i<size; i++ )
			image_ptr[i] = tmp_image_ptr[i] / 255.0f;
		delete tmp_image_ptr;
	}
	return bpp;
}

/**
	Write Image to file on hard disk

	This function write image raw data to file. The raw data
	could be LUMINANCE, RGB, RGBA of unsigned char type.

	\param	image_file_name	file name of the image
	\param	image_ptr		pointer to the raw data image storage
	\param	width			width of the raw data
	\param	height			height of the raw data
	\param	bpp				bits per pixel of the raw data, 8/24/32
	\return					whether write succeed
*/
inline int writeImageToFile(const char* image_file_name, unsigned char* image_ptr, int width, int height, int bpp){
	if( bpp!=8 && bpp!=24 && bpp!=32 ){
		std::cout << "writeImageToFile only accept LUMINANCE, RGB, RGBA format with unsigned char type, " 
			<< "write file failed" << std::endl;
		return 0;
	}

	//	check file extension
	FREE_IMAGE_FORMAT fif = FreeImage_GetFileType(image_file_name, 0);
	if(fif == FIF_UNKNOWN)
		fif = FreeImage_GetFIFFromFilename(image_file_name);
	if(fif == FIF_UNKNOWN){
		#ifndef BE_QUIET
			std::cout << "FreeImage doesn't support export ." << utils::getFileExtensionFromString(image_file_name)
				<< " image, write file failed" << std::endl;
		#endif
		return 0;
	}

	//	create FIBITMAP
	FIBITMAP* img = FreeImage_ConvertFromRawBits(image_ptr, width, height, width*bpp/8, bpp, 
		FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

	//	swap R & B
	int bytes_pp = bpp / 8;
	if( bytes_pp == 3 || bytes_pp == 4 )
		for(int j=0; j<height; j++ ){
			unsigned char* ptr = FreeImage_GetScanLine(img, j);
			for( int i=0; i<width; i++ ){
				unsigned char tmp = ptr[i*bytes_pp  ];
				ptr[i*bytes_pp  ] = ptr[i*bytes_pp+2];
				ptr[i*bytes_pp+2] = tmp;
			}
		}

	if( !FreeImage_Save(fif, img, image_file_name) ){
		std::cout << "write " << image_file_name << " failed" << std::endl;
		return 0;
	}

	FreeImage_Unload(img);
	return 1;
}

/**
	Write Image to file on hard disk

	This function write image raw data to file. The raw data
	could be LUMINANCE, RGB, RGBA of normalzied float type.

	\param	image_file_name	file name of the image
	\param	image_ptr		pointer to the raw data image storage
	\param	width			width of the raw data
	\param	height			height of the raw data
	\param	bpp				bits per pixel of the raw data, 32/96/128
	\return					whether write succeed
*/
inline int writeImageToFile(const char* image_file_name, const float* image_ptr, int width, int height, int bpp){
	if( bpp!=32 && bpp!=96 && bpp!=128 )	//	not LUM, RGB, RGBA type
		return 0;

	int size = width * height * bpp/32;
	unsigned char* tmp_image_ptr = new unsigned char[size];
	for( int i=0; i<size; i++ )
		tmp_image_ptr[i] = roundToClosestInteger(image_ptr[i] * 255);
	int success = writeImageToFile(image_file_name, tmp_image_ptr, width, height, bpp/4);
	delete tmp_image_ptr;
	return success;
}
///@}

}}	//	end namespace yz::image

#endif	//	__YZ_IMAGE_RW_FREEIMAGE_H__