/***********************************************************/
/**	\file
	\brief		Read & Write Image
	\details	All image rw related files has been moved to yz::image,
				we just leave interface here for compatibility
	\author		Yizhong Zhang
	\date		6/22/2016
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RW_UTILS_H__
#define __YZ_IMAGE_RW_UTILS_H__

#include "yzLib/yz_image/yz_image_rw.h"

namespace yz{ namespace utils{

//	========================================
///@{
/**	@name Read / Write Image
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
#ifndef YZ_FreeImage_h
	return yz::image::readBmpFromFile(image_file_name, image_ptr, width, height);
#else					//	read using API provided by FreeImage
	return yz::image::readImageFromFile(image_file_name, image_ptr, width, height, target_bpp);
#endif
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
		delete[] tmp_image_ptr;
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
#ifndef YZ_FreeImage_h	//	without FreeImage, we can only read .bmp file
	return yz::image::writeBmpToFile(image_file_name, image_ptr, width, height, bpp);
#else
	return yz::image::writeImageToFile(image_file_name, image_ptr, width, height, bpp);
#endif
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
	delete[] tmp_image_ptr;
	return success;
}
///@}

//	========================================
///@{
/**	@name Read / Write .BMP file
*/
//	========================================

/**
	Read .bmp file from hard disk

	Currently, only LUMINANCE, RGB, RGBA can be read

	If read image failed, don't touch old data

	\param	bmp_file_name	name of file, file extension must be .bmp
	\param	image_ptr		pointer of the image. new space will be allocated. 
							Pixels are stored in RGB sequence.
							Storage from top-left to bottom-right without padding
	\param	width			width of image in pixels
	\param	height			height of image in pixels
	\return					bits per pixel, 0: read image failed
*/
inline int readBmpFromFile(const char* bmp_file_name, unsigned char* &image_ptr, int& width, int& height){
	return image::readBmpFromFile(bmp_file_name, image_ptr, width, height);
}

/**
	Write .bmp file to hard disk

	Only LUMINANCE, RGB, RGBA file can be written using this function

	\param	bmp_file_name	name of file, file extension must be .bmp
	\param	image_ptr		pointer of the image. new space will be allocated. 
							Pixels are stored in RGBA sequence.
							Storage from top-left to bottom-right without padding
	\param	width			width of image in pixels
	\param	height			height of image in pixels
	\param	bits_per_pixel	bits per pixel, typically 8/24/32
	\return					write success flag, 1: succeed,  0: failed
*/
inline int writeBmpToFile(const char* bmp_file_name, const unsigned char* image_ptr, int width, int height, int bits_per_pixel){
	return image::writeBmpToFile(bmp_file_name, image_ptr, width, height, bits_per_pixel);
}
///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_IMAGE_RW_UTILS_H__