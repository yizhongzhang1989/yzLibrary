/***********************************************************/
/**	\file
	\brief		Read & Write BMP Image
	\details	Read .BMP files
	\author		Yizhong Zhang
	\date		6/24/2012
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RW_BMP_H__
#define __YZ_IMAGE_RW_BMP_H__

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

//	========================================
///@{
/**	@name Read / Write .BMP file
*/
//	========================================
#if (defined(_WIN32) || defined(__WIN32__))
#	pragma pack(push, 2)
#else
#	pragma pack(2)
#endif
/**
	BMP file header
*/
class BitmapFileHeader{
public:
	short	signature;	///<	inverse "BM", because the machine is little endian
	int		file_size;	///<	size of the whole file
	short	reserved1;	///<	not used
	short	reserved2;	///<	not used
	int		offset;		///<	starting position (from file start) of pixel data

	BitmapFileHeader() : signature(0x4D42) {}
};

/**
	BMP file info header
*/
class BitmapInfoHeader{
public:
	int		header_size;		///<	size of this header
	int		width;				///<	image width in pixels
	int		height;				///<	image height in pixels
	short	planes;				///<	must be 1
	short	bits_per_pixel;		///<	typical 1/4/8/16/24/32
	int		compression;		///<	compression type, 0: none
	int		image_size;			///<	size of image in bytes, with padding
	int		pixels_per_meter_x;
	int		pixels_per_meter_y;
	int		colors_used;
	int		colors_important;

	BitmapInfoHeader(){
		header_size		= sizeof(BitmapInfoHeader);
		width			= 0;
		height			= 0;
		planes			= 1;
		bits_per_pixel	= 24;
		compression		= 0;
		image_size		= 0;
		pixels_per_meter_x	= 0;
		pixels_per_meter_y	= 0;
		colors_used			= 0;
		colors_important	= 0;
	}
};
#if (defined(_WIN32) || defined(__WIN32__))
#	pragma pack(pop)
#else
#	pragma pack()
#endif

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
	std::ifstream bmp( bmp_file_name, std::ios::binary );
	if(!bmp.is_open()){
		#ifndef BE_QUIET
			std::cout << "cannot open " << bmp_file_name << ", read bmp failed" << std::endl;
		#endif
		return 0;
	}

	BitmapFileHeader file_header;
	BitmapInfoHeader info_header;
	bmp.read((char*)&file_header, sizeof(file_header));
	bmp.read((char*)&info_header, sizeof(info_header));

	//	file check
	if( file_header.signature != 0x4D42 ){
		#ifndef BE_QUIET
			std::cout << bmp_file_name << " signature error" << std:: endl;
		#endif
		return 0;
	}
	else if( info_header.width == 0 || info_header.height == 0 ){
		#ifndef BE_QUIET
			std::cout << bmp_file_name << " width or height equal zero" << std:: endl;
		#endif
		return 0;
	}
	else if( info_header.compression ){
		#ifndef BE_QUIET
			std::cout << bmp_file_name << " compressed" << std:: endl;
		#endif
		return 0;
	}
	else if( info_header.bits_per_pixel != 8 && //LUM
			info_header.bits_per_pixel != 24 && //RGB
			info_header.bits_per_pixel != 32 ){	//RGBA
		#ifndef BE_QUIET
			std::cout << bmp_file_name << " only support LUM, RGB, RGBA format" << std:: endl;
		#endif
		return 0;
	}

	width	= info_header.width;
	height	= info_header.height;
	int bytes_per_pixel = info_header.bits_per_pixel / 8;

	//	read to tmp array
	int row_size = (bytes_per_pixel * width + 3)/4 * 4;
	int image_size = row_size * height;
	unsigned char* tmp_image_ptr = new unsigned char[image_size];
	bmp.seekg(file_header.offset, std::ios::beg);
	bmp.read((char*)tmp_image_ptr, image_size);

	//	copy and flip from tmp to image
	image_ptr = new unsigned char[width * height * bytes_per_pixel];	//	image_ptr is new created
	for( int j=0; j<height; j++ ){
		memcpy(	image_ptr + j * width * bytes_per_pixel, 
				tmp_image_ptr + (height-1-j) * row_size,
				sizeof(unsigned char) * width * bytes_per_pixel );
	}

	if( bytes_per_pixel == 3 || bytes_per_pixel == 4 ){
		for(int j=0; j<height; j++){
			unsigned char* row = image_ptr + width * bytes_per_pixel * j;
			for(int i=0; i<width; i++ ){
				unsigned char tmp = row[i*bytes_per_pixel];
				row[i*bytes_per_pixel] = row[i*bytes_per_pixel+2];
				row[i*bytes_per_pixel+2] = tmp;
			}
		}
	}

	delete[] tmp_image_ptr;
	bmp.close();
	return info_header.bits_per_pixel;
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
	//	image check
	if( width<=0 || height <=0 ){
		#ifndef BE_QUIET
			std::cout << "bmp size error, fail to write bmp to disk" << std::endl;
		#endif
		return 0;
	}
	if( bits_per_pixel!=8 && bits_per_pixel!=24 && bits_per_pixel!=32 ){
		#ifndef BE_QUIET
			std::cout << "unsupported bits per pixel (" << bits_per_pixel << "), fail to write bmp to disk" << std::endl;
		#endif
		return 0;
	}

	std::ofstream bmp( bmp_file_name, std::ios::binary );
	if(!bmp.is_open()){
		#ifndef BE_QUIET
			std::cout << "cannot open " << bmp_file_name << ", write bmp failed" << std::endl;
		#endif
		return 0;
	}

	//	setup file header and info header
	int bytes_per_pixel = bits_per_pixel / 8;
	int row_size = (bytes_per_pixel * width + 3)/4 * 4;

	BitmapFileHeader file_header;
	BitmapInfoHeader info_header;
	info_header.width			= width;
	info_header.height			= height;
	info_header.bits_per_pixel	= bits_per_pixel;
	info_header.image_size		= row_size * height;
	file_header.offset		= sizeof(file_header) + sizeof(info_header);
	file_header.file_size	= file_header.offset + info_header.image_size;

	//	create tmp array
	unsigned char* tmp_image_ptr = new unsigned char[info_header.image_size];

	for( int j=0; j<height; j++ ){
		memcpy(	tmp_image_ptr + (height-1-j) * row_size,
				image_ptr + j * width * bytes_per_pixel, 
				sizeof(unsigned char) * width * bytes_per_pixel );
	}

	if( bytes_per_pixel == 3 || bytes_per_pixel == 4 ){
		for(int j=0; j<height; j++){
			unsigned char* row = tmp_image_ptr + row_size * j;
			for(int i=0; i<width; i++ ){
				unsigned char tmp = row[i*bytes_per_pixel];
				row[i*bytes_per_pixel] = row[i*bytes_per_pixel+2];
				row[i*bytes_per_pixel+2] = tmp;
			}
		}
	}

	//	write to file
	bmp.write((char*)&file_header, sizeof(file_header));
	bmp.write((char*)&info_header, sizeof(info_header));
	bmp.write((char*)tmp_image_ptr, info_header.image_size);

	delete[] tmp_image_ptr;
	bmp.close();
	return 1;
}

///@}
}}	//	end namespace yz::image

#endif	//	__YZ_IMAGE_RW_BMP_H__