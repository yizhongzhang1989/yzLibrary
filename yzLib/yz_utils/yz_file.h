/***********************************************************/
/**	\file
	\brief		File Operation
	\author		Yizhong Zhang
	\date		9/18/2012
*/
/***********************************************************/
#ifndef __YZ_FILE_H__
#define __YZ_FILE_H__

#ifdef _MSC_VER
#	include <io.h>
#else
#	include <sys/io.h>
#endif

#include <fstream>

#ifndef _MSC_VER 
#	include <unistd.h>
#else
#	define access _access
#endif


namespace yz{ namespace utils{

//	========================================
///@{
/**	@name File Operations
*/
//	========================================

/**
	check whether a file exist
*/
inline int isFileExist(const char* file_name){
	return access(file_name, 0) == 0;
}

/**
	get size of file in bytes

	if file doesn't exist, return -1
*/
inline long getFileSizeInBytes(std::ifstream& file_stream){
	file_stream.seekg(0, std::ifstream::end);
	long size = file_stream.tellg(); 
	file_stream.seekg(0, std::ifstream::beg);
	return size;
}

/**
	get size of file in bytes

	if file doesn't exist, return -1
*/
inline long getFileSizeInBytes(const char* file_name){
	if( !isFileExist(file_name) )
		return -1;

	std::ifstream in(file_name, std::ifstream::in | std::ifstream::binary);
	long size = getFileSizeInBytes(in);
	in.close();
	return size;
}

///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_FILE_H__