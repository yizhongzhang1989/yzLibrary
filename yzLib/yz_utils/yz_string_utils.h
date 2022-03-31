/***********************************************************/
/**	\file
	\brief		Utilities to Process String
	\author		Yizhong Zhang
	\date		6/8/2012
*/
/***********************************************************/
#ifndef __YZ_STRING_UTILS_H__
#define __YZ_STRING_UTILS_H__

#include <vector>
#include <list>
#include <string>
#include <string.h>

namespace yz{ namespace utils{

/**
	A poll to hold a lot of string

	We can add string to the pool. If the number of string
	has exceeded limit, old string will be deleted

	printStringList() can be used to display StringPool
	on OpenGL window
*/
class StringPool{
public:
	int pool_depth;					///< how many string can the pool hold
	std::list<std::string> pool;	///< the pool to hold string

public:
	/**
		Default pool depth: 20
	*/
	StringPool() : pool_depth(20){};
	/**
		Set the depth of the string pool (string capacity)
	*/
	inline void SetPoolDepth(int depth){
		if( depth <= 1 )	depth = 1;
		pool_depth = depth;
		RemoveOverFlow();
	}

	/**
		Add string to the pool.
	*/
	inline void AddString(std::string str){
		pool.push_back(str);
		RemoveOverFlow();
	}

protected:
	/**
		remove string that more than the pool can sustain
	*/
	inline void RemoveOverFlow(){
		while( pool.size() > pool_depth ){
			pool.pop_front();
		}
	}
};

//	========================================
///@{
/**	@name Bland Character in string
*/
//	========================================

/**
	check whether a character is a blank character

	blank character include: space(' '), tab('\\t'), enter('\\n')
*/
inline int isBlankCharacter(char c){
	if( c==' ' || c=='\t' || c=='\n' )
		return 1;
	return 0;
}

/**
	remove prefix blank characters from a string
*/
inline void removePrefixBlankCharacters(char* str){
	char* str2 = str;
	while( isBlankCharacter(*str2) && (*str2) != '\0' )	//	find the first non-blank character
		str2 ++;
	if( str == str2 )	return;
	while(*str2 != '\0'){
		*str = *str2;
		str ++;
		str2 ++;
	}
	*str = *str2;
}

/**
	remove prefix blank characters from a string
*/
inline void removePrefixBlankCharacters(std::string& str){
	std::string::iterator iter;
	for( iter = str.begin(); iter != str.end(); iter++ )
		if( !isBlankCharacter(*iter) )	break;
	if( iter == str.begin() )	return;
	str.erase(str.begin(), iter);
}

/**
	remove suffix blank characters from a string
*/
inline void removeSuffixBlankCharacters(char* str){
	char* str2 = str;
	while( (*str2) != '\0' )	//	find the end of the string
		str2 ++;
	str2 --;
	while( str2 != str-1 ){
		if( !isBlankCharacter(*str2) )
			break;
		*str2 = '\0';
		str2 --;
	}
}

/**
	remove suffix blank characters from a string
*/
inline void removeSuffixBlankCharacters(std::string& str){
	std::string::reverse_iterator iter;
	int count = 0;
	for( iter = str.rbegin(); iter != str.rend(); iter++ ){
		if( !isBlankCharacter(*iter) )	break;
		count ++;
	}
	if( !count )	return;
	str.erase(str.end()-count, str.end());
}

/**
	remove prefix and suffix blank characters from a string
*/
inline void removePrefixSuffixBlankCharacters(char* str){
	removePrefixBlankCharacters(str);
	removeSuffixBlankCharacters(str);
}

/**
	remove prefix and suffix blank characters from a string
*/
inline void removePrefixSuffixBlankCharacters(std::string& str){
	removePrefixBlankCharacters(str);
	removeSuffixBlankCharacters(str);
}

///@}

//	========================================
///@{
/**	@name File Name Operations
*/
//	========================================

/**
	Combine directory and file name

	\param	directory	directory of the file, without sufix '/'
	\param	file_name	file name, without prefix '/'
	\return				file name string
*/
inline std::string getFileNameCombineDirentory(const char* directory, const char* file_name){
	std::string real_directory = directory;
	for(int i=0; i<real_directory.size(); i++){	//	check whether directory contian only blank characters
		if( !isBlankCharacter( real_directory[i] ) )
			break;
		if( i == real_directory.size() )	//	all the characters are blank, we just return file name
			return file_name;
	}

	//	the directory is not blank
	if( real_directory.size() ){
		return real_directory + '/' + file_name;
	}
	else
		return file_name;
}

/**
	Get file name without directory

	if the file name contain directory in the string, then this function
	removes the directory. 

	\param	file_name	file name
	\return				file name string, 
						if the string don't have directory, we just return the original string
*/
inline std::string getFileNameWithoutDirectory(const char* file_name){
	if( !file_name || *file_name == '\0' )
		return "";

	const char* p = file_name + strlen(file_name) - 1;
	while( p!=file_name && *p!='\\' && *p!='/' )
		p --;
	if(p == file_name)	//	the whole string don't contain directory
		return file_name;

	std::string dir(p+1, file_name + strlen(file_name));
	return dir;
}

/**
	Get directory from a string

	if the file name contain directory in the string, then this function
	return the directory. The endding \ or / is not included.

	\param	file_name	file name
	\return				file directory string, 
						if the string don't have directory, we just return empty string
*/
inline std::string getDirectoryFromString(const char* file_name){
	if( !file_name || *file_name == '\0' )
		return "";

	const char* p = file_name + strlen(file_name) - 1;
	while( p!=file_name && *p!='\\' && *p!='/' )
		p --;
	if(p == file_name)	//	the whole string don't contain directory
		return "";

	std::string dir(file_name, p);
	return dir;
}

/**
	Get file extension from a string.

	file extension is the last non-empty sub string that don't 
	include '.' but has a '.' in front of it in the string

	If the last non-blank charactor of the string is '.', then
	the file extension is empty. 

	The return string don't include '.'

	\param	file_name	file name
	\return				file extension string, 
						if file_name don't have extension, return empty string
*/
inline std::string getFileExtensionFromString(const char* file_name){
	if( !file_name || *file_name == '\0' )
		return "";

	//	parse the string from end to front
	const char* p = file_name + strlen(file_name) - 1;
	while( p!=file_name && (*p==' ' || *p=='\t' || *p=='\n') )	//	skip blank charactors
		p --;
	if(p == file_name)	//	the whole string is made up of at most one non-blank charactor
		return "";
	else if(*p == '.')	//	the last non-blank charactor is '.'
		return "";

	const char* ext_end = p + 1;
	while(p!=file_name && *p!='.')
		p --;
	if(*p != '.')
		return "";

	std::string extension(p+1, ext_end);
	return extension;
}

/**
	append a string to the file name before file extension

	for example, getAppendedFileName("file.txt", "_01") returns "file_01.txt"

	\param	file_name		name of the file
	\param	append_str		string that append to the end of the file
	\return					the resulting file name
*/
inline std::string getFileNameAppendString(const char* file_name, const char* append_str){
	std::string extension = getFileExtensionFromString(file_name);
	std::string new_file_name(file_name);
	int insert_pos = new_file_name.size() - extension.size() - 1;
	if (insert_pos < 0)
		insert_pos = 0;
	new_file_name.insert(insert_pos, append_str);
	return new_file_name;
}

/**
	change file extension

	\param	file_name		original name of the file
	\param	new_extension	new file extension
	\return					return the changed extension file name
*/
inline std::string getFileNameChangeExtension(const char* file_name, const char* new_extension){
	std::string extension = getFileExtensionFromString(file_name);
	std::string new_file_name(file_name);
	new_file_name.erase(new_file_name.size() - extension.size(), extension.size());
	new_file_name.append(new_extension);
	return new_file_name;
}

///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_STRING_UTILS_H__