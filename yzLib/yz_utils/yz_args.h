/***********************************************************/
/**	\file
	\brief		Args
	\author		Yizhong Zhang
	\date		6/20/2012
*/
/***********************************************************/
#ifndef __YZ_ARGS_H__
#define __YZ_ARGS_H__

#include <string.h>

namespace yz{ namespace utils{

/**
	This class provide functions to process args in main(int argc, char* argv[]){}
*/
class ArgUtils{
public:
	/**
		This function creates arg according to string

		new space will be allocated in this function, they are:
		1, argv; 2, argv[0] (if argc != 0)
		
		\param	argc	main(int argc, char* argv[]), call Str2Arg(&argc, ...)
		\param	argv	main(int argc, char* argv[]), call Str2Arg(..., &argv, ...)
		\param	str		the input string
		\return			argc
	*/
	static int Str2Arg(int* argc, char*** argv, const char* str){
		int str_len = strlen(str);
		char* arg_str = new char[str_len+1];	//	!!! this memory will not be released
		memcpy(arg_str, str, sizeof(char)*(str_len+1));		//	we need to append a '\0' at the end of the string

		//	cut the string into several pieces
		int count = 0;
		char* p = arg_str;
		while(p){
			count ++;

			while(' '!=*p && '\t'!=*p && '\n'!=*p && '\0'!=*p)	//	skip meaningfull characters
				p++;

			if('\0'==*p)	
				break;

			while(' '==*p || '\t'==*p || '\n'==*p)				//	skip meaningfull characters
				*(p++) = '\0';	
		}

		//	create space for string pointers
		*argc = count;
		if( 0 == count ){
			*argv = NULL;
			return 0;
		}
		*argv = new char*[count];				//	!!! this memory will not be released

		//	set string pointers
		p = arg_str;
		for(int i=0; i<count; i++){
			(*argv)[i] = p;

			if(count-1 != i){	//	search for the next start
				while('\0' != *p)	p++;
				while('\0' == *p)	p++;
			}
		}

		return count;
	};

	/**
		This function append args in the string into existing args

		there is memory leak since we don't release existing memory of argv

		new space will be allocated in this function, they are:
		1, argv; 2, argv[argc] (the input argc, if argc induced by str is not 0)
		
		\param	argc	main(int argc, char* argv[]), call Str2Arg(&argc, ...)
		\param	argv	main(int argc, char* argv[]), call Str2Arg(..., &argv, ...)
		\param	str		the input string
		\return			argc
	*/
	static int AppendStr2Arg(int* argc, char*** argv, const char* str){
		int app_argc;
		char** app_argv;
		Str2Arg(&app_argc, &app_argv, str);

		if(0 == app_argc)	//	no arg in the str
			return *argc;

		//	set the new argv list
		char** new_argv = new char*[*argc+app_argc];
		for(int i=0; i<*argc; i++)
			new_argv[i] = (*argv)[i];
		for(int i=0; i<app_argc; i++)
			new_argv[*argc+i] = app_argv[i];

		//	set the data
		*argc += app_argc;
		*argv = new_argv;

		//	delete old data
		delete app_argv;

		return *argc;
	};

	/**
		print the args
	*/
	static void Print(int argc, const char* argv[]){
		for(int i=0; i<argc; i++)
			printf("%d: %s\n", i, argv[i]);
	}

};


}}	//	end namespace yz::utils

#endif	//	__YZ_ARGS_H__