/***********************************************************/
/**	\file
	\brief		Threads
	\author		Yizhong Zhang
	\date		6/19/2012
*/
/***********************************************************/
#ifndef __YZ_THREAD_H__
#define __YZ_THREAD_H__

#ifdef _WINDOWS
#	include <process.h>
#endif

namespace yz{

/**
	start a simple thread

	The thread function should return void and have no arguments
*/
inline void startSimpleThread(void (*func)()){
#ifdef _WINDOWS
	_beginthread((void(*)(void*))func, 0, NULL);
#else
	std::cout << "error: yz::startSimpleThread(), this functions is only implemented in windows" << std::endl;
#endif
}


}	//	namespace yz

#endif	//	__YZ_THREAD_H__