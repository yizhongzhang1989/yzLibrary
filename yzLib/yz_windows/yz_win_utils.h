/***********************************************************/
/**	\file
	\brief		utilities of windows
	\author		Yizhong Zhang
	\date		6/20/2012
*/
/***********************************************************/
#ifndef __YZ_WIN_UTILS_H__
#define __YZ_WIN_UTILS_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_windows_h
#	error yz_win_process.h must be included after windows.h
#endif

#include <iostream>
#include "yzLib/yz_utils/yz_timer.h"

namespace yz{

/**
	namespace contain windows related classes and functions
*/
namespace windows{
/**
	Get last error code, convert to string and print
*/
inline void PrintLastError(const char* info = NULL){
#ifndef BE_QUIET
	LPSTR lpMsgBuf;
	DWORD dw = GetLastError(); 

	FormatMessageA(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		dw,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPSTR)&lpMsgBuf,
		0, 
		NULL );

	if( info )
		std::cout << info << ", Last Error " << dw << ": " << lpMsgBuf << std::endl;
	else
		std::cout << "Last Error " << dw << ": " << lpMsgBuf << std::endl;

	LocalFree(lpMsgBuf);
#endif
}


/**
	A precise timer based on windows API
*/
class WinTimer : public utils::BaseTimer{
public:
	WinTimer() : BaseTimer() {
		t_start.QuadPart	= 0;
		t_end.QuadPart		= 0;
		if( !QueryPerformanceFrequency(&frequency) )	PrintLastError();
		frequency.QuadPart /= 1000;	//	change the unit to counter/ms
	}
	inline void Start(){
		if( status == TIMER_STATUS_STOPED ){
			if( !QueryPerformanceCounter(&t_start) )	PrintLastError();
			status = TIMER_STATUS_TICKING;
		}
	}
	inline float Stop(){
		if( status == TIMER_STATUS_TICKING ){
			if( !QueryPerformanceCounter(&t_end) )	PrintLastError();
			status = TIMER_STATUS_STOPED;
			return float(t_end.QuadPart - t_start.QuadPart) / frequency.QuadPart;
		}
		return 0;
	}
	inline float Elapsed(){
		if( status == TIMER_STATUS_TICKING ){
			LARGE_INTEGER t_now;
			if( !QueryPerformanceCounter(&t_now) )	PrintLastError();
			return float(t_now.QuadPart - t_start.QuadPart) / frequency.QuadPart;
		}
		else	//	status == TIMER_STATUS_STOPED
			return float(t_end.QuadPart - t_start.QuadPart) / frequency.QuadPart;
	}

protected:
	LARGE_INTEGER t_start;		///<	start time
	LARGE_INTEGER t_end;		///<	end time
	LARGE_INTEGER frequency;	///<	frequency of the counter
};


/**
	A real time fps calculator using WinTimer

	to use the fps calculator, you should create a static class
	instance of FPSCalculator, then call GetFPS() every frame.

	CAUTION!!! it is very important to create a static FPSCalculator

	fps is calculated using the time interval of each frame
*/
class RealTimeFPSCalculator{
public:
	/**
		Get the fps calculated real time
	*/
	float GetFPS(){
		if( timer.status == timer.TIMER_STATUS_STOPED ){
			timer.Start();
			return 0;
		}

		timer.Stop();
		float fps = 1000.0f / timer.Elapsed();
		timer.Start();
		return fps;
	}
private:
	WinTimer timer;
};

}	//	namespace windows
}	//	namespace yz

#endif	//	__YZ_WIN_UTILS_H__