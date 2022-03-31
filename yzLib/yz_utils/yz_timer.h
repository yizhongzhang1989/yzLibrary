/***********************************************************/
/**	\file
	\brief		Timer based on c++ time.h
	\author		Yizhong Zhang
	\date		6/20/2012
*/
/***********************************************************/
#ifndef __YZ_TIMER_H__
#define __YZ_TIMER_H__

#include <time.h>

namespace yz{ namespace utils{

/**
	Base Timer, abstract class

	This abstract class is used as base for implementation
*/
class BaseTimer{
public:
	int status;		///<	status of the timer, one of the value in the following enum
	enum{
		TIMER_STATUS_STOPED		= 0x00,
		TIMER_STATUS_TICKING	= 0x01
	};

public:
	BaseTimer() : status(TIMER_STATUS_STOPED){}

	/**
		Restart the timer

		Stop() then Start
	*/
	inline void Restart(){
		Stop();
		Start();
	}

	/**
		Start the timer

		work only when status is stoped
	*/
	virtual void Start()		= 0;

	/**
		Stop the timer

		work only when status is ticking

		\return		elapsed time from start in ms
	*/
	virtual float Stop()		= 0;

	/**
		Get elapsed time

		\return		if status is ticking, return elapsed time in ms from start to current\n
					if status is stoped, return elapsed time in ms from start to stoped
	*/
	virtual float Elapsed()	= 0;
};

/**
	A timer implementation based on basic C++ time.h
*/
class Timer : public BaseTimer{
public:
	Timer() : BaseTimer(), t_start(0), t_end(0){}

	virtual void Start(){
		if( status == TIMER_STATUS_STOPED ){
			t_start = clock();
			status = TIMER_STATUS_TICKING;
		}
	}

	virtual float Stop(){
		if( status == TIMER_STATUS_TICKING ){
			t_end = clock();
			status = TIMER_STATUS_STOPED;
			return float(t_end-t_start) * 1000 / CLOCKS_PER_SEC;
		}
		return 0;
	}

	virtual float Elapsed(){
		if( status == TIMER_STATUS_TICKING )
			return float(clock()-t_start) * 1000 / CLOCKS_PER_SEC;
		else	//	 status == TIMER_STATUS_STOPED
			return float(t_end-t_start) * 1000 / CLOCKS_PER_SEC;
	}

protected:
	clock_t t_start;	///<	ticks when start counting
	clock_t	t_end;		///<	ticks when end counting
};

/**
	FPS calculator

	to use the fps calculator, you should create a static class
	instance of FPSCalculator, then call GetFPS() every frame.

	CAUTION!!! it is very important to create a static FPSCalculator

	fps value is updated every second
*/
class FPSCalculator{
public:
	/**
		constructor
	*/
	FPSCalculator(){
		counter	= 0;
		fps		= 0;
	}

	/**
		Each frame call this function once, and return the fps value
	*/
	float GetFPS(){
		if( timer.status == timer.TIMER_STATUS_STOPED ){	//	first frame
			timer.Start();
			return 0;
		}

		counter ++;
		if( timer.Elapsed() > 1000 ){	//	update fps value every second
			timer.Stop();
			fps = float(counter) * 1000.0f / timer.Elapsed();
			counter = 0;
			timer.Start();
		}
		
		return fps;
	}

private:
	Timer	timer;
	int		counter;
	float	fps;
};


}}	//	end namespace yz::utils

#endif	//	__YZ_TIMER_H__