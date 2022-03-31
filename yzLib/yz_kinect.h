/***********************************************************/
/**	\file
	\brief		Microsoft Kinect Related Classes and Functions
	\author		Yizhong Zhang
	\date		11/9/2012
*/
/***********************************************************/
#ifndef __YZ_KINECT_H__
#define __YZ_KINECT_H__

//	setting
#include "yzLib/yz_setting.h"


//	the base of kinect, just virtual class
#include "yzLib/yz_kinect/yz_kinect_base.h"

//	implimentation of kinect classes are driver dependent
#ifdef	YZ_NuiApi_h		//	using microsoft kinect sdk
#include "yzLib/yz_kinect/yz_microsoft_kinect.h"
#endif


#endif	//	__YZ_KINECT_H__