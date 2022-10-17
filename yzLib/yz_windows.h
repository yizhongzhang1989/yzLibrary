/***********************************************************/
/**	\file
	\brief		something only used in windows
	\author		Yizhong Zhang
	\date		5/28/2012
*/
/***********************************************************/
#ifndef __YZ_WINDOWS_H__
#define __YZ_WINDOWS_H__

#pragma warning(push)
#pragma warning(disable: 4267)	//	disable warning of size_t int conversion

//	setting
#include "yzLib_config.h"
#include "yzLib/yz_setting.h"

#include "yzLib/yz_windows/yz_win_utils.h"
#include "yzLib/yz_windows/yz_win_process.h"

#pragma warning(pop)

#endif	//	__YZ_WINDOWS_H__