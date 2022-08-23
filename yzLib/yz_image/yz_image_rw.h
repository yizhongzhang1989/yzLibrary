/***********************************************************/
/**	\file
	\brief		Read & Write Image
	\details	This file allow 3rd party library to be used 
				for image read and write. \n
				A potential bug exist in QT moc when FreeImage 
				is enabled, so I seperate these files to hack the bug
	\author		Yizhong Zhang
	\date		6/21/2016
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RW_H__
#define __YZ_IMAGE_RW_H__

#include "yzLib_config.h"
#include "yzLib/yz_image/yz_image_rw_bmp.h"
#include "yzLib/yz_image/yz_image_rw_pfm.h"

#ifdef yzLib_ENABLE_FreeImage
#include "yzLib/yz_image/yz_image_rw_freeimage.h"
#endif


#endif	//	__YZ_IMAGE_RW_H__