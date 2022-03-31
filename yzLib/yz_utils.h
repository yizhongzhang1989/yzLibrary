/***********************************************************/
/**	\file
	\brief		a collection of util functions to make programming easier
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_UTILS_H__
#define __YZ_UTILS_H__

//	setting
#include "yzLib/yz_setting.h"

#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_math/yz_filter.h"

#include "yzLib/yz_utils/yz_args.h"
#include "yzLib/yz_utils/yz_color_bar.h"
#include "yzLib/yz_utils/yz_mem_pattern.h"
#include "yzLib/yz_utils/yz_string_utils.h"
#include "yzLib/yz_utils/yz_timer.h"
#include "yzLib/yz_utils/yz_file.h"
#include "yzLib/yz_utils/yz_reorder.h"
#include "yzLib/yz_utils/yz_array_check.h"

#include "yzLib/yz_utils/yz_image_rw_utils.h"


#ifdef YZ_imdebug_h	//	visualize sparse matrix
#	include "yzLib/yz_utils/yz_sparse_visualize.h"
#endif


#endif	//	__YZ_UTILS_H__