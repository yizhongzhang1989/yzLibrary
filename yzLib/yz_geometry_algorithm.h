/***********************************************************/
/**	\file
	\brief		Geometry Algorithm
	\author		Yizhong Zhang
	\date		6/13/2012
*/
/***********************************************************/
#ifndef __YZ_GEOMETRY_ALGORITHM_H__
#define __YZ_GEOMETRY_ALGORITHM_H__

//	setting
#include "yzLib/yz_setting.h"

//	clipping
#include "yzLib/yz_geometry/yz_clipping.h"

//	intersection test of geometry elements
#include "yzLib/yz_geometry/yz_intersection_test.h"
#include "yzLib/yz_geometry/yz_mesh_intersection_test.h"
#include "yzLib/yz_geometry/yz_mesh_nn.h"

//	calculate matrix for mesh
#include "yzLib/yz_geometry/yz_mesh_matrix.h"

//	floodfill
#include "yzLib/yz_geometry/yz_mesh_floodfill.h"
#include "yzLib/yz_geometry/yz_mesh_seperate.h"

//	matching
#include "yzLib/yz_geometry/yz_mesh_matching.h"

//	remeshing
#include "yzLib/yz_geometry/yz_remeshing.h"

//	deformation
#include "yzLib/yz_geometry/yz_mesh_deform.h"

//	mesh generation
#include "yzLib/yz_geometry/yz_mesh_generation.h"

//	UV
#include "yzLib/yz_geometry/yz_mesh_uv.h"

//	ICP
#include "yzLib/yz_geometry/yz_icp.h"

//	fitting
#include "yzLib/yz_geometry/yz_fitting.h"

#endif	//	__YZ_GEOMETRY_ALGORITHM_H__