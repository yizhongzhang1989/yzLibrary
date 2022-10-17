/***********************************************************/
/**	\file
	\brief		Geometry
	\author		Yizhong Zhang
	\date		6/1/2012
*/
/***********************************************************/
#ifndef __YZ_GEOMETRY_H__
#define __YZ_GEOMETRY_H__

#pragma warning(push)
#pragma warning(disable: 4267)	//	disable warning of size_t int conversion

//	setting
#include "yzLib/yz_setting.h"

//	topology
#include "yzLib/yz_geometry_topology.h"

//	mesh
#include "yzLib/yz_geometry_mesh.h"

//	algorithm
#include "yzLib/yz_geometry_algorithm.h"

//	polygon
#include "yzLib/yz_geometry/yz_polygon.h"

//	field
#include "yzLib/yz_geometry/yz_field.h"

//	tetrahedron
#include "yzLib/yz_geometry/yz_tetrahedron.h"

#pragma warning(pop)

#endif	//	__YZ_GEOMETRY_H__