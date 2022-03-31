/***********************************************************/
/**	\file
	\brief		all kinds of mesh and mesh tools
	\author		Yizhong Zhang
	\date		6/1/2012
*/
/***********************************************************/
#ifndef __YZ_GEOMETRY_MESH_H__
#define __YZ_GEOMETRY_MESH_H__

//	setting
#include "yzLib/yz_setting.h"

//	mesh utils (functions)
#include "yzLib/yz_geometry/yz_create_topology.h"	//	create connectivity information
#include "yzLib/yz_geometry/yz_mesh_subd.h"			//	subdivide
#include "yzLib/yz_geometry/yz_mesh_rw.h"			//	mesh read / write file
#include "yzLib/yz_geometry/yz_mesh_normal.h"		//	calculate normal of mesh
#include "yzLib/yz_geometry/yz_mesh_smoothing.h"	//	smooth functions
#include "yzLib/yz_geometry/yz_mesh_transform.h"	//	transform of mesh
#include "yzLib/yz_geometry/yz_mesh_repair.h"		//	repair a mesh
#include "yzLib/yz_geometry/yz_mesh_area.h"			//	calculate area
#include "yzLib/yz_geometry/yz_mesh_utils.h"		//	not grouped utility functions

#include "yzLib/yz_geometry/yz_mesh_curvature.h"	//	calculate curvature
#include "yzLib/yz_geometry/yz_bonding_box.h"		//	calculate bonding box
#include "yzLib/yz_geometry/yz_aabb_tree.h"			//	calculate AABB tree
#include "yzLib/yz_geometry/yz_sphere.h"			//	calculate sphere
#include "yzLib/yz_geometry/yz_sphere_tree.h"		//	calculate sphere tree

//	mesh structures (classes & functions)
#include "yzLib/yz_geometry/yz_tri_mesh.h"
#include "yzLib/yz_geometry/yz_quad_mesh.h"
#include "yzLib/yz_geometry/yz_tri_mesh_2d.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"
#include "yzLib/yz_geometry/yz_offset_tri_mesh.h"
#include "yzLib/yz_geometry/yz_mesh_texture.h"

#if (defined(YZ_glut_h) || defined(YZ_freeglut_h))
#	include "yzLib/yz_geometry/yz_texture_tri_mesh.h"
#endif

#endif	//	__YZ_GEOMETRY_MESH_H__