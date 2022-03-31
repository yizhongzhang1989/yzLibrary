/***********************************************************/
/**	\file
	\brief		Interpolation functions on mesh
	\author		Yizhong Zhang
	\date		1/12/2013
*/
/***********************************************************/
#ifndef __YZ_MESH_INTERPOLATION_H__
#define __YZ_MESH_INTERPOLATION_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_geometry/yz_aabb_tree.h"
#include "yzLib/yz_geometry/yz_mesh_intersection_test.h"

namespace yz{  namespace geometry{  namespace remeshing{

//	========================================
///@{
/**	@name Cubic Bezier Interpolation
*/
//	========================================

/**
	calculate control points of cubic bezier triangle

	this algorithm comes from paper:
	Curved PN Triangles

	\param	control_point	control points arranged in sequence
							b300, b030, b003, b210, b120, b021, b012, b102, b201, b111
	\param	v0				corner vertex 0
	\param	v1				corner vertex 1
	\param	v2				corner vertex 2
	\param	normal0			normal corner vertex 0
	\param	normal1			normal corner vertex 1
	\param	normal2			normal corner vertex 2
*/
template<typename T>
void getCubicBezierTriangleControlPoints(std::vector<yz::Vec3<T>>&	control_point, 
										 const yz::Vec3<T>&			v0, 
										 const yz::Vec3<T>&			v1, 
										 const yz::Vec3<T>&			v2, 
										 const yz::Vec3<T>&			normal0, 
										 const yz::Vec3<T>&			normal1, 
										 const yz::Vec3<T>&			normal2){
	yz::Vec3<T> n0 = normal0.Normalize();
	yz::Vec3<T> n1 = normal1.Normalize();
	yz::Vec3<T> n2 = normal2.Normalize();

	control_point.resize(10);
	control_point[0] = v0;
	control_point[1] = v1;
	control_point[2] = v2;
	control_point[3] = (v0*2+v1 - dot(v1-v0, n0) * n0) / 3.0;
	control_point[4] = (v1*2+v0 - dot(v0-v1, n1) * n1) / 3.0;
	control_point[5] = (v1*2+v2 - dot(v2-v1, n1) * n1) / 3.0;
	control_point[6] = (v2*2+v1 - dot(v1-v2, n2) * n2) / 3.0;
	control_point[7] = (v2*2+v0 - dot(v0-v2, n2) * n2) / 3.0;
	control_point[8] = (v0*2+v2 - dot(v2-v0, n0) * n0) / 3.0;
	control_point[9] = (control_point[3] + control_point[4] + control_point[5] +
		control_point[6] + control_point[7] + control_point[8]) / 4.0 - (v0+v1+v2) / 6.0;
}

///@}

}}}	//	namespace yz::geometry::remeshing



#endif	//	__YZ_MESH_INTERPOLATION_H__