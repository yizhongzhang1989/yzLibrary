/***********************************************************/
/**	\file
	\brief		a collection of physics functions
	\author		Yizhong Zhang
	\date		10/17/2012
*/
/***********************************************************/
#ifndef __YZ_PHYSICS_H__
#define __YZ_PHYSICS_H__

#pragma warning(push)
#pragma warning(disable: 4267)	//	disable warning of size_t int conversion

//	setting
#include "yzLib/yz_setting.h"

//	mass pring
#include "yzLib/yz_physics/yz_mass_spring.h"

//	FEM
#include "yzLib/yz_physics/yz_fem/yz_fem_tri_mesh.h"

//	position based dynamics
#include "yzLib/yz_physics/yz_position_based_dynamics.h"
#include "yzLib/yz_physics/yz_wrinkle_mesh.h"

//	folds and wrinkles
#include "yzLib/yz_physics/yz_bending_force.h"

#pragma warning(pop)

#endif	//	__YZ_PHYSICS_H__