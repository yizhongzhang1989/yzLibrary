/***********************************************************/
/**	\file
	\brief		Mesh Deformer
	\details	
	\author		Yizhong Zhang
	\date		5/2/2018
*/
/***********************************************************/
#ifndef __YZ_MESH_DEFORM_H__
#define __YZ_MESH_DEFORM_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{  namespace geometry{

/**
	Deform polygon to target shape by calculating numerical derivative.
	Target shape is defined as target length of each edge and target angle of each vertex.
	This class doesn't guarantee the resulting polygon to be self-intersection free!

	How to use:	\n
	1, call SetTarget() to set target length and angle	\n
	2, call Deform() to deform to target	\n
	3, Draw() can be used to display with error	\n

	How to change energy function:	\n
	overload Energy(), CalcVertexEnergy(int vid), CalcCurrVertexEnergy(int vid), 
	to change energy function and energy derivatives.
*/
template <class T>
class MeshPatchDeformer2D : 
	public TriMesh2D<T>,
	public TriMeshEdge 
{
public:
	using TriMesh2D<T>::vertex;
	using TriMesh2D<T>::face;

public:
	/**
	Set the target length of a certain edge
	*/
	void SetTargetEdgeLength(int2 edge, T target_len) {
		if (edge.y < edge.x)
			mySwap(edge.x, edge.y);

		auto iter = target_edge_length.find(edge);
		if (iter == target_edge_length.end())
			target_edge_length.insert(std::pair<int2, T>(edge, target_len));
		else
			iter->second = target_len;
	}

	/**
	Deform the polygon to match target

	\param	max_iterations				max iterations for optimize
	\param	iteration_terminate_thre	if max movement of each vertices is small enough (with respect to min target len), iteration terminate
	\return								actual iterations, -1 if input error
	*/
	int Deform(int max_iterations = 1000, double iteration_terminate_thre = 1e-3) {
		//	check whether edge exist
		if (face.empty())
			return -1;
		if (edge.empty())
			CreateEdge(face);

		//	calculate step size
		double avg_len = 0;
		for (int eid = 0; eid < edge.size(); eid++) {
			avg_len += (vertex[edge[eid].x] - vertex[edge[eid].y]).Length();
		}
		avg_len /= edge.size();
		double step = avg_len * 0.1;

		//	deform step by step
		int iteration_count = 0;
		for (int i = 0; i < max_iterations; i++) {
			double max_move = DeformStep(step);
			iteration_count++;

			if (max_move < avg_len*iteration_terminate_thre)
				break;
		}

		return iteration_count;
	}

	/**
	Draw the polygon, with color indicating error
	*/
	void Draw() {
#ifdef YZ_gl_h
		glColor3f(0, 0, 1);
		opengl::drawMeshEdgeFromFace2D(vertex, face);
#else
		std::cout << "gl.h has to be included in order to use Draw() in MeshPatchDeformer2D" << std::endl;
		return;
#endif
	}

protected:
	/**
	Deform the polygon step.

	This function first calculate derivative of each vertex numerically,
	then go one step toward the derivative.

	\param	delta		delta to calculate jacobian, recommended: min_length / 1000
	\param	step		iteration step: recommended: min_length / 20
	\param	return		biggest step length of all vertices
	*/
	double DeformStep(double step) {
		CalculateJacobian();

		//	go step in jacobian direction
		double max_square_len = 0;
		for (int i = 0; i < vertex.size(); i++) {
			vertex[i] -= jacob[i] * step;

			double len = jacob[i].SquareLength();
			if (len > max_square_len)
				max_square_len = len;
		}

		return sqrt(max_square_len) * step;
	}

	/**
	Calculate Jacobian of the Energy

	Energy is defined as square error sum of all edge length with target
	*/
	void CalculateJacobian() {
		//	clear old jacobian
		jacob.clear();
		jacob.resize(vertex.size(), yz::Vec2d(0, 0));

		for (int eid = 0; eid < edge.size(); eid++) {
			//	skip edges that doesn't exist
			auto iter = target_edge_length.find(edge[eid]);
			if (iter == target_edge_length.end())
				continue;

			int vid0 = edge[eid].x;
			int vid1 = edge[eid].y;
			Vec2<T> v0 = vertex[vid0];
			Vec2<T> v1 = vertex[vid1];
			Vec2<T> r = v1 - v0;
			T edge_len = r.Length();

			// jacobian is calculated mathematically using the following equation
			T coef = 2.0 * (1.0 - iter->second / edge_len);
			jacob[vid0] -= coef * r;
			jacob[vid1] += coef * r;
		}
	}

protected:
	std::unordered_map<int2, T, utils::BitwiseHasher<int2>>		target_edge_length;		///< target length of edge

	std::vector<Vec2d>		jacob;				///< jacobian of total energy of each vertex
};


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_DEFORM_H__