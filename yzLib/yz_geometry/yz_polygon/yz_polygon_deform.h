/***********************************************************/
/**	\file
	\brief		Polygon Deformer
	\details	
	\author		Yizhong Zhang
	\date		3/29/2018
*/
/***********************************************************/
#ifndef __YZ_POLYGON_DEFORM_H__
#define __YZ_POLYGON_DEFORM_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_color_bar.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display mesh
#endif

namespace yz{  namespace geometry{  namespace polygon{

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
class PolygonDeformer {
public:
	std::vector<Vec2<T>>	v;		///< vertex list of the polygon

public:
	/**
		Setup target length and angle

		\param	target_edge_length			target length of each edge
		\param	target_vertex_angle_rad		target angle (in radius) of each vertex
	*/
	void SetTarget(
		const std::vector<T>& target_edge_length,
		const std::vector<T>& target_vertex_angle_rad
	) {
		//	check input
		if (target_edge_length.size() != target_vertex_angle_rad.size()) {
			std::cout << "error: PolygonDeformer::SetupTarget, target size don't match" << std::endl;
			return;
		}
		if (target_edge_length.size() < 3) {
			std::cout << "error: PolygonDeformer::SetupTarget, less than 3 vertices" << std::endl;
			return;
		}

		//	setup
		target_length = target_edge_length;
		target_angle = target_vertex_angle_rad;

		//	initialize the polygon to a circle
		InitPolygon();
	}

	/**
		Deform the polygon to match target

		\param	max_iterations				max iterations for optimize
		\param	iteration_terminate_thre	if max movement of each vertices is small enough (with respect to min target len), iteration terminate
		\return								actual iterations, -1 if input error
	*/
	int Deform(int max_iterations = 1000, double iteration_terminate_thre = 1e-3) {
		//	check input
		if (target_length.size() != target_angle.size()) {
			std::cout << "error: PolygonDeformer::Deform, target not properly set" << std::endl;
			return -1;
		}
		if (target_length.size() < 3) {
			std::cout << "error: PolygonDeformer::Deform, target not properly set" << std::endl;
			return -1;
		}

		//	calculate step size
		double min_len = target_length[0];
		for (int i = 1; i < target_length.size(); i++)
			if (target_length[i] < min_len)
				min_len = target_length[i];
		double delta = min_len * 0.001;
		double step = min_len * 0.05;

		//	deform step by step
		int iteration_count = 0;
		for (int i = 0; i < max_iterations; i++) {			
			double max_move = DeformStep(delta, step);
			iteration_count++;

			//	check whether terminate
			if (curr_energy < initial_energy*1e-5)
				break;
			if (max_move < min_len*iteration_terminate_thre)
				break;
		}

		return iteration_count;
	}

	/**
		Draw the polygon, with color indicating error
	*/
	void Draw() {
#ifdef YZ_gl_h
		//	get min length
		double min_len = target_length[0];
		for (int i = 1; i < target_length.size(); i++)
			if (target_length[i] < min_len)
				min_len = target_length[i];

		//	draw vertex
		for (int i = 0; i < v.size(); i++) {
			float color[3];
			utils::convertToColorSimple(color, curr_angle[i] - target_angle[i], -YZ_PI / 36, YZ_PI / 36);	//	5 degree
			glColor3fv(color);

			yz::opengl::drawPointAsCircle(v[i], min_len*0.05);
		}

		//	draw edge
		for (int i = 0; i < v.size(); i++) {
			float color[3];
			utils::convertToColorSimple(color, curr_length[i] - target_length[i], min_len*0.9, min_len*1.1);	//	10% min_len
			glColor3fv(color);

			yz::opengl::drawLineSegment(v[i], v[(i + 1) % v.size()]);
		}
#else
		std::cout << "gl.h has to be included in order to use Draw() in PolygonDeformer" << std::endl;
		return;
#endif
	}

protected:
	/**
		initialize the polygon to a circle
	*/
	void InitPolygon() {
		int v_num = target_length.size();

		//	calculate radius, we assume the length of the circle is the same as the polygon
		double radius = 0;
		for (int i = 0; i < v_num; i++)
			radius += target_length[i];
		radius /= (2 * YZ_PI);

		//	set the polygon as a uniformly distributed circle
		v.resize(v_num);
		for (int i = 0; i < v_num; i++) {
			double theta = 2 * YZ_PI * i / v_num;
			v[i].x = radius * cos(theta);
			v[i].y = radius * sin(theta);
		}

		//	calculate current length and angle
		curr_length.resize(target_length.size());
		curr_angle.resize(target_angle.size());
		CalculateCurr();

		//	calculate energy
		jacob.resize(v_num);
		initial_energy = Energy();
		curr_energy = initial_energy;
	}

	/**
		Deform the polygon step.

		This function first calculate derivative of each vertex numerically,
		then go one step toward the derivative.

		\param	delta		delta to calculate jacobian, recommended: min_length / 1000
		\param	step		iteration step: recommended: min_length / 20
		\param	return		biggest step length of all vertices
	*/
	double DeformStep(double delta, double step) {
		int v_num = v.size();

		//	calculate jacobian numerically
		for (int i = 0; i < v_num; i++) {
			double org_energy = CalcCurrVertexEnergy(i);
			Vec2d org_v = v[i];

			v[i].x += delta;
			double dx_energy = (CalcVertexEnergy(i) - org_energy) / delta;
			v[i] = org_v;

			v[i].y += delta;
			double dy_energy = (CalcVertexEnergy(i) - org_energy) / delta;
			v[i] = org_v;

			jacob[i].x = dx_energy;
			jacob[i].y = dy_energy;
		}

		//	normalize jacobian
		Vec2d offset(0, 0);
		for (int i = 0; i < v_num; i++) {
			offset += jacob[i];
		}
		offset /= v_num;
		for (int i = 0; i < v_num; i++) {
			jacob[i] -= offset;
		}

		//	go step in jacobian direction
		double max_square_len = 0;
		for (int i = 0; i < v_num; i++) {
			v[i] -= jacob[i] * step;

			double len = jacob[i].SquareLength();
			if (len > max_square_len)
				max_square_len = len;
		}

		//	calculate current
		CalculateCurr();

		return sqrt(max_square_len) * step;
	}

	/**
		calculate current length, angle and energy
	*/
	void CalculateCurr() {
		int v_num = v.size();

		for (int i = 0; i < v_num; i++) {
			curr_length[i] = (v[(i + 1) % v_num] - v[i]).Length();
		}

		for (int i = 0; i < v_num; i++) {
			curr_angle[i] = YZ_PI - angleRadBetweenVectors(v[i] - v[(i + v_num - 1) % v_num], v[(i + v_num + 1) % v_num] - v[i]);
		}

		curr_energy = Energy();
	}

	/**
		Energy of the polygon, can be overloaded

		we use energy defined in the paper: Designing Inflatable Structure
	*/
	virtual double Energy() {
		int v_num = v.size();

		double energy_sum = 0;
		for (int i = 0; i < v_num; i++) {
			double diff = curr_length[i] - target_length[i];
			energy_sum += diff * diff / target_length[i];
			diff = curr_angle[i] - target_angle[i];
			energy_sum += diff * diff * (target_length[i] + target_length[(i + v_num - 1) % v_num]) * 0.5;
		}

		return energy_sum;
	}

	/**
		Calculate energy associated with this vertex

		This function calculate energy using v only
	*/
	virtual inline double CalcVertexEnergy(int vid) {
		int v_num = v.size();

		int i = vid;
		int i_prev = (i + v_num - 1) % v_num;
		int i_next = (i + 1) % v_num;

		//	record original length and angle
		double org_curr_len_prev = curr_length[i_prev];
		double org_curr_len = curr_length[i];
		double org_curr_angle_prev = curr_angle[i_prev];
		double org_curr_angle = curr_angle[i];
		double org_curr_angle_next = curr_angle[i_next];

		//	calculate current length and angle
		curr_length[i_prev] = (v[(i_prev + 1) % v_num] - v[i_prev]).Length();
		curr_length[i] = (v[(i + 1) % v_num] - v[i]).Length();
		curr_angle[i_prev] = YZ_PI - angleRadBetweenVectors(v[i_prev] - v[(i_prev + v_num - 1) % v_num], v[(i_prev + v_num + 1) % v_num] - v[i_prev]);
		curr_angle[i] = YZ_PI - angleRadBetweenVectors(v[i] - v[(i + v_num - 1) % v_num], v[(i + v_num + 1) % v_num] - v[i]);
		curr_angle[i_next] = YZ_PI - angleRadBetweenVectors(v[i_next] - v[(i_next + v_num - 1) % v_num], v[(i_next + v_num + 1) % v_num] - v[i_next]);

		//	calculate energy
		double energy = CalcCurrVertexEnergy(vid);

		//	recover
		curr_length[i_prev] = org_curr_len_prev;
		curr_length[i] = org_curr_len;
		curr_angle[i_prev] = org_curr_angle_prev;
		curr_angle[i] = org_curr_angle;
		curr_angle[i_next] = org_curr_angle_next;

		return energy;
	}

	/**
		Calculate energy associated with this vertex.

		This function assumes that curr_length and curr_angle are up to date with v
	*/
	virtual inline double CalcCurrVertexEnergy(int vid) {
		int v_num = v.size();

		int i = vid;
		int i_prev = (i + v_num - 1) % v_num;
		int i_next = (i + 1) % v_num;

		//	calculate length energy
		double diff = curr_length[i] - target_length[i];
		double len_energy = diff * diff / target_length[i];
		diff = curr_length[i_prev] - target_length[i_prev];
		len_energy += diff * diff / target_length[i_prev];

		//	calculate angle energy
		diff = curr_angle[i] - target_angle[i];
		double angle_energy = diff * diff * (target_length[i] + target_length[i_prev]) * 0.5;
		diff = curr_angle[i_prev] - target_angle[i_prev];
		angle_energy += diff * diff * (target_length[i_prev] + target_length[(i + v_num - 2) % v_num]) * 0.5;
		diff = curr_angle[i_next] - target_angle[i_next];
		angle_energy += diff * diff * (target_length[i] + target_length[i_next]) * 0.5;

		return len_energy + angle_energy;
	}

protected:
	std::vector<T>			target_length;		///< target length of each edge (corresponding vertex and next vertex)
	std::vector<T>			target_angle;		///< target angle of each vertex (the angle inside the polygon)

	std::vector<T>			curr_length;		///< current length
	std::vector<T>			curr_angle;			///< current angle
	double					curr_energy;		///< current total energy

	std::vector<Vec2d>		jacob;				///< jacobian of total energy of each vertex
	double					initial_energy;		///< initial total energy

};


}}}	//	namespace yz::geometry::polygon

#endif	//	__YZ_POLYGON_DEFORM_H__