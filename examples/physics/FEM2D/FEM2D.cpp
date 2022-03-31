#include <iostream>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>
#include "data_path.h"


//	========================================================
namespace yz{
	namespace physics{

/**
	TriMesh calculated using FEM
*/
template<class T>
class FEMTriMesh : public geometry::TriMesh2D<T>{
public:
	void Draw(){
		glColor3f(0, 0, 1);
		for (int i = 0; i < f.size(); i++){
			yz::opengl::drawLineSegment(vertex[i], vertex[i] + f[i]);
		}

		glColor3f(0, 0, 0);
		Display();

		glColor3f(1, 0, 0);
		for (int i = 0; i < vertex.size(); i++){
			yz::opengl::drawPointAsBall(vertex[i], 0.05);
		}
	}

	void PickingDraw(){
		for (int i = 0; i < vertex.size(); i++){
			yz::opengl::setPickingIndex(i);
			yz::opengl::drawPointAsBall(vertex[i], 0.05);
		}
	}

public:
	void SetupPhysicalValues(){
		f.resize(vertex.size());

		Dm.resize(face.size());
		Bm.resize(face.size());
		W.resize(face.size());
	}

	void PreComputation(){
		k_youngs = 1;	//	Young's modulus
		nu_poisson = 0.4;	//	Poisson ratio

		for (int i = 0; i < face.size(); i++){
			Dm[i] = yz::Matrix2x2<T>(
				vertex[face[i][1]] - vertex[face[i][0]],
				vertex[face[i][2]] - vertex[face[i][0]]
				);
			Bm[i] = Dm[i].Inverse();
			W[i] = fabs(Dm[i].Det()) / 2.0;
		}
	}

	void ComputeElasticForce(){
		T k = k_youngs;
		T nu = nu_poisson;
		T mu = k / (2.0 * (1.0 + nu));
		T lambda = (k*nu) / ((1.0 + nu)*(1.0 - 2.0*nu));
		yz::Matrix2x2<T> I;
		I.SetIdentity();

		for (int i = 0; i < f.size(); i++)
			f[i][0] = f[i][1] = 0;

		for (int i = 0; i < face.size(); i++){
			yz::Matrix2x2<T> Ds(
				vertex[face[i][1]] - vertex[face[i][0]],
				vertex[face[i][2]] - vertex[face[i][0]]
				);
			yz::Matrix2x2<T> F = Ds * Bm[i];
			yz::Matrix2x2<T> E = (F.Transpose() * F - I) * 0.5;
			yz::Matrix2x2<T> P = F * (E * 2 * mu + I * lambda*E.Trace());
			//yz::Matrix2x2<T> P = (F + F.Transpose() - I * 2) * mu + I * (F - I).Trace() * lambda;
			yz::Matrix2x2<T> H = P * Bm[i].Transpose() * (-W[i]);

			yz::Vec2<T> f1(H[0][0], H[1][0]);
			yz::Vec2<T> f2(H[0][1], H[1][1]);

			f[face[i][1]] += f1;
			f[face[i][2]] += f2;
			f[face[i][0]] += (-f1-f2);
		}
	}

	void ComputeForceDifferentials(std::vector<yz::Vec2<T>>& df, const std::vector<yz::Vec2<T>>& dx){
		assert(vertex.size() == dx.size());
		df.resize(vertex.size());

		T k = k_youngs;
		T nu = nu_poisson;
		T mu = k / (2.0 * (1.0 + nu));
		T lambda = (k*nu) / ((1.0 + nu)*(1.0 - 2.0*nu));
		yz::Matrix2x2<T> I;
		I.SetIdentity();

		for (int i = 0; i < df.size(); i++)
			df[i] = yz::Vec2<T>(0, 0);

		//	calculate df
		for (int i = 0; i < face.size(); i++){
			yz::Matrix2x2<T> Ds(
				vertex[face[i][1]] - vertex[face[i][0]],
				vertex[face[i][2]] - vertex[face[i][0]]
				);
			yz::Matrix2x2<T> dDs(
				dx[face[i][1]] - dx[face[i][0]],
				dx[face[i][2]] - dx[face[i][0]]
				);
			yz::Matrix2x2<T> F = Ds * Bm[i];
			yz::Matrix2x2<T> dF = dDs * Bm[i];
			yz::Matrix2x2<T> E = (F.Transpose() * F - I) * T(0.5);
			yz::Matrix2x2<T> dE = (dF.Transpose() * F + F.Transpose() * dF) * T(0.5);
			yz::Matrix2x2<T> dP = dF * (E * 2 * mu + I * lambda*E.Trace()) + F * (dE * 2 * mu + I * lambda*dE.Trace());
			yz::Matrix2x2<T> dH = dP * Bm[i].Transpose() * (-W[i]);

			yz::Vec2<T> df1(dH[0][0], dH[1][0]);
			yz::Vec2<T> df2(dH[0][1], dH[1][1]);
			df[face[i][1]] += df1;
			df[face[i][2]] += df2;
			df[face[i][0]] += -df1 - df2;
		}
	}

	void NewtonStep(){
		ComputeElasticForce();

		std::vector<yz::Vec2<T>>	x;
		x.resize(vertex.size(), yz::Vec2<T>(0, 0));

		std::vector<yz::Vec2<T>>	r = f;
		std::vector<yz::Vec2<T>>	p = r;
		T rTr = XTY(r, r);

		for (int k = 0; k < 10; k++){
			std::vector<yz::Vec2<T>>	Ap;
			ComputeForceDifferentials(Ap, p);
			T pTAp = XTY(p, Ap);
			T alpha = rTr / pTAp;

			for (int i = 0; i < x.size(); i++)
				x[i] += alpha * p[i];
			for (int i = 0; i < r.size(); i++)
				r[i] -= alpha * Ap[i];

			T old_rTr = rTr;
			rTr = XTY(r, r);

			//	exit loop if r is small enough
			std::cout << rTr << std::endl;

			T beta = old_rTr / rTr;
			for (int i = 0; i < p.size(); i++)
				p[i] = r[i] + beta * p[i];
		}

		for (int i = 0; i < vertex.size(); i++)
			vertex[i] -= x[i];
	}

	T XTY(const std::vector<yz::Vec2<T>>& X, const std::vector<yz::Vec2<T>>& Y){
		assert(X.size() == Y.size());

		T sum = 0;
		for (int i = 0; i < X.size(); i++){
			sum += X[i].x*Y[i].x + X[i].y*Y[i].y;
		}

		return sum;
	}


public:
	T								k_youngs;	//	Young's Module
	T								nu_poisson;	//	Poisson ratio

	std::vector<yz::Vec2<T>>		f;			//	force on each vertex

	std::vector<yz::Matrix2x2<T>>	Dm;			//	vector matrix of reference mesh
	std::vector<yz::Matrix2x2<T>>	Bm;			//	inverse of Dm
	std::vector<T>					W;			//	
};


}}	//	yz::physics
//	========================================================

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow2D		win2d;

yz::physics::FEMTriMesh<double>	mesh;

int			picked_vertex_id = -1;
yz::Vec2d	picked_point;

void draw(){
	mesh.Draw();
}

void picking_draw(){
	mesh.PickingDraw();
}

void process_picking(int idx){
	GLdouble z;
	yz::opengl::getWorldCoordinate(picked_point[0], picked_point[1], z, win2d.old_x, win2d.old_y, 0);

	if (idx >= 0 && idx < mesh.vertex.size())
		picked_vertex_id = idx;
	else
		picked_vertex_id = -1;
}

void mouse_func(int button, int state, int x, int y){
	win2d.mouse_state[button] = state;
	win2d.old_x = x;
	win2d.old_y = y;
	win2d.modifier = glutGetModifiers();

	//	if picking_displayFunc is set and Ctrl is pressed, then perform color picking
	if (win2d.picking_displayFunc && (win2d.modifier&GLUT_ACTIVE_CTRL) && win2d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN){
		//	check whether picking is legal
		int picking_legal = 0;

		//	picking_displayFunc has been changed, then we know the user has specifically written the function
		if (win2d.picking_displayFunc != win2d.DefaultPickingDisplayFunc)
			picking_legal = 1;
		//	picking_displayFunc is not changed, but picking_draw is changed, legal
		else if (picking_draw && picking_draw != win2d.DefaultPickingDraw)
			picking_legal = 1;
		//	picking is default function, then only display functions are default, can picking be legal
		else if (win2d.displayFunc == win2d.DefaultDisplayFunc && draw == win2d.DefaultDraw)
			picking_legal = 1;

		if (picking_legal){
			(*win2d.picking_displayFunc)();
			int index = win2d.PickingIndex(x, y);
			if (win2d.process_picking)	(*win2d.process_picking)(index);
		}
	}
	else if (win2d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_UP){
		picked_vertex_id = -1;
	}


	glutPostRedisplay();
}

void motion_func(int x, int y){
	if (win2d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN){			//	move
		if (picked_vertex_id >= 0){
			//	move selected vertex
			double u, v, d;
			yz::opengl::getWindowCoordinate(u, v, d, picked_point.x, picked_point.y, 0);
			yz::Vec2d xy;
			double z;
			yz::opengl::getWorldCoordinate(xy.x, xy.y, z, x, y, d);

			mesh.vertex[picked_vertex_id] += xy - picked_point;
			picked_point = xy;

			mesh.ComputeElasticForce();
		}
		else{
			GLdouble width = win2d.right - win2d.left;
			GLdouble height = win2d.top - win2d.bottom;
			GLdouble move_x = -width * GLdouble(x - win2d.old_x) / win2d.win_width;
			GLdouble move_y = height * GLdouble(y - win2d.old_y) / win2d.win_height;
			win2d.left += move_x;
			win2d.right += move_x;
			win2d.bottom += move_y;
			win2d.top += move_y;
			win2d.DefaultReshapeFunc(win2d.win_width, win2d.win_height);
		}
	}
	else if (win2d.mouse_state[GLUT_MIDDLE_BUTTON] == GLUT_DOWN){	//	move
		GLdouble width = win2d.right - win2d.left;
		GLdouble height = win2d.top - win2d.bottom;
		GLdouble move_x = -width * GLdouble(x - win2d.old_x) / win2d.win_width;
		GLdouble move_y = height * GLdouble(y - win2d.old_y) / win2d.win_height;
		win2d.left += move_x;
		win2d.right += move_x;
		win2d.bottom += move_y;
		win2d.top += move_y;
		win2d.DefaultReshapeFunc(win2d.win_width, win2d.win_height);
	}
	else if (win2d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN){		//	scale
		GLdouble width = win2d.right - win2d.left;
		GLdouble height = win2d.top - win2d.bottom;
		GLdouble center_x = (win2d.left + win2d.right) / 2;
		GLdouble center_y = (win2d.bottom + win2d.top) / 2;
		GLdouble zoom = 1 + GLdouble(y - win2d.old_y) / 200;
		width *= zoom;
		height *= zoom;
		win2d.left = center_x - width / 2;
		win2d.right = center_x + width / 2;
		win2d.bottom = center_y - height / 2;
		win2d.top = center_y + height / 2;
		win2d.DefaultReshapeFunc(win2d.win_width, win2d.win_height);
	}

	win2d.old_x = x;
	win2d.old_y = y;
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y){
	switch (key){
	case 27:
		exit(0);
	case ' ':
		mesh.NewtonStep();
		break;
	}

	glutPostRedisplay();
}

int main(){
	char obj_filename[1024];
	sprintf(obj_filename, "%s/sheet.obj", data_path);

	mesh.ReadMeshFromFile(obj_filename, 'z');
	mesh.SetupPhysicalValues();
	mesh.PreComputation();
	mesh.ComputeElasticForce();

	////	validate force differential
	//{
	//	std::vector<yz::Vec2d>	ref_vertex = mesh.vertex;
	//	std::vector<yz::Vec2d>	df;

	//	std::vector<yz::Vec2d> dx;
	//	dx.resize(mesh.vertex.size());
	//	dx[0] = yz::Vec2d(0, 0.001);

	//	for (int i = 0; i < mesh.vertex.size(); i++)
	//		mesh.vertex[i] += dx[i];
	//	mesh.ComputeElasticForce();
	//	for (int i = 0; i < mesh.f.size(); i++)
	//		std::cout << mesh.f[i];
	//	std::cout << std::endl;

	//	mesh.vertex = ref_vertex;
	//	mesh.ComputeForceDifferentials(df, dx);
	//	for (int i = 0; i < df.size(); i++)
	//		std::cout << df[i];
	//	std::cout << std::endl;
	//}

	win2d.SetDraw(draw);
	win2d.SetPickingDraw(picking_draw);
	win2d.SetProcessPicking(process_picking);
	win2d.mouseFunc = mouse_func;
	win2d.motionFunc = motion_func;
	win2d.keyboardFunc = keyboard;
	win2d.CreateGLUTWindow();

	manager.EnterMainLoop();
}
