/***********************************************************/
/**	\file
	\brief		Example of ICP using Eigen
	\details	ICP requires SVD. In this example, we use the 
				SVD module of Eigen to calculate the transform.
	\author		Yizhong Zhang
	\date		3/29/2022
*/
/***********************************************************/
#include <iostream>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <Eigen/Dense>	//	include Eigen before yzLib, so that icp will use Eigen
#include "yzLib/yz_lib.h"

yz::opengl::DemoWindowManager	manager;
yz::opengl::GLUTWindow3D<0>		win3d;

std::vector<yz::Vec3d>			point_source;
std::vector<yz::Vec3d>			point_target;
yz::Matrix4x4d					trans_s2t;

void draw_3d() {
	yz::opengl::drawXYZAxis();

	glColor3f(1, 0, 0);
	for (int i = 0; i < point_target.size(); i++)
		yz::opengl::drawPointAsCube(point_target[i], 0.05);

	glColor3f(0, 1, 0);
	for (int i = 0; i < point_source.size(); i++)
		yz::opengl::drawPointAsWireSphere(point_source[i], 0.05);

	glColor3f(1, 1, 0);
	yz::Matrix4x4d	s2t = trans_s2t;
	s2t.SetTranspose();
	glPushMatrix();
	glMultMatrixd(s2t.data[0]);
	for (int i = 0; i < point_source.size(); i++)
		yz::opengl::drawPointAsWireSphere(point_source[i], 0.05);
	glPopMatrix();
}

void print_3d() {
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "target points");

	glColor3f(0, 1, 0);
	yz::opengl::printInfo(0, 30, "source points");

	glColor3f(1, 1, 0);
	yz::opengl::printInfo(0, 60, "transformed source points");
}

int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	//	prepare transform matrix
	yz::Matrix4x4d trans;
	yz::Matrix4x4d rot_x, rot_y, rot_z, trans_xyz;
	rot_x.SetRotationDeg(yz::Vec3d(1, 0, 0), 10);
	rot_y.SetRotationDeg(yz::Vec3d(0, 1, 0), 10);
	rot_z.SetRotationDeg(yz::Vec3d(0, 0, 1), 10);
	trans_xyz.SetTranslation(0.1, 0.2, 0.3);
	trans = rot_x * rot_y * rot_z * trans_xyz;
	double noise_level = 1e-2;

	//	create point array
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 3; j++) {
			yz::Vec3d p_t((i - 2) * 0.3, (j - 1) * 0.3, 0);
			yz::Vec3d p_s = trans.TransformVector(p_t)
				+ yz::Vec3d(
					yz::randNumber(-noise_level, noise_level), 
					yz::randNumber(-noise_level, noise_level), 
					yz::randNumber(-noise_level, noise_level));
			point_target.push_back(p_t);
			point_source.push_back(p_s);
		}
	}

	//	calculate transform
	yz::geometry::icp::calculateTransformPointToPoint(
		trans_s2t, point_source, point_target);

	std::cout << "transform source to target" << std::endl;
	std::cout << trans_s2t << std::endl;

	//	opengl context
	win3d.use_arcball_flag = 1;
	win3d.SetDraw(draw_3d);
	win3d.SetDrawAppend(print_3d);
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();

	return 0;
}