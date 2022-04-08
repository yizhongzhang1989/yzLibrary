#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include "yzLib/yz_lib.h"


yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow2D		win2d;

std::vector<yz::Vec2d>			point;
yz::geometry::AABBTree<yz::Vec2d> aabb;

int depth = 0;
int center_id = 0;

void draw() {
	yz::opengl::drawXYZAxis();

	glPointSize(5);
	glColor3f(0, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 0; i < point.size(); i++) {
		glVertex2d(point[i].x, point[i].y);
	}
	glEnd();
	glPointSize(1);

	//		
	std::vector<yz::Vec2d>	np;
	std::vector<int>		npid;
	yz::geometry::getKNearestPointsOnVertices(np, npid, 10, point[center_id], point, aabb);

	glColor3f(1, 0, 0);
	yz::opengl::drawCircle(point[center_id], 0.02);

	glColor3f(0, 1, 0);
	for (int i = 0; i < npid.size(); i++)
		yz::opengl::drawCircle(point[npid[i]], 0.03 + 0.01*i);

}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case 'q':
		center_id++;
		break;
	case 'a':
		center_id--;
		break;
	}

	glutPostRedisplay();
}


int main() {

	yz::geometry::TriMesh<double>	mesh;
	mesh.ReadMeshFromFile("non_ground.obj");

	point.resize(mesh.vertex.size());
	for (int i = 0; i < point.size(); i++) {
		point[i] = yz::Vec2d(mesh.vertex[i].x, mesh.vertex[i].y);
	}


	//point.resize(100);
	//for (int i = 0; i < point.size(); i++) {
	//	point[i].x = yz::randFloatingPointNumber(-1, 1);
	//	point[i].y = yz::randFloatingPointNumber(-1, 1);
	//}

	aabb.BuildVertexAABBTree(point);


	win2d.back_ground_red = 1;
	win2d.back_ground_green = 1;
	win2d.back_ground_blue = 1;
	win2d.keyboardFunc = keyboard;
	win2d.SetDraw(draw);
	win2d.CreateGLUTWindow();

	manager.EnterMainLoop();
}
