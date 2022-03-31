#include <iostream>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>
using namespace std;

typedef float REAL;

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow2D		win2d;

yz::geometry::TriMesh2D<REAL>	mesh;
std::vector<yz::Vec2<REAL>>		points;

void draw(){
	glColor3f(0, 0, 1);
	mesh.Display();

	glColor3f(1, 0, 0);
	for(int i=0; i<points.size(); i++){
		yz::opengl::drawPointAsBall(points[i], 0.003);
	}
}

int main(){
	mesh.ReadMeshFromFile("mesh.obj");

	yz::geometry::remeshing::samplingPoissonDiskOn2DMesh(points, mesh.vertex, mesh.face, REAL(0.01), "robertbridson");
	//yz::geometry::remeshing::samplingPoissonDisk(points, 0.02f, yz::Vec2f(0.5, 0.5), yz::Vec2f(1, 1) );

	//yz::geometry::remeshing::samplingPointsOn2DMeshStaggered(points, mesh.vertex, mesh.face, REAL(0.02));

	win2d.SetDraw(draw);
	win2d.CreateGLUTWindow();

	manager.AddIdleFunc(win2d.idleFunc);
	manager.EnterMainLoop();
}