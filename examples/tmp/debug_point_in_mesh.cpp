#include <iostream>
#include <string>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::TriMesh<double>	mesh;

std::vector<yz::Vec3d>			point;
std::vector<int>				inside_flag;

void draw() {
	for (int i = 0; i < point.size(); i++) {
		if (inside_flag[i])
			glColor3f(1, 0, 0);
		else
			glColor3f(1, 1, 0);
		yz::opengl::drawPointAsCube(point[i], 0.01);
	}
}

int main() {
	mesh.ReadMeshFromFile("icosphere.obj");
	yz::geometry::applyNoiseToVertices(mesh.vertex, 0.0001);
	mesh.CalculateNormals();

	point = mesh.vertex;
	for (int i = 0; i < point.size(); i++) {
		point[i] -= mesh.vertex_normal[i] * 0.1;
	}

	for (int i = 0; i < point.size(); i++) {
		point[i] *= 0.000001;
		mesh.vertex[i] *= 0.000001;
	}

	yz::geometry::AABBTree3D<double>	aabb;
	aabb.BuildTriangleAABBTree(mesh.vertex, mesh.face);

	inside_flag.resize(point.size(), 0);
	for (int i = 0; i < point.size(); i++) {
		inside_flag[i] = yz::geometry::isPointInsideClosedMesh(point[i], mesh.vertex, mesh.face, aabb);
	}


	for (int i = 0; i < point.size(); i++) {
		point[i] /= 0.000001;
		mesh.vertex[i] /= 0.000001;
	}

	win3d.use_arcball_flag = 1;
	win3d.SetDraw(draw);
	win3d.CreateGLUTWindow();

	manager.EnterMainLoop();
}

