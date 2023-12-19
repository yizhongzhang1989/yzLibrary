#include <iostream>
#include <vector>
#include <GL/glew.h>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>
#include "data_path.h"

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::SingleDisplayTextureTriMesh<float>	tex_mesh;

void print() {
	static yz::utils::FPSCalculator fps_calc;
	float fps = fps_calc.GetFPS();
	glColor3f(1, 1, 1);
	yz::opengl::printInfo(5, 0, "fps: %f", fps);
}

void draw(){
	yz::opengl::drawXYZAxis();

	glDisable(GL_LIGHTING);
	glColor3f(0, 1, 1);

	tex_mesh.Display();

	glEnable(GL_LIGHTING);
}

int main(){
	win3d.CreateGLUTWindow();
	win3d.SetDraw(draw);
	win3d.SetDrawAppend(print);

	std::string model_path = std::string(data_path) + "/human_sitting/human_sitting.obj";

	//	read mesh
	tex_mesh.ReadMeshFromFile(model_path.c_str());

	//	create texture
	tex_mesh.map_Ka.ReadTexImage();
	tex_mesh.map_Kd.ReadTexImage();
	tex_mesh.map_Ks.ReadTexImage();

	tex_mesh.CreateTexture();

	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}