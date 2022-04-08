#include <iostream>
#include <vector>
#include <GL/glew.h>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;
yz::opengl::AsciiDisplayer		ascii_disp;

void draw(){

	yz::opengl::drawXYZAxis();

	glDisable(GL_LIGHTING);
	glColor3f(0, 1, 1);

	for(int j=0; j<16; j++)
		for(int i=0; i<16; i++){
			ascii_disp.Draw(j*16+i, i, j, i+0.5, j+1);
		}

	glEnable(GL_LIGHTING);
}

int main(){
	win3d.CreateGLUTWindow();
	win3d.SetDraw(draw);

	ascii_disp.Setup8x16();


	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}