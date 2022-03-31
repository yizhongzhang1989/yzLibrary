/***********************************************************/
/**	\file
	\brief		Example of glut window and fbo usage
	\details	In this demo, three glut windows of different 
				kinds are created. A FBO and a cuda enabled FBO
				are created in context of win1. It uses display
				function of win0 as its own texture to display,
				then create fbo of itself and transfer to win2.

				When there are multi-windows, the context of
				OpenGL may be very complex. You must be very 
				careful with context.
	\author		Yizhong Zhang
	\date		6/27/2012
*/
/***********************************************************/

//	This example shows how to use GLUTWindow and FBO

#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
#include <FreeImage.h>
#include <yzLib/yz_lib.h>
using namespace std;

//	3 windows and a manager, ID of each window must be unique
yz::opengl::DemoWindowManager		manager;
yz::opengl::GLUTWindow3D<0>			win0(0,	0, 800, 600);
yz::opengl::GLUTWindow2D<1>			win1(820, 0, 400, 300);
yz::opengl::GLUTTextureWindow2D<2>	win2(820, 350, 400, 300);

//	picking variables
int current_select = -1;
yz::Vec3f picking_point;
yz::Vec3f cube_center(1, 0, 0);
yz::Vec3f teapot_center(-1, 0, 0);

//	other variables
float				fps;
yz::opengl::FBO		fbo;
yz::cuda::CudaFBO	cuda_fbo;
float*				tex_ptr;

//	over-ridden functions
void calculateFps(){
	static yz::utils::FPSCalculator fps_calc;
	fps = fps_calc.GetFPS();
}

void renderToFBO1(){
//!	[Render to FBO]
	//	If we call display function of another context, we have to 
	//	save all attributes since we don't know what the function did
	glutSetWindow(win1.win_id);			//	context must be set correctly
	fbo.BeginRender();

	yz::opengl::pushAllAttributesAndMatrices();	//	push all attributes and matrices
	win0.DefaultReshapeFunc(win0.win_width, win0.win_height);	//	set projection matrix
	win0.auto_swap_buffers = 0;					//	disable glutSwapBuffers() in win0.displayFunc()
	win0.DefaultDisplayFunc();					//	call display function of another context
	win0.auto_swap_buffers = 1;					//	enable glutSwapBuffers() in win0.displayFunc()
	yz::opengl::popAllAttributesAndMatrices();	//	pop all attributes and matrices

	fbo.EndRender();
	//	Still some variables are possibly changed, such as the matrix.
	//	Special attention must be paid
//!	[Render to FBO]
}

void renderToFBO2(){
	glutSetWindow(win1.win_id);
	cuda_fbo.BeginRender();

	win1.DefaultDisplayFunc();

	cuda_fbo.EndRender();
	cuda_fbo.GetFBOData();	//	get the texture to GPU memory

	cpyd2h(tex_ptr, cuda_fbo.tex_ptr_d, sizeof(float)*cuda_fbo.tex_width*cuda_fbo.tex_height*4);
	glutSetWindow(win2.win_id);
	win2.texture.LoadPtrToTexture();
}

void win0_showFps(){
	glColor3f(1, 0, 1);
	yz::opengl::printInfo(0, 0, "fps: %f\nCtrl + Left click can pick object", fps);
}

void win0_draw(){
	//	draw picking point
	if( current_select != -1 ){
		glColor3f(1, 0, 1);
		yz::opengl::drawPointAsSphere(picking_point, 0.02);
	}

	glPushMatrix();
	glTranslatef(cube_center.x, cube_center.y, cube_center.z);
	glColor3f(1, current_select==0, 0);
	glutSolidCube(0.5);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(teapot_center.x, teapot_center.y, teapot_center.z);
	glColor3f(1, current_select==1, 0);
	glutSolidTeapot(0.5);
	glPopMatrix();
}

void win0_picking_draw(){
	glPushMatrix();
	glTranslatef(cube_center.x, cube_center.y, cube_center.z);
	yz::opengl::setPickingIndex(0);
	glutSolidCube(0.5);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(teapot_center.x, teapot_center.y, teapot_center.z);
	yz::opengl::setPickingIndex(1);
	glutSolidTeapot(0.5);
	glPopMatrix();
}

void win0_process_picking(int idx){
	current_select = idx;

	//	calculate world coordinate of the picking point
	float depth = win0.PickingDepth(win0.old_x, win0.old_y);
	yz::opengl::getWorldCoordinate(picking_point[0], picking_point[1], picking_point[2], win0.old_x, win0.old_y, depth);
}

void win0_mouse_func(int button, int state, int x, int y){
	win0.mouse_state[button] = state;
	win0.old_x = x;
	win0.old_y = y;
	win0.modifier = glutGetModifiers();

	//	if picking_displayFunc is set and Ctrl is pressed, then perform color picking
	if( win0.picking_displayFunc && (win0.modifier&GLUT_ACTIVE_CTRL) && win0.mouse_state[GLUT_LEFT_BUTTON]==GLUT_DOWN ){
		(*win0.picking_displayFunc)();
		int index = win0.PickingIndex(x, y);
		(*win0.process_picking)(index);
	}
	else if( win0.mouse_state[GLUT_LEFT_BUTTON]==GLUT_UP ){
		current_select = -1;
	}

	glutPostRedisplay();
}

void win0_motion_func(int x, int y){
	if( win0.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN ){
		if( current_select >= 0 ){
			//	move the selected objects
			double u, v, d;
			yz::opengl::getWindowCoordinate( u, v, d, picking_point.x, picking_point.y, picking_point.z );
			yz::Vec3f xyz;
			yz::opengl::getWorldCoordinate(xyz.x, xyz.y, xyz.z, x, y, d);

			if( current_select == 0 )
				cube_center += xyz - picking_point;
			else if( current_select == 1 )
				teapot_center += xyz - picking_point;
			picking_point = xyz;
		}
		else{
			//	rotate the screen
			if( win0.modifier & GLUT_ACTIVE_SHIFT ){
				win0.eye_x -= float( x - win0.old_x ) / 20;
				win0.eye_y += float( y - win0.old_y ) / 20;
			}
			else{
				win0.spin_x += float( x - win0.old_x );
				win0.spin_y += float( y - win0.old_y );
			}
		}
	}
	else if( win0.mouse_state[GLUT_MIDDLE_BUTTON] == GLUT_DOWN ){
		win0.eye_x -= float( x - win0.old_x ) / 20;
		win0.eye_y += float( y - win0.old_y ) / 20;
	}
	else if( win0.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN ){
		win0.eye_z += float( y - win0.old_y ) / 10;
	}

	win0.old_x = x;
	win0.old_y = y;
	glutPostRedisplay();
}

void win1_showInfo(){
	glColor3f(0, 0, 1);
	yz::opengl::printInfo(0, 30, "texture got from FBO");
}

void win2_showInfo(){
	glColor3f(0, 1, 0);
	yz::opengl::printInfo(0, 60, "press 'w' to write the image");
}

void win1_draw(){
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, fbo.tex_id);
	glColor3f(1, 1, 1);
	yz::opengl::drawWholeTexture(0, 0, 1, 1);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void win2_keyboard(unsigned char key, int x, int y){
	switch( key ){
		case 'r':
			win2.left	= 0;
			win2.right	= 1;
			win2.bottom	= 0;
			win2.top	= 1;
			win2.DefaultReshapeFunc(win2.win_width, win2.win_height);	//	reset projection matrix
			break;
		case 'w':
			cout << "write texture to image";
			if( yz::utils::writeImageToFile("teapot.png", tex_ptr, cuda_fbo.tex_width, cuda_fbo.tex_height, 32*4) )
				cout << " succeed" << endl;
			else
				cout << " failed" << endl;
			break;
		case 27:
			exit(0);
	}
}

int main(){
	yz::cuda::InitCudaGL();

//!	[Context awareness]
	//	OpenGL is context aware, resources are valid only inside the context. 
	//	So if you need to create FBO or texture, do it right after create window
	//	create window 0
	win0.SetDrawAppend(win0_showFps);				//	set draw append callback
	win0.SetDraw(win0_draw);						//	set draw callback
	win0.SetBackgroundColor(0, 0.4, 0.5, 0);		//	set background color
	win0.SetPickingDraw(win0_picking_draw);			//	set picking draw callback
	win0.SetProcessPicking(win0_process_picking);	//	set process picking callback
	win0.mouseFunc	= win0_mouse_func;				//	mouse function
	win0.motionFunc	= win0_motion_func;				//	motion function
	win0.CreateGLUTWindow();						//	create window 0

	//	create window 1, then setup fbo
	win1.SetDraw(win1_draw);							//	set draw callback
	win1.SetDrawAppend(win1_showInfo);					//	set draw append callback
	win1.SetBackgroundColor(0, 1, 1, 1);				//	set background color
	win1.CreateGLUTWindow();							//	create window 1
	fbo.InitFBO(win0.win_width, win0.win_height);		//	init FBO
	cuda_fbo.InitFBO(win1.win_width, win1.win_height);	//	init cuda FBO

	//	create window 2, then setup texture
	win2.keyboardFunc = win2_keyboard;			//	set keyboard callback
	win2.SetDrawAppend(win2_showInfo);			//	set draw append callback
	win2.CreateGLUTWindow();					//	create window 2, texture created inside
	tex_ptr = new float[cuda_fbo.tex_width * cuda_fbo.tex_height * 4];	//	alloc texture memory
	win2.SetupTexturePtr(cuda_fbo.tex_width, cuda_fbo.tex_height, tex_ptr, GL_RGBA, GL_RGBA, GL_FLOAT);	//	setup texture
//!	[Context awareness]

//!	[Window Manager Setup]
	manager.AddIdleFunc(calculateFps);	//	add each idle function in sequence
	manager.AddIdleFunc(win0.idleFunc);
	manager.AddIdleFunc(renderToFBO1);
	manager.AddIdleFunc(win1.idleFunc);
	manager.AddIdleFunc(renderToFBO2);
	manager.AddIdleFunc(win2.idleFunc);

	manager.EnterMainLoop();			//	then enter main loop
//!	[Window Manager Setup]
}