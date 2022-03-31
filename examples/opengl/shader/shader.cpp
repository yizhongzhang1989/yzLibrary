/***********************************************************/
/**	\file
	\brief		Example of shader
	\details	In this demo, three glut windows of different 
				kinds are created. A FBO and a cuda enabled FBO
				are created in context of win1. It uses display
				function of win0 as its own texture to display,
				then create fbo of itself and transfer to win2.

				When there are multi-windows, the context of
				OpenGL may be very complex. You must be very 
				careful with context.
	\author		Yizhong Zhang
	\date		5/28/2020
*/
/***********************************************************/
#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>

//	Two windows and a manager, ID of each window must be unique
yz::opengl::DemoWindowManager		manager;
yz::opengl::GLUTWindow3D<0>			win0(0,	0, 800, 600);
yz::opengl::GLUTWindow2D<1>			win1(820, 0, 400, 300);

//	picking variables
int current_select = -1;
yz::Vec3f picking_point;
yz::Vec3f cube_center(1, 0, 0);
yz::Vec3f teapot_center(-1, 0, 0);

//	other variables
float				fps;
yz::opengl::FBO		fbo;
yz::opengl::Shader	shader;
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

	shader.UseMe();

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

	shader.UseDefault();
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

void win1_draw(){
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, fbo.tex_id);
	glColor3f(1, 1, 1);
	yz::opengl::drawWholeTexture(0, 0, 1, 1);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void init_shader() {
	static const char vertex_shader[] =
		"\
			#version 330 compatibility\n \
			#extension GL_EXT_gpu_shader4 : enable\n \
			out vec4 point;\
			void main()\
			{\
				gl_Position = ftransform();\
				point = gl_ModelViewMatrix * gl_Vertex;\
			}\
		";

	static const char fragment_shader[] =
		"\
			#version 330 compatibility\n \
			#extension GL_EXT_gpu_shader4 : enable\n \
			in vec4 point;\
			void main()\
			{\
				gl_FragColor[0] = point.x;\
				gl_FragColor[1] = point.y;\
				gl_FragColor[2] = point.z;\
				gl_FragColor[3] = 1.0;\
			}\
		";

	shader.Setup(vertex_shader, fragment_shader);
}

int main(){
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

	init_shader();


	//	create window 1, then setup fbo
	win1.SetDraw(win1_draw);							//	set draw callback
	win1.SetDrawAppend(win1_showInfo);					//	set draw append callback
	win1.SetBackgroundColor(0, 1, 1, 1);				//	set background color
	win1.CreateGLUTWindow();							//	create window 1
	fbo.InitFBO(win0.win_width, win0.win_height);		//	init FBO

//!	[Window Manager Setup]
	manager.AddIdleFunc(calculateFps);	//	add each idle function in sequence
	manager.AddIdleFunc(renderToFBO1);
	manager.AddIdleFunc(win0.idleFunc);
	manager.AddIdleFunc(win1.idleFunc);

	manager.EnterMainLoop();			//	then enter main loop
//!	[Window Manager Setup]
}