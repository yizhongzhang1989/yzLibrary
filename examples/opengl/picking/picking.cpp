/***********************************************************/
/**	\file
	\brief		Example of picking
	\details	this example performs picking on 3D and 2D.
				3D is shown by picking vertices of a mesh.
				2D is shown by picking vertices of a graph.
				The picking use color coding framework provided
				in yzLib. \n
				This file is a template of picking programs.
				If you need a 3D or 2D picking program, just 
				copy all the code, remove the window that you don't 
				need. Variables and functions of 3D and 2D are seperated.
	\author		Yizhong Zhang
	\date		10/27/2012
*/
/***********************************************************/

//	This example shows how to use implement picking

#include <iostream>
#include <string>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <yzLib/yz_lib.h>
using namespace std;
#include "data_path.h"


#define N 10
typedef double REAL;

yz::opengl::DemoWindowManager	manager;

//	3d values
yz::opengl::GLUTWindow3D<0>		win3d;
yz::geometry::TriMesh<REAL>		mesh;
int			picking_flag_3d = 0;
int			picked_vertex_id_3d = 0;
GLdouble	picking_depth_3d;
yz::Vec3d	picking_point_3d;
yz::Vec3d	picking_offset_3d;
int			draw_face_flag_3d = 1;
int			draw_edge_flag_3d = 1;
REAL		cube_size_3d = 0.01;

//	2d values
yz::opengl::GLUTWindow2D<1>		win2d(800, 0, 800, 600);
std::vector<yz::Vec2<REAL>>		point;
int			picking_flag_2d = 0;
int			picked_vertex_id_2d = 0;
yz::Vec2d	picking_point_2d;
yz::Vec2d	picking_offset_2d;
REAL		point_radius = 0.01;

//	functions for 3d window
void showInfo_3d(){
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "ctrl + mouse left buttom to select and move vertex");
	yz::opengl::printInfo(0, 30, "up / dowm key to change size of vertex cubes");
}

void draw_3d(){
	yz::opengl::drawXYZAxis();

	if( draw_face_flag_3d ){
		glColor3f(1, 1, 1);
		mesh.Display();
	}

	if( draw_edge_flag_3d ){
		glColor3f(0, 0, 1);
		yz::opengl::drawMeshEdgeFromFace(mesh.vertex, mesh.face, mesh.vertex_normal);
	}

	//	draw_3d vertices
	glColor3f(0, 1, 0);
	for(int i=0; i<mesh.vertex.size(); i++){
		if( picking_flag_3d && picked_vertex_id_3d==i )
			glColor3f(1, 0, 0);
		yz::opengl::drawPointAsCube(mesh.vertex[i], cube_size_3d);
		if( picking_flag_3d && picked_vertex_id_3d==i )
			glColor3f(0, 1, 0);
	}
}

void picking_draw_3d(){
	if( draw_face_flag_3d ){
		yz::opengl::setPickingIndex(mesh.vertex.size());	//	set the index to be bigger than any vertex
		mesh.Display();
	}

	//	draw_3d vertices
	for(int i=0; i<mesh.vertex.size(); i++){
		yz::opengl::setPickingIndex(i);
		yz::opengl::drawPointAsCube(mesh.vertex[i], cube_size_3d);
	}
}

void process_picking_3d(int idx){
	float depth = win3d.PickingDepth(win3d.old_x, win3d.old_y);
	if( idx != 0xffffffff && idx < mesh.vertex.size() ){	//	vertex selected
		picking_flag_3d = 1;
		picked_vertex_id_3d = idx;
		picking_depth_3d = depth;
		yz::opengl::getWorldCoordinate(picking_point_3d[0], picking_point_3d[1], picking_point_3d[2], 
			win3d.old_x, win3d.old_y, depth);
		picking_offset_3d = picking_point_3d - mesh.vertex[idx];
	}
}

void mouse_3d(int button, int state, int x, int y){
	win3d.mouse_state[button] = state;
	win3d.old_x = x;
	win3d.old_y = y;
	win3d.modifier = glutGetModifiers();

	if( (win3d.modifier&GLUT_ACTIVE_CTRL) && win3d.mouse_state[GLUT_LEFT_BUTTON]==GLUT_DOWN ){
		(*win3d.picking_displayFunc)();
		int index = win3d.PickingIndex(x, y);
		process_picking_3d(index);
	}
	else if(picking_flag_3d && win3d.mouse_state[GLUT_LEFT_BUTTON]==GLUT_UP){
		picking_flag_3d = 0;
	}

	glutPostRedisplay();
}

void motion_3d(int x, int y){
	if( picking_flag_3d ){
		yz::opengl::getWorldCoordinate(picking_point_3d[0], picking_point_3d[1], picking_point_3d[2], 
			x, y, picking_depth_3d);

		mesh.vertex[picked_vertex_id_3d] = picking_point_3d - picking_offset_3d;
		mesh.CalculateNormals();

		win3d.old_x = x;
		win3d.old_y = y;
		glutPostRedisplay();
	}
	else
		win3d.DefaultMotionFunc(x, y);
}

void keyboard_3d(unsigned char key, int x, int y){
	switch(key){
		case 27:
			exit(0);
		case 'e':
			draw_edge_flag_3d = !draw_edge_flag_3d;
			break;
		case 'f':
			draw_face_flag_3d = !draw_face_flag_3d;
			break;
	}
}

void special_keys_3d(int key, int x, int y){
	switch(key){
		case GLUT_KEY_UP:
			cube_size_3d *= 2;
			break;
		case GLUT_KEY_DOWN:
			cube_size_3d *= 0.5;
			break;
	}
}

//	functions for 2d window
void showInfo_2d(){
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "ctrl + mouse left buttom to select and move vertex");
	yz::opengl::printInfo(0, 30, "up / dowm key to change size of points");
}

void draw_2d(){
	//	draw curve
	glColor3f(0, 0, 1);
	yz::opengl::drawCurve(point);

	//	draw points
	glColor3f(0, 1, 0);
	for(int i=0; i<point.size(); i++){
		if( picking_flag_2d && picked_vertex_id_2d==i )
			glColor3f(1, 0, 0);
		yz::opengl::drawPointAsBall(point[i], point_radius);
		if( picking_flag_2d && picked_vertex_id_2d==i )
			glColor3f(0, 1, 0);
	}
}

void picking_draw_2d(){
	for(int i=0; i<point.size(); i++){
		yz::opengl::setPickingIndex(i);
		yz::opengl::drawPointAsBall(point[i], point_radius);
	}
}

void process_picking_2d(int idx){
	if( idx != 0xffffffff && idx < point.size() ){	//	vertex selected
		picking_flag_2d = 1;
		picked_vertex_id_2d = idx;
		win2d.CalculateCoordinate(picking_point_2d[0], picking_point_2d[1], win2d.old_x, win2d.old_y);
		picking_offset_2d = picking_point_2d - point[idx];
	}
}

void mouse_2d(int button, int state, int x, int y){
	win2d.mouse_state[button] = state;
	win2d.old_x = x;
	win2d.old_y = y;
	win2d.modifier = glutGetModifiers();

	if( (win2d.modifier&GLUT_ACTIVE_CTRL) && win2d.mouse_state[GLUT_LEFT_BUTTON]==GLUT_DOWN ){
		(*win2d.picking_displayFunc)();
		int index = win2d.PickingIndex(x, y);
		process_picking_2d(index);
	}
	else if(picking_flag_2d && win2d.mouse_state[GLUT_LEFT_BUTTON]==GLUT_UP){
		picking_flag_2d = 0;
	}

	glutPostRedisplay();
}

void motion_2d(int x, int y){
	if( picking_flag_2d ){
		win2d.CalculateCoordinate(picking_point_2d[0], picking_point_2d[1], x, y);

		point[picked_vertex_id_2d] = picking_point_2d - picking_offset_2d;

		win2d.old_x = x;
		win2d.old_y = y;
		glutPostRedisplay();
	}
	else
		win2d.DefaultMotionFunc(x, y);
}

void keyboard_2d(unsigned char key, int x, int y){
	switch(key){
		case 27:
			exit(0);
	}
}

void special_keys_2d(int key, int x, int y){
	switch(key){
		case GLUT_KEY_UP:
			point_radius *= 2;
			break;
		case GLUT_KEY_DOWN:
			point_radius *= 0.5;
			break;
	}
}

//	main function
int main(int argc, char* argv[]){
	glutInit(&argc, argv);
	
	char obj_filename[1024];
	sprintf(obj_filename, "%s/unit_icosphere.obj", data_path);

	//	setup 3d mesh
	string mesh_name = obj_filename;
	if( argc == 2 ){	//	parse command line
		string extension = yz::utils::getFileExtensionFromString(argv[1]);
		if( extension == "obj" ){
			mesh_name = argv[1];
		}
	}
	mesh.ReadMeshFromFile( mesh_name.c_str() );

	//	setup 2d points
	point.resize(N);
	for(int i=0; i<point.size(); i++){
		point[i].x = yz::rand0to1f();
		point[i].y = yz::rand0to1f();
	}

	//	setup 3d window
	win3d.SetDraw(draw_3d);
	win3d.SetDrawAppend(showInfo_3d);
	win3d.SetPickingDraw(picking_draw_3d);
	win3d.mouseFunc		= mouse_3d;
	win3d.motionFunc	= motion_3d;
	win3d.keyboardFunc	= keyboard_3d;
	win3d.specialFunc	= special_keys_3d;
	win3d.CreateGLUTWindow();

	//	setup 2d window
	win2d.SetDraw(draw_2d);
	win2d.SetDrawAppend(showInfo_2d);
	win2d.SetPickingDraw(picking_draw_2d);
	win2d.mouseFunc		= mouse_2d;
	win2d.motionFunc	= motion_2d;
	win2d.keyboardFunc	= keyboard_2d;
	win2d.specialFunc	= special_keys_2d;
	win2d.CreateGLUTWindow();

	//	setup window manager
	manager.AddIdleFunc(win3d.idleFunc);
	manager.AddIdleFunc(win2d.idleFunc);
	manager.EnterMainLoop();
}