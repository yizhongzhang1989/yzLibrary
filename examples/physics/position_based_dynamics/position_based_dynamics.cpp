/***********************************************************/
/**	\file
	\brief		Example of Position Based Dynamics
	\details	This file illustrates the use of PBD in yzLib
	\author		Yizhong Zhang
	\date		12/18/2015
	*/
/***********************************************************/
#include <iostream>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <yzLib/yz_lib.h>
using namespace std;
#include "data_path.h"

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::physics::PositionBasedDynamicTriMesh<double>	sheet;

yz::geometry::TriMesh<double>						rigid;
yz::geometry::AABBTree3D<double>					rigid_aabb;

int			picked_type = -1;			//	0:sheet, 1:rigid
int			picked_vertex_id = -1;
yz::Vec3d	picked_point;


void print(){
	static yz::utils::FPSCalculator	fps_calc;
	float fps = fps_calc.GetFPS();

	glColor3f(1, 1, 0);
	yz::opengl::printInfo(0, 0, "Fps: %f", fps);
}

void draw() {
	yz::opengl::drawXYZAxis();

	//	draw sheet face
	glColor3f(1, 1, 1);
	sheet.Display();

	//	draw sheet edge
	glColor3f(0, 0, 1);
	yz::opengl::drawMeshEdge(sheet.vertex, sheet.edge, sheet.vertex_normal);

	//	draw sheet vertices
	glColor3f(1, 1, 0);
	for (int i = 0; i < sheet.vertex.size(); i++)
		yz::opengl::drawPointAsCube(sheet.vertex[i], 0.02);

	//	draw sheet constraint
	glColor3f(1, 0, 0);
	for (int i = 0; i < sheet.vertex.size(); i++) {
		if (sheet.attach_flag[i])
			yz::opengl::drawPointAsCube(sheet.vertex[i], 0.04);
	}

	//	draw rigid
	glColor3f(0, 1, 1);
	rigid.DisplayFlat();
}

void picking_draw(){
	//	draw sheet face
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < sheet.face.size(); i++){
		yz::opengl::setPickingIndex(i);
		yz::Vec3d v0 = sheet.vertex[sheet.face[i].x];
		yz::Vec3d v1 = sheet.vertex[sheet.face[i].y];
		yz::Vec3d v2 = sheet.vertex[sheet.face[i].z];
		glVertex3d(v0.x, v0.y, v0.z);
		glVertex3d(v1.x, v1.y, v1.z);
		glVertex3d(v2.x, v2.y, v2.z);
	}
	glEnd();

	//	draw constraint
	for (int i = 0; i < sheet.vertex.size(); i++) {
		yz::opengl::setPickingIndex(i + sheet.face.size());
		if (sheet.attach_flag[i])
			yz::opengl::drawPointAsCube(sheet.vertex[i], 0.04);
	}

	//	draw rigid
	yz::opengl::setPickingIndex(sheet.face.size() + sheet.vertex.size());
	rigid.DisplayFlat();
}

void process_picking(int idx){
	double depth = win3d.PickingDepth(win3d.old_x, win3d.old_y);
	yz::opengl::getWorldCoordinate(picked_point[0], picked_point[1], picked_point[2], win3d.old_x, win3d.old_y, depth);

	if (idx >= 0 && idx < sheet.face.size()){	//	picked a face, choose one of the three neighboring vertices that is closest to the picked point
		double min_squ_dist = 1e10;
		int min_dist_vid = -1;
		for (int i = 0; i < 3; i++){
			int vid = sheet.face[idx][i];
			yz::Vec3d v = sheet.vertex[vid];
			if ((v - picked_point).SquareLength() < min_squ_dist){
				min_dist_vid = vid;
				min_squ_dist = (v - picked_point).SquareLength();
			}
		}
		picked_type = 0;
		picked_vertex_id = min_dist_vid;
	}
	else if (idx < sheet.face.size() + sheet.vertex.size()){		//	picked a constraint
		picked_type = 0;
		picked_vertex_id = idx - sheet.face.size();
	}
	else if (idx < sheet.face.size() + sheet.vertex.size() + 1){	//	picked rigid
		picked_type = 1;
	}
	else{	//	picked nothing
		picked_type = -1;
		picked_vertex_id = -1;
	}
}

void keyboard_func(unsigned char key, int x, int y){
	switch (key){
	case 27:
		exit(0);
	}
}

void mouse_func(int button, int state, int x, int y){
	win3d.mouse_state[button] = state;
	win3d.old_x = x;
	win3d.old_y = y;
	win3d.modifier = glutGetModifiers();

	//	if picking_displayFunc is set and Ctrl is pressed, then perform color picking
	if (win3d.picking_displayFunc && (win3d.modifier&GLUT_ACTIVE_CTRL) && state == GLUT_DOWN){
		(*win3d.picking_displayFunc)();
		int index = win3d.PickingIndex(x, y);
		(*win3d.process_picking)(index);

		if (picked_type == 0 && picked_vertex_id >= 0){		//	picked sheet
			if (win3d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN)
				sheet.AddAttachConstraint(picked_vertex_id);
			else if (win3d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN)
				sheet.RemoveAttachConstraint(picked_vertex_id);
		}
		else if (picked_type == 1){		//	picked rigid

		}
	}
	else if (win3d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_UP){
		picked_type = -1;
		picked_vertex_id = -1;
	}

	glutPostRedisplay();
}

void motion_func(int x, int y){
	if (win3d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN){
		if (picked_type == 0 && picked_vertex_id >= 0){
			//	move selected vertex
			double u, v, d;
			yz::opengl::getWindowCoordinate(u, v, d, picked_point.x, picked_point.y, picked_point.z);
			yz::Vec3d xyz;
			yz::opengl::getWorldCoordinate(xyz.x, xyz.y, xyz.z, x, y, d);

			yz::Vec3d attach_point;
			sheet.GetAttachedPosition(attach_point, picked_vertex_id);
			attach_point += xyz - picked_point;
			picked_point = xyz;
			sheet.AddAttachConstraint(picked_vertex_id, attach_point);
		}
		else if (picked_type == 1){
			//	move rigid
			double u, v, d;
			yz::opengl::getWindowCoordinate(u, v, d, picked_point.x, picked_point.y, picked_point.z);
			yz::Vec3d xyz;
			yz::opengl::getWorldCoordinate(xyz.x, xyz.y, xyz.z, x, y, d);

			yz::geometry::translateVertices(rigid.vertex, xyz - picked_point);
			picked_point = xyz;

			rigid_aabb.UpdateTriangleAABBTree(rigid.vertex, rigid.face, 0.001);
		}
		else{
			//	rotate the screen
			if (win3d.modifier & GLUT_ACTIVE_SHIFT){
				win3d.eye_x -= float(x - win3d.old_x) / 20;
				win3d.eye_y += float(y - win3d.old_y) / 20;
			}
			else{
				win3d.spin_x += float(x - win3d.old_x);
				win3d.spin_y += float(y - win3d.old_y);
			}
		}
	}
	else if (win3d.mouse_state[GLUT_MIDDLE_BUTTON] == GLUT_DOWN){
		win3d.eye_x -= float(x - win3d.old_x) / 20;
		win3d.eye_y += float(y - win3d.old_y) / 20;
	}
	else if (win3d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN){
		win3d.eye_z += float(y - win3d.old_y) / 10;
	}

	win3d.old_x = x;
	win3d.old_y = y;
	glutPostRedisplay();
}

void idle(){
//!	[Position Based Dynamics Simulation]
	sheet.SimulateStepPart1();
	sheet.AddCollisionConstraintsWithClosedTriMesh(rigid.vertex, rigid.face, rigid_aabb);
	sheet.SimulateStepPart2();
	sheet.SimulateStepPart3();
//!	[Position Based Dynamics Simulation]
}

//	main function
int main(int argc, char* argv[]){
	char sheet_filename[1024], rigid_filename[1024];
	sprintf(sheet_filename, "%s/sheet.obj", data_path);
	sprintf(rigid_filename, "%s/sheet_rigid.obj", data_path);

	//	read models
	sheet.k_bending = 0.01;
	sheet.ReadMeshFromFile(sheet_filename);
	sheet.AddAttachConstraint(0);
	sheet.AddAttachConstraint(1);

	rigid.ReadMeshFromFile(rigid_filename);
	rigid_aabb.BuildTriangleAABBTree(rigid.vertex, rigid.face, 0.001);

	//	setup OpenGL context
	win3d.SetDraw(draw);
	win3d.SetDrawAppend(print);
	win3d.SetPickingDraw(picking_draw);
	win3d.SetProcessPicking(process_picking);
	win3d.keyboardFunc = keyboard_func;
	win3d.mouseFunc = mouse_func;
	win3d.motionFunc = motion_func;
	win3d.CreateGLUTWindow();

	//	setup window manager
	manager.AddIdleFunc(win3d.idleFunc);
	manager.AddIdleFunc(idle);
	manager.EnterMainLoop();
}