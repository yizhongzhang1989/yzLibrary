#include <iostream>
#include <string>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <yzLib/yz_lib.h>
using namespace std;

typedef double REAL;

yz::opengl::DemoWindowManager	manager;

//	3d values
yz::opengl::GLUTWindow3D<0>		win3d;

yz::geometry::SmoothShadingTriMesh<REAL>	mesh;

int			picking_flag_3d = 0;
int			picking_normal_flag_3d = 0;
int			picked_vertex_id_3d = 0;
int			picked_normal_id_3d = 0;
GLdouble	picking_depth_3d;
yz::Vec3d	picking_point_3d;
yz::Vec3d	picking_offset_3d;
int			draw_face_flag_3d = 1;
int			draw_edge_flag_3d = 1;
REAL		cube_size_3d = 0.01;
REAL		normal_scale = 0.1;

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

	//	draw vertex normal
	glColor3f(1, 0, 0);
	yz::opengl::drawVertexNormal(mesh.vertex, mesh.vertex_normal, normal_scale);
	glColor3f(1, 1, 0);
	for(int i=0; i<mesh.vertex.size(); i++){
		yz::Vec3f normal_point = mesh.vertex[i] + mesh.vertex_normal[i] * normal_scale;
		yz::opengl::drawPointAsSphere(normal_point, cube_size_3d);
	}	

	//	draw cubic bezier control points
	for(int i=0; i<mesh.face.size(); i++){
		glColor3f(0,0,1);
		int v0id = mesh.face[i].x;
		int v1id = mesh.face[i].y;
		int v2id = mesh.face[i].z;
		std::vector<yz::Vec3<REAL>> control_points;
		yz::geometry::remeshing::getCubicBezierTriangleControlPoints(control_points,
			mesh.vertex[v0id], mesh.vertex[v1id], mesh.vertex[v2id], 
			mesh.vertex_normal[v0id], mesh.vertex_normal[v1id], mesh.vertex_normal[v2id]);

		for(int j=0; j<control_points.size(); j++){
			yz::opengl::drawPointAsSphere(control_points[j], cube_size_3d);
		}

		//	draw a lot of points on bezier surface
		glColor3f(0,1,1);
		for(float u=0; u<=1; u+=0.01){
			for( float v=0; v<=1-u; v+=0.01 ){
				yz::Vec3f p = yz::interpCubicBezierTriangle(&control_points[0], u, v);
				yz::opengl::drawPointAsCube(p, cube_size_3d);
			}
		}
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
		yz::opengl::drawPointAsSphere(mesh.vertex[i], cube_size_3d*2);
	}

	for(int i=0; i<mesh.vertex.size(); i++){
		yz::opengl::setPickingIndex(mesh.vertex.size() + i);
		yz::Vec3f normal_point = mesh.vertex[i] + mesh.vertex_normal[i] * normal_scale;
		yz::opengl::drawPointAsSphere(normal_point, cube_size_3d);
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
	else if( idx != 0xffffffff && idx>=mesh.vertex.size() && idx<mesh.vertex.size()*2 ){	//	normal selected
		picking_normal_flag_3d = 1;
		picked_normal_id_3d = idx - mesh.vertex.size();
		picking_depth_3d = depth;
		yz::opengl::getWorldCoordinate(picking_point_3d[0], picking_point_3d[1], picking_point_3d[2], 
			win3d.old_x, win3d.old_y, depth);
		picking_offset_3d = picking_point_3d - (mesh.vertex[picked_normal_id_3d] + mesh.vertex_normal[picked_normal_id_3d]*normal_scale);
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
	else if(picking_normal_flag_3d && win3d.mouse_state[GLUT_LEFT_BUTTON]==GLUT_UP){
		picking_normal_flag_3d = 0;
	}

	glutPostRedisplay();
}

void motion_3d(int x, int y){
	if( picking_flag_3d ){
		yz::opengl::getWorldCoordinate(picking_point_3d[0], picking_point_3d[1], picking_point_3d[2], 
			x, y, picking_depth_3d);

		mesh.vertex[picked_vertex_id_3d] = picking_point_3d - picking_offset_3d;

		win3d.old_x = x;
		win3d.old_y = y;
		glutPostRedisplay();
	}
	else if( picking_normal_flag_3d ){
		yz::opengl::getWorldCoordinate(picking_point_3d[0], picking_point_3d[1], picking_point_3d[2], 
			x, y, picking_depth_3d);

		yz::Vec3f normal_point = picking_point_3d - picking_offset_3d;
		yz::Vec3f nor = normal_point - mesh.vertex[picked_normal_id_3d];
		nor.SetNormalize();
		mesh.vertex_normal[picked_normal_id_3d] = nor;

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

//	main function
int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	//	setup 3d mesh
	mesh.vertex.resize(3);
	mesh.vertex[0] = yz::Vec3f(1, 0, 0);
	mesh.vertex[1] = yz::Vec3f(0, 1, 0);
	mesh.vertex[2] = yz::Vec3f(0, 0, 1);
	mesh.face.resize(1, yz::int3(0,1,2));
	mesh.CalculateVertexNormal();

	//	setup 3d window
	win3d.SetDraw(draw_3d);
	win3d.SetDrawAppend(showInfo_3d);
	win3d.SetPickingDraw(picking_draw_3d);
	win3d.mouseFunc		= mouse_3d;
	win3d.motionFunc	= motion_3d;
	win3d.keyboardFunc	= keyboard_3d;
	win3d.specialFunc	= special_keys_3d;
	win3d.CreateGLUTWindow();

	//	setup window manager
	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}