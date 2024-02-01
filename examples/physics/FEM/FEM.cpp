#include <iostream>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <yzLib/yz_lib.h>
#include "data_path.h"


template<class T>
class Balloon :
	public yz::physics::fem::FEMTriMesh<T>,
	public yz::geometry::MeshTextureCoordinate<T>
{
public:
	int ReadMeshFromFile(const char* file_name) {
		if (!yz::physics::fem::FEMTriMesh<T>::ReadMeshFromFile(file_name))
			return 0;

		pressure = 0.1;
		face_area.resize(this->face.size());
		CalculateFaceArea();

		return 1;
	}

	virtual inline int AddExternalForce() {
		CalculateFaceArea();

		for (int fid = 0; fid < this->face.size(); fid++) {
			yz::Vec3<T> pres = this->face_normal[fid] * pressure * face_area[fid];
			pres /= 3;
			this->f[this->face[fid].x] += pres;
			this->f[this->face[fid].y] += pres;
			this->f[this->face[fid].z] += pres;
		}

		return 1;
	}

	void CalculateFaceArea() {
		for (int i = 0; i < face_area.size(); i++) {
			yz::Vec3<T> r1 = this->vertex[this->face[i].y] - this->vertex[this->face[i].x];
			yz::Vec3<T> r2 = this->vertex[this->face[i].z] - this->vertex[this->face[i].x];
			face_area[i] = yz::cross(r1, r2).Length() * 0.5;
		}
	}

public:
	T				pressure;
	std::vector<T>	face_area;
};


yz::opengl::DemoWindowManager	manager;
yz::opengl::GLUTWindow3D<0>		win3d;

Balloon<double>	mesh;

int			picked_vertex_id = -1;
yz::Vec3d	picked_point;

int			simulate_flag = 0;

void draw() {
	yz::opengl::drawXYZAxis();

	mesh.Draw();
}

void picking_draw() {
	mesh.PickingDraw();
}

void process_picking(int idx) {
	double depth = win3d.PickingDepth(win3d.old_x, win3d.old_y);
	yz::opengl::getWorldCoordinate(picked_point[0], picked_point[1], picked_point[2], win3d.old_x, win3d.old_y, depth);

	if (idx >= 0 && idx < mesh.vertex.size()) {
		picked_vertex_id = idx;
	}
	else {
		picked_vertex_id = -1;
	}
}

void mouse_func(int button, int state, int x, int y) {
	win3d.mouse_state[button] = state;
	win3d.old_x = x;
	win3d.old_y = y;
	win3d.modifier = glutGetModifiers();

	//	if picking_displayFunc is set and Ctrl is pressed, then perform color picking
	if (win3d.picking_displayFunc && (win3d.modifier&GLUT_ACTIVE_CTRL) && state == GLUT_DOWN) {
		(*win3d.picking_displayFunc)();
		int index = win3d.PickingIndex(x, y);
		(*win3d.process_picking)(index);
	}
	else if (win3d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_UP) {
		picked_vertex_id = -1;
	}

	glutPostRedisplay();
}

void motion_func(int x, int y) {
	if (win3d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN) {
		if (picked_vertex_id >= 0) {
			//	move selected vertex
			double u, v, d;
			yz::opengl::getWindowCoordinate(u, v, d, picked_point.x, picked_point.y, picked_point.z);
			yz::Vec3d xyz;
			yz::opengl::getWorldCoordinate(xyz.x, xyz.y, xyz.z, x, y, d);

			mesh.vertex[picked_vertex_id] += xyz - picked_point;
			picked_point = xyz;

			mesh.CalculateNormals();
			mesh.ComputeStrain();
			mesh.ClearForce();
			mesh.AddElasticForce();
		}
		else {
			//	rotate the screen
			if (win3d.modifier & GLUT_ACTIVE_SHIFT) {
				win3d.eye_x -= float(x - win3d.old_x) / 20;
				win3d.eye_y += float(y - win3d.old_y) / 20;
			}
			else {
				win3d.spin_x += float(x - win3d.old_x);
				win3d.spin_y += float(y - win3d.old_y);
			}
		}
	}
	else if (win3d.mouse_state[GLUT_MIDDLE_BUTTON] == GLUT_DOWN) {
		win3d.eye_x -= float(x - win3d.old_x) / 20;
		win3d.eye_y += float(y - win3d.old_y) / 20;
	}
	else if (win3d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN) {
		win3d.eye_z += float(y - win3d.old_y) / 10;
	}

	win3d.old_x = x;
	win3d.old_y = y;
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case 's':
		simulate_flag = !simulate_flag;
		break;
	case ' ':
		std::cout << mesh.SimulateStep(0.01) << std::endl;
		break;
	}

	glutPostRedisplay();
}

void special(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP:
		break;
	case GLUT_KEY_DOWN:
		break;
	}
	glutPostRedisplay();
}

void idle() {
	if (!simulate_flag)
		return;

	static yz::utils::Timer timer;
	timer.Start();

	int steps = mesh.SimulateStep(-1);
	std::cout << steps << std::endl;
}

int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	char obj_filename[1024];
	sprintf(obj_filename, "%s/cube.obj", data_path);

	mesh.ReadMeshFromFile(obj_filename);
	mesh.AddAttachConstraint(0);

	win3d.SetDraw(draw);
	win3d.SetPickingDraw(picking_draw);
	win3d.SetProcessPicking(process_picking);
	win3d.mouseFunc = mouse_func;
	win3d.motionFunc = motion_func;
	win3d.keyboardFunc = keyboard;
	win3d.specialFunc = special;
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(idle);
	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}
