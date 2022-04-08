/***********************************************************/
/**	\file
	\brief		Example of UV unwrap
	\author		Yizhong Zhang
	\date		4/17/2018
*/
/***********************************************************/

#include <iostream>
#include <queue>
#include <GL/glew.h>
#include <GL/glut.h>
#include "yzLib/yz_lib.h"
#include "data_path.h"

/**
	randomly create seam on mesh 
*/
class RandomSeamCreator {
public:
	//	mesh and topology
	yz::geometry::TriMesh<double>	mesh;
	std::set<yz::int2>				e_map;
	std::multimap<yz::int2, int>	ef_multimap;
	std::vector<int>				ff;
	std::vector<int>				ff_start;

	//	seam
	std::set<yz::int2>				seam_edge;				///< all seam edges
	std::multimap<int, int>			seam_vv;				///< vertex-vertex on boundary
	std::vector<int>				seam_joint_vertex;		///< vertex connecting different segments
	std::vector<int>				seam_end_vertex;		///< vertex that is the end of a segment, exclude joint vertex
	std::vector<int>				seam_normal_vertex;		///< boundary vertex exclude joint vertex and end vertex

	//	patch
	std::vector<int>				face_label;
	std::vector<yz::int3>			face_tex_coor;
	std::vector<yz::Vec2d>			tex_coor;

public:
	void ReadMesh(std::vector<yz::Vec3d>& v, std::vector<yz::int3>& f) {
		mesh.vertex = v;
		mesh.face = f;
		mesh.CalculateNormals();

		//	create topology
		e_map.clear();
		ef_multimap.clear();
		ff.clear();
		ff.reserve(mesh.face.size() * 3);
		ff_start.clear();
		ff_start.reserve(mesh.face.size() + 1);
		ff_start.push_back(0);

		if (mesh.face.empty())
			return;

		//	first parse, create edge-face multimap		
		for (unsigned int f = 0; f != mesh.face.size(); f++) {	//	for each face
			for (int i = 0; i < 3; i++) {						//	for each edge of the face
				yz::int2 edge(mesh.face[f][i], mesh.face[f][(i + 1) % 3]);
				if (edge.x > edge.y)
					yz::mySwap(edge.x, edge.y);

				e_map.insert(edge);

				//	insert into hash table
				ef_multimap.insert(std::pair<yz::int2, int>(edge, f));
			}
		}

		//	second parse, extract face-face
		for (unsigned int f = 0; f != mesh.face.size(); f++) {	//	for each face
			for (int i = 0; i < 3; i++) {						//	for each edge of the face
				yz::int2 edge(mesh.face[f][i], mesh.face[f][(i + 1) % 3]);
				if (edge.x > edge.y)
					yz::mySwap(edge.x, edge.y);

				//	for each neighbor of this edge
				auto range = ef_multimap.equal_range(edge);
				for (auto it = range.first; it != range.second; ++it) {
					if (it->second == f)	//	if it is not f itself, it is neighbor
						continue;
					ff.push_back(it->second);	//	a neighbor face
				}
			}

			//	all neighbor faces detected, count the number
			ff_start.push_back(ff.size());
		}
	}

	void ReadMesh(const char* file_name) {
		if (!mesh.ReadMeshFromFile(file_name)) {
			std::cout << "read mesh failed" << std::endl;
			return;
		}

		ReadMesh(mesh.vertex, mesh.face);
	}

	void RandLabel(int label_num) {
		face_label.clear();
		face_label.resize(mesh.face.size(), -1);

		//	random pick several faces, set as initial seed
		std::queue<int> face_queue;
		for (int i = 0; i < label_num; i++) {
			int fid = rand() % mesh.face.size();
			while (face_label[fid] != -1)
				fid = rand() % mesh.face.size();
			face_label[fid] = i;
			face_queue.push(fid);
		}

		//	broad first search to floodfill
		while (!face_queue.empty()) {
			int fid = face_queue.front();
			face_queue.pop();
			for (int i = ff_start[fid]; i < ff_start[fid + 1]; i++) {
				int nfid = ff[i];
				if (face_label[nfid] != -1)
					continue;
				face_label[nfid] = face_label[fid];
				face_queue.push(nfid);
			}
		}

		MarkSeamFromLabel();
		CreatePatch();
	}

	void MarkSeamFromLabel() {
		seam_edge.clear();
		seam_vv.clear();
		seam_joint_vertex.clear();
		seam_end_vertex.clear();
		seam_normal_vertex.clear();

		//	create seam edge
		for (auto iter = e_map.begin(); iter != e_map.end(); iter++) {	//	for each edge
			auto range = ef_multimap.equal_range(*iter);
			int fid = -1;
			for (auto it = range.first; it != range.second; ++it) {		//	for each neighbor face of the edge
				if (fid < 0) {				//	the first face
					fid = it->second;
				}
				else {						//	other neighbor faces
					int nfid = it->second;
					if (face_label[nfid] != face_label[fid])	//	if face label is different, this edge is seam edge
						seam_edge.insert(*iter);
				}
			}
		}

		//	create vertex-vertex
		for (auto iter = seam_edge.begin(); iter != seam_edge.end(); iter++) {
			seam_vv.insert(std::pair<int, int>(iter->x, iter->y));
			seam_vv.insert(std::pair<int, int>(iter->y, iter->x));
		}

		//	parse vv, count number of neighbors of each vertex
		std::map<int, int> vv_count;
		for (auto iter = seam_vv.begin(); iter != seam_vv.end(); iter++) {
			auto vv_iter = vv_count.find(iter->first);
			if (vv_iter == vv_count.end())
				vv_count.insert(std::pair<int, int>(iter->first, 1));
			else
				vv_iter->second++;
		}

		//	create seam vertices of different type
		for (auto iter = vv_count.begin(); iter != vv_count.end(); iter++) {
			if (iter->second == 1)		// only 1 neighbor seam vertex
				seam_end_vertex.push_back(iter->first);
			else if (iter->second > 2)	// 3 and more neighbor seam vertices
				seam_joint_vertex.push_back(iter->first);
			else 						// 2 neighbor seam vertices
				seam_normal_vertex.push_back(iter->first);
		}

	}

	void CreatePatch() {
		//	detect all labels
		std::set<int> label;
		for (int i = 0; i < face_label.size(); i++)
			label.insert(face_label[i]);

		face_tex_coor.resize(mesh.face.size());

		//	for each existing label
		for (auto iter = label.begin(); iter != label.end(); iter++) {
			int label = *iter;

			//	mark all vertices shared by this label
			std::vector<int> vertex_reindex;
			vertex_reindex.resize(mesh.vertex.size(), -1);
			for (int fid = 0; fid < mesh.face.size(); fid++) {
				if (face_label[fid] == label) {
					for (int i = 0; i < 3; i++) {
						int vid = mesh.face[fid][i];
						vertex_reindex[vid] = 0;
					}
				}
			}

			//	reindex these vertices
			int patch_vertex_count = 0;
			for (int vid = 0; vid < vertex_reindex.size(); vid++) {
				if (vertex_reindex[vid] == 0)
					vertex_reindex[vid] = tex_coor.size() + patch_vertex_count++;
			}

			//	record the face tex_coord
			for (int fid = 0; fid < mesh.face.size(); fid++) {
				if (face_label[fid] == label) {
					for (int i = 0; i < 3; i++) {
						int vid = mesh.face[fid][i];
						face_tex_coor[fid][i] = vertex_reindex[vid];
					}
				}
			}

			tex_coor.resize(tex_coor.size() + patch_vertex_count);
		}


	}

	void Draw3D() {
		std::vector<yz::uchar3> face_color;
		face_color.resize(face_label.size());
		for (int i = 0; i < face_color.size(); i++) {
			yz::opengl::getSequentialDisplayColor(&face_color[i].x, face_label[i]);
		}

		//	draw mesh
		glDisable(GL_LIGHTING);
		yz::opengl::drawFlatColorTriMesh(mesh.vertex, mesh.face, face_color);
		glEnable(GL_LIGHTING);

		//	draw mesh edge
		glColor3f(1, 1, 1);
		yz::opengl::drawMeshEdgeFromFace(mesh.vertex, mesh.face, mesh.vertex_normal);

		//	draw seam
		glColor3f(1, 0, 1);
		for (auto iter = seam_edge.begin(); iter != seam_edge.end(); iter++) {
			yz::opengl::drawCylinder(mesh.vertex[iter->x], mesh.vertex[iter->y], 0.01);
		}

		//	draw joint vertex
		glColor3f(1, 0, 0);
		for (int i = 0; i < seam_joint_vertex.size(); i++) {
			int vid = seam_joint_vertex[i];
			yz::opengl::drawPointAsCube(mesh.vertex[vid], 0.05);
		}

		//	draw non-manifold boundary end vertex
		glColor3f(0, 1, 0);
		for (int i = 0; i < seam_end_vertex.size(); i++) {
			int vid = seam_end_vertex[i];
			yz::opengl::drawPointAsCube(mesh.vertex[vid], 0.05);
		}

	}
};

/**
	UV generator, with random seam
*/
template <class T>
class UVGenerator : public yz::geometry::UVunwrapper<T> {
public:
	void DebugCreateSeam(int patches) {
		RandomSeamCreator seam_creator;
		seam_creator.ReadMesh(this->vertex, this->face);
		seam_creator.RandLabel(patches);

		std::vector<yz::int2> tmp_seam;
		for (auto iter = seam_creator.seam_edge.begin(); iter != seam_creator.seam_edge.end(); iter++)
			tmp_seam.push_back(*iter);

		this->SetSeam(tmp_seam);
	}
};


yz::opengl::DemoWindowManager	manager;
yz::opengl::GLUTWindow3D<0>		win3d(0, 0, 600, 500);
yz::opengl::GLUTWindow2D<1>		win2d(600, 0, 600, 500);

UVGenerator<double> unwrapper;
int			picked_patch_id = -1;	//	picked patch index
int			rotate_flag = 0;		//	0: translate, 1: rotate
yz::Vec2d	picked_point;			//	the picked point
yz::int2	picked_coord;			//	the picked coordinate

void draw_3d() {
	unwrapper.Draw3D();
}

void showInfo_2d() {
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "ctrl + mouse left buttom to move patch");
	yz::opengl::printInfo(0, 30, "ctrl + mouse right buttom to rotate patch");
	yz::opengl::printInfo(0, 60, "d to dump uv.obj");
}

void draw_2d() {
	yz::opengl::drawXYZAxis();

	unwrapper.Draw2D();
}

void picking_draw_2d() {
	unwrapper.PickingDraw2D();
}

void process_picking_2d(int idx) {
	if (idx != 0xffffffff && idx < unwrapper.PatchNumber()) {
		picked_patch_id = idx;
		win2d.CalculateCoordinate(picked_point[0], picked_point[1], win2d.old_x, win2d.old_y);
		picked_coord = yz::int2(win2d.old_x, win2d.old_y);
	}
}

void mouse_2d(int button, int state, int x, int y) {
	win2d.mouse_state[button] = state;
	win2d.old_x = x;
	win2d.old_y = y;
	win2d.modifier = glutGetModifiers();

	if (picked_patch_id >= 0) {
		if (rotate_flag == 0 && win2d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_UP)
			picked_patch_id = -1;
		else if (rotate_flag == 1 && win2d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_UP) {
			picked_patch_id = -1;
		}
	}
	else {
		if ((win2d.modifier&GLUT_ACTIVE_CTRL) && win2d.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN) {
			(*win2d.picking_displayFunc)();
			int index = win2d.PickingIndex(x, y);
			process_picking_2d(index);
			rotate_flag = 0;
		}
		else if ((win2d.modifier&GLUT_ACTIVE_CTRL) && win2d.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN) {
			(*win2d.picking_displayFunc)();
			int index = win2d.PickingIndex(x, y);
			process_picking_2d(index);
			rotate_flag = 1;
		}
	}

	glutPostRedisplay();
}

void motion_2d(int x, int y) {
	if (picked_patch_id >= 0) {
		yz::Vec2d curr_point;
		win2d.CalculateCoordinate(curr_point[0], curr_point[1], x, y);

		if (rotate_flag == 0)
			unwrapper.TranslatePatchUV(picked_patch_id, curr_point - picked_point);
		else
			unwrapper.RotatePatchUV(picked_patch_id, (y - picked_coord.y) * 0.2);

		picked_point = curr_point;
		picked_coord = yz::int2(x, y);

		win2d.old_x = x;
		win2d.old_y = y;
		glutPostRedisplay();
	}
	else
		win2d.DefaultMotionFunc(x, y);
}

void keyboard_2d(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case 'd':
		unwrapper.WriteMeshToFile("uv.obj");
		break;
	}
}

int main() {
	char obj_filename[1024];
	sprintf(obj_filename, "%s/unit_icosphere.obj", data_path);

	unwrapper.ReadMeshFromFile(obj_filename);
	unwrapper.DebugCreateSeam(3);
	unwrapper.Unwrap();

	win3d.use_arcball_flag = 1;
	win3d.SetDraw(draw_3d);
	win3d.CreateGLUTWindow();

	win2d.keyboardFunc = keyboard_2d;
	win2d.mouseFunc = mouse_2d;
	win2d.motionFunc = motion_2d;
	win2d.CreateGLUTWindow();
	win2d.SetDrawAppend(showInfo_2d);
	win2d.SetDraw(draw_2d);
	win2d.SetPickingDraw(picking_draw_2d);
	win2d.SetProcessPicking(process_picking_2d);

	manager.AddIdleFunc(win3d.idleFunc);
	manager.AddIdleFunc(win2d.idleFunc);
	manager.EnterMainLoop();

	return 0;
}