#include <iostream>
#include <set>
#include <unordered_map>
#include <GL/glut.h>
#include "yzLib/yz_lib.h"

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::field::MultiMaterialCubeMarcher		marcher;
yz::geometry::TriMeshd mesh1, mesh2, mesh3;
yz::geometry::AABBTree3Dd aabb1, aabb2, aabb3;

int idx = 0;

class NonManifoldMesh : public yz::geometry::TriMeshd {
public:
	std::vector<yz::int2>	face_label;
	std::vector<yz::uchar3> face_color;

	std::vector<std::vector<yz::int3>>		seperate_face;
	std::vector<std::vector<yz::uchar3>>	seperate_face_color;

	void CalculateFaceColor() {
		int count = 0;
		std::unordered_map<yz::Vec2i, int, yz::utils::BitwiseHasher<yz::Vec2i>> label;
		for (int i = 0; i < face_label.size(); i++) {
			if (label.find(face_label[i]) == label.end()) {
				label.insert(std::pair<yz::Vec2i, int>(face_label[i], count));
				count++;
			}
		}

		face_color.resize(face_label.size());
		for (int i = 0; i < face_color.size(); i++) {
			int idx = label.find(face_label[i])->second;
			yz::opengl::getSequentialDisplayColor(&face_color[i].x, idx);
		}
	}

	void CreateSeperateFace() {
		for (int i = 0; i < face.size(); i++) {
			for (int j = 0; j < 2; j++) {
				int idx = face_label[i][j];
				if (seperate_face.size() <= idx) {
					seperate_face.resize(idx + 1);
					seperate_face_color.resize(idx + 1);
				}

				seperate_face[idx].push_back(face[i]);
				seperate_face_color[idx].push_back(face_color[i]);
			}
		}
	}

	void Disp() {
		glDisable(GL_LIGHTING);
		yz::opengl::drawFlatColorTriMesh(
			vertex,
			face,
			face_color
		);
		glEnable(GL_LIGHTING);

		glColor3f(0, 0, 0);
		yz::opengl::drawMeshEdgeFromFace(
			vertex,
			face
		);
	}

	void DispSeperate(int idx) {
		if (idx < 0 || idx >= seperate_face.size()) {
			Disp();
			return;
		}

		glDisable(GL_LIGHTING);
		yz::opengl::drawFlatColorTriMesh(
			vertex,
			seperate_face[idx],
			seperate_face_color[idx]
		);
		glEnable(GL_LIGHTING);

		glColor3f(0, 0, 0);
		yz::opengl::drawMeshEdgeFromFace(
			vertex,
			seperate_face[idx],
			vertex_normal
		);

		yz::opengl::drawMeshEdgeFromFace(
			vertex,
			seperate_face[idx],
			vertex_normal,
			-0.001
		);
	}

}mesh;

void print() {
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "%d", idx);
}

void draw() {
	//mesh.Disp();
	mesh.DispSeperate(idx);

	//glColor3f(1, 0, 0);
	//yz::opengl::drawMeshEdgeFromFace(mesh1.vertex, mesh1.face, mesh1.vertex_normal);
	//glColor3f(0, 1, 0);
	//yz::opengl::drawMeshEdgeFromFace(mesh2.vertex, mesh2.face, mesh2.vertex_normal);
	//glColor3f(0, 0, 1);
	//yz::opengl::drawMeshEdgeFromFace(mesh3.vertex, mesh3.face, mesh3.vertex_normal);

	//glPushMatrix();
	//glTranslatef(-(marcher.dim.x - 1.), -(marcher.dim.y - 1.), -(marcher.dim.z - 1.));

	//for (unsigned int k = 0; k < marcher.dim.z; k++) {
	//	for (unsigned int j = 0; j < marcher.dim.y; j++) {
	//		for (unsigned int i = 0; i < marcher.dim.x; i++) {
	//			yz::opengl::setSequentialDisplayColor(marcher.GetData(i, j, k));
	//			yz::opengl::drawNumber(marcher.GetData(i, j, k), yz::Vec3i(i, j, k) * 2);
	//		}
	//	}
	//}

	//glColor3f(1, 1, 1);
	//std::vector<yz::Vec3d> face_normal;
	//yz::geometry::calculateFaceNormal(face_normal, marcher.vertex, marcher.face);
	//yz::opengl::drawFlatShadingTriMesh(marcher.vertex, marcher.face, face_normal);
	//glColor3f(0, 0, 1);
	//yz::opengl::drawMeshEdgeFromFace(marcher.vertex, marcher.face);

	//for (int i = 0; i < 8; i++) {
	//	yz::opengl::setSequentialDisplayColor(cube.label[i]);
	//	//yz::opengl::drawPointAsSphere(vertex[i], 0.3);

	//	//glColor3f(0, 0, 0);
	//	yz::opengl::drawNumber(cube.label[i], vertex[i]);
	//}

	//glColor3f(1, 1, 1);
	//for (int i = 0; i < cube.edge.size(); i++) {
	//	yz::opengl::drawArrow(
	//		cube.edge[i].first,
	//		cube.edge[i].second,
	//		0.05
	//	);
	//}

	//for (int i = 0; i < cube.edge_loop.size(); i++) {
	//	yz::opengl::setSequentialDisplayColor(i);
	//	
	//	yz::Vec3d center;
	//	for (int j = 0; j < cube.edge_loop[i].size(); j++) {
	//		center += cube.edge_loop[i][j];
	//	}
	//	center /= cube.edge_loop[i].size();

	//	for (int j = 0; j < cube.edge_loop[i].size(); j++) {
	//		yz::Vec3d v1 = cube.edge_loop[i][j];
	//		yz::Vec3d v2 = cube.edge_loop[i][(j + 1) % cube.edge_loop[i].size()];
	//		yz::opengl::drawArrow(
	//			v1 + (center - v1) * 0.1,
	//			v2 + (center - v2) * 0.1,
	//			0.05
	//		);
	//	}
	//}


	//for (int i = 0; i < marcher.edge_loop.size(); i++) {
	//	yz::opengl::setSequentialDisplayColor(i);

	//	yz::Vec3d center;
	//	for (int j = 0; j < marcher.edge_loop[i].size(); j++) {
	//		center += marcher.edge_loop[i][j];
	//	}
	//	center /= marcher.edge_loop[i].size();

	//	for (int j = 0; j < marcher.edge_loop[i].size(); j++) {
	//		yz::Vec3d v1 = marcher.edge_loop[i][j];
	//		yz::Vec3d v2 = marcher.edge_loop[i][(j + 1) % marcher.edge_loop[i].size()];
	//		yz::opengl::drawArrow(
	//			v1 + (center - v1) * 0.1,
	//			v2 + (center - v2) * 0.1,
	//			0.05
	//		);
	//	}
	//}


	//glPopMatrix();
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case ' ':
		for (int i = 0; i < marcher.data.size(); i++)
			marcher.data[i] = rand() % 2;
		marcher.Marcher(2, yz::int3(-1, -1, -1));
		break;
	case 'd':
		marcher.WriteTriMeshToObj("marcher_mesh.obj");
		break;
	}
}

void special(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_LEFT:
		idx--;
		break;
	case GLUT_KEY_RIGHT:
		idx++;
		break;
	}
}

int main() {
	mesh1.ReadMeshFromFile("mesh1.obj");
	mesh2.ReadMeshFromFile("mesh2.obj");
	mesh3.ReadMeshFromFile("mesh3.obj");
	aabb1.BuildTriangleAABBTree(mesh1.vertex, mesh1.face);
	aabb2.BuildTriangleAABBTree(mesh2.vertex, mesh2.face);
	aabb3.BuildTriangleAABBTree(mesh3.vertex, mesh3.face);

	marcher.SetupVolume(40, 30, 40);
	for (int k = 0; k < marcher.dim.z; k++)
		for (int j = 0; j < marcher.dim.y; j++)
			for (int i = 0; i < marcher.dim.x; i++) {
				yz::Vec3d xyz(i*0.1, j*0.1, k*0.1);

				int inside1 = yz::geometry::isPointInsideClosedMesh(
					xyz, mesh1.vertex, mesh1.face, aabb1);
				int inside2 = yz::geometry::isPointInsideClosedMesh(
					xyz, mesh2.vertex, mesh2.face, aabb2);
				int inside3 = yz::geometry::isPointInsideClosedMesh(
					xyz, mesh3.vertex, mesh3.face, aabb3);

				int code = (inside1 << 2) | (inside2 << 1) | inside3;
				//int code = 0;
				//if (inside1)
				//	code = 1;
				//if (inside2)
				//	code = 2;
				//if (inside3)
				//	code = 3;

				marcher.data[marcher.GetVoxelID(i, j, k)] = code;
			}

	marcher.Marcher(0.1, yz::int3(0, 0, 0));

	mesh.vertex = marcher.vertex;
	mesh.face = marcher.face;
	mesh.face_label = marcher.face_diff_label;
	mesh.CalculateNormals();
	mesh.CalculateFaceColor();
	mesh.CreateSeperateFace();

	win3d.keyboardFunc = keyboard;
	win3d.specialFunc = special;
	win3d.CreateGLUTWindow();
	win3d.SetDraw(draw);
	win3d.SetDrawAppend(print);

	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}

template <typename T>
void SmoothPreserveVolume(
	T							volume,
	std::vector<yz::Vec3<T>>&	vertex,
	std::vector<yz::Vec3<T>>&	vertex_normal,
	std::vector<yz::int3>&		face,
	std::vector<yz::Vec3<T>>&	face_normal,
	T							coef
) {
	std::vector<yz::Vec3<T>>	tmp_vertex;
	std::vector<int>			count;
	tmp_vertex.resize(vertex.size());
	count.resize(vertex.size());
	for (int i = 0; i < face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			count[face[i][j]] += 2;
			tmp_vertex[face[i][j]] += vertex[face[i][(j + 1) % 3]];
			tmp_vertex[face[i][j]] += vertex[face[i][(j + 2) % 3]];
		}
	}

	for (int i = 0; i < vertex.size(); i++) {
		vertex[i] = tmp_vertex[i] / count[i] * coef + vertex[i] * (1.0 - coef);
	}

	yz::geometry::calculateFaceNormal(face_normal, vertex, face);
	yz::geometry::calculateVertexNormal(vertex_normal, vertex, face);

	yz::Vec2<T> mesh_data = CalculateMeshVolumeSurfaceArea(vertex, face, face_normal);
	T delta_vol = mesh_volume - mesh_data.y;
	T surface_area = mesh_data.x;
	T offset = delta_vol / surface_area;

	for (int i = 0; i < vertex.size(); i++)
		vertex[i] += vertex_normal[i] * offset;
}

template<typename T>
yz::Vec2<T> CalculateMeshVolumeSurfaceArea(
	const std::vector<yz::Vec3<T>>&		vertex,
	const std::vector<yz::int3>&		face,
	const std::vector<yz::Vec3<T>>&		face_normal
) {
	T surface_area_sum = 0;
	T volume_sum = 0;
	for (int i = 0; i < face.size(); i++) {
		yz::Vec3<T> r1 = vertex[face[i].y] - vertex[face[i].x];
		yz::Vec3<T> r2 = vertex[face[i].z] - vertex[face[i].x];
		T tri_area = yz::cross(r1, r2).Length();
		T tet_vol = tri_area * yz::dot(face_normal[i], vertex[face[i].x]);
		surface_area_sum += tri_area;
		volume_sum += tet_vol;
	}
	return yz::Vec2<T>(surface_area_sum / 2, volume_sum / 6);
}

int IsManifold(yz::geometry::TriMeshd& mesh) {
	std::unordered_set<yz::int2, yz::utils::BitwiseHasher<yz::int2>> hash1, hash2, hash3;
	for (int i = 0; i < mesh.face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int v0 = mesh.face[i][j];
			int v1 = mesh.face[i][(j + 1) % 3];
			if (v0 > v1)
				yz::mySwap(v0, v1);
			yz::int2 edge(v0, v1);
			if (hash3.find(edge) != hash3.end()) {
				//	do nothing
			}
			else if (hash2.find(edge) != hash2.end()) {
				hash2.erase(hash2.find(edge));
				hash3.insert(edge);
			}
			else if (hash1.find(edge) != hash1.end()) {
				hash1.erase(hash1.find(edge));
				hash2.insert(edge);
			}
			else {
				hash1.insert(edge);
			}
		}
	}

	if (hash1.empty() && hash3.empty())
		return 1;

	return 0;
}
