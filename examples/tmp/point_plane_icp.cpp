#include <iostream>
#include <string>
#include <GL/glut.h>
#include <Eigen/Dense>
#include <yzLib/yz_lib.h>


yz::opengl::DemoWindowManager	manager;
yz::opengl::GLUTWindow3D<0>		win3d;

yz::geometry::TriMeshf			mesh_s, mesh_d;
yz::geometry::AABBTree3Df		mesh_d_aabb;

namespace yz {  namespace geometry {  namespace icp {

inline int solveTransFrom6x6(
	Matrix4x4d&	trans,
	double		ATA[][6],
	double		ATb[6]
) {
	trans.SetIdentity();

#ifdef YZ_eigen_dense_h		//	if Eigen/Dense is included
	//	solve the matrix by Eigen
	Eigen::Matrix<double, 6, 6, Eigen::RowMajor> A;
	Eigen::Matrix<double, 6, 1> b;

	for (int i = 0; i < 36; i++)
		A.data()[i] = ATA[0][i];
	for (int i = 0; i < 6; i++)
		b.data()[i] = ATb[i];

	double det = A.determinant();
	if (fabs(det) < 1e-15 || _isnan(det)) {		//	cannot fit
		std::cout << "error: solveTransFrom6x6, det error" << std::endl;
		return 0;
	}

	Eigen::Matrix<double, 6, 1> result = A.llt().solve(b).cast<double>();

	double alpha = result(0);
	double beta = result(1);
	double gamma = result(2);

	Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Rinc =
		(Eigen::Matrix<double, 3, 3, Eigen::RowMajor>)
		Eigen::AngleAxisd(gamma, Eigen::Vector3d::UnitZ()) *
		Eigen::AngleAxisd(beta, Eigen::Vector3d::UnitY()) *
		Eigen::AngleAxisd(alpha, Eigen::Vector3d::UnitX());
	Eigen::Vector3d tinc = result.tail<3>();

	//	write the result
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++)
			trans[j][i] = Rinc.data()[j * 3 + i];
		trans[j][3] = tinc(j);
	}

	return 1;

#else
	std::cout << "<Eigen/Dense> is required to call yz::geometry::icp::solveTransFrom6x6()" << std::endl;

	return 0;
#endif // YZ_eigen_dense_h


}

template <typename T>
int calculateTransformPointToPlane(
	Matrix4x4<T>&					trans,
	const std::vector<Vec3<T>>&		source,
	const std::vector<Vec3<T>>&		target,
	const std::vector<Vec3<T>>&		target_normal,
	const std::vector<T>&			weight
) {
	trans.SetIdentity();

	if (source.empty())
		return 0;
	if (source.size() != target.size() ||
		source.size() != target_normal.size() ||
		source.size() != weight.size())
		return 0;

	//	data matrix
	double ATA[6][6], ATb[6];
	for (int i = 0; i < 36; i++)
		ATA[0][i] = 0.0;
	for (int i = 0; i < 6; i++)
		ATb[i] = 0.0;

	//	fill the matrix
	for (int k = 0; k < source.size(); k++) {
		//	for each vertex, add to the matrix
		yz::Vec3d s = source[k];
		yz::Vec3d t = target[k];
		yz::Vec3d n = target_normal[k];
		double w	= weight[k];

		double row[7];
		*(yz::Vec3d*)&row[0] = w * yz::cross(s, n);
		*(yz::Vec3d*)&row[3] = w * n;
		row[6] = w * yz::dot(n, t - s);

		//	add row to matrix
		for (int j = 0; j < 6; j++) {
			for (int i = 0; i < 6; i++) {
				ATA[j][i] += row[j] * row[i];
			}
			ATb[j] += row[j] * row[6];
		}
	}

	//	solve
	Matrix4x4d trans_d;
	int res = solveTransFrom6x6(trans_d, ATA, ATb);

	if (!res)
		std::cout << "error: calculateTransformPointToPlane, solve error" << std::endl;

	trans = trans_d;

	return res;
}

template <typename T>
int calculateTransformPointToPlane(
	Matrix4x4<T>&					trans,
	const std::vector<Vec3<T>>&		source,
	const std::vector<Vec3<T>>&		target,
	const std::vector<Vec3<T>>&		target_normal
) {
	trans.SetIdentity();

	if (source.empty())
		return 0;
	if (source.size() != target.size() ||
		source.size() != target_normal.size())
		return 0;

	std::vector<T> weight;
	weight.resize(source.size(), 1);

	return calculateTransformPointToPlane(trans, source, target, target_normal, weight);
}


}}}

void fitStepPointPlane() {
	std::vector<yz::Vec3f> nearest_point, nearest_normal;
	nearest_point.resize(mesh_s.vertex.size());
	nearest_normal.resize(mesh_s.vertex.size());

	for (int k = 0; k<mesh_s.vertex.size(); k++) {
		//	for each vertex of mesh_s, find the nearest point on mesh_d
		yz::Vec3f	np, nn;
		int			nearest_face_id;

		yz::geometry::getNearestPointOnMesh(np, nearest_face_id,
			mesh_s.vertex[k], mesh_d.vertex, mesh_d.face, mesh_d_aabb);
		nn = mesh_d.face_normal[nearest_face_id];

		nearest_point[k] = np;
		nearest_normal[k] = nn;
	}

	//	calculate transform
	yz::Matrix4x4f R;
	yz::geometry::icp::calculateTransformPointToPlane(R, mesh_s.vertex, nearest_point, nearest_normal);

	yz::geometry::transformVertices(mesh_s.vertex, R);
	mesh_s.CalculateNormals();

}

void fitStepPointPoint() {
	std::vector<yz::Vec3f> nearest_point;
	nearest_point.resize(mesh_s.vertex.size());

	for (int k = 0; k<mesh_s.vertex.size(); k++) {
		//	for each vertex of mesh_s, find the nearest point on mesh_d
		yz::Vec3f	np;
		int			nearest_face_id;

		yz::geometry::getNearestPointOnMesh(np, nearest_face_id,
			mesh_s.vertex[k], mesh_d.vertex, mesh_d.face, mesh_d_aabb);

		nearest_point[k] = np;
	}

	//	calculate transform
	yz::Matrix4x4f R;
	yz::geometry::icp::calculateTransformPointToPoint(R, mesh_s.vertex, nearest_point);

	yz::geometry::transformVertices(mesh_s.vertex, R);
	mesh_s.CalculateNormals();

}


void keyboard(unsigned char key, int x, int y){
	switch(key){
	case 27:
		exit(0);
	case 'p':
		fitStepPointPoint();
		break;
	case 'f':
		fitStepPointPlane();
		break;
	}
}

void print(){
	glColor3f(1,1,1);
	yz::opengl::printInfo(0,0,"press 'f' to fit step");
}

void draw(){
	glColor3f(1, 0, 0);
	mesh_s.Display(GL_SMOOTH);

	glColor3f(0, 1, 0);
	mesh_d.Display(GL_SMOOTH);
}

int main(int argc, char* argv[]){
	yz::Matrix4x4d M, Rx, Ry, Rz, T;
	Rx.SetRotationDeg(yz::int3(1, 0, 0), 30);
	Ry.SetRotationDeg(yz::int3(0, 1, 0), 20);
	Rz.SetRotationDeg(yz::int3(0, 0, 1), 10);
	T.SetTranslation(0.1, 0, 0);
	M = T * Rx * Ry * Rz;

	//	point test
	std::vector<yz::Vec3d> point_source, point_target;
	point_source.resize(30);
	point_target.resize(point_source.size());
	for (int i = 0; i < point_source.size(); i++) {
		point_source[i] = yz::Vec3d(yz::rand0to1d(), yz::rand0to1d(), yz::rand0to1d());
		point_target[i] = M * point_source[i];
		//point_target[i] = M * point_source[i] + yz::Vec3d(yz::rand0to1d(), yz::rand0to1d(), yz::rand0to1d()) * 0.1;
		//point_target[i] = point_source[i];
	}

	yz::Matrix4x4d trans;
	yz::geometry::icp::calculateTransformPointToPoint(trans, point_source, point_target);

	std::cout << M << std::endl;
	std::cout << trans << std::endl;

	std::cout << trans.Det() << std::endl;

	return 0;

	//	setup 3d mesh
	//mesh_s.ReadMeshFromFile("bunny_s.stl");
	//mesh_d.ReadMeshFromFile("bunny_d.stl");
	mesh_s.ReadMeshFromFile("bunny_s.obj");
	mesh_d.ReadMeshFromFile("bunny_d.obj");
	mesh_d_aabb.BuildTriangleAABBTree( mesh_d.vertex, mesh_d.face, 0.001f );


	yz::geometry::transformVertices(mesh_s.vertex, M);
	mesh_s.CalculateNormals();
	mesh_d.CalculateNormals();

	//yz::geometry::writeTriMeshToFile(
	//	"bunny_d.stl", 
	//	&mesh_d.vertex, 
	//	&mesh_d.face, 
	//	&mesh_d.vertex_normal, 
	//	&mesh_d.face_normal);

	//	setup 3d window
	win3d.keyboardFunc = keyboard;
	win3d.SetDraw(draw);
	win3d.SetDrawAppend(print);
	win3d.CreateGLUTWindow();

	//	setup window manager
	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}

