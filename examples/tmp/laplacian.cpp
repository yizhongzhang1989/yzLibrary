#include <iostream>
#include <GL/glut.h>
#include <mkl.h>
#include <mkl_spblas.h>
#include <yzLib/yz_lib.h>

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::CurvTriMeshf					curv_mesh;

yz::geometry::SmoothShadingTriMesh<float>	mesh;
yz::SparseMatrixCSR<float>					lap_mat;
yz::DenseVector<float>						value;

std::vector<yz::float3>						vertex_color;


void draw(){
	glColor3f(1, 1, 1);

	yz::opengl::drawSmoothColorTriMesh(
		curv_mesh.vertex, 
		curv_mesh.face,
		vertex_color );

	//mesh.Display();

	glColor3f(1, 0, 0);
	yz::opengl::drawMeshEdgeFromFace(
		mesh.vertex, 
		mesh.face,
		mesh.vertex_normal );

	
}

int main(){
	//	setup data
	mesh.ReadMeshFromFile("mesh.obj");
	yz::geometry::createLaplacianMatrixForMesh(lap_mat, mesh.vertex, mesh.face);

	curv_mesh.vertex = mesh.vertex;
	curv_mesh.face = mesh.face;
	curv_mesh.CalculateMeanCurvature();

	vertex_color.resize( curv_mesh.vertex.size() );
	yz::utils::convertToColorArrayJet( 
		(float*)&vertex_color[0], 
		&curv_mesh.mixed_area[0],
		curv_mesh.mixed_area.size(),
		0, 
		0.001);

	win3d.SetDraw(draw);
	win3d.CreateGLUTWindow();


	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}
