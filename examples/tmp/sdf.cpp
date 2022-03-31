#include <iostream>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::field::SignedDistanceField<double>	sdf;
yz::geometry::TriMesh<double>	mesh;

std::vector<yz::Vec3d> xyz;
std::vector<double> dist;

namespace yz {
	namespace geometry {
		namespace field {

			/**
			Marching Cube Algorithm
			*/
			template<class T>
			class MarchingCube {
			public:
				/**
				*/
				MarchingCube(T* sdf, uint3 dim, T voxel_size, Vec3<T> xyz0) {
					sdf_ptr = sdf;
					this->dim = dim;
					this->voxel_size = voxel_size;
					this->xyz0 = xyz0;
				}

				void ExtractIsoSurface(std::vector<Vec3<T>>& vertex, std::vector<int3>& face) {
				}

			protected:
				Vec3ui	dim;		///<	dimension of the volume
				T*		sdf_ptr;	///<	distance to the surface of each voxel center to the surface

				T voxel_size;		///<	the size of each voxel
				Vec3<T>	xyz0;		///<	the left most corner of the volume

			};


		}
	}
}	//	namespace yz::geometry::field

int dim_id = 0;
int slice_id = 0;

void special(int key, int x, int y){
	switch (key){
	case GLUT_KEY_LEFT:
		slice_id--;
		break;
	case GLUT_KEY_RIGHT:
		slice_id++;
		break;
	}
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case '0':
	case 'x':
	case 'X':
		dim_id = 0;
		break;
	case '1':
	case 'y':
	case 'Y':
		dim_id = 1;
		break;
	case '2':
	case 'z':
	case 'Z':
		dim_id = 2;
		break;
	}
}

void print(){
	glColor3f(1, 0, 0);
	yz::opengl::printInfo(0, 0, "slice_id: %d", slice_id);
}

void draw(){
	yz::opengl::drawXYZAxis();

	glColor3f(1, 1, 1);
	yz::opengl::drawMeshEdgeFromFace(mesh.vertex, mesh.face, mesh.vertex_normal, 0);
	//mesh.DisplayFlat();

	//sdf.DrawSliceNormal(slice_id);
	sdf.DrawSlice(dim_id, slice_id, 8, 8);

	glDisable(GL_LIGHTING);
	glColor3f(1, 1, 1);
	sdf.DrawVolume();
	glEnable(GL_LIGHTING);
}

int main(){
	mesh.ReadMeshFromFile("../FEM/torus.obj");
	yz::geometry::applyNoiseToVertices(mesh.vertex, 1e-5);
	//sdf.SetupVolume(20, 30, 40, 0.1, yz::Vec3d(-1, -1, -1));
	//sdf.CreateSDFFromClosedMesh(mesh.vertex, mesh.face);
	//sdf.WriteSDFToBinaryFile("torus.sdf");
	sdf.ReadSDFFromBinaryFile("torus.sdf");

	//xyz.resize(101);
	//dist.resize(101);
	//for (int i = 0; i <= 100; i++){
	//	xyz[i] = yz::Vec3d(i / 50.0 - 1.0, 0.125, 0);
	//	dist[i] = sdf.GetSDFInterpTriLinear(xyz[i]);
	//}
	//std::ofstream file("dump.txt");
	//for (int i = 0; i < dist.size(); i++)
	//	file << dist[i] << std::endl;

	yz::geometry::field::MarchingCube<double> marching_cube(&sdf.data[0], sdf.dim.x, sdf.dim.y, sdf.dim.z);

	win3d.SetDrawAppend(print);
	win3d.SetDraw(draw);
	win3d.specialFunc = special;
	win3d.keyboardFunc = keyboard;
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(win3d.idleFunc);
	manager.EnterMainLoop();
}