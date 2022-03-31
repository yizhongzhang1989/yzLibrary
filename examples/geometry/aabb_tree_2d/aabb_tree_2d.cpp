#include <iostream>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>
#include "data_path.h"


typedef float REAL;

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow2D		win2d;

yz::geometry::AABBTree2D<REAL>	aabb;
yz::geometry::TriMesh2D<REAL>	mesh;

std::vector<yz::float3>			point_color;
std::vector<yz::Vec2<REAL>>		sample_point;

int disp_depth = 0;

template<class T>
class Polygon : public yz::geometry::polygon::Polygon<T>{
public:

	yz::geometry::AABBTree2D<T>	edge_aabb;

public:
	void CreatePolygonFromMeshEdge(const yz::geometry::TriMesh2D<REAL>& mesh){
		//vertex = mesh.vertex;
		////yz::geometry::applyNoiseToVertices(vertex, 1e-4);

		//yz::geometry::TriMeshEdge mesh_edge;
		//mesh_edge.CreateEdge(mesh.face);

		//yz::geometry::getBoundaryEdges(edge, mesh_edge.edge, mesh.face);
		yz::geometry::polygon::createPolygonAs2DMeshBoundary(this->vertex, this->edge, mesh.vertex, mesh.face);

		edge_aabb.BuildEdgeAABBTree(this->vertex, this->edge, 1e-4);
	}

	void Draw(){
		this->Display();
	}
};

Polygon<REAL> polygon;

void draw(){
	//glColor3f(0, 0, 1);
	//mesh.Display();

	//glColor3f(1, 0, 0);
	//aabb.Display(disp_depth);
	////for(int i=0; i<aabb.leaf_number; i++)
	////	yz::opengl::drawAABBWire(aabb.node[i].bb_min, aabb.node[i].bb_max);


	for(int i=0; i<sample_point.size(); i++){
		glColor3fv((float*)&point_color[i]);
		yz::opengl::drawPointAsBall(sample_point[i], 0.003);

		////	draw ray
		//yz::opengl::drawRay(sample_point[i], sample_ray_next[i], 1);

		////	draw inter p
		//for(int j=0; j<sample_ray_inter_p[i].size(); j++){
		//	yz::opengl::drawPointAsBall(sample_ray_inter_p[i][j], 0.003);
		//}
	}


	glColor3f(0, 0, 0);
	polygon.Draw();

}

void special(int key, int x, int y){
	switch(key){
		case GLUT_KEY_UP:
			disp_depth++;
			break;
		case GLUT_KEY_DOWN:
			disp_depth--;
			break;
	}
}

int main(){
	char obj_filename[1024];
	sprintf(obj_filename, "%s/sheet.obj", data_path);

	mesh.ReadMeshFromFile(obj_filename, 'z');
	aabb.BuildTriangleAABBTree(mesh.vertex, mesh.face, 0);

	yz::utils::Timer timer;
	timer.Start();

	yz::geometry::remeshing::samplingPoissonDiskRobertBridson(sample_point, 
		REAL(0.03), aabb.node[aabb.leaf_number].bb_min, aabb.node[aabb.leaf_number].bb_max);

	timer.Stop();
	std::cout << "sampling: " << timer.Elapsed() << std::endl;
	std::cout << "points: " << sample_point.size() << std::endl;

	//	polygon
	polygon.CreatePolygonFromMeshEdge(mesh);

	timer.Start();

	//	calculate point color
	point_color.resize(sample_point.size());
	for(int i=0; i<sample_point.size(); i++){
		yz::Vec2<REAL> p;
		yz::geometry::polygon::getNearestPointOnPolygon(p, sample_point[i], polygon.vertex, polygon.edge, polygon.edge_aabb);
		REAL length = (p - sample_point[i]).Length();

		if( ! yz::geometry::polygon::isPointInsidePolygon(sample_point[i], polygon.vertex, polygon.edge, polygon.edge_aabb) )
			length = -length;

		yz::utils::convertToColorJet((float*)&point_color[i], length, -0.2, 0.2);
	}	

	timer.Stop();
	std::cout << "nearest: " << timer.Elapsed() << std::endl;

	win2d.SetDraw(draw);
	win2d.specialFunc = special;
	win2d.CreateGLUTWindow();

	manager.AddIdleFunc(win2d.idleFunc);
	manager.EnterMainLoop();
}