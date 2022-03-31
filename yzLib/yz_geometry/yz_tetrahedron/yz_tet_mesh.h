/***********************************************************/
/**	\file
	\brief		Tetrahedral Mesh
	\details	Data structure needed to represent a tetrahedral mesh
	\author		Yizhong Zhang
	\date		12/6/2016
*/
/***********************************************************/
#ifndef __YZ_TET_MESH_H__
#define __YZ_TET_MESH_H__

namespace yz {  namespace geometry {  namespace tetrahedron {


/**
	base tetrahedral mesh using vector, just vertex and tetrahedron
*/
template<class T>
class BaseTetMesh {
public:
	std::vector<Vec3<T>>	vertex;
	std::vector<int4>		tetrahedron;

public:
	//	constructor
	BaseTetMesh() {};

	//	functions same as TriMesh
	/**
		Reset all members, clear all memory
	*/
	virtual void Reset() {
		vertex.clear();
		tetrahedron.clear();
	}

	/**
		read TetMesh from file

		There is no standard tet mesh file format, so this function is left for user to implement
	*/
	virtual int ReadMeshFromFile(const char* file_name) {
		std::cout << "BaseTetMesh::ReadMeshFromFile() is not implemented" << std::endl;
		return	0;
	}

	virtual int Display(int display_mode = 0) const {
		std::cout << "BaseTetMesh con't have Display() function implemented" << std::endl;
		return 0;
	}

	virtual void AppendTetMesh(const BaseTetMesh<T>& mesh) {
		int old_vertex_number = vertex.size();
		int old_tet_number = tetrahedron.size();
		vertex.insert(vertex.end(), mesh.vertex.begin(), mesh.vertex.end());
		tetrahedron.insert(tetrahedron.end(), mesh.tetrahedron.begin(), mesh.tetrahedron.end());
		for (int t = old_tet_number; t < tetrahedron.size(); t++) {
			tetrahedron[t].x += old_vertex_number;
			tetrahedron[t].y += old_vertex_number;
			tetrahedron[t].z += old_vertex_number;
			tetrahedron[t].w += old_vertex_number;
		}
	}
};

/**
	The surface of a tet_mesh
*/
template<class T>
class TetMeshSurface {
public:
	std::vector<int3>		surface_face;
	std::vector<Vec3<T>>	surface_vertex_normal;
	std::vector<Vec3<T>>	surface_face_normal;

public:
	/**
		create the surface mesh of the tet_mesh
	*/
	void CreateSurfaceFace(
		const std::vector<Vec3<T>>&	vertex,
		const std::vector<int4>&	tetrahedron) 
	{
		surface_face.clear();

		//	we use a vector to record all faces and corresponding tetrahedron index
		std::vector<std::pair<int3, int>> triangle;
		triangle.resize(tetrahedron.size() * 4);

		//	collect triangles of each tatrahedron
		for (int i = 0; i < tetrahedron.size(); i++) {
			int tet[4] = { tetrahedron[i][0], tetrahedron[i][1], tetrahedron[i][2], tetrahedron[i][3] };
			std::sort(tet, tet + 4);
			triangle[(i << 2)] = std::pair<int3, int>(int3(tet[0], tet[1], tet[2]), i);
			triangle[(i << 2) | 0x01] = std::pair<int3, int>(int3(tet[0], tet[1], tet[3]), i);
			triangle[(i << 2) | 0x02] = std::pair<int3, int>(int3(tet[0], tet[2], tet[3]), i);
			triangle[(i << 2) | 0x03] = std::pair<int3, int>(int3(tet[1], tet[2], tet[3]), i);
		}

		//	sort the triangles to find duplicates
		std::sort(triangle.begin(), triangle.end());

		for (int i = 0; i < triangle.size(); ) {
			int j = i + 1;
			while (j < triangle.size() && triangle[j].first == triangle[i].first)
				j++;
			if (j == i + 1) {	//	triangle i is unique, it is a surface triangle
				int3 face(triangle[i].first);
				int vid = -1;
				for (int k = 0; k < 4; k++) {
					vid = tetrahedron[triangle[i].second][k];
					if (vid != face.x && vid != face.y && vid != face.z)
						break;
				}

				//	check the normal direction of this face
				Vec3<T> nor = cross(vertex[face.y] - vertex[face.x], vertex[face.z] - vertex[face.x]);
				Vec3<T> r = vertex[vid] - vertex[face.x];
				if (dot(nor, r) > 0)
					mySwap(face.y, face.z);

				//	record this face
				surface_face.push_back(face);
			}
			else {	//	triangle i is shared by more than one tetrahedron
				if (j > i + 2) {
					std::cout << "error: TetMeshSurface::CreateSurfaceFace, triangle shared by more than two tetrahedrons" << std::endl;
				}
			}

			i = j;
		}
	}

	/**
	*/
	void CalculateSurfaceNormals(
		const std::vector<Vec3<T>>& vertex) 
	{
		//	calculate face normal
		calculateFaceNormal(surface_face_normal, vertex, surface_face);

		//	calculate vertex normal
		surface_vertex_normal.clear();
		surface_vertex_normal.resize(vertex.size());
		for (int i = 0; i < surface_face.size(); i++) {
			Vec3<T> nor = surface_face_normal[i];
			surface_vertex_normal[surface_face[i].x] += nor;
			surface_vertex_normal[surface_face[i].y] += nor;
			surface_vertex_normal[surface_face[i].z] += nor;
		}
		for (int i = 0; i < surface_vertex_normal.size(); i++) {
			if (surface_vertex_normal[i].SquareLength() > 1e-6)
				surface_vertex_normal[i].SetNormalize();
		}
	}
};

/**
	tetrahedron mesh, and corresponding surface
*/
template<class T>
class TetMesh : 
	public BaseTetMesh<T>, 
	public TetMeshSurface<T>
{
public:
	/**
		find the surface triangles of this tet mesh
	*/
	void CreateSurface() {
		CreateSurfaceFace(this->vertex, this->tetrahedron);
		CalculateSurfaceNormals(this->vertex);
	}

};



}}}	//	namespace yz::geometry::tetrahedron


#endif	//	__YZ_TET_MESH_H__