/***********************************************************/
/**	\file
	\brief		Marching Cube
	\author		Yizhong Zhang
	\date		12/18/2017
*/
/***********************************************************/
#ifndef __YZ_MARCHING_CUBE_H__
#define __YZ_MARCHING_CUBE_H__

#include <iostream>
#include <vector>
#include <map>
#include "yzLib/yz_math/yz_vector.h"


namespace yz {	namespace geometry { namespace field {


/**
	Single cube multi-material marching cube algorithm

	Given the label of 8 corners of a (0,0,0)-(2,2,2) cube, we create non-manifold faces
	(represented using edge loop) using multi-material marching cube algorithm. The vertices
	of edge loop are all integer vertices inside the cube. The direction of the edge loop is
	pointing to smaller label region (counter-clockwise normal). To create triangle face,
	you should insert a center point for each edge loop, and triangulate the edge loop.
	The label seperated by each edge loop can be read from edge_loop_label, represented using int2,
	while x is the label under face normal, y is the label in front of face normal

	How to use: \n
	1, setup labels for each corner. coordinate xyz arranged in : 000, 100, 010, 110, 001, 101, 011, 111	\n
	2, call Marcher(), then edge_loop will be generated		\n
	3, read the data from edge_loop and edge_loop_label \n
*/
class MultiMaterialSingleCubeMarcher {
public:
	int									label[8];			///< input: coordinate xyz arranged in : 000, 100, 010, 110, 001, 101, 011, 111
	std::vector<std::vector<Vec3i>>		edge_loop;			///< output: edge loop of each face in (0,0,0)-(2,2,2) cube
	std::vector<int2>					edge_loop_label;	///< output: label on different side of the edge loop, (x, y) - (under face normal, up face normal)

	void Marcher() {
		edge.clear();
		edge_loop.clear();
		edge_loop_label.clear();

		CreateEdge();

		CreateLoop();
	}

private:
	struct EdgeLabel{
		Vec3i	v0;			///< start vertex of the edge
		Vec3i	v1;			///< end vertex of the edge
		int2	label_r_l;	///< label on each side of the edge, seeing in the direction of the edge, label on (right, left)
	};

	std::vector<EdgeLabel>	edge;			///< edges extracted on each square

private:
	/**
		Create all edges

		Create edge for each square face
	*/
	void CreateEdge() {
		//	XY plane
		SquareMarcher(int3(0, 0, 0), int3(2, 0, 0), int3(0, 2, 0), int3(2, 2, 0), label[0], label[1], label[2], label[3]);
		SquareMarcher(int3(0, 2, 2), int3(2, 2, 2), int3(0, 0, 2), int3(2, 0, 2), label[6], label[7], label[4], label[5]);
		//	XZ plane
		SquareMarcher(int3(0, 0, 2), int3(2, 0, 2), int3(0, 0, 0), int3(2, 0, 0), label[4], label[5], label[0], label[1]);
		SquareMarcher(int3(0, 2, 0), int3(2, 2, 0), int3(0, 2, 2), int3(2, 2, 2), label[2], label[3], label[6], label[7]);
		//	YZ plane
		SquareMarcher(int3(0, 0, 0), int3(0, 2, 0), int3(0, 0, 2), int3(0, 2, 2), label[0], label[2], label[4], label[6]);
		SquareMarcher(int3(2, 0, 2), int3(2, 2, 2), int3(2, 0, 0), int3(2, 2, 0), label[5], label[7], label[1], label[3]);
	}

	/**
		Marching a single square, and insert edges into the edge list

		Input vertices are arranged in the following order	\n

		v00---v10	\n
		|		|	\n
		|		|	\n
		v01---v11	\n
	*/
	void SquareMarcher(Vec3i v00, Vec3i v10, Vec3i v01, Vec3i v11, int l00, int l10, int l01, int l11) {
		//	 count number of distinct labels
		int label_counter[4] = { l00, l10, l01, l11 };
		std::sort(label_counter, label_counter + 4);
		int distinct_label_count = 1;
		for (int i = 1; i < 4; i++) {
			if (label_counter[i] != label_counter[i - 1])
				distinct_label_count++;
		}

		//	if the square only contain one label, do nothing
		if (distinct_label_count == 1)
			return;

		//	in the following, this square contain at least 2 labels, so edges must be inserted

		//	rotate the square, so that the first vertex is min
		int min_label = label_counter[0];
		Vec3i square_coord[4] = { v00, v10, v11, v01 };
		int square_label[4] = { l00, l10, l11, l01 };
		while (square_label[0] != min_label || square_label[3] == min_label) {
			int tmp_label = square_label[0];
			Vec3i tmp_coord = square_coord[0];
			for (int i = 0; i < 3; i++) {
				square_label[i] = square_label[i + 1];
				square_coord[i] = square_coord[i + 1];
			}
			square_label[3] = tmp_label;
			square_coord[3] = tmp_coord;
		}

		//	insert edges according to label
		if (distinct_label_count == 2) {
			if (label_counter[1] != label_counter[0]) {			//	only 1 min value
				EdgeLabel e;
				e.v0 = (square_coord[0] + square_coord[3]) / 2;
				e.v1 = (square_coord[0] + square_coord[1]) / 2;
				e.label_r_l = int2(square_label[3], square_label[0]);
				edge.push_back(e);
			}
			else if (label_counter[2] != label_counter[1]) {	//	2 min values
				if (square_label[0] == square_label[1]) {	//	parallel to square edge
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[3]) / 2;
					e.v1 = (square_coord[1] + square_coord[2]) / 2;
					e.label_r_l = int2(square_label[3], square_label[0]);
					edge.push_back(e);
				}
				else {										//	diagonal
					EdgeLabel e;
					e.v0 = (square_coord[1] + square_coord[2]) / 2;
					e.v1 = (square_coord[0] + square_coord[1]) / 2;
					e.label_r_l = int2(square_label[1], square_label[2]);
					edge.push_back(e);
					e.v0 = (square_coord[0] + square_coord[3]) / 2;
					e.v1 = (square_coord[2] + square_coord[3]) / 2;
					e.label_r_l = int2(square_label[3], square_label[0]);
					edge.push_back(e);
				}
			}
			else if (label_counter[3] != label_counter[2]) {	//	3 min values
				EdgeLabel e;
				e.v0 = (square_coord[0] + square_coord[3]) / 2;
				e.v1 = (square_coord[2] + square_coord[3]) / 2;
				e.label_r_l = int2(square_label[3], square_label[0]);
				edge.push_back(e);
			}
			else {
				std::cout << "error: MultiMaterialSingleCubeMarcher::SquareMarcher, if (distinct_label_count == 2) {" << std::endl;
			}
		}
		else if (distinct_label_count == 3) {
			if (square_label[0] == square_label[2] || square_label[1] == square_label[3]) {		//	diagonal
				if (square_label[1] == square_label[3]) {	//	rotate the square to main diagonal
					int tmp_label = square_label[0];
					Vec3i tmp_coord = square_coord[0];
					for (int i = 0; i < 3; i++) {
						square_label[i] = square_label[i + 1];
						square_coord[i] = square_coord[i + 1];
					}
					square_label[3] = tmp_label;
					square_coord[3] = tmp_coord;
				}
				if (square_label[0] > square_label[1]) {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[1]) / 2;
					e.v1 = (square_coord[1] + square_coord[2]) / 2;
					e.label_r_l = int2(square_label[0], square_label[1]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[1] + square_coord[2]) / 2;
					e.v1 = (square_coord[0] + square_coord[1]) / 2;
					e.label_r_l = int2(square_label[1], square_label[0]);
					edge.push_back(e);
				}
				if (square_label[0] > square_label[3]) {
					EdgeLabel e;
					e.v0 = (square_coord[2] + square_coord[3]) / 2;
					e.v1 = (square_coord[0] + square_coord[3]) / 2;
					e.label_r_l = int2(square_label[0], square_label[3]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[3]) / 2;
					e.v1 = (square_coord[2] + square_coord[3]) / 2;
					e.label_r_l = int2(square_label[3], square_label[0]);
					edge.push_back(e);
				}
			}
			else {
				while (square_label[0] != square_label[1]) {	//	rotate the square so that 0 and 1 have the same label
					int tmp_label = square_label[0];
					Vec3i tmp_coord = square_coord[0];
					for (int i = 0; i < 3; i++) {
						square_label[i] = square_label[i + 1];
						square_coord[i] = square_coord[i + 1];
					}
					square_label[3] = tmp_label;
					square_coord[3] = tmp_coord;
				}
				if (square_label[1] > square_label[2]) {
					EdgeLabel e;
					e.v0 = (square_coord[1] + square_coord[2]) / 2;
					e.v1 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.label_r_l = int2(square_label[1], square_label[2]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.v1 = (square_coord[1] + square_coord[2]) / 2;
					e.label_r_l = int2(square_label[2], square_label[1]);
					edge.push_back(e);
				}
				if (square_label[3] > square_label[2]) {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.v1 = (square_coord[2] + square_coord[3]) / 2;
					e.label_r_l = int2(square_label[3], square_label[2]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[2] + square_coord[3]) / 2;
					e.v1 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.label_r_l = int2(square_label[2], square_label[3]);
					edge.push_back(e);
				}
				if (square_label[0] > square_label[3]) {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.v1 = (square_coord[0] + square_coord[3]) / 2;
					e.label_r_l = int2(square_label[0], square_label[3]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[3]) / 2;
					e.v1 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.label_r_l = int2(square_label[3], square_label[0]);
					edge.push_back(e);
				}
			}
		}
		else {	//	4 distinct labels
			for (int i = 0; i < 4; i++) {
				int l = i, r = (i + 1) % 4;
				if (square_label[l] > square_label[r]) {
					EdgeLabel e;
					e.v0 = (square_coord[l] + square_coord[r]) / 2;
					e.v1 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.label_r_l = int2(square_label[l], square_label[r]);
					edge.push_back(e);
				}
				else {
					EdgeLabel e;
					e.v0 = (square_coord[0] + square_coord[1] + square_coord[2] + square_coord[3]) / 4;
					e.v1 = (square_coord[l] + square_coord[r]) / 2;
					e.label_r_l = int2(square_label[r], square_label[l]);
					edge.push_back(e);
				}
			}
		}
	}

	/**
		Create edge loop from existing edges.
	*/
	void CreateLoop() {
		//	count number of face center points
		int face_center_test[6] = { 0,0,0,0,0,0 };
		for (int i = 0; i < edge.size(); i++) {
			if (FaceCenterPointIndex(edge[i].v0) >= 0)
				face_center_test[FaceCenterPointIndex(edge[i].v0)] ++;
			if (FaceCenterPointIndex(edge[i].v1) >= 0)
				face_center_test[FaceCenterPointIndex(edge[i].v1)] ++;
		}
		int face_center_count = 0;	//	number of faces that contain node
		int face_center_node_count = 0;			//	number of edge nodes of face center
		for (int i = 0; i < 6; i++) {
			face_center_count += face_center_test[i] ? 1 : 0;
			face_center_node_count += face_center_test[i];
		}

		//	find all loops that start and stop at face center
		while (face_center_node_count) {
			//	find an edge that start from face center
			for (std::vector<EdgeLabel>::iterator iter = edge.begin(); iter != edge.end(); iter++) {
				if (FaceCenterPointIndex(iter->v0) >= 0) {
					edge_loop.resize(edge_loop.size() + 1);
					edge_loop_label.push_back(iter->label_r_l);
					edge_loop.back().push_back(iter->v0);
					edge_loop.back().push_back(iter->v1);
					edge.erase(iter);
					break;
				}
			}

			//	find new edges until the loop closes
			while (FaceCenterPointIndex(edge_loop.back().back()) < 0) {
				for (std::vector<EdgeLabel>::iterator iter = edge.begin(); iter != edge.end(); iter++) {
					if (iter->v0 == edge_loop.back().back()) {	//	this edge can extend the loop
						edge_loop.back().push_back(iter->v1);
						edge.erase(iter);
						break;
					}
				}
			}

			face_center_node_count -= 2;

			//	more than three center faces, this loop should go through volume center
			if (face_center_count > 2)
				edge_loop.back().push_back(int3(1, 1, 1));
		}

		//	find all close loops
		while (!edge.empty()) {
			//	create an initial point
			edge_loop.resize(edge_loop.size() + 1);
			edge_loop_label.push_back(edge.back().label_r_l);
			edge_loop.back().push_back(edge.back().v0);
			edge_loop.back().push_back(edge.back().v1);
			edge.pop_back();

			//	continue until the loop closed
			while (edge_loop.back().front() != edge_loop.back().back()) {
				for (std::vector<EdgeLabel>::iterator iter = edge.begin(); iter != edge.end(); iter++) {
					if (iter->v0 == edge_loop.back().back()) {	//	this edge can extend the loop
						edge_loop.back().push_back(iter->v1);
						edge.erase(iter);
						break;
					}
				}
			}

			//	remove the duplicate node
			edge_loop.back().pop_back();
		}
	}

	/**
		Give each face center point a unique index
	*/
	inline int FaceCenterPointIndex(Vec3i p) {
		if (p.y == 1 && p.z == 1) {
			return p.x == 0 ? 0 : 1;
		}
		if (p.x == 1 && p.z == 1) {
			return p.y == 0 ? 2 : 3;
		}
		if (p.x == 1 && p.y == 1) {
			return p.z == 0 ? 4 : 5;
		}
		return -1;
	}

};

/**
	Perform multi-material marching cube algorithm

	Given a labeled volume, extract a non-manifold triangle mesh that seperate different labels.

	How to use: \n
	1,	setup the volume label by call SetupVolume() or setup label directly \n
	2,	call Marcher() to create the mesh, size of the volume can be set as parameters \n
	3,	read the mesh from vertex and face, normal of each face points to smaller label region. 
		face_diff_label is the label that this face seperated represented using int2, x : under the face, y: in front of the face \n
*/
class MultiMaterialCubeMarcher : public LabelField {
public:
	std::vector<Vec3d>	vertex;				///< output: vertex of generated mesh
	std::vector<int3>	face;				///< output: face of generated mesh
	std::vector<int2>	face_diff_label;	///< output: label of each face, (x, y) - (label under face normal, label in front of face normal)

	/**
		Setup the volume

		\param	dim_x		dimension in x direction
		\param	dim_y		dimension in y direction
		\param	dim_z		dimension in z direction
		\param	label_ptr	label array of each voxel
	*/
	void SetupVolume(
		unsigned int	dim_x,
		unsigned int	dim_y,
		unsigned int	dim_z,
		const int*		label_ptr = NULL)
	{
		dim = uint3(dim_x, dim_y, dim_z);
		data.resize(dim.x * dim.y * dim.z);

		if (label_ptr)
			memcpy(&data[0], label_ptr, sizeof(int)*data.size());
	}

	/**
		Marching a labeled volume

		After calling this function, vertex face and face_diff_label will be built.

		\param	voxel_size		voxel size of the volume
		\param	xyz0			min coordinate of the volume
	*/
	void Marcher(double voxel_size = 1.0, Vec3d xyz0 = Vec3d(0, 0, 0)) {
		CreateEdgeLoops();
		CreateTriangles();

		//	scale the mesh
		for (unsigned int i = 0; i < vertex.size(); i++) {
			vertex[i] = xyz0 + vertex[i] * 0.5 * voxel_size;
		}
	}

	/**
		Read obj file to the mesh that written by WriteTriMeshFromFile()
	*/
	int ReadTriMeshFromObj(const char* file_name) {
		std::ifstream obj(file_name);
		if (!obj.is_open()) {
			std::cout << "error: MultiMaterialCubeMarcher::ReadTriMeshFromFile, cannot open " << file_name << std::endl;
			return 0;
		}

		std::string line;
		std::getline(obj, line);
		if (line != "# yz::geometry::field::MultiMaterialCubeMarcher") {
			obj.close();
			std::cout << "error: MultiMaterialCubeMarcher::ReadTriMeshFromFile, " << file_name << " is not legal MultiMaterialCubeMarcher mesh" << std::endl;
			return 0;
		}

		vertex.clear();
		face.clear();
		face_diff_label.clear();

		while (std::getline(obj, line)) {
			std::stringstream ss;
			std::string cmd;
			ss << line;
			ss >> cmd;

			if (cmd == "v") {			//	got a vertex, insert into vertex
				Vec3d xyz;
				ss >> xyz[0] >> xyz[1] >> xyz[2];
				vertex.push_back(xyz);
			}
			else if (cmd == "f") {		//	got a face, insert into face
				int3 xyz;
				ss >> xyz[0] >> xyz[1] >> xyz[2];
				face.push_back(int3(xyz[0] - 1, xyz[1] - 1, xyz[2] - 1));
			}
			else if (cmd == "fl") {		//	got a face label, insert into face_diff_label
				int2 fl;
				ss >> fl[0] >> fl[1];
				face_diff_label.push_back(fl);
			}
		}

		obj.close();

		//	check whether face size match face_diff_label
		if (face.size() != face_diff_label.size()) {
			std::cout << "error: MultiMaterialCubeMarcher::ReadTriMeshFromFile, " << file_name << " face.size() != face_labe.size()" << std::endl;
			return 0;
		}

		return 1;
	}

	/**
		Write the mesh to obj file

		embed flag inside this file
	*/
	int WriteTriMeshToObj(const char* file_name) {
		std::ofstream obj(file_name);
		if (!obj.is_open()) {
			std::cout << "error: MultiMaterialCubeMarcher::WriteTriMeshFromFile, cannot open " << file_name << std::endl;
			return 0;
		}

		// the first line serve as a flag, which will be checked when reading this file
		obj << "# yz::geometry::field::MultiMaterialCubeMarcher" << std::endl;	
		obj << "# obj file generated in yzLib" << std::endl;

		for (unsigned int i = 0; i != vertex.size(); i++)
			obj << "v " << vertex[i].x << ' ' << vertex[i].y << ' ' << vertex[i].z << std::endl;

		for (unsigned int i = 0; i != face.size(); i++)
			obj << "f " << face[i].x + 1 << ' ' << face[i].y + 1 << ' ' << face[i].z + 1 << std::endl;

		for (unsigned int i = 0; i != face_diff_label.size(); i++)
			obj << "fl " << face_diff_label[i].x << ' ' << face_diff_label[i].y << std::endl;

		obj.close();

		return 1;
	}

private:
	/**
		Create edge loops for the whole volume
	*/
	void CreateEdgeLoops() {
		edge_loop.clear();
		edge_loop_label.clear();

		MultiMaterialSingleCubeMarcher cube;
		for (unsigned int k = 1; k < dim.z; k++) {
			for (unsigned int j = 1; j < dim.y; j++) {
				for (unsigned int i = 1; i < dim.x; i++) {
					//	march each cube
					cube.label[0] = GetData(i - 1, j - 1, k - 1);
					cube.label[1] = GetData(i, j - 1, k - 1);
					cube.label[2] = GetData(i - 1, j, k - 1);
					cube.label[3] = GetData(i, j, k - 1);
					cube.label[4] = GetData(i - 1, j - 1, k);
					cube.label[5] = GetData(i, j - 1, k);
					cube.label[6] = GetData(i - 1, j, k);
					cube.label[7] = GetData(i, j, k);
					cube.Marcher();

					//	offset the cube
					Vec3i offset(i * 2 - 2, j * 2 - 2, k * 2 - 2);
					for (unsigned int m = 0; m != cube.edge_loop.size(); m++) {
						for (unsigned int n = 0; n != cube.edge_loop[m].size(); n++) {
							cube.edge_loop[m][n] += offset;
						}
					}

					//	record edge loops
					edge_loop.insert(edge_loop.end(), cube.edge_loop.begin(), cube.edge_loop.end());
					edge_loop_label.insert(edge_loop_label.end(), cube.edge_loop_label.begin(), cube.edge_loop_label.end());
				}
			}
		}
	}

	/**
		Create triangles from edge loops
	*/
	void CreateTriangles() {
		vertex.clear();
		face.clear();
		face_diff_label.clear();
		vertex_hash.clear();

		for (unsigned int j = 0; j != edge_loop.size(); j++) {
			if (edge_loop[j].size() < 3) {
				std::cout << "error: MultiMaterialCubeMarcher::CreateTriangles, unclosed edge loop detected" << std::endl;
				continue;
			}

			//	if the edge loop only contain 3 vertices, create the triangle directly
			if (edge_loop[j].size() == 3) {
				face.push_back(int3(
					SetVertexIndex(edge_loop[j][0]),
					SetVertexIndex(edge_loop[j][1]),
					SetVertexIndex(edge_loop[j][2])
				));
				face_diff_label.push_back(edge_loop_label[j]);
				continue;
			}

			//	if the edge loop contian more than 3 vertices, insert a center point
			yz::Vec3d center;
			for (unsigned int i = 0; i != edge_loop[j].size(); i++) {
				center += edge_loop[j][i];
			}
			center /= edge_loop[j].size();

			unsigned int center_vertex_index = vertex.size();
			vertex.push_back(center);

			for (int i = 0; i != edge_loop[j].size(); i++) {
				face.push_back(int3(
					center_vertex_index,
					SetVertexIndex(edge_loop[j][i]),
					SetVertexIndex(edge_loop[j][(i + 1) % edge_loop[j].size()])
				));
				face_diff_label.push_back(edge_loop_label[j]);
			}
		}
	}

	/**
		Set vertex.

		For an existing vertex, return the index directly.
		For an un-existing vertex, create the vertex then return the index
	*/
	inline unsigned int SetVertexIndex(Vec3i v) {
		std::unordered_map<Vec3i, unsigned int, utils::BitwiseHasher<Vec3i>>::iterator iter = vertex_hash.find(v);
		if (iter == vertex_hash.end()) {	//	this vertex not exist yet, create the vertex
			int nv_index = vertex.size();
			vertex_hash.insert(std::pair<Vec3i, unsigned int>(v, nv_index));
			vertex.push_back(v);
			return nv_index;
		}
		else {	//	this vertex already exist, return its index
			return iter->second;
		}
	}

private:
	std::vector<std::vector<Vec3i>>		edge_loop;
	std::vector<int2>					edge_loop_label;
	std::unordered_map<Vec3i, unsigned int, utils::BitwiseHasher<Vec3i>> vertex_hash;
};


}}}	//	namespace yz::geometry::field



#endif	//	__YZ_MARCHING_CUBE_H__