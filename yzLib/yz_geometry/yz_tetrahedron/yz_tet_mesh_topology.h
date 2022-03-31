/***********************************************************/
/**	\file
	\brief		Topology of Tetrahedral Mesh
	\details	
	\author		Yizhong Zhang
	\date		12/6/2016
*/
/***********************************************************/
#ifndef __YZ_TET_MESH_TOPOLOGY_H__
#define __YZ_TET_MESH_TOPOLOGY_H__

#include <vector>
#include <algorithm>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_numerical_utils.h"

namespace yz {  namespace geometry {  namespace tetrahedron {

/**
	Topology of tetrahedron mesh
*/
class TetMeshEdge {
public:
	std::vector<int2>	edge;

public:
	void CreateEdge(const std::vector<int4>& tetrahedron) {
		//	gather all potential edges
		edge.resize(tetrahedron.size() * 6);	//	6 edges each tetrahedron
		for (int i = 0; i < tetrahedron.size(); i++) {
			edge[i * 6] = int2(tetrahedron[i][0], tetrahedron[i][1]);
			edge[i * 6 + 1] = int2(tetrahedron[i][0], tetrahedron[i][2]);
			edge[i * 6 + 2] = int2(tetrahedron[i][0], tetrahedron[i][3]);
			edge[i * 6 + 3] = int2(tetrahedron[i][1], tetrahedron[i][2]);
			edge[i * 6 + 4] = int2(tetrahedron[i][1], tetrahedron[i][3]);
			edge[i * 6 + 5] = int2(tetrahedron[i][2], tetrahedron[i][3]);
			for (int j = 0; j < 6; j++) {
				if (edge[i * 6 + j].x > edge[i * 6 + j].y)
					mySwap(edge[i * 6 + j].x, edge[i * 6 + j].y);
			}
		}

		//	sort edge
		std::sort(edge.begin(), edge.end());

		//	remove duplicates
		int edge_count = 0, scan_idx = 0;
		while (scan_idx < edge.size()) {
			int2 curr_edge = edge[scan_idx];
			scan_idx++;
			while (scan_idx < edge.size() && edge[scan_idx] == curr_edge)
				scan_idx++;
			edge[edge_count++] = curr_edge;
		}

		//	resize edge
		edge.resize(edge_count);
	}

	void Reset() {
		edge.clear();
	}
};

/**
	Vertex connectivity of tetrahedron mesh

	since this structure cannot stand alone without edge, we don't provide creation functino inside
*/
class TetMeshVVE {
public:
	std::vector<int>	vv, ve, vve_start;

public:
	void Reset() {
		vv.clear();
		ve.clear();
		vve_start.clear();
	}
};

/**
	Edge and vertex connectivity of tetrahedron mesh
*/
class TetMeshEdgeVVE : public TetMeshEdge, public TetMeshVVE {
public:
	
public:
	void CreateVVEFromEdge(int vertex_number) {
		//	edge must be created ahead
		if (edge.empty()) {
			vv.clear();
			ve.clear();
			vve_start.clear();
			vve_start.resize(vertex_number + 1, 0);
			return;
		}

		//	create a connectivity list
		std::vector<int3> connectivity;	//	v1, v2, edge_idx
		connectivity.resize(edge.size() * 2);
		for (int i = 0; i < edge.size(); i++) {
			connectivity[i << 1] = int3(edge[i].x, edge[i].y, i);			//	v1,v2,edge
			connectivity[(i << 1) | 0x01] = int3(edge[i].y, edge[i].x, i);	//	v2,v1,edge
		}

		//	sort according to vertex order
		std::sort(connectivity.begin(), connectivity.end());

		//	create vv, ve and vve_start
		vv.resize(connectivity.size());
		ve.resize(connectivity.size());
		vve_start.resize(vertex_number + 1);

		for (int i = 0; i<connectivity.size(); i++) {
			vv[i] = connectivity[i].y;
			ve[i] = connectivity[i].z;
		}

		vve_start[0] = 0;
		for (int i = 0, j = 0; i<vertex_number; i++) {
			while (j<connectivity.size() && connectivity[j].x == i)
				j++;
			vve_start[i + 1] = j;
		}
	}

	void CreateEdgeVVE(int vertex_number, const std::vector<int4>& tetrahedron) {
		CreateEdge(tetrahedron);
		CreateVVEFromEdge(vertex_number);
	}
};

/**
	edge, edge-vertex, edge-tet connectivity
*/
class TetMeshEdgeVVETEET : public TetMeshEdgeVVE {
public:
	std::vector<int>	te;		//	each tetrahedron has 6 adjacent edges
	std::vector<int>	et, et_start;

public:
	void CreateEdgeVVETEET(int vertex_number, const std::vector<int4>& tetrahedron) {
		std::vector<int3>	tmp_edge;	//	v1_idx, v2_idx, tet_idx, edge_idx

		//	gather all potential edges
		tmp_edge.resize(tetrahedron.size() * 6);	//	6 edges each tetrahedron
		for (int i = 0; i < tetrahedron.size(); i++) {
			tmp_edge[i * 6] = int3(tetrahedron[i][0], tetrahedron[i][1], i);
			tmp_edge[i * 6 + 1] = int3(tetrahedron[i][0], tetrahedron[i][2], i);
			tmp_edge[i * 6 + 2] = int3(tetrahedron[i][0], tetrahedron[i][3], i);
			tmp_edge[i * 6 + 3] = int3(tetrahedron[i][1], tetrahedron[i][2], i);
			tmp_edge[i * 6 + 4] = int3(tetrahedron[i][1], tetrahedron[i][3], i);
			tmp_edge[i * 6 + 5] = int3(tetrahedron[i][2], tetrahedron[i][3], i);
			for (int j = 0; j < 6; j++) {
				if (tmp_edge[i * 6 + j].x > tmp_edge[i * 6 + j].y)
					mySwap(tmp_edge[i * 6 + j].x, tmp_edge[i * 6 + j].y);
			}
		}

		//	sort edge
		std::sort(tmp_edge.begin(), tmp_edge.end());

		//	create edge, et, te
		edge.clear();
		et.clear();
		et_start.clear();
		te.clear();
		edge.reserve(tetrahedron.size() * 1.5);		//	approximate the number of edges
		et.reserve(tetrahedron.size() * 6);
		et_start.reserve(tetrahedron.size() * 1.5);
		te.resize(tetrahedron.size() * 6, -1);
		for (int i = 0; i < tmp_edge.size(); i++) {
			//	create new edge
			if (edge.empty() || edge.back() != int2(tmp_edge[i].x, tmp_edge[i].y)) {
				edge.push_back(int2(tmp_edge[i].x, tmp_edge[i].y));
				et_start.push_back(et.size());
				et.push_back(tmp_edge[i].z);
			}

			//	record te
			int tid = tmp_edge[i].z;
			for (int j = 0; j < 6; j++) {
				if (te[tid * 6 + j] < 0) {
					te[tid * 6 + j] = edge.size() - 1;
					break;
				}
			}

			//	record et
			et.push_back(tid);
		}
		et_start.push_back(et.size());

		//	create vv, ve
		vv.resize(edge.size() * 2);
		ve.resize(edge.size() * 2);
		vve_start.clear();
		vve_start.resize(vertex_number + 1, 0);
		for (int i = 0; i < edge.size(); i++) {
			vve_start[edge[i].x + 1] ++;
			vve_start[edge[i].y + 1] ++;
		}
		for (int i = 1; i < vve_start.size(); i++) {
			vve_start[i] += vve_start[i - 1];
		}
		for (int i = 0; i < edge.size(); i++) {
			int v1 = edge[i].x, v2 = edge[i].y;
			ve[vve_start[v1]] = i;
			ve[vve_start[v2]] = i;
			vv[vve_start[v1]++] = v2;
			vv[vve_start[v2]++] = v1;
		}
		for (int i = vve_start.size() - 1; i > 0; i--)
			vve_start[i] = vve_start[i - 1];
		vve_start[0] = 0;
	}
};

/**
	Triangle faces of each tet
*/
class TetMeshFace {
public:
	std::vector<int3>	face;

public:
	void CreateFace(const std::vector<int4>& tetrahedron) {
		face.clear();
		face.resize(tetrahedron.size() * 4);
		if (face.empty())
			return;

		//	collect triangles of each tatrahedron
		for (int i = 0; i < tetrahedron.size(); i++) {
			int tet[4] = { tetrahedron[i][0], tetrahedron[i][1], tetrahedron[i][2], tetrahedron[i][3] };
			std::sort(tet, tet + 4);
			face[(i << 2)] = int3(tet[0], tet[1], tet[2]);
			face[(i << 2) | 0x01] = int3(tet[0], tet[1], tet[3]);
			face[(i << 2) | 0x02] = int3(tet[0], tet[2], tet[3]);
			face[(i << 2) | 0x03] = int3(tet[1], tet[2], tet[3]);
		}

		//	sort the triangles to find duplicates
		std::sort(face.begin(), face.end());

		//	remove duplicates
		int count = 1;
		for (int i = 1; i < face.size(); i++) {
			if (face[i] != face[count - 1]) {
				face[count++] = face[i];
			}
		}
		face.resize(count);
	}

	void Reset() {
		face.clear();
	}
};

/**
	tetrahedrom-face connectivity
*/
class TetMeshTFFT{
public:
	std::vector<int4>	tf;
	std::vector<int2>	ft;

public:
	inline void Reset() {
		tf.clear();
		ft.clear();
	}
};

/**
	tet mesh face tet-face
*/
class TetMeshFaceTFFT : public TetMeshFace, public TetMeshTFFT {
public:
	void CreateFaceTFFT(const std::vector<int4>& tetrahedron) {
		face.clear();
		tf.clear();
		ft.clear();
		if (tetrahedron.empty())
			return;

		face.resize(tetrahedron.size() * 4);
		tf.resize(tetrahedron.size(), int4(-1, -1, -1, -1));
		ft.resize(tetrahedron.size() * 4, int2(-1, -1));

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

		//	create topology
		int face_count = 0;
		for (int i = 0; i < triangle.size(); i++) {
			int tid = triangle[i].second;

			//	record face
			if (!face_count || face[face_count - 1] != triangle[i].first) {
				face[face_count++] = triangle[i].first;
			}

			//	record t-f
			int j = 0;
			for (; j < 4; j++) {
				if (tf[tid][j] == -1) {
					tf[tid][j] = face_count - 1;
					break;
				}
			}
			if (j >= 4) {
				std::cout << "error: TetMeshFaceTFFT::CreateFaceTFFT record t-f" << std::endl;
			}

			//	record f-t
			j = 0;
			for (; j < 2; j++) {
				if (ft[face_count - 1][j] == -1) {
					ft[face_count - 1][j] = tid;
					break;
				}
			}
			if (j >= 2) {
				std::cout << "error: TetMeshFaceTFFT::CreateFaceTFFT record f-t" << std::endl;
			}
		}

		face.resize(face_count);
		ft.resize(face_count);
	}
};

/**
	vertex-tet connectivity
*/
class TetMeshVT {
public:
	std::vector<int>	vt, vt_start;

public:
	void CreateVT(int vertex_number, const std::vector<int4>& tetrahedron) {
		vt.clear();
		vt_start.clear();
		vt.resize(tetrahedron.size() * 4);
		vt_start.resize(vertex_number + 1, 0);
		if (vt.empty())
			return;

		std::vector<int2> tmp_vt;
		tmp_vt.reserve(tetrahedron.size() * 4);
		for (int i = 0; i<tetrahedron.size(); i++) {
			tmp_vt.push_back(int2(tetrahedron[i][0], i));
			tmp_vt.push_back(int2(tetrahedron[i][1], i));
			tmp_vt.push_back(int2(tetrahedron[i][2], i));
			tmp_vt.push_back(int2(tetrahedron[i][3], i));
		}

		std::sort(tmp_vt.begin(), tmp_vt.end());

		for (int i = 0; i<tmp_vt.size(); i++)
			vt[i] = tmp_vt[i].y;

		vt_start[0] = 0;
		for (int i = 0, j = 0; i<vertex_number; i++) {
			while (j<tmp_vt.size() && tmp_vt[j].x == i)
				j++;
			vt_start[i + 1] = j;
		}
	}

	inline void Reset() {
		vt.clear();
		vt_start.clear();
	}
};

/**
	tet-tet connectivity
*/
class TetMeshTT {
public:
	std::vector<int4>	tt;

public:
	void CreateTT(const std::vector<int4>& tf, const std::vector<int2>& ft) {
		tt.resize(tf.size(), int4(-1, -1, -1, -1));
		for (int i = 0; i < tf.size(); i++) {
			for (int j = 0; j < 4; j++) {
				int fid = tf[i][j];
				int nt = ft[fid][0] == i ? ft[fid][1] : ft[fid][0];
				tt[i][j] = nt;
			}
		}
	}

	inline void Reset() {
		tt.clear();
	}
};

/**
	topology of tet mesh
*/
class TetMeshTopology : 
	public TetMeshEdgeVVETEET, 
	public TetMeshFaceTFFT, 
	public TetMeshVT, 
	public TetMeshTT 
{
public:
	void CreateTopology(int vertex_number, const std::vector<int4>& tetrahedron) {
		CreateEdgeVVETEET(vertex_number, tetrahedron);
		//CreateEdgeVVE(vertex_number, tetrahedron);
		CreateFaceTFFT(tetrahedron);
		CreateVT(vertex_number, tetrahedron);
		CreateTT(tf, ft);
	}
};

/**
	mark whether a vertex is a boundary (surface) vertex
*/
class TetMeshBoundaryVertex {
public:
	std::vector<int>	boundary_vertex_flag;	///<	1: boundary vertex, 0: not boundary vertex

public:
	void MarkBoundaryVertexFromSurface(
		int		vertex_number, 
		const	std::vector<int3>& surface_face) 
	{
		boundary_vertex_flag.clear();
		boundary_vertex_flag.resize(vertex_number, 0);

		for (int i = 0; i < surface_face.size(); i++) {
			boundary_vertex_flag[surface_face[i][0]] = 1;
			boundary_vertex_flag[surface_face[i][1]] = 1;
			boundary_vertex_flag[surface_face[i][2]] = 1;
		}
	}
};


/**
*/
class TetMeshBoundary :
	public TetMeshBoundaryVertex
{
public:
};



}}}	//	namespace yz::geometry::tetrahedron


#endif	//	__YZ_TET_MESH_TOPOLOGY_H__