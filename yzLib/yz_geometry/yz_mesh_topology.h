/***********************************************************/
/**	\file
	\brief		Mesh Topology
	\details	Different modules are implemented so that
				we can create different modules as we need
	\author		Yizhong Zhang
	\date		6/6/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_TOPOLOGY_H__
#define __YZ_MESH_TOPOLOGY_H__

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_create_topology.h"

namespace yz{	namespace geometry{

/**
	edge of triangle mesh
*/
class PtrTriMeshEdge{
public:
	int*	edge_ptr;

public:
	PtrTriMeshEdge() : edge_ptr(NULL){}
};

/**
	edge of triangle mesh
*/
class TriMeshEdge{
public:
	std::vector<int2>	edge;

public:
	inline void CreateEdge(const std::vector<int3>& face){	//	create edge according to face of triangle mesh
		createEdgeFromFace(edge, face);
	}

	inline void Reset(){
		edge.clear();
	}
};

/**
	vertex-face connectivity
*/
class PtrTriMeshVF{
public:
	int*	vf_ptr;
	int*	vf_start_ptr;

public:
	PtrTriMeshVF() : vf_ptr(NULL), vf_start_ptr(NULL){}
};

/**
	vertex-face connectivity
*/
class TriMeshVF{
public:
	std::vector<int>	vf, vf_start;

public:
	inline void CreateVF(int vertex_number, const std::vector<int3>& face){
		createVFFromFace(vf, vf_start, vertex_number, face);
	}

	inline void Reset(){
		vf.clear();
		vf_start.clear();
	}
};

/**
	edge-face connectivity
*/
class PtrTriMeshEFFE{
public:
	int*	ef_ptr;
	int*	fe_ptr;

public:
	PtrTriMeshEFFE() : ef_ptr(NULL), fe_ptr(NULL){}
};

/**
	edge-face connectivity

	if edge i is the bondary of face j, and the vertices 
	sequence of face j is the same as edge i, then ef[i].x = j;
	if the sequence is not the same, then ef[i].y = j;

	if edge i is the bondary of face j, and is in the opposite
	position of vertex face[j].x, then fe[j].x = i; the same for 
	fe[j].y and fe[j].z;
*/
class TriMeshEFFE{
public:
	std::vector<int2>	ef;
	std::vector<int3>	fe;

public:
	inline void Reset(){
		ef.clear();
		fe.clear();
	}
};

/**
	vertex-vertex vertex-edge connectivity
*/
class PtrTriMeshVVE{
public:
	int*	vv_ptr;
	int*	ve_ptr;
	int*	vve_start_ptr;

public:
	PtrTriMeshVVE() : vv_ptr(NULL), ve_ptr(NULL), vve_start_ptr(NULL){}
};

/**
	vertex-vertex vertex-edge connectivity
*/
class TriMeshVVE{
public:
	std::vector<int>	vv, ve, vve_start;

public:
	inline void Reset(){
		vv.clear();
		ve.clear();
		vve_start.clear();
	}
};

/**
	edge vv ve ef fe
*/
class PtrTriMeshVEF : public PtrTriMeshEdge, public PtrTriMeshEFFE, public PtrTriMeshVVE{
public:
	PtrTriMeshVEF() : PtrTriMeshEdge(), PtrTriMeshEFFE(), PtrTriMeshVVE(){}
};

/**
	edge vv ve ef fe
*/
class TriMeshVEF : public TriMeshEdge, public TriMeshEFFE, public TriMeshVVE{
public:
	inline void CreateTopology(int vertex_number, const std::vector<int3>& face){
		//	create edge, ef, fe
		createEdgeEFFEFromFace(edge, ef, fe, face);
		//	create vv, ve, vve_start
		createVVEFromEdge(vv, ve, vve_start, vertex_number, edge);
	}

	inline void Reset(){
		TriMeshEdge::Reset();
		TriMeshEFFE::Reset();
		TriMeshVVE::Reset();
	}
};

/**
	topology of triangle mesh
*/
class PtrTriMeshTopology : public PtrTriMeshVEF, public PtrTriMeshVF{
public:
	PtrTriMeshTopology() : PtrTriMeshVEF(), PtrTriMeshVF(){}
};

/**
	topology of triangle mesh
*/
class TriMeshTopology : public TriMeshVEF, public TriMeshVF{
public:
	inline void CreateTopology(int vertex_number, const std::vector<int3>& face){
		//	create edge, vv, ef, fe, vv, ve
		TriMeshVEF::CreateTopology(vertex_number, face);
		//	create vf
		TriMeshVF::CreateVF(vertex_number, face);
	}

	inline void Reset(){
		TriMeshVEF::Reset();
		TriMeshVF::Reset();
	}

};

/**
	edge and edge-face
*/
class PtrTriMeshEdgeEFFE : public PtrTriMeshEdge, public PtrTriMeshEFFE{
public:
	PtrTriMeshEdgeEFFE() : PtrTriMeshEdge(), PtrTriMeshEFFE(){}
};

/**
	edge and edge-face
*/
class TriMeshEdgeEFFE : public TriMeshEdge, public TriMeshEFFE{
public:
	inline void CreateTopology(int vertex_number, const std::vector<int3>& face){
		createEdgeEFFEFromFace(edge, ef, fe, face);
	}

	inline void Reset(){
		TriMeshEdge::Reset();
		TriMeshEFFE::Reset();
	}
};

/**
	edge and vv ve
*/
class PtrTriMeshEdgeVVE : public PtrTriMeshEdge, public PtrTriMeshVVE{
public:
	PtrTriMeshEdgeVVE() : PtrTriMeshEdge(), PtrTriMeshVVE(){}
};

/**
	edge and vv ve
*/
class TriMeshEdgeVVE : public TriMeshEdge, public TriMeshVVE{
public:
	inline void CreateTopology(int vertex_number, const std::vector<int3>& face){
		TriMeshEdge::CreateEdge(face);
		createVVEFromEdge(vv, ve, vve_start, vertex_number, edge);
	}

	inline void Reset(){
		TriMeshEdge::Reset();
		TriMeshVVE::Reset();
	}
};

/**
	mark whether vertex is boundary vertex
*/
class PtrTriMeshBoundaryVertex{
public:
	int*	boundary_vertex_flag_ptr;	///<	1: boundary vertex, 0: not boundary vertex

public:
	PtrTriMeshBoundaryVertex() : boundary_vertex_flag_ptr(NULL) {}
};

/**
	mark whether vertex is boundary vertex
*/
class TriMeshBoundaryVertex{
public:
	std::vector<int>	boundary_vertex_flag;	///<	1: boundary vertex, 0: not boundary vertex

public:
	inline int MarkBoundaryVertexFromEF(int vertex_number, const std::vector<int2>& edge, const std::vector<int2>& ef){
		return markBoundaryVertexFromEF(boundary_vertex_flag, vertex_number, edge, ef);
	}

	inline int MarkBoundaryVertex(int vertex_number, const std::vector<int3>& face){
		return markBoundaryVertex(boundary_vertex_flag, vertex_number, face);
	}

};
/**
	mark whether edge is boundary edge
*/
class PtrTriMeshBoundaryEdge{
public:
	int*	boundary_edge_flag_ptr;		///<	1: boundary edge, 0: not boundary edge

public:
	PtrTriMeshBoundaryEdge() : boundary_edge_flag_ptr(NULL) {}
};

/**
	mark whether edge is boundary edge
*/
class TriMeshBoundaryEdge{
public:
	std::vector<int>	boundary_edge_flag;		///<	1: boundary edge, 0: not boundary edge

public:
	inline int MarkBoundaryEdgeFromEF(const std::vector<int2>& edge, const std::vector<int2>& ef){
		return markBoundaryEdgeFromEF(boundary_edge_flag, edge, ef);
	}
	inline int MarkBoundaryEdge(const std::vector<int2>& edge, const std::vector<int3>&face){
		return markBoundaryEdge(boundary_edge_flag, edge, face);
	}
};

/**
	mesh boundary vertices and edges
*/
class PtrTriMeshBoundary : public PtrTriMeshBoundaryVertex, public PtrTriMeshBoundaryEdge{
public:
	PtrTriMeshBoundary() : PtrTriMeshBoundaryVertex(), PtrTriMeshBoundaryEdge() {}
};

/**
	mesh boundary vertices and edges
*/
class TriMeshBoundary : public TriMeshBoundaryVertex, public TriMeshBoundaryEdge{
public:
	inline void MarkBoundaryFromEF(int vertex_number, const std::vector<int2>& edge, const std::vector<int2>& ef){
		MarkBoundaryVertexFromEF(vertex_number, edge, ef);
		MarkBoundaryEdgeFromEF(edge, ef);
	}

	inline void MarkBoundary(int vertex_number, const std::vector<int2>& edge, const std::vector<int3>& face){
		MarkBoundaryVertex(vertex_number, face);
		MarkBoundaryEdge(edge, face);
	}
};

//	========================================
///@{
/**	@name Topology Utitity Functions
*/
//	========================================

/**
	get vertex-vertex hash (represented using std::unordered_multimap)

	\param	vv_hash			output, the vertex-vertex hash table
	\param	face_ptr		triangle face array
	\param	face_number		number of triangle faces
*/
inline void createVVHashFromTriFace(
	std::unordered_multimap<int, int>&	vv_hash,
	const int*							face_ptr,
	unsigned int						face_number
) {
	if (!face_ptr)
		return;

	vv_hash.clear();
	vv_hash.reserve(face_number * 6);

	const int* f_end_ptr = face_ptr + face_number * 3;
	for (const int* f_ptr = face_ptr; f_ptr != f_end_ptr; f_ptr += 3) {
		vv_hash.insert(std::pair<int, int>(f_ptr[0], f_ptr[1]));
		vv_hash.insert(std::pair<int, int>(f_ptr[1], f_ptr[0]));
		vv_hash.insert(std::pair<int, int>(f_ptr[0], f_ptr[2]));
		vv_hash.insert(std::pair<int, int>(f_ptr[2], f_ptr[0]));
		vv_hash.insert(std::pair<int, int>(f_ptr[1], f_ptr[2]));
		vv_hash.insert(std::pair<int, int>(f_ptr[2], f_ptr[1]));
	}
}

/**
	get vertex-vertex hash (represented using std::unordered_multimap)

	\param	vv_hash			output, the vertex-vertex hash table
	\param	face			triangle face array
*/
inline void createVVHashFromTriFace(
	std::unordered_multimap<int, int>&	vv_hash,
	const std::vector<int3>&			face
) {
	if (face.empty())
		return;

	createVVHashFromTriFace(vv_hash, &face[0].x, face.size());
}

/**
	Given a face and two vertices, return the third vertex
*/
inline int getThirdVertexOfFace(int3 face, int v0, int v1){
	if( face.x == v0 ){
		if( face.y == v1 )
			return face.z;
		else
			return face.y;
	}
	else if( face.y == v0 ){
		if( face.x == v1 )
			return face.z;
		else
			return face.x;
	}
	else{
		if( face.x == v1 )
			return face.y;
		else
			return face.x;
	}
}


///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_TOPOLOGY_H__