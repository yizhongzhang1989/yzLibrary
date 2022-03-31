/***********************************************************/
/**	\file
	\brief		create topology of a given structure
	\author		Yizhong Zhang
	\date		6/3/2012
*/
/***********************************************************/
#ifndef __YZ_CREATE_TOPOLOGY_H__
#define __YZ_CREATE_TOPOLOGY_H__

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "yzLib/yz_math/yz_vector.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name compare functions used in sort
*/
//	========================================
inline bool ctIsSmallerInt2(int2 v1, int2 v2){
	return v1.x<v2.x ? true : ( v1.x>v2.x ? false : (v1.y<v2.y) );
}
inline bool ctIsSmallerInt2x(int2 v1, int2 v2){
	return v1.x < v2.x;
}
inline bool ctIsSmallerInt3x(int3 v1, int3 v2){
	return v1.x < v2.x;
}
inline bool ctIsSmallerInt4xy(int4 v1, int4 v2){
	return v1.x<v2.x ? true : ( v1.x>v2.x ? false : (v1.y<v2.y) );
}
///@}

//	========================================
///@{
/**	@name Create Connectivity Information of Mesh
*/
//	========================================

/**
	create edge from given face list

	\param	edge			array to hold result, arranged in v0v1_v0v1_, minimal safe size = sizeof(int)*(face_number*6)
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\param	face_number		number of faces
	\return					number of edges created
*/
inline int createEdgeFromFace(int*			edge, 
							  const int*	face, 
							  int			face_number){
	//	memory check
	if( edge==NULL || face==NULL || face_number==0 )
		return 0;

	std::vector<int2> tmp_edge;

	//	insert all potential edge to tmp_edge
	for( int i=0; i<face_number; i++ ){
		int2 e[3] = {int2(face[i*3], face[i*3+1]), int2(face[i*3+1], face[i*3+2]), int2(face[i*3+2], face[i*3])};
		for( int j=0; j<3; j++ ){
			if( e[j].x > e[j].y )
				e[j] = int2(e[j].y, e[j].x);
			tmp_edge.push_back( e[j] );
		}
	}

	//	sort edge
	std::sort(tmp_edge.begin(), tmp_edge.end(), ctIsSmallerInt2);

	//	add non-duplicate to edge
	int edge_number = 0;
	if( !tmp_edge.empty() ){
		edge[edge_number*2  ] = tmp_edge[0].x;
		edge[edge_number*2+1] = tmp_edge[0].y;
		edge_number ++;
	}
	for( int i=1; i<tmp_edge.size(); i++ )
		if( (tmp_edge[i].x != tmp_edge[i-1].x) || (tmp_edge[i].y != tmp_edge[i-1].y) ){
			edge[edge_number*2  ] = tmp_edge[i].x;
			edge[edge_number*2+1] = tmp_edge[i].y;
			edge_number ++;
		}

	return edge_number;
}

/**
	create edge from given face list using vector

	\param	edge	array to hold result
	\param	face	list of input face
	\return			number of edges created
*/
inline int createEdgeFromFace(std::vector<int2>&		edge, 
							  const std::vector<int3>&	face){
	if( face.empty() )
		return 0;
	edge.resize(face.size()*3);

	int edge_number = createEdgeFromFace((int*)&edge[0], (int*)&face[0], face.size());

	edge.resize(edge_number);

	return edge_number;
}

/**
	create edge, edge-face, face-edge connectivity

	If we don't want certain list, just set the input pointer to NULL, input must have value

	\param	edge			array to hold result, arranged in v0v1_v0v1_, size unknown, minimal safe size = sizeof(int)*(face_number*6)
	\param	ef				array of ef connectivity, size equal to edge size, minimal safe size = sizeof(int)*(face_number*6)
	\param	fe				array of fe connectivity, size = sizeof(int)*(face_number*3)
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\param	face_number		number of faces
	\return					number of edges created
*/
inline int createEdgeEFFEFromFace(int*			edge, 
								  int*			ef, 
								  int*			fe, 
								  const int*	face, 
								  int			face_number ){
	if( face==NULL || face_number==0 )			//	don't have input
		return 0;
	if( edge==NULL )				//	cannot miss edge
		return -1;
	else if( ef==NULL && fe==NULL )	//	don't have ef and fe output, simply create edge
		return createEdgeFromFace(edge, face, face_number);

	//	create edge
	std::vector<int4> tmp_edge;	//	vertex 1, vertex 2, face index, edge index
	tmp_edge.reserve(face_number * 3);

	for( int i=0; i<face_number; i++ ){
		int4 e[3] = {int4(face[i*3], face[i*3+1], i, -1), int4(face[i*3+1], face[i*3+2], i, -1), int4(face[i*3+2], face[i*3], i, -1)};
		for( int j=0; j<3; j++ ){
			if( e[j].x > e[j].y )
				e[j] = int4(e[j].y, e[j].x, e[j].z, e[j].w);
			tmp_edge.push_back( e[j] );
		}
	}

	std::sort(tmp_edge.begin(), tmp_edge.end(), ctIsSmallerInt4xy);

	int edge_number = 0;
	if( !tmp_edge.empty() ){
		edge[edge_number*2  ] = tmp_edge[0].x;
		edge[edge_number*2+1] = tmp_edge[0].y;
		tmp_edge[0].w = 0;
		edge_number ++;
	}
	for( int i=1; i<tmp_edge.size(); i++ ){
		if( (tmp_edge[i].x != tmp_edge[i-1].x) || (tmp_edge[i].y != tmp_edge[i-1].y) ){
			edge[edge_number*2  ] = tmp_edge[i].x;
			edge[edge_number*2+1] = tmp_edge[i].y;
			edge_number ++;
		}
		tmp_edge[i].w = edge_number-1;
	}

	//	create ef, fe
	if(ef)	std::fill(&ef[0], &ef[edge_number*2], -1);//memset(&ef[0], -1, sizeof(int)*edge_number*2);
	if(fe)	std::fill(&fe[0], &fe[face_number*3], -1);//memset(&fe[0], -1, sizeof(int)*face_number*3);
	for( int i=0; i<tmp_edge.size(); i++ ){
		int e = tmp_edge[i].w;
		int f = tmp_edge[i].z;
		int v[2] = {edge[e*2], edge[e*2+1]};

		//	set ef
		if( ef ){
			if(		(face[f*3]==v[0] && face[f*3+1]==v[1]) ||	
					(face[f*3+1]==v[0] && face[f*3+2]==v[1]) || 
					(face[f*3+2]==v[0] && face[f*3]==v[1]) ){
				ef[e*2] = f;	//	if the direction of edge is the same as face, put it in front
			}
			else
				ef[e*2+1] = f;	//	or, put it back
		}

		//	set fe
		if( fe ){
			int j = 0;
			for(; j<3; j++ ){
				if( face[f*3+j]!=v[0] && face[f*3+j]!=v[1] ){
					fe[f*3+j] = e;	//	put the edge in the opposite position of the vertex
					break;
				}
			}
		}
	}

	return edge_number;
}

/**
	create edge, edge-face, face-edge connectivity using vector

	\param	edge			array to hold result
	\param	ef				array of ef connectivity, 
							ef[i].x = j means edge i is a boundary of face j, 
							and their vertices sequence are the same; 
							ef[i].y = j means edge i is a boundary of face j, 
							but their vertices sequence are the opposite; 
	\param	fe				array of fe connectivity,
							fe[i].x = j means edge j is a boundary of face i,
							and x vertex of face i is opposite to edge j
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\return					number of edges created
*/
inline int createEdgeEFFEFromFace(std::vector<int2>&		edge, 
								  std::vector<int2>&		ef, 
								  std::vector<int3>&		fe, 
								  const std::vector<int3>&	face){
	if( face.empty() )
		return 0;

	edge.resize(face.size()*3);
	ef.resize(face.size()*3);
	fe.resize(face.size());

	int edge_number = createEdgeEFFEFromFace((int*)&edge[0], (int*)&ef[0], (int*)&fe[0], (int*)&face[0], face.size());

	edge.resize(edge_number);
	ef.resize(edge_number);

	return edge_number;
}

/**
	create edge, edge-face, connectivity using vector

	\param	edge			array to hold result
	\param	ef				array of ef connectivity, 
							ef[i].x = j means edge i is a boundary of face j, 
							and their vertices sequence are the same; 
							ef[i].y = j means edge i is a boundary of face j, 
							but their vertices sequence are the opposite; 
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\return					number of edges created
*/
inline int createEdgeEFFromFace(std::vector<int2>&			edge, 
								std::vector<int2>&			ef, 
								const std::vector<int3>&	face){
	if( face.empty() )
		return 0;
	edge.resize(face.size()*3);
	ef.resize(face.size()*3);

	int edge_number = createEdgeEFFEFromFace((int*)&edge[0], (int*)&ef[0], NULL, (int*)&face[0], face.size());

	edge.resize(edge_number);
	ef.resize(edge_number);

	return edge_number;
}

/**
	create edge, face-edge connectivity using vector

	\param	edge			array to hold result
	\param	fe				array of fe connectivity,
							fe[i].x = j means edge j is a boundary of face i,
							and x vertex of face i is opposite to edge j
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\return					number of edges created
*/
inline int createEdgeFEFromFace(std::vector<int2>&			edge, 
								std::vector<int3>&			fe, 
								const std::vector<int3>&	face){
	if( face.empty() )
		return 0;

	edge.resize(face.size()*3);
	fe.resize(face.size());

	int edge_number = createEdgeEFFEFromFace((int*)&edge[0], NULL, (int*)&fe[0], (int*)&face[0], face.size());

	edge.resize(edge_number);

	return edge_number;
}

/**
	create vertex - face connectivity
	\param	vf				list to hold connectivity information, size = sizeof(int)*(face_number*3)
	\param	vf_start		list to hold vf start information, size = sizeof(int)*(vertex_number+1)
	\param	vertex_number	vertex number
	\param	face			list of input face, arraged in v0v1v2_v0v1v2_
	\param	face_number		number of faces
*/
inline void createVFFromFace(int*		vf, 
							 int*		vf_start, 
							 int		vertex_number, 
							 const int* face, 
							 int		face_number){
	if( vertex_number==0 || face==NULL || face_number==0 )	//	no input
		return;
	if( vf==NULL || vf_start==NULL )						//	no output
		return;
									
	std::vector<int2> tmp_vf;
	for( int i=0; i<face_number; i++ ){
		tmp_vf.push_back( int2(face[i*3  ], i) );
		tmp_vf.push_back( int2(face[i*3+1], i) );
		tmp_vf.push_back( int2(face[i*3+2], i) );
	}

	std::sort(tmp_vf.begin(), tmp_vf.end(), ctIsSmallerInt2x);

	for( int i=0; i<tmp_vf.size(); i++ )
		vf[i] = tmp_vf[i].y;

	vf_start[0] = 0;
	for( int i=0, j=0; i<vertex_number; i++ ){
		while(j<tmp_vf.size() && tmp_vf[j].x==i)
			j++;
		vf_start[i+1] = j;
	}
}

/**
	create vertex - face connectivity using vector
	\param	vf				list to hold connectivity information
	\param	vf_start		list to hold vf start information
	\param	vertex_number	vertex number
	\param	face			list of input face
*/
inline void createVFFromFace(std::vector<int>&			vf, 
							 std::vector<int>&			vf_start, 
							 int						vertex_number, 
							 const std::vector<int3>&	face){
	if( vertex_number<=0 || face.empty() )
		return;

	vf.resize( face.size()*3 );
	vf_start.resize( vertex_number+1 );

	createVFFromFace((int*)&vf[0], (int*)&vf_start[0], vertex_number, (int*)&face[0], face.size());
}

/**
	create vertex-vertex, vertex-edge connectivity

	If we don't want certain list, just set the input pointer to NULL, input must have value

	\param	vv				list to hold vertex-vertex connectivity, size = sizeof(int)*(edge_number*2)
	\param	ve				list to hold vertex-edge connectivity, size = sizeof(int)*(edge_number*2)
	\param	vve_start		list to hold vv,ve start information, size = sizeof(int)*(vertex_number+1)
	\param	vertex_number	vertex numbe
	\param	edge			list of input edge
	\param	edge_number		edge number
*/
inline void createVVEFromEdge(int*			vv, 
							  int*			ve, 
							  int*			vve_start, 
							  int			vertex_number, 
							  const int*	edge, 
							  int			edge_number ){
	if( vertex_number==0 || edge==NULL || edge_number==0 )	//	no input
		return;
	if( vv==NULL && ve==NULL && vve_start==NULL )			//	no output
		return;

	std::vector<int3> tmp_vve;
	for( int i=0; i<edge_number; i++ ){
		tmp_vve.push_back( int3(edge[i*2  ], edge[i*2+1], i) );
		tmp_vve.push_back( int3(edge[i*2+1], edge[i*2  ], i) );
	}

	std::sort(tmp_vve.begin(), tmp_vve.end(), ctIsSmallerInt3x);

	for( int i=0; i<tmp_vve.size(); i++ ){
		if(vv)	vv[i] = tmp_vve[i].y;
		if(ve)	ve[i] = tmp_vve[i].z;
	}

	if(vve_start){
		vve_start[0] = 0;
		for( int i=0, j=0; i<vertex_number; i++ ){
			while(j<tmp_vve.size() && tmp_vve[j].x==i)
				j++;
			vve_start[i+1] = j;
		}
	}
}

/**
	create vertex-vertex, vertex-edge connectivity using vector

	\param	vv				list to hold vertex-vertex connectivity
	\param	ve				list to hold vertex-edge connectivity
	\param	vve_start		list to hold vv,ve start information
	\param	vertex_number	vertex number
	\param	edge			list of input edge
*/
inline void createVVEFromEdge(std::vector<int>&			vv, 
							  std::vector<int>&			ve, 
							  std::vector<int>&			vve_start, 
							  int						vertex_number, 
							  const std::vector<int2>&	edge){
	if( vertex_number<=0|| edge.empty() )
		return;
	
	vv.resize( edge.size()*2 );
	ve.resize( edge.size()*2 );
	vve_start.resize( vertex_number+1 );
	createVVEFromEdge((int*)&vv[0], (int*)&ve[0], (int*)&vve_start[0], vertex_number, (int*)&edge[0], edge.size());
}

/**
	create vertex-vertex, connectivity using vector

	\param	vv				list to hold vertex-vertex connectivity
	\param	vv_start		list to hold vv start information
	\param	vertex_number	vertex number
	\param	edge			list of input edge
*/
inline void createVVFromEdge(std::vector<int>&			vv, 
							 std::vector<int>&			vv_start, 
							 int						vertex_number, 
							 const std::vector<int2>&	edge){
	if( vertex_number<=0 || edge.empty() )
		return;
	vv.resize( edge.size()*2 );
	vv_start.resize( vertex_number+1 );
	createVVEFromEdge((int*)&vv[0], NULL, (int*)&vv_start[0], vertex_number, (int*)&edge[0], edge.size());
}

/**
	create vertex-edge connectivity using vector

	\param	ve				list to hold vertex-edge connectivity
	\param	ve_start		list to hold ve start information
	\param	vertex_number	vertex number
	\param	edge			list of input edge
*/
inline void createVEFromEdge(std::vector<int>&			ve, 
							 std::vector<int>&			ve_start, 
							 int						vertex_number, 
							 const std::vector<int2>&	edge){
	if( vertex_number<=0 || edge.empty() )
		return;
	ve.resize( edge.size()*2 );
	ve_start.resize( vertex_number+1 );
	createVVEFromEdge(NULL, (int*)&ve[0], (int*)&ve_start[0], vertex_number, (int*)&edge[0], edge.size());
}

/**
	Create face-face (via edge) connectivity.

	In this function, we first create a hash table for edge-face lookup.
	Then find all neighbor faces of each face. Non-manifold is also appliable.

	\param	ff				face-face connectivity
	\param	vf_start		start of each face in ff
	\param	face			input faces
*/
inline void createFFFromFace(
	std::vector<int>&			ff,
	std::vector<int>&			ff_start,
	const std::vector<int3>&	face
) {
	ff.clear();
	ff.reserve(face.size() * 3);
	ff_start.clear();
	ff_start.reserve(face.size() + 1);
	ff_start.push_back(0);

	if (face.empty())
		return;

	//	first parse, create edge-face hash
	std::unordered_multimap<int2, int, utils::BitwiseHasher<int2>> ef_hash;
	ef_hash.reserve(face.size() * 3);
	for (unsigned int f = 0; f != face.size(); f++) {	//	for each face
		for (int i = 0; i < 3; i++) {					//	for each edge of the face
			int2 edge(face[f][i], face[f][(i + 1) % 3]);
			if (edge.x > edge.y)
				mySwap(edge.x, edge.y);

			//	insert into hash table
			ef_hash.insert(std::pair<int2, int>(edge, f));
		}
	}

	//	second parse, extract face-face
	for (unsigned int f = 0; f != face.size(); f++) {	//	for each face
		for (int i = 0; i < 3; i++) {					//	for each edge of the face
			int2 edge(face[f][i], face[f][(i + 1) % 3]);
			if (edge.x > edge.y)
				mySwap(edge.x, edge.y);

			//	for each neighbor of this edge
			auto range = ef_hash.equal_range(edge);
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

///@}

//	========================================
///@{
/**	@name Mark Boundary on Mesh
*/
//	========================================

/**
	mark boundary vertex of mesh, if the vertex is on mesh boundary, set it to 1, otherwise, set to 0

	isolated vertices are not marked in this function 

	\param	boundary_vertex_flag	return array, 1 is boundary;		0 is not boundary
	\param	vertex_number			vertex number
	\param	edge					edge list
	\param	edge_number				edge number
	\param	ef						edge-face connectivity
	\return							number of boundary vertices
*/
inline int markBoundaryVertexFromEF(int*		boundary_vertex_flag, 
									int			vertex_number, 
									const int*	edge, 
									int			edge_number, 
									const int*	ef){
	std::fill(&boundary_vertex_flag[0], &boundary_vertex_flag[vertex_number], 0); //memset(boundary_vertex_flag, 0, sizeof(int)*vertex_number);
	int boundary_vertex_count = 0;
	for( int e=0; e<edge_number; e++ ){
		if( ef[e*2]==-1 || ef[e*2+1]==-1 ){
			if( boundary_vertex_flag[edge[e*2  ]] == 0 ){
				boundary_vertex_flag[edge[e*2  ]] = 1;
				boundary_vertex_count ++;
			}
			if( boundary_vertex_flag[edge[e*2+1]] == 0 ){
				boundary_vertex_flag[edge[e*2+1]] = 1;
				boundary_vertex_count ++;
			}
		}
	}
	return boundary_vertex_count;
}

/**
	mark boundary vertex of mesh using vector

	if the vertex is on mesh boundary, set it to 1, otherwise, set to 0

	isolated vertices are not marked in this function 

	\param	boundary_vertex_flag	return array, 1 is boundary;		0 is not boundary
	\param	vertex_number			vertex number
	\param	edge					edge list
	\param	ef						edge-face connectivity
	\return							number of boundary vertices
*/
inline int markBoundaryVertexFromEF(std::vector<int>&			boundary_vertex_flag,
									int							vertex_number,
									const std::vector<int2>&	edge,
									const std::vector<int2>&	ef	){
	if( vertex_number<=0 || edge.empty() || ef.empty() )
		return 0;

	boundary_vertex_flag.resize(vertex_number);
	int boundary_vertex_count = markBoundaryVertexFromEF((int*)&boundary_vertex_flag[0], vertex_number, (int*)&edge[0], edge.size(), (int*)&ef[0]);

	if( !boundary_vertex_count )	boundary_vertex_flag.clear();	//	if the mesh don't have boundary, we delete the space
	return boundary_vertex_count;
}

/**
	mark boundary vertex of mesh, if the vertex is on mesh boundary, set it to 1, otherwise, set to 0

	isolated vertices are not marked in this function 

	\param	boundary_vertex_flag	return array, 1 is boundary;		0 is not boundary
	\param	vertex_number			vertex number
	\param	face					face list arranged in v0v1v2_v0v1v2_ 
	\param	face_number				number of faces
	\return							number of boundary vertices
*/
inline int markBoundaryVertex(int*			boundary_vertex_flag, 
							  int			vertex_number, 
							  const int*	face,
							  int			face_number){
	if( vertex_number<=0 || !face_number )
		return 0;

	int* edge	= new int[face_number*6];
	int* ef		= new int[face_number*6];
	int edge_number = createEdgeEFFEFromFace(edge, ef, NULL, face, face_number);

	int boundary_vertex_number = markBoundaryVertexFromEF(boundary_vertex_flag,
		vertex_number, edge, edge_number, ef);

	delete[] edge;
	delete[] ef;
	return boundary_vertex_number;
}

/**
	mark boundary vertex of mesh using vector

	if the vertex is on mesh boundary, set it to 1, otherwise, set to 0

	isolated vertices are not marked in this function 

	\param	boundary_vertex_flag	return array, 1 is boundary;		0 is not boundary
	\param	vertex_number			vertex number
	\param	face					list of input face, arraged in v0v1v2_v0v1v2_
	\return							number of boundary vertices
*/
inline int markBoundaryVertex(std::vector<int>&			boundary_vertex_flag,
							  int						vertex_number,
							  const std::vector<int3>&	face ){
	if(vertex_number<=0 || face.empty() )
		return 0;
	boundary_vertex_flag.resize(vertex_number);
	return markBoundaryVertex(&boundary_vertex_flag[0], vertex_number, (int*)&face[0], face.size());
}

/**
	mark the edge whether the edge is on the boundary of a mesh

	edge, ef must be created using yzLib from face, or the calculation
	may not be correct

	\param	boundary_edge_flag	return whether on boundary, 1: on boundary; 0: not on boundary
	\param	edge				edge list of the mesh, created from face
	\param	edge_number			number of edges
	\param	ef					edge-face connectivity, created from face and the above edge
	\return						how many boundary edges are found
*/
inline int markBoundaryEdgeFromEF(int*			boundary_edge_flag, 
								  const int*	edge, 
								  int			edge_number, 
								  const int*	ef){
	std::fill(&boundary_edge_flag[0], &boundary_edge_flag[edge_number], 0);//memset(boundary_edge_flag, 0, sizeof(int)*edge_number);
	int boundary_edge_count = 0;
	for( int e=0; e<edge_number; e++ ){
		if( ef[e*2]==-1 || ef[e*2+1]==-1 ){
			boundary_edge_flag[e] = 1;
			boundary_edge_count ++;
		}
	}
	return boundary_edge_count;
}

/**
	mark the edge whether the edge is on the boundary of a mesh

	edge, ef must be created using yzLib from face, or the calculation
	may not be correct

	\param	boundary_edge_flag	return whether on boundary, 1: on boundary; 0: not on boundary
	\param	edge				edge list of the mesh, created from face
	\param	ef					edge-face connectivity, created from face and the above edge
	\return						how many boundary edges are found
*/
inline int markBoundaryEdgeFromEF(std::vector<int>&			boundary_edge_flag,
								  const std::vector<int2>&	edge,
								  const std::vector<int2>&	ef	){
	if(edge.empty() || ef.empty())
		return 0;
	boundary_edge_flag.resize(edge.size());
	return markBoundaryEdgeFromEF(&boundary_edge_flag[0], (int*)&edge[0], edge.size(), (int*)&ef[0]);
}

/**
	mark the edge whether the edge is on the boundary of a mesh

	edge, ef must be created using yzLib from face, or the calculation
	may not be correct

	\param	boundary_edge_flag	return whether on boundary, 1: on boundary; 0: not on boundary
	\param	edge				edge list of the mesh, created from face
	\param	edge_number			number of edges
	\param	face				face list of the mesh
	\param	face_number			number of faces
	\return						how many boundary edges are found
*/
inline int markBoundaryEdge(int*		boundary_edge_flag, 
							const int*	edge, 
							int			edge_number, 
							const int*	face,
							int			face_number){
	int* edge2	= new int[face_number*6];	//	edge2 should be identical to edge
	int* ef		= new int[face_number*6];
	int edge_number2 = createEdgeEFFEFromFace(edge2, ef, NULL, face, face_number);
	if( edge_number != edge_number2 ){
		#ifndef	BE_QUIET
			std::cout << "error: markBoundaryEdge, edge of input is not created by face" << std::endl;
		#endif
		return 0;
	}

	return markBoundaryEdgeFromEF(boundary_edge_flag, edge2, edge_number, ef);
}

/**
	mark the edge whether the edge is on the boundary of a mesh

	edge, ef must be created using yzLib from face, or the calculation
	may not be correct

	\param	boundary_edge_flag	return whether on boundary, 1: on boundary; 0: not on boundary
	\param	edge				edge list of the mesh, created from face
	\param	face				face of the mesh
	\return						how many boundary edges are found
*/
inline int markBoundaryEdge(std::vector<int>&			boundary_edge_flag,
							const std::vector<int2>&	edge,
							const std::vector<int3>&	face){
	if(edge.empty() || face.empty())
		return 0;
	boundary_edge_flag.resize(edge.size());
	return markBoundaryEdge(&boundary_edge_flag[0], (int*)&edge[0], edge.size(), (int*)&face[0], face.size());
}
///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_CREATE_TOPOLOGY_H__