/***********************************************************/
/**	\file
	\brief		Seperate the Mesh
	\author		Yizhong Zhang
	\date		9/22/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_SEPERATE_H__
#define __YZ_MESH_SEPERATE_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include "yzLib/yz_geometry/yz_mesh_topology.h"
#include "yzLib/yz_geometry/yz_mesh_floodfill.h"
#include "yzLib/yz_geometry/yz_mesh_repair.h"
#include "yzLib/yz_utils/yz_reorder.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Mark Mesh Loose Parts
*/
//	========================================

/**
	mark each loose part on the mesh according to vertex connectivity

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark		vertex mark to be colored, start from 0, space must be allocated
	\param	vertex_number	number of vertex
	\param	vv				vertex-vertex connectivity
	\param	vv_start		vv start index
	\return					how many loose parts on the mesh
*/
inline int markMeshVertexLoosePartsByVV(int*		vertex_mark,
										int			vertex_number,
										const int*	vv,
										const int*	vv_start ){
	int count = 0;

	memset(vertex_mark, -1, sizeof(int)*vertex_number);	//	it's safe to set to -1

	for(int v_id=0; v_id<vertex_number; v_id++){
		if( vertex_mark[v_id] == -1 ){	//	this vertex is not flooded yet
			floodfillMeshVertexByVV(vertex_mark, count, v_id, vertex_number, vv, vv_start);
			count ++;
		}
	}

	return count;
}

/**
	mark each loose part on the mesh according to vertex connectivity

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark		vertex mark to be colored, start from 0
	\param	vertex_number	number of vertex
	\param	vv				vertex-vertex connectivity
	\param	vv_start		vv start index
	\return					how many loose parts on the mesh
*/
inline int markMeshVertexLoosePartsByVV(std::vector<int>&		vertex_mark,
										int						vertex_number,
										const std::vector<int>&	vv,
										const std::vector<int>&	vv_start ){
	if( vertex_number<=0 || vertex_number != vv_start.size()-1 || vv.empty() || vv_start.empty() ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshVertexLooseParts, array size error" << std::endl;
		#endif
		return -1;
	}

	vertex_mark.resize(vertex_number);

	return markMeshVertexLoosePartsByVV((int*)&vertex_mark[0],
		vertex_number, (int*)&vv[0], (int*)&vv_start[0]);
}

/**
	mark each loose part on the mesh according to vertex connectivity

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark			vertex mark to be colored, start from 0
	\param	count_of_color		how many vertices are flooded with given color
	\param	vertex_number		number of vertex
	\param	vv					vertex-vertex connectivity
	\param	vv_start			vv start index
	\return						how many loose parts on the mesh, same as count_of_color.size()
*/
inline int markMeshVertexLoosePartsByVV(std::vector<int>&		vertex_mark,
										std::vector<int>&		count_of_color,
										int						vertex_number,
										const std::vector<int>&	vv,
										const std::vector<int>&	vv_start ){
	if( vertex_number <= 0 )
		return 0;
	if( vertex_number != vv_start.size()-1 ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshVertexLooseParts, vertex number vv_start size don't match" << std::endl;
		#endif
		return -1;
	}

	vertex_mark.clear();	//	to make all labeled -1
	vertex_mark.resize(vertex_number, -1);
	count_of_color.clear();

	int count = 0;
	for(int v_id=0; v_id<vertex_number; v_id++){
		if( vertex_mark[v_id] == -1 ){	//	this vertex is not flooded yet
			int this_color_count = floodfillMeshVertexByVV(vertex_mark, count, v_id, vertex_number, vv, vv_start);
			count ++;
			count_of_color.push_back(this_color_count);
		}
	}

	return count;
}

/**
	mark each loose part vertex on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark		vertex mark to be colored, start from 0, space must be allocated
	\param	vertex_number	number of vertices
	\param	face			face list
	\param	face_number		number of faces
	\return					how many loose parts on the mesh
*/
inline int markMeshVertexLooseParts(int*		vertex_mark,
									int			vertex_number,
									const int*	face,
									int			face_number){
	if( vertex_number == 0 || face_number == 0 )
		return 0;

	int* edge = new int[face_number*6];
	int edge_number = createEdgeFromFace(edge, face, face_number);

	int* vv			= new int[edge_number*2];
	int* vv_start	= new int[vertex_number+1];
	createVVEFromEdge(vv, NULL, vv_start, vertex_number, edge, edge_number);

	int count = markMeshVertexLoosePartsByVV(vertex_mark, vertex_number, vv, vv_start);

	delete[] vv_start;
	delete[] vv;
	delete[] edge;
	
	return count;
}

/**
	mark each loose part on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark		vertex mark to be colored, start from 0
	\param	vertex_number	number of vertex
	\param	face			face list
	\return					how many loose parts on the mesh
*/
inline int markMeshVertexLooseParts(std::vector<int>&			vertex_mark,
									int							vertex_number,
									const std::vector<int3>&	face){
	if( vertex_number == 0 || face.empty() )
		return 0;
	vertex_mark.resize(vertex_number);
	return markMeshVertexLooseParts((int*)&vertex_mark[0], vertex_number, (int*)&face[0], face.size());
}
/**
	mark each loose part face on the mesh, loose parts number is the same as vertex loose parts

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	vertex_number		number of vertices
	\param	face				face list
	\param	face_number			number of faces
	\param	vertex_loose_part	the already marked loose part vertex of the mesh
*/
inline void markMeshFaceLoosePartsByVertexLooseParts(int*		face_mark,
													 int		vertex_number,
													 const int*	face,
													 int		face_number,
													 const int*	vertex_loose_part){
	if( vertex_number == 0 || face_number == 0 )
		return;
	for(int i=0; i<face_number; i++){
		face_mark[i] = vertex_loose_part[ face[i*3] ];
	}
}

/**
	mark each loose part face on the mesh, loose parts number is the same as vertex loose parts

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	vertex_number		number of vertices
	\param	face				face list
	\param	vertex_loose_part	the already marked loose part vertex of the mesh
*/
inline void markMeshFaceLoosePartsByVertexLooseParts(std::vector<int>&			face_mark,
													 int						vertex_number,
													 const std::vector<int3>&	face,
													 const std::vector<int>&	vertex_loose_part){
	if( vertex_number == 0 || face.empty() )
		return;
	face_mark.resize( face.size() );
	markMeshFaceLoosePartsByVertexLooseParts((int*)&face_mark[0], 
		vertex_number, (int*)&face[0], face.size(), (int*)&vertex_loose_part[0]);
}

/**
	mark each loose part face on the mesh, loose parts number is the same as vertex loose parts

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	count_of_color		count of each color
	\param	vertex_number		number of vertices
	\param	face				face list
	\param	vertex_loose_part	the already marked loose part vertex of the mesh
*/
inline void markMeshFaceLoosePartsByVertexLooseParts(std::vector<int>&			face_mark,
													 std::vector<int>&			count_of_color,
													 int						vertex_number,
													 const std::vector<int3>&	face,
													 const std::vector<int>&	vertex_loose_part){
	if( vertex_number == 0 || face.empty() )
		return;
	int face_number = face.size();
	face_mark.resize( face_number );
	count_of_color.clear();

	//	set face mark to the same as adjacent vertex
	for(int i=0; i<face_number; i++){
		int color = vertex_loose_part[ face[i][0] ];
		face_mark[i] = color;

		if( color >= count_of_color.size() ){
			count_of_color.resize(color+1, 0);
		}
		count_of_color[color] ++;
	}	
}

/**
	mark each loose part face on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	vertex_number		number of vertices
	\param	face				face list
	\param	face_number			number of faces
	\return						the number of loose parts
*/
inline int markMeshFaceLooseParts(int*			face_mark,
								  int			vertex_number,
								  const int*	face,
								  int			face_number){
	if( vertex_number == 0 || face_number == 0 )
		return 0;
	int* vertex_mark = new int[vertex_number];
	int pieces = markMeshVertexLooseParts(vertex_mark, vertex_number, face, face_number);
	markMeshFaceLoosePartsByVertexLooseParts(face_mark, vertex_number, face, face_number, vertex_mark);
	delete[] vertex_mark;
	return pieces;
}

/**
	mark each loose part on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0
	\param	vertex_number	number of vertex
	\param	face			face list
	\return					how many loose parts on the mesh
*/
inline int markMeshFaceLooseParts(std::vector<int>&			face_mark,
								  int						vertex_number,
								  const std::vector<int3>&	face){
	if( vertex_number == 0 || face.empty() )
		return 0;
	face_mark.resize(face.size());
	return markMeshFaceLooseParts((int*)&face_mark[0], vertex_number, (int*)&face[0], face.size());
}
/**
	mark each loose part on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	we use pointer to vector, so if we don't need certain list, just set to NULL

	\param	vertex_mark_ptr				pointer to return vertex mark to be colored, start from 0
	\param	vertex_count_of_color_ptr	pointer to return how many vertices are flooded with given color
	\param	face_mark_ptr				pointer to return face mark to be colored, start from 0
	\param	face_count_of_color_ptr		pointer to return how many faces are flooded with given color
	\param	vertex_number				number of vertex
	\param	face						face list
*/
inline void markMeshLooseParts(std::vector<int>*		vertex_mark_ptr,
							   std::vector<int>*		vertex_count_of_color_ptr,
							   std::vector<int>*		face_mark_ptr,
							   std::vector<int>*		face_count_of_color_ptr,
							   int						vertex_number,
							   const std::vector<int3>&	face){
	if( !vertex_mark_ptr && !face_mark_ptr && !vertex_count_of_color_ptr && !face_count_of_color_ptr )
		return;		//	nothing to output

	//	create topology
	std::vector<int2>	edge;
	std::vector<int>	vv;
	std::vector<int>	vv_start;
	createEdgeFromFace(edge, face);
	createVVFromEdge(vv, vv_start, vertex_number, edge);

	//	calculate vertex loose parts
	std::vector<int>	vertex_mark;
	std::vector<int>	vertex_count_of_color;
	markMeshVertexLoosePartsByVV(vertex_mark, vertex_count_of_color, vertex_number, vv, vv_start);

	//	calculate face loose parts conditionally
	std::vector<int>	face_mark;
	std::vector<int>	face_count_of_color;
	if( face_mark_ptr || face_count_of_color_ptr )
		markMeshFaceLoosePartsByVertexLooseParts(face_mark, face_count_of_color, vertex_number, face, vertex_mark);

	//	write data
	if( vertex_mark_ptr )	vertex_mark_ptr->swap(vertex_mark);
	if( face_mark_ptr )		face_mark_ptr->swap(face_mark);
	if( vertex_count_of_color_ptr )	vertex_count_of_color_ptr->swap(vertex_count_of_color);
	if( face_count_of_color_ptr )	face_count_of_color_ptr->swap(face_count_of_color);
}
/**
	mark each loose part on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part vertex of the mesh are labeled
	sequentially starting from 0.

	\param	vertex_mark		vertex mark to be colored, start from 0
	\param	count_of_color	how many vertices are flooded with given color
	\param	vertex_number	number of vertex
	\param	face			face list
	\return					how many loose parts on the mesh
*/
inline int markMeshVertexLooseParts(std::vector<int>&			vertex_mark,
									std::vector<int>&			count_of_color,
									int							vertex_number,
									const std::vector<int3>&	face){
	markMeshLooseParts(&vertex_mark, &count_of_color, NULL, NULL, vertex_number, face);
	return count_of_color.size();
}
/**
	mark each loose part on the mesh

	if two complete mesh pieces don't have common vertex, then they are loose part.

	after calling this function, each loose part face of the mesh are labeled
	sequentially starting from 0.

	\param	face_mark		vertex mark to be colored, start from 0
	\param	count_of_color	how many faces are flooded with given color
	\param	vertex_number	number of vertex
	\param	face			face list
	\return					how many loose parts on the mesh
*/
inline int markMeshFaceLooseParts(std::vector<int>&			face_mark,
								  std::vector<int>&			count_of_color,
								  int						vertex_number,
								  const std::vector<int3>&	face){
	markMeshLooseParts(NULL, NULL, &face_mark, &count_of_color, vertex_number, face);
	return count_of_color.size();
}
///@}

//	========================================
///@{
/**	@name Mark Mesh Faces Connected Parts
*/
//	========================================

/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	face_number		number of faces
	\param	ef				edge - neighbor face
	\param	fe				face - surrounding edge
	\return					how many face connected parts on the mesh
*/
inline int markMeshFaceConnectedByEFFE(int*			face_mark,
									   int			face_number,
									   const int*	ef,
									   const int*	fe ){
	int count = 0;
	memset(face_mark, -1, sizeof(int)*face_number);	//	it's safe to set to -1

	for(int f_id=0; f_id<face_number; f_id++){
		if( face_mark[f_id] == -1 ){	//	this face is not flooded yet
			floodfillMeshFaceByEFFE(face_mark, count, f_id, face_number, ef, fe);
			count ++;
		}
	}

	return count;
}



/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	face_number		number of faces
	\param	ef				edge - neighbor face
	\param	fe				face - surrounding edge
	\return					how many face connected parts on the mesh
*/
inline int markMeshFaceConnectedByEFFE(std::vector<int>&		face_mark,
									   int						face_number,
									   const std::vector<int2>&	ef,
									   const std::vector<int3>&	fe){
	if( face_number<=0 || fe.size()!=face_number || ef.empty() ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshFaceConnectedByEFFE, array size error" << std::endl;
		#endif
		return 0;
	}
	face_mark.resize(face_number);

	return markMeshFaceConnectedByEFFE((int*)&face_mark[0], face_number, (int*)&ef[0], (int*)&fe[0]);
}

/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	count_of_color	how many faces are flooded with given color
	\param	face_number		number of faces
	\param	ef				edge - neighbor face
	\param	fe				face - surrounding edge
	\return					how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnectedByEFFE(std::vector<int>&		face_mark,
									   std::vector<int>&		count_of_color,
									   int						face_number,
									   const std::vector<int2>&	ef,
									   const std::vector<int3>&	fe){
	if( face_number<=0 || fe.size()!=face_number || ef.empty() ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshFaceConnectedByEFFE, array size error" << std::endl;
		#endif
		return -1;
	}

	face_mark.clear();	//	to make all labeled -1
	face_mark.resize(face_number, -1);
	count_of_color.clear();

	int count = 0;
	for(int f_id=0; f_id<face_number; f_id++){
		if( face_mark[f_id] == -1 ){	//	this face is not flooded yet
			int this_color_count = floodfillMeshFaceByEFFE(face_mark, count, f_id, face_number, ef, fe);
			count ++;
			count_of_color.push_back(this_color_count);
		}
	}

	return count;
}
/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	face_number		number of faces
	\param	face			face list
	\return					how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnected(int*		face_mark,
								 int		face_number,
								 const int*	face ){
	if( face_number == 0 )
		return 0;

	int* edge	= new int[face_number*6];
	int* ef		= new int[face_number*6];
	int* fe		= new int[face_number*3];
	int edge_number = createEdgeEFFEFromFace(edge, ef, fe, face, face_number);

	int count = markMeshFaceConnectedByEFFE(face_mark, face_number, ef, fe);

	delete[] fe;
	delete[] ef;
	delete[] edge;
	
	return count;
}

/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	face			face list
	\return					how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnected(std::vector<int>&			face_mark,
								 const std::vector<int3>&	face ){
	if( face.empty() )
		return 0;
	face_mark.resize(face.size());
	return markMeshFaceConnected((int*)&face_mark[0], face.size(), (int*)&face[0]);
}

/**
	mark each face connected parts on the mesh

	There is no common edges between face connected parts. Different
	from loose parts, just sharing a vertex is not treated as face connected.

	after calling this function, each face of face connected parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark		face mark to be colored, start from 0, space must be allocated
	\param	count_of_color	how many faces are flooded with given color
	\param	face			face list
	\return					how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnected(std::vector<int>&			face_mark,
								 std::vector<int>&			count_of_color,
								 const std::vector<int3>&	face ){
	if( face.empty() )
		return 0;

	std::vector<int2>	edge;
	std::vector<int2>	ef;
	std::vector<int3>	fe;
	createEdgeEFFEFromFace(edge, ef, fe, face);

	return markMeshFaceConnectedByEFFE(face_mark, count_of_color, face.size(), ef, fe);
}


///@}

//	========================================
///@{
/**	@name Mark Mesh Faces With Edge Boundary
*/
//	========================================

/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	face_number			number of faces
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many face connected parts on the mesh
*/
inline int markMeshFaceConnectedWithEdgeBoundaryByEFFE(int*			face_mark,
													   int			face_number,
													   const int*	boundary_edge_flag,
													   const int*	ef,
													   const int*	fe ){
	int count = 0;
	memset(face_mark, -1, sizeof(int)*face_number);	//	it's safe to set to -1

	for(int f_id=0; f_id<face_number; f_id++){
		if( face_mark[f_id] == -1 ){	//	this face is not flooded yet
			floodfillMeshFaceWithEdgeBoundaryByEFFE(face_mark, count, f_id, face_number, boundary_edge_flag, ef, fe);
			count ++;
		}
	}

	return count;
}

/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	face_number			number of faces
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many face connected parts on the mesh
*/
inline int markMeshFaceConnectedWithEdgeBoundaryByEFFE(std::vector<int>&		face_mark,
													   int						face_number,
													   const std::vector<int>&	boundary_edge_flag,
													   const std::vector<int2>&	ef,
													   const std::vector<int3>&	fe){
	if( face_number<=0 || boundary_edge_flag.size()!=ef.size() || fe.size()!=face_number || ef.empty() ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshFaceConnectedByEFFE, array size error" << std::endl;
		#endif
		return 0;
	}
	face_mark.resize(face_number);

	return markMeshFaceConnectedWithEdgeBoundaryByEFFE((int*)&face_mark[0], 
		face_number, (int*)&boundary_edge_flag[0], (int*)&ef[0], (int*)&fe[0]);
}

/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	count_of_color		how many faces are flooded with given color
	\param	face_number			number of faces
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\param	ef					edge - neighbor face
	\param	fe					face - surrounding edge
	\return						how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnectedWithEdgeBoundaryByEFFE(std::vector<int>&		face_mark,
													   std::vector<int>&		count_of_color,
													   int						face_number,
													   const std::vector<int>&	boundary_edge_flag,
													   const std::vector<int2>&	ef,
													   const std::vector<int3>&	fe){
	if( face_number<=0 || boundary_edge_flag.size()!=ef.size() || fe.size()!=face_number || ef.empty() ){
		#ifndef	BE_QUIET
			std::cout << "error: markMeshFaceConnectedByEFFE, array size error" << std::endl;
		#endif
		return -1;
	}

	face_mark.clear();	//	to make all labeled -1
	face_mark.resize(face_number, -1);
	count_of_color.clear();

	int count = 0;
	for(int f_id=0; f_id<face_number; f_id++){
		if( face_mark[f_id] == -1 ){	//	this face is not flooded yet
			int this_color_count = floodfillMeshFaceWithEdgeBoundaryByEFFE(face_mark, 
				count, f_id, face_number, boundary_edge_flag, ef, fe);
			count ++;
			count_of_color.push_back(this_color_count);
		}
	}

	return count;
}
/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	face_number			number of faces
	\param	face				face list
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\return						how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnectedWithEdgeBoundary(int*		face_mark,
												 int		face_number,
												 const int*	face,
												 const int*	boundary_edge_flag){
	if( face_number == 0 )
		return 0;

	int* edge	= new int[face_number*6];
	int* ef		= new int[face_number*6];
	int* fe		= new int[face_number*3];
	int edge_number = createEdgeEFFEFromFace(edge, ef, fe, face, face_number);

	int count = markMeshFaceConnectedWithEdgeBoundaryByEFFE(face_mark, face_number, boundary_edge_flag, ef, fe);

	delete[] fe;
	delete[] ef;
	delete[] edge;
	
	return count;
}

/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	face				face list
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\return						how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnectedWithEdgeBoundary(std::vector<int>&			face_mark,
												 const std::vector<int3>&	face,
												 const std::vector<int>&	boundary_edge_flag){
	if( face.empty() )
		return 0;
	face_mark.resize(face.size());
	return markMeshFaceConnectedWithEdgeBoundary((int*)&face_mark[0], 
		face.size(), (int*)&face[0], (int*)&boundary_edge_flag[0]);
}

/**
	mark each face with edge boundary on the mesh

	after calling this function, each face of face not seperated parts of the mesh 
	are labeled sequentially starting from 0.

	\param	face_mark			face mark to be colored, start from 0, space must be allocated
	\param	count_of_color		how many faces are flooded with given color
	\param	face				face list
	\param	boundary_edge_flag	1: if the edge is boundary; 0: the edge is not boundary.	\n
								boundary edge are considered boundary implicitly
	\return						how many face connected parts on the mesh, -1 on error
*/
inline int markMeshFaceConnectedWithEdgeBoundary(std::vector<int>&			face_mark,
												 std::vector<int>&			count_of_color,
												 const std::vector<int3>&	face,
												 const std::vector<int>&	boundary_edge_flag){
	if( face.empty() )
		return 0;

	std::vector<int2>	edge;
	std::vector<int2>	ef;
	std::vector<int3>	fe;
	createEdgeEFFEFromFace(edge, ef, fe, face);

	return markMeshFaceConnectedWithEdgeBoundaryByEFFE(face_mark, 
		count_of_color, face.size(), boundary_edge_flag, ef, fe);
}


///@}

//	========================================
///@{
/**	@name Mesh Seperate
*/
//	========================================

/**
	get mapping of duplicated vertices so that it is not shared by different mesh pieces.

	this function doesn't really change vertex, it only calculate the mapping
	from after to original

	When the mesh is seperated into several pieces, each piece has a unique 
	color, and stored in face_color. For the boundary of a mesh piece, it is
	possible that one vertex is shared by different pieces. By calling this
	function, we will duplicate such vertices to make each piece of mesh has 
	a unique copy of it.

	\param	mapping			map of vertices after duplication to original vertices.
							its size will be the same as vertices after this function
	\param	vertex_number	number of vertices
	\param	face			face of mesh, will be modified
	\param	face_color		color of each face
	\return					how many vertices are added. 
							if return 0, old data didn't change
*/
inline int getMappingDuplicateSharedVerticesFromFaceColor(std::vector<int>&			mapping,
														  int						vertex_number,
														  std::vector<int3>&		face,
														  const std::vector<int>&	face_color){
	int added_count = 0;

	//	create v-f
	std::vector<int>	vf;
	std::vector<int>	vf_start;
	createVFFromFace(vf, vf_start, vertex_number, face);

	//	setup mapping
	mapping.resize(vertex_number);
	for(int i=0; i<mapping.size(); i++)
		mapping[i] = i;

	//	add new vertices
	int old_vertex_number = vertex_number;
	for(int v=0; v<old_vertex_number; v++){
		std::vector<int> color;
		for(int j=vf_start[v]; j<vf_start[v+1]; j++)
			color.push_back( face_color[vf[j]] );
		if(color.size() < 2)	//	if vertex is connected to less than two faces, 
			continue;			//	no need for further check

		std::sort(color.begin(), color.end());
		if( color.front() == color.back() )	//	just one color
			continue;

		for(int i=1; i<color.size(); i++){
			if( color[i] != color[i-1] ){	//	a different color
				int new_vertex_id = mapping.size();
				mapping.push_back(mapping[v]);	//	add new vertex
				added_count ++;
				//	for all faces with this color, change the vertex index
				for(int j=vf_start[v]; j<vf_start[v+1]; j++)
					if( face_color[vf[j]] == color[i] ){
						if( face[vf[j]].x == v )
							face[vf[j]].x = new_vertex_id;
						else if( face[vf[j]].y == v )
							face[vf[j]].y = new_vertex_id;
						else
							face[vf[j]].z = new_vertex_id;
					}
			}
		}
	}

	return added_count;
}

/**
	after mapping vertices, reorder the vertices according to mapping

	this function is specifically designed for duplicate shared vertices,
	only mapping calculated by getMappingDuplicateSharedVerticesFromFaceColor()
	is allowed to be used. 

	\param	vertex		vertex of mesh, will be modified
	\param	mapping		mapping of vertices calculated by getMappingDuplicateSharedVerticesFromFaceColor()
	\return				how many vertices are added. if return 0, old data didn't change
*/
template<typename T>
int duplicateSharedVerticesFromMapping(std::vector<Vec2<T>>&	vertex,
									   const std::vector<int>&	mapping ){
	int added = mapping.size() - vertex.size();
	assert(added >= 0);
	if( added == 0 )
		return 0;

	int old_vertex_number = vertex.size();
	vertex.resize(mapping.size());
	for(int i=old_vertex_number; i<vertex.size(); i++)
		vertex[i] = vertex[ mapping[i] ];

	return added;
}

/**
	after mapping vertices, reorder the vertices according to mapping

	this function is specifically designed for duplicate shared vertices,
	only mapping calculated by getMappingDuplicateSharedVerticesFromFaceColor()
	is allowed to be used. 

	\param	vertex		vertex of mesh, will be modified
	\param	mapping		mapping of vertices calculated by getMappingDuplicateSharedVerticesFromFaceColor()
	\return				how many vertices are added. if return 0, old data didn't change
*/
template<typename T>
int duplicateSharedVerticesFromMapping(std::vector<Vec3<T>>&	vertex,
									   const std::vector<int>&	mapping ){
	int added = mapping.size() - vertex.size();
	assert(added >= 0);
	if( added == 0 )
		return 0;

	int old_vertex_number = vertex.size();
	vertex.resize(mapping.size());
	for(int i=old_vertex_number; i<vertex.size(); i++)
		vertex[i] = vertex[ mapping[i] ];

	return added;
}

/**
	duplicate vertices so that it is not shared by different mesh pieces

	When the mesh is seperated into several pieces, each piece has a unique 
	color, and stored in face_color. For the boundary of a mesh piece, it is
	possible that one vertex is shared by different pieces. By calling this
	function, we will duplicate such vertices to make each piece of mesh has 
	a unique copy of it.

	\param	vertex			vertex of mesh, will be modified
	\param	face			face of mesh, will be modified
	\param	face_color		color of each face
	\return					how many vertices are added. 
							if return 0, old data didn't change
*/
template<typename T>
int duplicateSharedVerticesFromFaceColor(std::vector<Vec2<T>>&		vertex,
										 std::vector<int3>&			face,
										 const std::vector<int>&	face_color){
	std::vector<int> mapping;
	getMappingDuplicateSharedVerticesFromFaceColor(mapping, vertex.size(), face, face_color);
	return duplicateSharedVerticesFromMapping(vertex, mapping);
}

/**
	duplicate vertices so that it is not shared by different mesh pieces

	When the mesh is seperated into several pieces, each piece has a unique 
	color, and stored in face_color. For the boundary of a mesh piece, it is
	possible that one vertex is shared by different pieces. By calling this
	function, we will duplicate such vertices to make each piece of mesh has 
	a unique copy of it.

	\param	vertex			vertex of mesh, will be modified
	\param	face			face of mesh, will be modified
	\param	face_color		color of each face
	\return					how many vertices are added. 
							if return 0, old data didn't change
*/
template<typename T>
int duplicateSharedVerticesFromFaceColor(std::vector<Vec3<T>>&		vertex,
										 std::vector<int3>&			face,
										 const std::vector<int>&	face_color){
	std::vector<int> mapping;
	getMappingDuplicateSharedVerticesFromFaceColor(mapping, vertex.size(), face, face_color);
	return duplicateSharedVerticesFromMapping(vertex, mapping);
}

/**
	Seperate loose part of the mesh

	Loose part is identified by vertex-vertex connectivity, that means if singular vertex
	exist, the two parts are considered to be one mesh. If isolated vertex exist, it will be
	removed and doesn't be treated as one mesh segment.

	After calling this function, the order of vertex and face of original mesh will be changed.	

	\param	loose_part_vertex_start		return start index of each start vertx
	\param	loose_part_face_start		return start index of each 
	\param	vertex						vertex of the original mesh, will be changed
	\param	face						face of the original mesh, will be changed
*/
template<typename T>
int seperateMeshLooseParts(std::vector<int>&		loose_part_vertex_start,
						   std::vector<int>&		loose_part_face_start,
						   std::vector<Vec3<T>>&	vertex,
						   std::vector<int3>&		face){
	//	1, create vv
	std::vector<int2>	edge;
	std::vector<int>	vv;
	std::vector<int>	vv_start;
	createEdgeFromFace(edge, face);
	createVVFromEdge(vv, vv_start, vertex.size(), edge);

	//	2, floodfill each loose parts by v-v connectivity
	int class_count = 0;
	std::vector<int> floodfill_mark;
	floodfill_mark.resize(vertex.size(), -1);
	loose_part_vertex_start.clear();
	loose_part_vertex_start.push_back(0);
	for(int v_id=0; v_id<vertex.size(); v_id++){
		if( floodfill_mark[v_id] == -1 ){	//	this vertex is not flooded yet
			int this_color_count = floodfillMeshVertexByVV(floodfill_mark, class_count, v_id, vertex.size(), vv, vv_start);
			if( this_color_count == 1 ){	//	this vertex is isolated, so it should be droped in the final mesh
				floodfill_mark[v_id] = -1;
			}
			else if( this_color_count == 2 ){	//	this should never happen, maybe bug exist in createVVFromEdge()
				std::cout << "error: int seperateMeshLooseParts, a segment with 2 vertices detected" << std::endl;
				return 0;
			}
			else{
				loose_part_vertex_start.push_back( this_color_count );
				class_count ++;
			}
		}
	}

	//	3, calculate vertex loose part start index
	for(int i=1; i<=class_count; i++)
		loose_part_vertex_start[i] += loose_part_vertex_start[i-1];

	//	4, calculate face loose part and reorder index
	std::vector<int> counter;
	counter.resize(class_count);
	std::vector<int> reorder_index;
	reorder_index.resize(face.size());
	for(int f=0; f<face.size(); f++){
		int cid = floodfill_mark[ face[f][0] ];
		counter[cid] ++;
	}
	loose_part_face_start.clear();
	loose_part_face_start.resize(class_count+1, 0);
	for(int i=1; i<=class_count; i++){
		loose_part_face_start[i] = loose_part_face_start[i-1] + counter[i-1];
		counter[i-1] = 0;
	}
	for(int f=0; f<face.size(); f++){
		int cid = floodfill_mark[ face[f][0] ];
		reorder_index[f] = loose_part_face_start[cid] + counter[cid];
		counter[cid] ++;
	}
	utils::reorderDataToTargetOrder(face, reorder_index);

	//	5, reorder the vertices
	reorder_index.resize(vertex.size());
	for(int i=0; i<class_count; i++)
		counter[i] = 0;
	for(int i=0; i<floodfill_mark.size(); i++){
		int cid = floodfill_mark[i];
		if(cid == -1)	//	this vertex should be removed
			reorder_index[i] = -1;
		else{
			reorder_index[i] = loose_part_vertex_start[cid] + counter[cid];
			counter[cid] ++;
		}
	}

	//	finally, reorder the vertex list
	reorderVertices(vertex, face, reorder_index);

	return class_count;
}


///@}


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_SEPERATE_H__