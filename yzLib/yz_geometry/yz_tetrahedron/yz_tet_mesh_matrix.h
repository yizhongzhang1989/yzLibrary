/***********************************************************/
/**	\file
	\brief		Matrix for Tetrahedral Mesh
	\details	Matrix used for tetrahedral mesh calculation
	\author		Yizhong Zhang
	\date		12/8/2016
*/
/***********************************************************/
#ifndef __YZ_TET_MESH_MATRIX_H__
#define __YZ_TET_MESH_MATRIX_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_geometry/yz_mesh_matrix.h"

namespace yz {  namespace geometry {  namespace tetrahedron {

/**
	create simple sparse laplacian matrix for a tetrahedral mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex.

	Mesh connectivity must satisfy the topology calculation rules.
	If the mesh contain stand alone vertices, the corresponding row
	in the matrix is just an element 0 on the diagonal

	In fact, simple laplacian for tet_mesh is the same as tri_mesh

	\param	spa_mat			the resulting sparse matrix
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful
*/
template<typename T>
int createSimpleLaplacianMatrixForTetMesh(
	SparseMatrixCSR<T>* spa_mat,
	const int*			vv,
	const int*			vv_start,
	int					vertex_number) 
{
	return(createSimpleLaplacianMatrixForMesh(spa_mat, vv, vv_start, vertex_number));
}

/**
	create sparse laplacian matrix for a tet mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex.

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat			the resulting sparse matrix
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForTetMesh(
	SparseMatrixCSR<T>&		spa_mat,
	const std::vector<int>&	vv,
	const std::vector<int>&	vv_start,
	int						vertex_number) 
{
	if (vertex_number <= 0 )
		return 0;
	if (vv.size() < vertex_number || vv_start.size() < vertex_number + 1) {
		std::cout << "error: createSimpleLaplacianMatrixForTetMesh, vv or vv_start size don't match vertex_number" << std::endl;
		return 0;
	}
	return createSimpleLaplacianMatrixForTetMesh(&spa_mat, &vv[0], &vv_start[0], vertex_number);
}


}}}	//	namespace yz::geometry::tetrahedron


#endif	//	__YZ_TET_MESH_MATRIX_H__