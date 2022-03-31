/***********************************************************/
/**	\file
	\brief		Matrix created for mesh
	\details	
	\author		Yizhong Zhang
	\date		9/7/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_MATRIX_H__
#define __YZ_MESH_MATRIX_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_sparse_matrix.h"
#include "yzLib/yz_geometry/yz_create_topology.h"
#include "yzLib/yz_geometry/yz_mesh_curvature.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Create Sparse Laplacian Matrix for Mesh
*/
//	========================================


/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.
	The matrix is not symmetric created in this way.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules. 
	If the mesh contain stand alone vertices, the corresponding row
	in the matrix is just an element 0 on the diagonal

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	ve				vertex-edge connectivity
	\param	vve_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\param	mixed_area		mixed area of each vertex
	\param	edge_weight		cotangent weight of each edge
	\return					1 on successful
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCoo<T>* spa_mat,
								 const int*	vv,
								 const int*	ve,
								 const int*	vve_start,
								 int		vertex_number,
								 const T*	mixed_area,
								 const T*	edge_weight ){
	//	setup the matrix
	if( spa_mat->IsSymmetric() ){
		#ifndef	BE_QUIET
			std::cout << "error: createLaplacianMatrixForMesh, matrix cannot be symmetric" << std::endl;
		#endif
		return 0;
	}
	spa_mat->Reset();
	spa_mat->SetDimension(vertex_number, vertex_number);
	spa_mat->SetIndexing(0);

	//	setup space for matrix COO
	int elements = vve_start[vertex_number] + vertex_number;
	spa_mat->value.resize(elements);
	spa_mat->row_id.resize(elements);
	spa_mat->col_id.resize(elements);

	//	set matrix value
	for( int v_id = 0; v_id < vertex_number; v_id ++ ){
		T	prefix_coef = 0.5;
		T	weight_sum	= 0;
		for( int i = vve_start[v_id]; i < vve_start[v_id+1]; i ++ ){
			int n_id = vv[i];
			int e_id = ve[i];
			T	weight = edge_weight[e_id];

			spa_mat->value[i+v_id] = -prefix_coef*weight;
			spa_mat->row_id[i+v_id] = v_id;
			spa_mat->col_id[i+v_id] = n_id;

			weight_sum += weight;
		}
		spa_mat->value[vve_start[v_id+1]+v_id] = prefix_coef*weight_sum;
		spa_mat->row_id[vve_start[v_id+1]+v_id] = v_id;
		spa_mat->col_id[vve_start[v_id+1]+v_id] = v_id;
	}

	//	reorder the sparse matrix
	spa_mat->Reorder();

	return 1;
}

/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.
	The matrix is not symmetric created in this way.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules. 
	If the mesh contain stand alone vertices, the corresponding row
	in the matrix is just an element 0 on the diagonal

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	ve				vertex-edge connectivity
	\param	vve_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\param	mixed_area		mixed area of each vertex
	\param	edge_weight		cotangent weight of each edge
	\return					1 on successful
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCSR<T>* spa_mat,
								 const int*	vv,
								 const int*	ve,
								 const int*	vve_start,
								 int		vertex_number,
								 const T*	mixed_area,
								 const T*	edge_weight ){
	//	setup the matrix
	if( spa_mat->IsSymmetric() ){
		#ifndef	BE_QUIET
			std::cout << "error: createLaplacianMatrixForMesh, matrix cannot be symmetric" << std::endl;
		#endif
		return 0;
	}
	spa_mat->Reset();
	spa_mat->SetDimension(vertex_number, vertex_number);
	spa_mat->SetIndexing(0);

	//	setup space for matrix CSR
	int elements = vve_start[vertex_number] + vertex_number;
	spa_mat->value.resize(elements);
	spa_mat->col_id.resize(elements);

	//	set matrix value
	for( int v_id = 0; v_id < vertex_number; v_id ++ ){
		T	prefix_coef = 0.5;
		T	weight_sum	= 0;
		for( int i = vve_start[v_id]; i < vve_start[v_id+1]; i ++ ){
			int n_id = vv[i];
			int e_id = ve[i];
			T	weight = edge_weight[e_id];

			spa_mat->value[i+v_id] = -prefix_coef*weight;
			spa_mat->col_id[i+v_id] = n_id;

			weight_sum += weight;
		}
		spa_mat->value[vve_start[v_id+1]+v_id] = prefix_coef*weight_sum;
		spa_mat->col_id[vve_start[v_id+1]+v_id] = v_id;

		spa_mat->row_start[v_id] = vve_start[v_id] + v_id;
	}
	spa_mat->row_start[vertex_number] = elements;

	//	reorder the sparse matrix
	spa_mat->Reorder();

	return 1;
}

/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	ve				vertex-edge connectivity
	\param	vve_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\param	mixed_area		mixed area of each vertex
	\param	edge_weight		cotangent weight of each edge
	\return					1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCoo<T>&		spa_mat,
								 const std::vector<int>&	vv,
								 const std::vector<int>&	ve,
								 const std::vector<int>&	vve_start,
								 int						vertex_number,
								 const std::vector<T>&		mixed_area,
								 const std::vector<T>&		edge_weight ){
	if( vertex_number <= 0 )
		return 0;

	return createLaplacianMatrixForMesh( &spa_mat, 
		&vv[0], &ve[0], &vve_start[0], vertex_number, &mixed_area[0], &edge_weight[0] );
}

/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	ve				vertex-edge connectivity
	\param	vve_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\param	mixed_area		mixed area of each vertex
	\param	edge_weight		cotangent weight of each edge
	\return					1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCSR<T>&		spa_mat,
								 const std::vector<int>&	vv,
								 const std::vector<int>&	ve,
								 const std::vector<int>&	vve_start,
								 int						vertex_number,
								 const std::vector<T>&		mixed_area,
								 const std::vector<T>&		edge_weight ){
	if( vertex_number <= 0 )
		return 0;

	return createLaplacianMatrixForMesh( &spa_mat, 
		&vv[0], &ve[0], &vve_start[0], vertex_number, &mixed_area[0], &edge_weight[0] );
}

/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCoo<T>&			spa_mat,
								 const std::vector<Vec3<T>>&	vertex,
								 const std::vector<int3>&		face ){
	//	step 1, create vv ve fe
	std::vector<int2>	edge;
	std::vector<int3>	fe;
	std::vector<int>	vv;
	std::vector<int>	ve;
	std::vector<int>	vve_start;
	createEdgeFEFromFace(edge, fe, face);
	createVVEFromEdge(vv, ve, vve_start, vertex.size(), edge);

	//	step 2, calculate edge weight, mixed area via PtrCurvTriMesh
	std::vector<T>		mixed_area;
	std::vector<T>		edge_weight;
	mixed_area.resize(vertex.size());
	edge_weight.resize(edge.size());

	PtrCurvTriMesh<T> tmp_ptr_mesh;
	tmp_ptr_mesh.SetupPtr((T*)&vertex[0], vertex.size(), 
		(int*)&face[0], face.size(), 
		(int*)&edge[0], edge.size(), 
		NULL, NULL, (int*)&fe[0], 
		&vv[0], &ve[0], &vve_start[0], 
		NULL, NULL, NULL, 
		&mixed_area[0], &edge_weight[0], NULL);

	if( !tmp_ptr_mesh.CalculateEdgeWeightMixedArea() )
		return 0;

	//	step 3, calculate laplacian matrix
	return createLaplacianMatrixForMesh(spa_mat, 
		vv, ve, vve_start, vertex.size(), mixed_area, edge_weight);

}

/**
	create sparse laplacian matrix for a mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCSR<T>&			spa_mat,
								 const std::vector<Vec3<T>>&	vertex,
								 const std::vector<int3>&		face ){
	if(vertex.empty() || face.empty()){
		std::cout << "error: createLaplacianMatrixForMesh, no input mesh" << std::endl;
		return 0;
	}

	//	step 1, create vv ve fe
	std::vector<int2>	edge;
	std::vector<int3>	fe;
	std::vector<int>	vv;
	std::vector<int>	ve;
	std::vector<int>	vve_start;
	createEdgeFEFromFace(edge, fe, face);
	createVVEFromEdge(vv, ve, vve_start, vertex.size(), edge);

	//	step 2, calculate edge weight, mixed area via PtrCurvTriMesh
	std::vector<T>		mixed_area;
	std::vector<T>		edge_weight;
	mixed_area.resize(vertex.size());
	edge_weight.resize(edge.size());

	PtrCurvTriMesh<T> tmp_ptr_mesh;
	tmp_ptr_mesh.SetupPtr((T*)&vertex[0], vertex.size(), 
		(int*)&face[0], face.size(), 
		(int*)&edge[0], edge.size(), 
		NULL, NULL, (int*)&fe[0], 
		&vv[0], &ve[0], &vve_start[0], 
		NULL, NULL, NULL, 
		&mixed_area[0], &edge_weight[0], NULL);

	if( !tmp_ptr_mesh.CalculateEdgeWeightMixedArea() )
		return 0;

	//	step 3, calculate laplacian matrix
	return createLaplacianMatrixForMesh(spa_mat, 
		vv, ve, vve_start, vertex.size(), mixed_area, edge_weight);

}


/**
	create sparse laplacian matrix for a 2D mesh

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the 2D vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCoo<T>&			spa_mat,
								 const std::vector<Vec2<T>>&	vertex,
								 const std::vector<int3>&		face ){
	std::vector<Vec3<T>> tmp_vertex_3d;
	tmp_vertex_3d.resize(vertex.size());
	for(int i=0; i<tmp_vertex_3d.size(); i++){
		tmp_vertex_3d[i].x = vertex[i].x;
		tmp_vertex_3d[i].y = vertex[i].y;
	}

	return createLaplacianMatrixForMesh(spa_mat, tmp_vertex_3d, face);
}

/**
	create sparse laplacian matrix for a mesh 2D

	matrix is created according to the paper:
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the 2D vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createLaplacianMatrixForMesh(SparseMatrixCSR<T>&			spa_mat,
								 const std::vector<Vec2<T>>&	vertex,
								 const std::vector<int3>&		face ){
	std::vector<Vec3<T>> tmp_vertex_3d;
	tmp_vertex_3d.resize(vertex.size());
	for(int i=0; i<tmp_vertex_3d.size(); i++){
		tmp_vertex_3d[i].x = vertex[i].x;
		tmp_vertex_3d[i].y = vertex[i].y;
	}

	return createLaplacianMatrixForMesh(spa_mat, tmp_vertex_3d, face);
}


///@}

//	========================================
///@{
/**	@name Create Simple Sparse Laplacian Matrix for Mesh
*/
//	========================================

/**
	create simple sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules. 
	If the mesh contain stand alone vertices, the corresponding row
	in the matrix is just an element 0 on the diagonal

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCoo<T>* spa_mat,
									   const int*	vv,
									   const int*	vv_start,
									   int			vertex_number){
	//	setup the matrix
	spa_mat->Reset();
	spa_mat->SetDimension(vertex_number, vertex_number);

	//	setup space for matrix COO
	int elements = (spa_mat->IsSymmetric() ? vv_start[vertex_number]/2 : vv_start[vertex_number]) + vertex_number;
	spa_mat->value.resize(elements);
	spa_mat->row_id.resize(elements);
	spa_mat->col_id.resize(elements);

	//	set matrix value
	int curr_id = 0;
	for( int v_id = 0; v_id < vertex_number; v_id ++ ){
		for( int i = vv_start[v_id]; i < vv_start[v_id+1]; i ++ ){
			int n_id = vv[i];
			//	if not symmetric, insert anyway
			//	if is symmetric, insert only is at upper triangular
			if( !spa_mat->IsSymmetric() || (spa_mat->IsSymmetric() && n_id>v_id) ){
				spa_mat->value[curr_id] = -1;
				spa_mat->row_id[curr_id] = v_id;
				spa_mat->col_id[curr_id] = n_id;
				curr_id ++;
			}
		}
		spa_mat->value[curr_id] = vv_start[v_id+1] - vv_start[v_id];
		spa_mat->row_id[curr_id] = v_id;
		spa_mat->col_id[curr_id] = v_id;
		curr_id ++;
	}

	//	reorder the sparse matrix
	spa_mat->Reorder();

	return 1;
}

/**
	create simple sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules. 
	If the mesh contain stand alone vertices, the corresponding row
	in the matrix is just an element 0 on the diagonal

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCSR<T>* spa_mat,
									   const int*	vv,
									   const int*	vv_start,
									   int			vertex_number){
	//	setup the matrix
	spa_mat->Reset();
	spa_mat->SetDimension(vertex_number, vertex_number);

	//	setup space for matrix CSR
	int elements = (spa_mat->IsSymmetric() ? vv_start[vertex_number]/2 : vv_start[vertex_number]) + vertex_number;
	spa_mat->value.resize(elements);
	spa_mat->col_id.resize(elements);

	//	set matrix value
	int curr_id = 0;
	for( int v_id = 0; v_id < vertex_number; v_id ++ ){
		spa_mat->row_start[v_id] = curr_id;
		for( int i = vv_start[v_id]; i < vv_start[v_id+1]; i ++ ){
			int n_id = vv[i];
			//	if not symmetric, insert anyway
			//	if is symmetric, insert only is at upper triangular
			if( !spa_mat->IsSymmetric() || (spa_mat->IsSymmetric() && n_id>v_id) ){
				spa_mat->value[curr_id] = -1;
				spa_mat->col_id[curr_id] = n_id;
				curr_id ++;
			}
		}
		spa_mat->value[curr_id] = vv_start[v_id+1] - vv_start[v_id];
		spa_mat->col_id[curr_id] = v_id;
		curr_id ++;
	}
	spa_mat->row_start[vertex_number] = elements;

	//	reorder the sparse matrix
	spa_mat->Reorder();

	return 1;
}

/**
	create sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCoo<T>&		spa_mat,
									   const std::vector<int>&	vv,
									   const std::vector<int>&	vv_start,
									   int						vertex_number){
	if( vertex_number <= 0 )
		return 0;
	if (vv.size() < vertex_number || vv_start.size() < vertex_number + 1) {
		std::cout << "error: createSimpleLaplacianMatrixForMesh, vv or vv_start size don't match vertex_number" << std::endl;
		return 0;
	}

	return createSimpleLaplacianMatrixForMesh( &spa_mat, &vv[0], &vv_start[0], vertex_number);
}

/**
	create sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat			the resulting sparse matrix 
	\param	vv				vertex-vertex connectivity
	\param	vv_start		the start and end position of neighbors of each vertex
	\param	vertex_number	vertex number
	\return					1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCSR<T>&		spa_mat,
									   const std::vector<int>&	vv,
									   const std::vector<int>&	vv_start,
									   int						vertex_number){
	if( vertex_number <= 0 )
		return 0;
	if (vv.size() < vertex_number || vv_start.size() < vertex_number + 1) {
		std::cout << "error: createSimpleLaplacianMatrixForMesh, vv or vv_start size don't match vertex_number" << std::endl;
		return 0;
	}

	return createSimpleLaplacianMatrixForMesh( &spa_mat, &vv[0], &vv_start[0], vertex_number);
}

/**
	create sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCoo<T>&			spa_mat,
									   const std::vector<Vec3<T>>&	vertex,
									   const std::vector<int3>&		face ){
	//	create vv
	std::vector<int2>	edge;
	std::vector<int>	vv;
	std::vector<int>	vv_start;
	createEdgeFromFace(edge, face);
	createVVFromEdge(vv, vv_start, vertex.size(), edge);

	//	calculate laplacian matrix
	return createSimpleLaplacianMatrixForMesh(spa_mat, vv, vv_start, vertex.size());
}

/**
	create sparse laplacian matrix for a mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCSR<T>&			spa_mat,
									   const std::vector<Vec3<T>>&	vertex,
									   const std::vector<int3>&		face ){
	//	create vv
	std::vector<int2>	edge;
	std::vector<int>	vv;
	std::vector<int>	vv_start;
	createEdgeFromFace(edge, face);
	createVVFromEdge(vv, vv_start, vertex.size(), edge);

	//	calculate laplacian matrix
	return createSimpleLaplacianMatrixForMesh(spa_mat, vv, vv_start, vertex.size());
}


/**
	create sparse laplacian matrix for a 2D mesh

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the 2D vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCoo<T>&			spa_mat,
									   const std::vector<Vec2<T>>&	vertex,
									   const std::vector<int3>&		face ){
	std::vector<Vec3<T>> tmp_vertex_3d;
	tmp_vertex_3d.resize(vertex.size());
	for(int i=0; i<tmp_vertex_3d.size(); i++){
		tmp_vertex_3d[i].x = vertex[i].x;
		tmp_vertex_3d[i].y = vertex[i].y;
	}

	return createSimpleLaplacianMatrixForMesh(spa_mat, tmp_vertex_3d, face);
}

/**
	create sparse laplacian matrix for a mesh 2D

	laplacian(x) = A*x, x is a vector of parameter on each vertex. 

	Mesh connectivity must satisfy the topology calculation rules

	\param	spa_mat		the resulting sparse matrix 
	\param	vertex		the 2D vertex array
	\param	face		the face array
	\return				1 on successful, otherwise 0
*/
template<typename T>
int createSimpleLaplacianMatrixForMesh(SparseMatrixCSR<T>&			spa_mat,
									   const std::vector<Vec2<T>>&	vertex,
									   const std::vector<int3>&		face ){
	std::vector<Vec3<T>> tmp_vertex_3d;
	tmp_vertex_3d.resize(vertex.size());
	for(int i=0; i<tmp_vertex_3d.size(); i++){
		tmp_vertex_3d[i].x = vertex[i].x;
		tmp_vertex_3d[i].y = vertex[i].y;
	}

	return createSimpleLaplacianMatrixForMesh(spa_mat, tmp_vertex_3d, face);
}

/**
create simple sparse laplacian matrix for a mesh

laplacian(x) = A*x, x is a vector of parameter on each vertex.

Mesh connectivity must satisfy the topology calculation rules.
If the mesh contain stand alone vertices, the corresponding row
in the matrix is just an element 0 on the diagonal

\param	spa_mat			the resulting sparse matrix
\param	vv				vertex-vertex connectivity
\param	vve_start		the start and end position of neighbors of each vertex
\param	edge_weight		weight of each edge
\param	vertex_number	vertex number
\return					1 on successful
*/
template<typename T>
int createEdgeWeightedLaplacianMatrixForMesh(
	SparseMatrixCSR<T>* spa_mat,
	const int*			vv,
	const int*			ve,
	const int*			vve_start,
	const T*			edge_weight,
	int					vertex_number) 
{
	//	setup the matrix
	if (spa_mat->IsSymmetric()) {
#ifndef	BE_QUIET
		std::cout << "error: createEdgeWeightedLaplacianMatrixForMesh, matrix cannot be symmetric" << std::endl;
#endif
		return 0;
	}
	spa_mat->Reset();
	spa_mat->SetDimension(vertex_number, vertex_number);

	//	setup space for matrix CSR
	int elements = vve_start[vertex_number] + vertex_number;
	spa_mat->value.resize(elements);
	spa_mat->col_id.resize(elements);

	//	set matrix value
	int curr_id = 0;
	for (int v_id = 0; v_id < vertex_number; v_id++) {
		spa_mat->row_start[v_id] = curr_id;
		T w_sum = 0;
		for (int i = vve_start[v_id]; i < vve_start[v_id + 1]; i++) {
			int n_id = vv[i];
			int e_id = ve[i];
			T w = edge_weight[e_id];
			w_sum += w;

			spa_mat->value[curr_id] = -w;
			spa_mat->col_id[curr_id] = n_id;
			curr_id++;
		}
		spa_mat->value[curr_id] = w_sum;
		spa_mat->col_id[curr_id] = v_id;
		curr_id++;
	}
	spa_mat->row_start[vertex_number] = elements;

	//	reorder the sparse matrix
	spa_mat->Reorder();

	return 1;
}

/**
create sparse laplacian matrix for a mesh

laplacian(x) = A*x, x is a vector of parameter on each vertex.

Mesh connectivity must satisfy the topology calculation rules

\param	spa_mat			the resulting sparse matrix
\param	vv				vertex-vertex connectivity
\param	vve_start		the start and end position of neighbors of each vertex
\param	edge_weight		weight of each edge
\param	vertex_number	vertex number
\return					1 on successful
*/
template<typename T>
int createEdgeWeightedLaplacianMatrixForMesh(
	SparseMatrixCSR<T>&		spa_mat,
	const std::vector<int>&	vv,
	const std::vector<int>&	ve,
	const std::vector<int>&	vv_start,
	const std::vector<T>&	edge_weight,
	int						vertex_number) 
{
	if (vertex_number <= 0)
		return 0;
	if (vv.size() < vertex_number || ve.size() < vertex_number || vv_start.size() < vertex_number + 1) {
		std::cout << "error: createEdgeWeightedLaplacianMatrixForMesh, vv, ve or vv_start size don't match vertex_number" << std::endl;
		return 0;
	}
	if (edge_weight.empty()) {
		std::cout << "error: createEdgeWeightedLaplacianMatrixForMesh, edge weight empty" << std::endl;
		return 0;
	}

	return createEdgeWeightedLaplacianMatrixForMesh(&spa_mat, &vv[0], &ve[0], &vv_start[0], &edge_weight[0], vertex_number);
}

///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_MATRIX_H__