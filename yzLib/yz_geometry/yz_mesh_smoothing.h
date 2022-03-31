/***********************************************************/
/**	\file
	\brief		Some Smooth Functions.
				Some functions in this file are not limited 
				to be used on mesh. 
	\author		Yizhong Zhang
	\date		6/10/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_SMOOTH_H__
#define __YZ_MESH_SMOOTH_H__

#include <assert.h>
#include <vector>
#include <unordered_map>
#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_math/yz_interpolation.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Laplacian Smooth on Mesh
*/
//	========================================
/**
	Laplacian according to edge
	
	Data to apply laplacian can be any format, if operator +=, *=double are defined.
	Pass the right type of pointer as argument.

	If you want to laplacian 3D vertex position, pass pointer as Vec3*

	\param	data_ptr		pointer to data to apply laplacian on, if the type of data is a structure,
							then the type of pointer must be type of the structure
	\param	data_number		number of data
	\param	edge_ptr		pointer to edge list
	\param	edge_number		number of edge
	\param	coef			coefficient of laplacian, 0: original data, 1: full laplacian of one ring neighbor.
							other values means linear interpolation of the 0 and 1
*/
template<typename T>
inline void setLaplacianFromEdge(T* data_ptr, int data_number, const int* edge_ptr, int edge_number, double coef=1){
	coef = clamp(coef, 0, 1);

	T* tmp_data_ptr = new T[data_number];
	int* nei_count = new int[data_number];
	std::fill(&tmp_data_ptr[0], &tmp_data_ptr[data_number], 0);//memset(tmp_data_ptr, 0, sizeof(T)*data_number);
	memset(nei_count, 0, sizeof(int)*data_number);

	for( int e=0; e<edge_number; e++ ){
		int v0 = edge_ptr[e*2];
		int v1 = edge_ptr[e*2+1];
		tmp_data_ptr[v0] += data_ptr[v1];
		tmp_data_ptr[v1] += data_ptr[v0];
		nei_count[v0] ++;
		nei_count[v1] ++;
	}

	for( int i=0; i<data_number; i++ ){
		if( nei_count[i] <= 1 )
			tmp_data_ptr[i] = data_ptr[i];
		else
			tmp_data_ptr[i] *= (1.0/nei_count[i]);
	}

	if( coef > (1 - 1e-6) )
		memcpy(data_ptr, tmp_data_ptr, sizeof(T)*data_number);
	else{
		for( int i=0; i<data_number; i++ )
			data_ptr[i] = interpLinear(data_ptr[i], tmp_data_ptr[i], coef);
	}

	delete[] nei_count;
	delete[] tmp_data_ptr;
}

/**
	Laplacian data in vectors according to edge in vector
	
	Data to apply laplacian can be any format, if operator +=, *=double are defined.
	Pass the right type of pointer as argument.

	\param	data			vector of vertices
	\param	edge			vector of edges
	\param	coef			coefficient of laplacian, 0: original data, 1: full laplacian of one ring neighbor.
							other values means linear interpolation of the 0 and 1
*/
template<typename T>
inline void setLaplacianFromEdge(std::vector<T>& data, const std::vector<int2>& edge, double coef=1){
	setLaplacianFromEdge(&data[0], data.size(), (int*)&edge[0], edge.size(), coef);
}

/**
	Laplacian according to vertex-vertex connectivity

	Data to apply laplacian can be any format, if operator +=, *=double are defined.
	Pass the right type of pointer as argument.

	If you want to laplacian 3D vertex position, pass pointer as Vec3*

	\param	data_ptr		pointer to data to apply laplacian on, if the type of data is a structure,
							then the type of pointer must be type of the structure
	\param	data_number		number of data
	\param	vv_ptr			pointer to vertex-vertex connectivity, stored in [v0 neighbors list | v1 neighbors list | ...]
	\param	vv_start_ptr	pointer to neighbor list start index of each vertex, reference to TriMeshVVE
	\param	coef			coefficient of laplacian, 0: original data, 1: full laplacian of one ring neighbor.
							other values means linear interpolation of the 0 and 1
*/
template<typename T>
inline void setLaplacianFromVV(T* data_ptr, int data_number, const int* vv_ptr, const int* vv_start_ptr, double coef=1){
	coef = clamp(coef, 0, 1);

	T* tmp_data_ptr = new T[data_number];
	memset(tmp_data_ptr, 0, sizeof(T)*data_number);

	for( int i=0; i<data_number; i++ ){
		for( int j=vv_start_ptr[i]; j<vv_start_ptr[i+1]; j++ )
			tmp_data_ptr[i] += data_ptr[vv_ptr[j]];
		int num = vv_start_ptr[i+1] - vv_start_ptr[i];
		if( num==0 )
			tmp_data_ptr[i] = data_ptr[i];
		else
			tmp_data_ptr[i] *= 1.0 / num;
	}

	if( coef > (1 - 1e-6) )
		memcpy(data_ptr, tmp_data_ptr, sizeof(T)*data_number);
	else{
		for( int i=0; i<data_number; i++ )
			data_ptr[i] = interpLinear(data_ptr[i], tmp_data_ptr[i], coef);
	}

	delete[] tmp_data_ptr;
}

/**
	Laplacian data in vectors according to vertex vertex connectivity
	
	Data to apply laplacian can be any format, if operator +=, *=double are defined.
	Pass the right type of pointer as argument.

	\param	data			vector of vertices
	\param	vv				vertex-vertex connectivity, stored in [v0 neighbors list | v1 neighbors list | ...]
	\param	vv_start		neighbor list start index of each vertex, reference to TriMeshVVE
	\param	coef			coefficient of laplacian, 0: original data, 1: full laplacian of one ring neighbor.
							other values means linear interpolation of the 0 and 1
*/
template<typename T>
inline void setLaplacianFromVV(std::vector<T>& data, const std::vector<int>& vv, const std::vector<int>& vv_start, double coef=1){
	assert(data.size() == vv_start.size()-1);
	setLaplacianFromVV(&data[0], data.size(), &vv[0], &vv_start[0], coef);
}

/**
	perform Laplacian smooth on triangle mesh

	\param	vertex_ptr		vertex arranged in xyz-xyz-xyz-
	\param	vertex_number	number of vertices
	\param	face_ptr		triangle face arranged in xyz-xyz-xyz-
	\param	face_number		number of faces
	\param	vv_hash			vertex-vertex hash table
	\param	coef			linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
*/
template <typename T>
void smoothLaplacianTriMesh(
	T*											vertex_ptr,
	unsigned int								vertex_number,
	const int*									face_ptr,
	unsigned int								face_number,
	const std::unordered_multimap<int, int>&	vv_hash,
	T											coef = 1.0
) {
	if (!vertex_ptr || !face_ptr || !vertex_number || !face_number)
		return;

	Vec3<T>* vertex = (Vec3<T>*)vertex_ptr;
	Vec3<T>* tmp_vertex = new Vec3<T>[vertex_number];

	//	perform laplacian smooth
	for (unsigned int vid = 0; vid != vertex_number; vid++) {
		//	for each vertex, calculate the average vertex position of its neighbor
		auto its = vv_hash.equal_range(vid);
		Vec3<T> v_sum(0, 0, 0);
		unsigned int count = 0;
		for (auto iter = its.first; iter != its.second; iter++) {
			v_sum += vertex[iter->second];
			count++;
		}
		if (count) {	//	normal vertex, interpolate between original and Laplacian
			v_sum /= count;
			tmp_vertex[vid] = (1.0 - coef) * vertex[vid] + coef * v_sum;
		}
		else {			//	isolated vertex, just original 
			tmp_vertex[vid] = vertex[vid];
		}
	}

	memcpy(vertex, tmp_vertex, sizeof(T)*vertex_number * 3);

	delete[] tmp_vertex;
}

/**
	perform Laplacian smooth on triangle mesh

	\param	vertex_ptr		vertex arranged in xyz-xyz-xyz-
	\param	vertex_number	number of vertices
	\param	face_ptr		triangle face arranged in xyz-xyz-xyz-
	\param	face_number		number of faces
	\param	coef			linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
*/
template <typename T>
void smoothLaplacianTriMesh(
	T*				vertex_ptr,
	unsigned int	vertex_number,
	const int*		face_ptr,
	unsigned int	face_number,
	T				coef = 1.0
) {
	if (!vertex_ptr || !face_ptr || !vertex_number || !face_number)
		return;

	//	create vv-hash
	std::unordered_multimap<int, int>	vv_hash;
	createVVHashFromTriFace(vv_hash, face_ptr, face_number);

	smoothLaplacianTriMesh(
		vertex_ptr,
		vertex_number,
		face_ptr,
		face_number,
		vv_hash,
		coef
	);
}

/**
	perform Laplacian smooth on triangle mesh

	\param	vertex		vertex array
	\param	face		triangle face array
	\param	vv_hash		vertex-vertex hash table
	\param	coef		linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
*/
template <typename T>
void smoothLaplacianTriMesh(
	std::vector<yz::Vec3<T>>&					vertex,
	const std::vector<yz::int3>&				face,
	const std::unordered_multimap<int, int>&	vv_hash,
	T											coef = 1.0
) {
	if (vertex.empty() || face.empty())
		return;

	smoothLaplacianTriMesh(&vertex[0].x, vertex.size(), &face[0].x, face.size(), vv_hash, coef);
}

/**
	perform Laplacian smooth on triangle mesh

	\param	vertex		vertex array
	\param	face		triangle face array
	\param	coef		linear interpolation between original and Laplacian, 0: total original; 1: total Laplacian
*/
template <typename T>
void smoothLaplacianTriMesh(
	std::vector<yz::Vec3<T>>&		vertex,
	const std::vector<yz::int3>&	face,
	T								coef = 1.0
) {
	if (vertex.empty() || face.empty())
		return;

	std::unordered_multimap<int, int>	vv_hash;
	createVVHashFromTriFace(vv_hash, face);

	smoothLaplacianTriMesh(vertex, face, vv_hash, coef);
}


///@}

//	========================================
///@{
/**	@name Taubin Smooth on Mesh
*/
//	========================================

/**
	Perform Taubin smooth on triangle mesh

	Taubin smooth is to laplacian forward, then backward, so that the volume is almost preserved.
	suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

	\param	vertex				vertex array
	\param	face				triangle face array
	\param	vv_hash				vertex-vertex hash table
	\param	iterations			number of smoothing iterations
	\param	coef_forward		forward coefficient of Taubin smooth
	\param	coef_backward		backward coefficient of Taubin smooth
*/
template <typename T>
void smoothTaubinTriMesh(
	std::vector<yz::Vec3<T>>&					vertex,
	const std::vector<yz::int3>&				face,
	const std::unordered_multimap<int, int>&	vv_hash,
	int											iterations = 10,
	T											coef_forward = 0.4507499669,
	T											coef_backward = -0.4720265626
) {
	if (vertex.empty() || face.empty())
		return;

	for (int i = 0; i < iterations; i++) {
		smoothLaplacianTriMesh(&vertex[0].x, vertex.size(), &face[0].x, face.size(), vv_hash, coef_forward);
		smoothLaplacianTriMesh(&vertex[0].x, vertex.size(), &face[0].x, face.size(), vv_hash, coef_backward);
	}
}

/**
	perform Taubin smooth on triangle mesh

	Taubin smooth is to laplacian forward, then backward, so that the volume is almost preserved.
	suggested coefficient pairs are : (0.33, -0.34), (0.4507499669, -0.4720265626)

	\param	vertex				vertex array
	\param	face				triangle face array
	\param	iterations			number of smoothing iterations
	\param	coef_forward		forward coefficient of Taubin smooth
	\param	coef_backward		backward coefficient of Taubin smooth
*/
template <typename T>
void smoothTaobinTriMesh(
	std::vector<yz::Vec3<T>>&		vertex,
	const std::vector<yz::int3>&	face,
	int								iterations = 10,
	T								coef_forward = 0.4507499669,
	T								coef_backward = -0.4720265626
) {
	if (vertex.empty() || face.empty())
		return;

	std::unordered_multimap<int, int>	vv_hash;
	createVVHashFromTriFace(vv_hash, face);

	smoothTaubinTriMesh(vertex, face, vv_hash, iterations, coef_forward, coef_backward);
}



///@}

}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_SMOOTH_H__