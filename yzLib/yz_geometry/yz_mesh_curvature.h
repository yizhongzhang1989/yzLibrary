/***********************************************************/
/**	\file
	\brief		calculate curvature of mesh
	\details	The interface to calculate curvature is PtrCurvTriMesh. 
				See the instruction
	\author		Yizhong Zhang
	\date		6/4/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_CURVATURE_H__
#define __YZ_MESH_CURVATURE_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_vector_utils.h"
#include "yzLib/yz_geometry/yz_tri_mesh.h"
#include "yzLib/yz_geometry/yz_mesh_topology.h"
#include "yzLib/yz_geometry/yz_create_topology.h"

namespace yz{	namespace geometry{
//	declaration
template<class T> class CurvTriMesh;

/**
	curvature relatad parameter container
*/
template<class T>
class PtrTriMeshCurvature{
public:
	T*	mean_curvature_ptr;
	T*	gaussian_curvature_ptr;
	T*	curvature_normal_ptr;
	T*	mixed_area_ptr;
	T*	edge_weight_ptr;

public:
	PtrTriMeshCurvature() : mean_curvature_ptr(NULL), gaussian_curvature_ptr(NULL), 
		curvature_normal_ptr(NULL), mixed_area_ptr(NULL), edge_weight_ptr(NULL){}
};

/**
	curvature container
*/
template<class T>
class TriMeshCurvature{
public:
	std::vector<T>			mean_curvature;
	std::vector<T>			gaussian_curvature;
	std::vector<Vec3<T>>	curvature_normal;
	std::vector<T>			mixed_area;
	std::vector<T>			edge_weight;
};

/**
	Interface mesh used to calculate curvature, 
	including mesh ptr, topology ptr and curvature ptr

	if you want to calculate curvature of a triangle mesh, you should

	1,	alloc memory for both input and output

	2,	set input data, including mesh, normal, topology

	3,	create PtrCurvTriMesh object, and call SetupPtr

	4,	call Calculation function according to your need
*/
template<class T>
class PtrCurvTriMesh :	
	public PtrTriMeshCurvature<T>,
	public PtrSmoothShadingTriMesh<T>,
	public PtrTriMeshVEF,
	public PtrTriMeshBoundaryVertex
{
public:
	using PtrTriMeshCurvature<T>::mean_curvature_ptr;
	using PtrTriMeshCurvature<T>::gaussian_curvature_ptr;
	using PtrTriMeshCurvature<T>::curvature_normal_ptr;
	using PtrTriMeshCurvature<T>::mixed_area_ptr;
	using PtrTriMeshCurvature<T>::edge_weight_ptr;
	using PtrSmoothShadingTriMesh<T>::vertex_ptr;
	using PtrSmoothShadingTriMesh<T>::face_ptr;
	using PtrSmoothShadingTriMesh<T>::vertex_normal_ptr;

	int		vertex_number, face_number, edge_number;

	PtrCurvTriMesh() : PtrTriMeshCurvature<T>(), 
		PtrSmoothShadingTriMesh<T>(), 
		PtrTriMeshVEF(), 
		PtrTriMeshBoundaryVertex(),
		vertex_number(0), face_number(0), edge_number(0),
		edge_weight_guard(0), mixed_area_guard(0){}

	/**
		Setup Data using CurvTriMesh
	*/
	inline void SetupPtr(const CurvTriMesh<T>& curv_tri_mesh){
		vertex_number	= curv_tri_mesh.vertex.size();
		face_number		= curv_tri_mesh.face.size();
		edge_number		= curv_tri_mesh.edge.size();
		if( !curv_tri_mesh.vertex.empty() )				vertex_ptr				= (T*)&curv_tri_mesh.vertex[0];
		if( !curv_tri_mesh.face.empty() )				face_ptr				= (int*)&curv_tri_mesh.face[0];
		if( !curv_tri_mesh.vertex_normal.empty() )		vertex_normal_ptr		= (T*)&curv_tri_mesh.vertex_normal[0];
		if( !curv_tri_mesh.edge.empty() )				edge_ptr				= (int*)&curv_tri_mesh.edge[0];
		if( !curv_tri_mesh.ef.empty() )					ef_ptr					= (int*)&curv_tri_mesh.ef[0];
		if( !curv_tri_mesh.fe.empty() )					fe_ptr					= (int*)&curv_tri_mesh.fe[0];
		if( !curv_tri_mesh.vv.empty() )					vv_ptr					= (int*)&curv_tri_mesh.vv[0];
		if( !curv_tri_mesh.ve.empty() )					ve_ptr					= (int*)&curv_tri_mesh.ve[0];
		if( !curv_tri_mesh.vve_start.empty() )			vve_start_ptr			= (int*)&curv_tri_mesh.vve_start[0];
		if( !curv_tri_mesh.mean_curvature.empty() )		mean_curvature_ptr		= (T*)&curv_tri_mesh.mean_curvature[0];
		if( !curv_tri_mesh.gaussian_curvature.empty() )	gaussian_curvature_ptr	= (T*)&curv_tri_mesh.gaussian_curvature[0];
		if( !curv_tri_mesh.curvature_normal.empty() )	curvature_normal_ptr	= (T*)&curv_tri_mesh.curvature_normal[0];
		if( !curv_tri_mesh.mixed_area.empty() )			mixed_area_ptr			= (T*)&curv_tri_mesh.mixed_area[0];
		if( !curv_tri_mesh.edge_weight.empty() )		edge_weight_ptr			= (T*)&curv_tri_mesh.edge_weight[0];
		if( !curv_tri_mesh.boundary_vertex_flag.empty())boundary_vertex_flag_ptr= (int*)&curv_tri_mesh.boundary_vertex_flag[0];

		ResetGuard();
	}

	/**
		Setup Data by pointers directly

		If you don't need certain pointer, set it to NULL

		vertx_number, face_number, edge_number are needed, can't be 0
	*/
	inline void SetupPtr(
		T*		vertex,
		int		vertex_number,
		int*	face,
		int		face_number,
		int*	edge,
		int		edge_number,
		T*		vertex_normal,
		int*	ef,
		int*	fe,
		int*	vv,
		int*	ve,
		int*	vve_start,
		T*		mean_curvature,
		T*		gaussian_curvature,
		T*		curvature_normal,
		T*		mixed_area,
		T*		edge_weight,
		int*	boundary_vertex_flag)
	{
		vertex_ptr = vertex;
		this->vertex_number = vertex_number;
		face_ptr = face;
		this->face_number = face_number;
		edge_ptr				= edge;
		this->edge_number		= edge_number;
		vertex_normal_ptr		= vertex_normal;
		ef_ptr					= ef;
		fe_ptr					= fe;
		vv_ptr					= vv;
		ve_ptr					= ve;
		vve_start_ptr			= vve_start;
		mean_curvature_ptr		= mean_curvature;
		gaussian_curvature_ptr	= gaussian_curvature;
		curvature_normal_ptr	= curvature_normal;
		mixed_area_ptr			= mixed_area;
		edge_weight_ptr			= edge_weight;
		boundary_vertex_flag_ptr= boundary_vertex_flag;
	}

	/**
		calculate mean curvature of mesh

		This function doesn't allocate any space, you must allocate space yourself before calling this function

		Reference: Discrete Differential-Geometry Operators for Triangulated 2-Manifolds\n
		Mark Meyer, Mathieu Desbrun, Peter Schroder, and Alan H. Barr

		Requested:	vertex_ptr,		\n
					face_ptr,	\n
					vertex_normal_ptr,	\n
					fe_ptr	\n
					boundary_vertex_flag_ptr, vv_ptr, ve_ptr, vve_start_ptr (If mesh is not closed)	\n
		Output:		mean_curvature_ptr,	\n
					curvature_normal_ptr,	\n
					edge_weight_ptr,	\n
					mixed_area_ptr

		\return		1: succeed;		0: failed
	*/
	inline int CalculateMeanCurvature(){
		ResetGuard();

		//	memory check
		if( mean_curvature_ptr==NULL || curvature_normal_ptr==NULL || edge_weight_ptr==NULL || mixed_area_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
				std::cout << "memory of mean_curvature_ptr / curvature_normal_ptr / edge_weight_ptr / mixed_area_ptr are not set" << std::endl;
			#endif
			return 0;
		}
		if( vertex_ptr==NULL || vertex_number==0 || face_ptr==NULL || face_number==0 || vertex_normal_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
				std::cout << "vertex/face/vertex_normal of input mesh is not set." << std::endl;
				std::cout << "vertex_number: " << vertex_number << ", face_number: " << face_number << std::endl;
			#endif
			return 0;
		}
		if( fe_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
				std::cout << "fe is not set" << std::endl;
			#endif
			return 0;
		}
		if( boundary_vertex_flag_ptr!=NULL && (vv_ptr==NULL || ve_ptr==NULL || vve_start_ptr==NULL)){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
				std::cout << "mesh is not closed, vv, ve, vve_start must be provided" << std::endl;
			#endif
			return 0;
		}

		//	calculate edge weight & mixed area
		if( !CalculateEdgeWeightMixedArea() ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
			#endif
			return 0;
		}

		//	calculate curvature normal
		if( !CalculateCurvatureNormal() ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMeanCurvature()" << std::endl;
			#endif
			return 0;
		}

		//	calculate curvature
		for( int v=0; v<vertex_number; v++ ){
			Vec3<T>	v_nor	= ((Vec3<T>*)vertex_normal_ptr)[v];
			Vec3<T> c_nor	= ((Vec3<T>*)curvature_normal_ptr)[v];
			mean_curvature_ptr[v] = c_nor.Length();
			if( dot(v_nor, c_nor) < 0 )
				mean_curvature_ptr[v] = -mean_curvature_ptr[v];
		}

		return 1;
	}

	/**
		calculate gaussian curvature of mesh

		This function doesn't allocate any space, you must allocate space yourself before calling this function

		Reference: Discrete Differential-Geometry Operators for Triangulated 2-Manifolds\n
		Mark Meyer, Mathieu Desbrun, Peter Schroder, and Alan H. Barr

		Requested:	vertex_ptr,	\n
					face_ptr,	\n
					vv_ptr,	\n
					ve_ptr,	\n
					vve_start_ptr,	\n
					ef_ptr,	\n
					boundary_vertex_flag_ptr (If mesh is not closed)	\n
		Output:		gaussian_curvature_ptr,	\n
					mixed_area_ptr

		\return		1: succeed;		0: failed
	*/
	inline int CalculateGaussianCurvature(){
		ResetGuard();

		//	memory check
		if( gaussian_curvature_ptr==NULL || mixed_area_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateGaussianCurvature()" << std::endl;
				std::cout << "memory of mean_curvature_ptr / mixed_area_ptr are not set" << std::endl;
			#endif
			return 0;
		}
		if( vertex_ptr==NULL || vertex_number==0 || face_ptr==NULL || face_number==0 ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateGaussianCurvature()" << std::endl;
				std::cout << "vertex/face of input mesh is not set." << std::endl;
				std::cout << "vertex_number: " << vertex_number << ", face_number: " << face_number << std::endl;
			#endif
			return 0;
		}
		if( vv_ptr==NULL || ve_ptr==NULL || vve_start_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateGaussianCurvature()" << std::endl;
				std::cout << "vv / ve / vve_start must be provided" << std::endl;
			#endif
			return 0;
		}
		if( ef_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateGaussianCurvature()" << std::endl;
				std::cout << "ef is not set" << std::endl;
			#endif
			return 0;
		}

		//	calculate edge weight & mixed area
		if( !CalculateMixedArea() ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateGaussianCurvature()" << std::endl;
			#endif
			return 0;
		}

		//	calculate caussian curvature
		for( int v=0; v<vertex_number; v++ ){
			T angle_rad = 0;
			for( int id=vve_start_ptr[v]; id<vve_start_ptr[v+1]; id++ ){	//	calculate one ring neighbor angle sum
				int nei_v = vv_ptr[id];
				int nei_e = ve_ptr[id];
				int nei_f = ((int2*)ef_ptr)[nei_e].y;
				if( v==((int2*)edge_ptr)[nei_e].x && nei_v==((int2*)edge_ptr)[nei_e].y )
					nei_f = ((int2*)ef_ptr)[nei_e].x;

				if( nei_f < 0 )
					continue;

				int nei_v2 = getThirdVertexOfFace(((int3*)face_ptr)[nei_f], v, nei_v);


				//	we have got the neighbor face, then calculate the angle
				Vec3<T> r0 = ((Vec3<T>*)vertex_ptr)[nei_v] - ((Vec3<T>*)vertex_ptr)[v];
				Vec3<T> r1 = ((Vec3<T>*)vertex_ptr)[nei_v2] - ((Vec3<T>*)vertex_ptr)[v];

				angle_rad += angleRadBetweenVectors(r0, r1);
			}

			gaussian_curvature_ptr[v] = ((2*YZ_PI) - angle_rad) / mixed_area_ptr[v];
		}

		//	smooth boundary
		if( boundary_vertex_flag_ptr )
			setTriMeshBoundaryDataAsNeighborAverage(gaussian_curvature_ptr, 
				vertex_number, boundary_vertex_flag_ptr, vv_ptr, vve_start_ptr);

		return 1;
	}

	/**
		Calculate Edge Weight and Mixed Area

		If one pointer is NULL, just calculate another one

		calculateMixedArea() copied code from here
		if any bug detected in this function, check that function as well.

		\return		1: succeed;		0: failed
	*/
	inline int CalculateEdgeWeightMixedArea(){
		ResetGuard();

		//	memory check
		if( edge_weight_ptr==NULL && mixed_area_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateEdgeWeightMixedArea()" << std::endl;
				std::cout << "memory of edge_weight_ptr & mixed_area_ptr are not set" << std::endl;
			#endif
			return 0;
		}
		if( vertex_ptr==NULL || vertex_number==0 || face_ptr==NULL || face_number==0 ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateEdgeWeightMixedArea()" << std::endl;
				std::cout << "vertex/face of input mesh is not set." << std::endl;
				std::cout << "vertex_number: " << vertex_number << ", face_number: " << face_number << std::endl;
			#endif
			return 0;
		}
		if( fe_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateEdgeWeightMixedArea()" << std::endl;
				std::cout << "fe is not set" << std::endl;
			#endif
			return 0;
		}

		//	clear old data, incremental add is used
		if (edge_weight_ptr && edge_number)
			std::fill(&edge_weight_ptr[0], &edge_weight_ptr[edge_number], 0);//memset(edge_weight_ptr, 0, sizeof(T)*edge_number);
		if(mixed_area_ptr && vertex_number)	
			std::fill(&mixed_area_ptr[0], &mixed_area_ptr[vertex_number], 0);//memset(mixed_area_ptr, 0, sizeof(T)*vertex_number);

		//	for each face, calculate new data
		for( int f=0; f<face_number; f++ ){
			int		v_idx[3] = { face_ptr[f*3], face_ptr[f*3+1], face_ptr[f*3+2] };
			Vec3<T> v[3] = {	Vec3<T>(vertex_ptr + v_idx[0]*3),
								Vec3<T>(vertex_ptr + v_idx[1]*3),
								Vec3<T>(vertex_ptr + v_idx[2]*3)	};

			//	edge weight
			T cot_val[3];
			for( int i=0; i<3; i++ ){
				cot_val[i] = cotBetweenVectors(v[(i+1)%3]-v[i], v[(i+2)%3]-v[i]);
				if(edge_weight_ptr)	edge_weight_ptr[fe_ptr[f*3+i]] += cot_val[i];	//	edge should be at the opposite position of vertex
			}

			//	mixed area
			if( mixed_area_ptr ){
				if( cot_val[0]>=0 && cot_val[1]>=0 && cot_val[2]>=0 ){	//	non-obtuse triangle
					for(int i=0; i<3; i++ ){
						mixed_area_ptr[v_idx[i]] += cot_val[(i+1)%3] * (v[(i+2)%3]-v[i]).SquareLength() / 8;
						mixed_area_ptr[v_idx[i]] += cot_val[(i+2)%3] * (v[(i+1)%3]-v[i]).SquareLength() / 8;
					}
				}
				else{													//	obtuse triangle
					T area = cross(v[1]-v[0], v[2]-v[0]).Length() / 2;
					for(int i=0; i<3; i++ )
						mixed_area_ptr[v_idx[i]] += area/(cot_val[i]<0 ? 2 : 4);	//	obtuse angle take half the area
				}
			}
		}

		if( edge_weight_ptr )	edge_weight_guard	= 1;
		if( mixed_area_ptr )	mixed_area_guard	= 1;
		return 1;
	}

	/**
		Calculate Edge Weight

		\return		1: succeed;		0: failed
	*/
	inline int CalculateEdgeWeight(){
		//	memory check
		if( edge_weight_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateEdgeWeight()" << std::endl;
				std::cout << "memory of edge_weight_ptr is not set" << std::endl;
			#endif
			return 0;
		}

		//	set mixed_area_ptr to NULL and call CalculateEdgeWeightMixedArea()
		int tmp_mixed_area_guard	= mixed_area_guard;
		T* tmp_mixed_area_ptr		= mixed_area_ptr;
		mixed_area_ptr = NULL;		//	set to NULL, so only edge weight will be calculated
		int ret = CalculateEdgeWeightMixedArea();
		mixed_area_guard	= tmp_mixed_area_guard;
		mixed_area_ptr		= tmp_mixed_area_ptr;

		return ret;
	}

	/**
		Calculate Mixed Area

		\return		1: succeed;		0: failed
	*/
	inline int CalculateMixedArea(){
		//	memory check
		if( mixed_area_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateMixedArea()" << std::endl;
				std::cout << "memory of mixed_area_ptr is not set" << std::endl;
			#endif
			return 0;
		}

		//	set edge_weight_ptr to NULL and call CalculateEdgeWeightMixedArea()
		int tmp_edge_weight_guard	= edge_weight_guard;
		T* tmp_edge_weight_ptr		= edge_weight_ptr;
		edge_weight_ptr = NULL;		//	set to NULL, so only mixed area will be calculated
		int ret = CalculateEdgeWeightMixedArea();
		edge_weight_guard	= tmp_edge_weight_guard;
		edge_weight_ptr		= tmp_edge_weight_ptr;

		return ret;
	}

protected:
	int edge_weight_guard, mixed_area_guard;

	/**
		Reset Guard
	*/
	inline void ResetGuard(){
		edge_weight_guard		= 0;
		mixed_area_guard		= 0;
	}

	/**
		Calculate Curvature Normal

		\return		1: succeed;		0: failed
	*/
	inline int CalculateCurvatureNormal(){
		//	memory check
		if( curvature_normal_ptr==NULL ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateCurvatureNormal()" << std::endl;
				std::cout << "memory of curvature_normal_ptr is not set" << std::endl;
			#endif
			return 0;
		}
		if( vertex_ptr==NULL || vertex_number==0 || face_ptr==NULL || face_number==0 ){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateCurvatureNormal()" << std::endl;
				std::cout << "vertex/face of input mesh is not set." << std::endl;
				std::cout << "vertex_number: " << vertex_number << ", face_number: " << face_number << std::endl;
			#endif
			return 0;
		}
		if( fe_ptr==NULL && (vv_ptr==NULL && ve_ptr==NULL && vve_start_ptr==NULL)){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateCurvatureNormal()" << std::endl;
				std::cout << "fe & vve are not set, at least one should be provided" << std::endl;
			#endif
			return 0;
		}
		if( boundary_vertex_flag_ptr!=NULL && (vv_ptr==NULL || vve_start_ptr==NULL)){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateCurvatureNormal()" << std::endl;
				std::cout << "mesh is not closed, vv, vve_start must be provided" << std::endl;
			#endif
			return 0;
		}
		if( edge_weight_ptr==NULL || !edge_weight_guard || mixed_area_ptr==NULL || !mixed_area_guard){
			#ifndef BE_QUIET
				std::cout << "Error: PtrCurvTriMesh.CalculateCurvatureNormal()" << std::endl;
				std::cout << "edge_weight & mixed_area must be calculate first" << std::endl;
			#endif
			return 0;
		}

		//	calculate according to face-edge connectivity
		if( fe_ptr ){
			std::fill(&curvature_normal_ptr[0], &curvature_normal_ptr[vertex_number * 3], 0);//memset(curvature_normal_ptr, 0, sizeof(T)*vertex_number*3);

			//	for each edge, add weighted edge to corresponding vertex
			for( int f=0; f<face_number; f++ ){
				int		v_idx[3] = { face_ptr[f*3], face_ptr[f*3+1], face_ptr[f*3+2] };
				Vec3<T> v[3] = {	Vec3<T>(vertex_ptr + v_idx[0]*3),
									Vec3<T>(vertex_ptr + v_idx[1]*3),
									Vec3<T>(vertex_ptr + v_idx[2]*3)	};

				for( int i=0; i<3; i++ ){
					T	e_w		= edge_weight_ptr[fe_ptr[f*3+i]];
					((Vec3<T>*)curvature_normal_ptr)[v_idx[(i+2)%3]] += e_w * (v[(i+2)%3] - v[(i+1)%3]);
				}
			}

			//	for each vertex, divide mixed area
			for( int v=0; v<vertex_number; v++ )
				((Vec3<T>*)curvature_normal_ptr)[v] /= (2*mixed_area_ptr[v]);
		}
		//	calculate according to vertex-vertex, vertex-edge connectivity
		else if( vv_ptr && ve_ptr && vve_start_ptr ){
			for( int v_idx=0; v_idx<vertex_number; v_idx++ ){
				Vec3<T> v(vertex_ptr + v_idx*3);
				((Vec3<T>*)curvature_normal_ptr)[v_idx] = Vec3<T>(0, 0, 0);

				//	for each one ring neighbor of this vertex, add weighted edge
				for( int vn_idx=vve_start_ptr[v_idx]; vn_idx<vve_start_ptr[v_idx+1]; vn_idx++ ){
					Vec3<T> vn(vertex_ptr + vv_ptr[vn_idx]*3);
					int e_idx	= ve_ptr[vn_idx];
					T	e_w		= edge_weight_ptr[e_idx];
					((Vec3<T>*)curvature_normal_ptr)[v_idx] += e_w * (v - vn);
				}

				((Vec3<T>*)curvature_normal_ptr)[v_idx] /= (2*mixed_area_ptr[v_idx]);
			}
		}

		//	if the mesh is not closed, curvature normal on the boundary must be recalculated
		if( boundary_vertex_flag_ptr ){
			setTriMeshBoundaryDataAsNeighborAverage((Vec3<T>*)curvature_normal_ptr, 
				vertex_number, boundary_vertex_flag_ptr, vv_ptr, vve_start_ptr);
		}

		return 1;
	}


};

typedef PtrCurvTriMesh<float>	CurvatureCalculatorInterfacef;
typedef PtrCurvTriMesh<double>	CurvatureCalculatorInterfaced;

/**
	Triangle Mesh With Mean Curvature
*/
template<class T>
class CurvTriMesh :
	public TriMeshCurvature<T>,
	public TriMesh<T>,
	public TriMeshVEF,
	public TriMeshBoundaryVertex
{
public:
	using TriMeshCurvature<T>::mean_curvature;
	using TriMeshCurvature<T>::gaussian_curvature;
	using TriMeshCurvature<T>::curvature_normal;
	using TriMeshCurvature<T>::mixed_area;
	using TriMeshCurvature<T>::edge_weight;
	using TriMesh<T>::vertex;
	using TriMesh<T>::face;
	using TriMesh<T>::vertex_normal;
	using TriMesh<T>::face_normal;


public:
	/**
		Create mesh topology
	*/
	inline void CreateTopology(){
		TriMeshVEF::CreateTopology(vertex.size(), face);
		MarkBoundaryVertexFromEF(vertex.size(), edge, ef);
	}

	/**
		Calculate mean curvature of the mesh

		\return		1: succeed;		0: failed
	*/
	inline int CalculateMeanCurvature(){
		if(vertex_normal.empty())	this->CalculateVertexNormal();
		if(edge.empty())			CreateTopology();

		edge_weight.resize( edge.size() );
		mixed_area.resize( vertex.size() );
		curvature_normal.resize( vertex.size() );
		mean_curvature.resize( vertex.size() );

		//	calculate using PtrCurvTriMesh
		PtrCurvTriMesh<T>	ptr_mesh;
		ptr_mesh.SetupPtr( *this );
		return ptr_mesh.CalculateMeanCurvature();
	}

	/**
		Calculate Gaussian Curvature of the mesh

		\return		1: succeed;		0: failed
	*/
	inline int CalculateGaussianCurvature(){
		if(edge.empty())	CreateTopology();

		mixed_area.resize( vertex.size() );
		gaussian_curvature.resize( vertex.size() );

		//	calculate using PtrCurvTriMesh
		PtrCurvTriMesh<T>	ptr_mesh;
		ptr_mesh.SetupPtr( *this );
		return ptr_mesh.CalculateGaussianCurvature();
	}
};

typedef CurvTriMesh<float>	CurvTriMeshf;
typedef CurvTriMesh<double>	CurvTriMeshd;

/**
	Calculate Mean Curvature of a Mesh

	\param	mean_curvature	mean_curvature list ptr, size of list must be at least sizeof(T) * vertex_number
	\param	vertex			vertex ptr
	\param	vertex_number	vertex number
	\param	face			face ptr
	\param	face_number		face number
	\return					1: succeed;		0: failed
*/
template<typename T>
inline int calculateMeanCurvatureFromMesh(
	T*			mean_curvature,
	const T*	vertex,
	int			vertex_number,
	const int*	face,
	int			face_number)
{
	CurvTriMesh<T> curv_tri_mesh;

	curv_tri_mesh.vertex.resize(vertex_number);
	curv_tri_mesh.face.resize(face_number);
	memcpy(&curv_tri_mesh.vertex[0], vertex, sizeof(T)*vertex_number*3);
	memcpy(&curv_tri_mesh.face[0], face, sizeof(int)*face_number*3);

	if( curv_tri_mesh.CalculateMeanCurvature() ){
		memcpy(mean_curvature, &curv_tri_mesh.mean_curvature[0], sizeof(T)*vertex_number);
		return 1;
	}
	else{
		std::cout << "calculate mean curvature from mesh failed" << std::endl;
		return 0;
	}
}

/**
	Calculate Mean Curvature of a Mesh using vector

	\param	mean_curvature	mean_curvature vector
	\param	vertex			vertex vector
	\param	face			face vector
	\return					1: succeed;		0: failed
*/
template<typename T>
inline int calculateMeanCurvatureFromMesh(
	std::vector<T>&					mean_curvature,
	const std::vector<Vec3<T>>&		vertex,
	const std::vector<int3>&		face)
{
	mean_curvature.resize( vertex.size() );
	return calculateMeanCurvatureFromMesh((T*)&mean_curvature[0], (T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}


/**
	Set Boundary Data of TriMesh as Average of neighboring legal data

	type T must support +=, * operators

	\param	data					the data defined on each vertex, boundary vertex data will be replaced
	\param	vertex_number			vertex number
	\param	boundary_vertex_flag	flag array to check whether each vertex is a boundary vertex
	\param	vv						vertex-vertex connectivity
	\param	vv_start				vertex neighbor start index of each vertex
	\return							succeed
*/
template<typename T>
inline int setTriMeshBoundaryDataAsNeighborAverage(
	T*			data,
	int			vertex_number,
	const int*	boundary_vertex_flag,
	const int*	vv,
	const int*	vv_start)
{
	assert(data && vertex_number>=0 && boundary_vertex_flag && vv && vv_start);

	std::vector<int>	vertices_to_do;		//	boundary vertices
	std::vector<int>	vertices_st;		//	unconnected corner vertices, need special treate

	//	create boundary vertices
	for( int v=0; v<vertex_number; v++ )
		if( boundary_vertex_flag[v] )	vertices_to_do.push_back(v);

	for( int i=0; i<vertices_to_do.size(); i++ ){
		int count = 0;
		T	sum;
		//	scan one ring neighbor of this vertex, calculate the avarage of all legal curvature normal
		for( int vn_idx=vv_start[vertices_to_do[i]]; vn_idx<vv_start[vertices_to_do[i]+1]; vn_idx++ ){
			if( !boundary_vertex_flag[vv[vn_idx]] ){
				count ++;
				if( count == 1 )
					sum = data[vv[vn_idx]];
				else
					sum += data[vv[vn_idx]];
			}
		}
		if( count )
			data[vertices_to_do[i]] = sum * (1.0/count);
		else
			vertices_st.push_back(vertices_to_do[i]);
	}

	//	then treate corner vertices, or other special vertices
	for( int i=0; i<vertices_st.size(); i++ ){
		int count = 0;
		T	sum;
		for( int vn_idx=vv_start[vertices_st[i]]+1; vn_idx<vv_start[vertices_st[i]+1]; vn_idx++ ){
			if( std::find(vertices_st.begin(), vertices_st.end(), vv[vn_idx]) == vertices_st.end() ){	//	didn't find
				count ++;
				if( count == 1 )
					sum = data[vv[vn_idx]];
				else
					sum += data[vv[vn_idx]];
			}
		}
		if( count )
			data[vertices_st[i]] = sum * (1.0/count);
		else
			data[vertices_st[i]] = sum * 0;	//	make it zero
	}
	return 1;
}


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_CURVATURE_H__