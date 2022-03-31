/***********************************************************/
/**	\file
	\brief		Signed Distance Field
	\author		Yizhong Zhang
	\date		11/22/2017
*/
/***********************************************************/
#ifndef __YZ_SDF_H__
#define __YZ_SDF_H__

#include <iostream>
#include <fstream>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_nn.h"
#include "yzLib/yz_geometry/yz_field/yz_field_base.h"

namespace yz{	namespace geometry{		namespace field{

/**
	3D Axis Aligned Signed Distance Field
*/
template<class T>
class SignedDistanceField : public FieldBase<T>{
public:
	using FieldBase<T>::dim;
	using FieldBase<T>::data;

public:
	SignedDistanceField(){
		dim.x = dim.y = dim.z = 0;
		voxel_size = 0;
	}

	/**
		setup data of the volume

		\param	dim_x		dimension in x direction
		\param	dim_y		dimension in y direction
		\param	dim_z		dimension in z direction
		\param	voxel_size	size of each voxel
		\param	xyz0		coordinate origin of the volume
	*/
	void SetupVolume(
		unsigned int	dim_x,
		unsigned int	dim_y,
		unsigned int	dim_z,
		T				voxel_size,
		Vec3<T>			xyz0)
	{
		dim = uint3(dim_x, dim_y, dim_z);
		data.resize(dim.x * dim.y * dim.z);

		this->voxel_size = voxel_size;
		this->xyz0 = xyz0;
	}

	/**
		create the sdf from a closed triangle mesh

		we assume that the input mesh is closed manifold and self-intersect free

		\param	vertex		vertex list of a mesh
		\param	face		face list of a mesh
	*/
	void CreateSDFFromClosedMesh(
		const std::vector<Vec3<T>>& vertex,
		const std::vector<int3>&	face)
	{
		if (data.empty()){
			std::cout << "error: SignedDistanceField::CreateSDFFromClosedMesh, data not created. SetupVolume() must be called first" << std::endl;
			return;
		}

		AABBTree3D<T> aabb;
		aabb.BuildTriangleAABBTree(vertex, face);

#pragma omp parallel	//	the calculation of each voxel is independent
		for (int k = 0; k < dim.z; k++){
			Vec3<T> p(xyz0.x, xyz0.y, xyz0.z + k*voxel_size);
			for (int j = 0; j < dim.y; j++){
				p.y = xyz0.y + j * voxel_size;
				for (int i = 0; i < dim.x; i++){
					p.x = xyz0.x + i * voxel_size;
					Vec3<T> np;
					int nfid = getNearestPointOnMesh(np, p, vertex, face, aabb);
					T dist = (p - np).Length();
					//	add noise and perform intersection test in 3 directions
					int test1 = isPointInsideClosedMesh(p + Vec3<T>(0, 1, 2)*voxel_size*1e-2, vertex, face, aabb, 'X');
					int test2 = isPointInsideClosedMesh(p + Vec3<T>(1, 0, 2)*voxel_size*1e-2, vertex, face, aabb, 'Y');
					int test3 = isPointInsideClosedMesh(p + Vec3<T>(2, 1, 0)*voxel_size*1e-2, vertex, face, aabb, 'Z');
					if (test1 + test2 + test3 >= 2)
						dist = -dist;
					data[this->GetVoxelID(i, j, k)] = dist;
				}
			}
		}
	}

	/**
		read the distance field from .sdf file in binary format

		the file is arranged in the following format:\n
		dim				3 int	\n
		voxel_size		1 double	\n
		xyz0			3 double	\n
		sdf				dim.x*dim.y*dim.z double	\n
		*/
	int ReadSDFFromBinaryFile(const char* file_name){
		std::ifstream file(file_name, std::ifstream::binary);
		if (!file.is_open()){
			std::cout << "error: SignedDistanceField::ReadSDFFromBinaryFile, open file failed" << std::endl;
			return 0;
		}

		//	read parameters
		uint3 i_val;
		file.read((char*)&i_val[0], sizeof(unsigned int) * 3);
		dim = i_val;

		double3 d_val;
		file.read((char*)&d_val[0], sizeof(double));
		voxel_size = d_val[0];
		file.read((char*)&d_val[0], sizeof(double) * 3);
		xyz0 = d_val;

		//	read data
		data.clear();
		data.reserve(dim.x*dim.y*dim.z);
		std::vector<double> d_sdf;
		d_sdf.resize(dim.x*dim.y);
		if (!d_sdf.empty()){
			for (int k = 0; k < dim.z; k++){
				file.read((char*)&d_sdf[0], sizeof(double)*d_sdf.size());
				for (int i = 0; i < d_sdf.size(); i++)
					data.push_back(d_sdf[i]);
			}
		}

		file.close();

		return 1;
	}

	/**
		write the distance field to .sdf file in binary format

		file format is the same as ReadSDFFromBinaryFile()
	*/
	int WriteSDFToBinaryFile(const char* file_name){
		std::ofstream file(file_name, std::ofstream::binary);
		if (!file.is_open()){
			std::cout << "error: SignedDistanceField::WriteSDFToBinaryFile, open file failed" << std::endl;
			return 0;
		}

		//	write parameters
		uint3 i_val = dim;
		file.write((char*)&i_val[0], sizeof(unsigned int) * 3);

		double d_voxel_size = voxel_size;
		double3 d_val = xyz0;
		file.write((char*)&d_voxel_size, sizeof(double));
		file.write((char*)&d_val[0], sizeof(double) * 3);

		//	write data
		std::vector<double> d_sdf;
		d_sdf.resize(dim.x*dim.y);
		if (!d_sdf.empty()){
			for (int k = 0; k < dim.z; k++){
				for (int i = 0; i < d_sdf.size(); i++)
					d_sdf[i] = data[k*dim.x*dim.y + i];
				file.write((char*)&d_sdf[0], sizeof(double)*d_sdf.size());
			}
		}

		file.close();

		return 1;
	}

	/**
		get the distance value of a given point

		linear interpolation
	*/
	inline T GetSDFInterpTriLinear(Vec3<T> p){
		assert(dim.x > 1 && dim.y > 1 && dim.z > 1);

		Vec3i	ijk = (p - xyz0) / voxel_size;

		//	check range and clamp, if the point is outside the range, move it to the nearest point of the bonding box
		for (int i = 0; i < 3; i++) {
			if (ijk[i] < 0)
				ijk[i] = 0, p[i] = xyz0[i];
			else if (ijk[i] >= dim[i] - 1)
				ijk[i] = dim[i] - 2, p[i] = xyz0[i] + (dim[i] - 1) * voxel_size;
		}

		//	calculate interpolation coefficient inside this voxel
		Vec3<T> t0 = (p - xyz0) / voxel_size - ijk;
		Vec3<T> t1 = Vec3<T>(1, 1, 1) - t0;

		//	interpolate
		return
			GetSDF(ijk.x + 0, ijk.y + 0, ijk.z + 0) * t1.x * t1.y * t1.z +
			GetSDF(ijk.x + 1, ijk.y + 0, ijk.z + 0) * t0.x * t1.y * t1.z +
			GetSDF(ijk.x + 0, ijk.y + 1, ijk.z + 0) * t1.x * t0.y * t1.z +
			GetSDF(ijk.x + 1, ijk.y + 1, ijk.z + 0) * t0.x * t0.y * t1.z +
			GetSDF(ijk.x + 0, ijk.y + 0, ijk.z + 1) * t1.x * t1.y * t0.z +
			GetSDF(ijk.x + 1, ijk.y + 0, ijk.z + 1) * t0.x * t1.y * t0.z +
			GetSDF(ijk.x + 0, ijk.y + 1, ijk.z + 1) * t1.x * t0.y * t0.z +
			GetSDF(ijk.x + 1, ijk.y + 1, ijk.z + 1) * t0.x * t0.y * t0.z;
	}

	/**
		get the normal of a given point
	*/
	inline Vec3<T> GetNormalInterpTriLinear(Vec3<T> p){
		Vec3<T> n(
			GetSDFInterpTriLinear(p + Vec3<T>(voxel_size, 0, 0)) - GetSDFInterpTriLinear(p - Vec3<T>(voxel_size, 0, 0)),
			GetSDFInterpTriLinear(p + Vec3<T>(0, voxel_size, 0)) - GetSDFInterpTriLinear(p - Vec3<T>(0, voxel_size, 0)),
			GetSDFInterpTriLinear(p + Vec3<T>(0, 0, voxel_size)) - GetSDFInterpTriLinear(p - Vec3<T>(0, 0, voxel_size))
			);

		T len = n.Length();
		if (len < 1e-5)
			n = Vec3<T>(0, 0, 0);
		else
			n /= len;

		return n;
	}

	/**
		Draw the distance field
	*/
	void DrawVolume(){
#ifdef YZ_gl_h
		if (dim.x && dim.y && dim.z) {	//	don't draw empty volume
			Vec3<T> xyz1 = xyz0 + (dim - Vec3ui(1, 1, 1)) * voxel_size;
			opengl::drawAABBWire(xyz0, xyz1);
		}
#else
		std::cout << "gl.h has to be included in order to use DrawVolume() in SignedDistanceField" << std::endl;
		return;
#endif
	}

	/**
		Draw a single slice

		\param	dim_id				which dimension to draw, 0 x, 1 y, 2 z
		\param	slice_id			which slice in this dimension to draw
		\param	min_color_scale		blue color value with respect to voxel_size
		\param	max_color_scale		red color value with respect to voxel_size
	*/
	void DrawSlice(int dim_id, int slice_id, T min_color_scale = 3, T max_color_scale = 3) {
#ifdef YZ_gl_h
		if (dim_id < 0 || dim_id > 2 || slice_id < 0 || slice_id >= dim[dim_id])	//	dimension and slice id must be legal
			return;
		if (dim[(dim_id + 1) % 3] < 2 || dim[(dim_id + 2) % 3] < 2)	//	at least 4 points are required to draw the slice
			return;

		int dim_id_1 = (dim_id + 1) % 3;
		int dim_id_2 = (dim_id + 2) % 3;
		if (dim_id_1 > dim_id_2)
			mySwap(dim_id_1, dim_id_2);

		//	create content to draw
		std::vector<Vec3f>	vertex;		
		vertex.reserve(dim[dim_id_1] * dim[dim_id_2]);
		Vec3f xyz0_f = xyz0;
		for (int j = 0; j < dim[dim_id_2]; j++) {
			for (int i = 0; i < dim[dim_id_1]; i++) {
				Vec3f offset;
				offset[dim_id_1] = i;
				offset[dim_id_2] = j;
				offset[dim_id] = slice_id;
				vertex.push_back(xyz0_f + offset * voxel_size);
			}
		}

		std::vector<uchar3>	color;
		color.reserve(dim[dim_id_1] * dim[dim_id_2]);
		for (int j = 0; j < dim[dim_id_2]; j++) {
			for (int i = 0; i < dim[dim_id_1]; i++) {
				Vec3i voxel_id;
				voxel_id[dim_id_1] = i;
				voxel_id[dim_id_2] = j;
				voxel_id[dim_id] = slice_id;
				T dist = data[this->GetVoxelID(voxel_id.x, voxel_id.y, voxel_id.z)];
				uchar3 sdf_color;
				utils::convertToColorSimple(&sdf_color.x, dist, - min_color_scale * voxel_size, max_color_scale * voxel_size);
				color.push_back(sdf_color);
			}
		}

		std::vector<uint4>	quad;
		quad.reserve((dim[dim_id_1] - 1) * (dim[dim_id_2] - 1));
		for (int j = 1; j < dim[dim_id_2]; j++) {
			for (int i = 1; i < dim[dim_id_1]; i++) {
				quad.push_back(int4(
					(j - 1)*dim[dim_id_1] + (i - 1), 
					(j - 1)*dim[dim_id_1] + i, 
					j*dim[dim_id_1] + i, 
					j*dim[dim_id_1] + (i - 1)));
			}
		}

		//	draw quads
		Vec3f nor;
		nor[dim_id] = 1;
		glNormal3f(nor.x, nor.y, nor.z);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &vertex[0].x);
		glEnableClientState(GL_COLOR_ARRAY);		
		glColorPointer(3, GL_UNSIGNED_BYTE, 0, &color[0].x);

		glDrawElements(GL_QUADS, quad.size() * 4, GL_UNSIGNED_INT, &quad[0].x);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
#else
		std::cout << "gl.h has to be included in order to use DrawSlice() in SignedDistanceField" << std::endl;
		return;
#endif
	}

public:

	/**
		Get the sdf of a voxel
	*/
	inline T GetSDF(int x, int y, int z){
		return this->GetData(x, y, z);
	}


public:
	T		voxel_size;		///<	the size of each voxel
	Vec3<T>	xyz0;			///<	the min corner of the volume
};


}}}	//	namespace yz::geometry::field

#endif	//	__YZ_SDF_H__