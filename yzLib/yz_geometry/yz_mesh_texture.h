/***********************************************************/
/**	\file
	\brief		Mesh texture
	\details	Texture of the mesh
	\author		Yizhong Zhang
	\date		10/26/2014
*/
/***********************************************************/
#ifndef __YZ_MESH_TEXTURE_H__
#define __YZ_MESH_TEXTURE_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_image_rw_utils.h"

namespace yz{	namespace geometry{

/**
	texture of a mesh
*/
template<class T> 
class MeshTextureCoordinate{
public:
	std::vector<Vec2<T>>	tex_coord;
};

class MeshTextureFace {
public:
	std::vector<int3>	tex_face;
};

/**
	maps that defined in material, currently not complete
*/
class MeshMaterialMap{
public:
	std::string				file_dir;
	std::string				file_name;

	int						tex_w, tex_h;
	std::vector<uchar3>		tex_image;

	//	options of map_Ka, map_Kd, map_Ks
	int			blendu;
	int			blendv;
	int			cc;
	int			clamp;
	double		mm_base, mm_gain;
	double3		o, s, t;
	double		texres;

	//	options specially for map_Ns, map_d, decal, disp
	char		imfchan;

	//	options specially for bump
	double		bm;

public:
	MeshMaterialMap(){
		tex_w = tex_h = 0;
	}

	int ReadTexImage(){
		if( file_name.empty() )
			return 1;

		std::string name = utils::getFileNameCombineDirentory( file_dir.c_str(), file_name.c_str() );

		unsigned char* img_ptr = NULL;
		int w, h;
		int ret = utils::readImageFromFile(name.c_str(), img_ptr, w, h, 24);
		if (!ret || w <= 0 || h <= 0){
			//	read image failed, we create a chessboard
			w = 128;
			h = 128;
			img_ptr = new unsigned char[w*h * 3];
			for (int j = 0; j < 8; j++){
				for (int i = 0; i < 8; i++){
					int white_flag = (i + j) % 2;
					for (int jj = 0; jj < 16; jj++){
						for (int ii = 0; ii < 16; ii++){
							int jjj = j * 16 + jj;
							int iii = i * 16 + ii;
							int pid = jjj * 128 + iii;
							img_ptr[pid * 3] = white_flag ? 255 : 0;
							img_ptr[pid * 3 + 1] = 255;
							img_ptr[pid * 3 + 2] = 255;
						}
					}
				}
			}
		}

		//	copy flipped image 
		tex_w = w;
		tex_h = h;
		tex_image.resize( w * h );
		for(int j=0; j<h; j++){
			memcpy( &tex_image[j*w], img_ptr+(h-j-1)*w*3, sizeof(uchar3)*w );
		}

		delete[] img_ptr;

		return 1;
	}

	int WriteTexImage(){
		if( file_name.empty() )
			return 0;

		std::string name = utils::getFileNameCombineDirentory( file_dir.c_str(), file_name.c_str() );

		unsigned char* img_ptr = new unsigned char[tex_w*tex_h*3];
		if( !img_ptr )
			return 0;
		for(int j=0; j<tex_h; j++){
			memcpy( img_ptr+(tex_h-j-1)*tex_w*3, &tex_image[j*tex_w], sizeof(uchar3)*tex_w );
		}

		int ret = utils::writeImageToFile(name.c_str(), img_ptr, tex_w, tex_h, 24);

		delete[] img_ptr;

		return ret;
	}
};

/**
	material of a mesh, currently not complete
*/
class MeshMaterial{
public:
	std::string				name;

	double3					Ka, Kd, Ks, Tf;

	int						illum;

	double					d_factor;
	int						halo;

	double					Ns;
	double					sharpness;
	double					Ni;

	MeshMaterialMap			map_Ka, map_Kd, map_Ks;

	MeshMaterial(){
		Ka = double3(0, 0, 0);
		Kd = double3(1, 1, 1);
		Ks = double3(0, 0, 0);
	}
};


inline int readMtlKadsFromString(double3& rgb, std::string line){
	std::stringstream ss;
	std::string cmd;
	ss << line;
	ss >> cmd;

	//	Ka support 3 formats, check which one is it
	ss >> cmd;
	if( cmd == "spectral" ){	//	Ka spectral file.rfl factor
		rgb = double3(1,1,1);		//	currently not supported
	}
	else if( cmd == "xyz" ){	//	Ka xyz x y z
		rgb = double3(1,1,1);		//	currently not supported
	}
	else{						//	Ka r g b
		ss << line;
		ss >> cmd;	//	eat the first command
		ss >> rgb.x;
		if( ! (ss >> rgb.y) )
			rgb.y = rgb.x;
		if( ! (ss >> rgb.z) )
			rgb.z = rgb.x;
	}

	return 1;
}

inline int readMtlMapFromString(std::string& map_file_name, std::string line){
	std::stringstream ss;
	std::string cmd;
	ss << line;
	ss >> cmd;

	//	read options or filename
	while( ss >> cmd ){
		if( cmd[0] == '-' ){	//	an option
			if( cmd == "-blendu" ){
				ss >> cmd;
			}
			else if( cmd == "-blendv" ){
				ss >> cmd;
			}
			else if( cmd == "-cc" ){
				ss >> cmd;
			}
			else if( cmd == "-clamp" ){
				ss >> cmd;
			}
			else if( cmd == "-mm" ){
				ss >> cmd >> cmd;
			}
			else if( cmd == "-o" ){
				ss >> cmd >> cmd >> cmd;
			}
			else if( cmd == "-s" ){
				ss >> cmd >> cmd >> cmd;
			}
			else if( cmd == "-t" ){
				ss >> cmd >> cmd >> cmd;
			}
			else if( cmd == "-texres" ){
				ss >> cmd;
			}
			else{
				std::cout << "error: readMtlMapFromString, undefined option : " << cmd << std::endl;
			}
		}
		else					//	filename
			map_file_name = cmd;
	}

	return 1;
}

/**
	Read a material library file, .mtl
*/
inline int readMaterialFromMtl(
	const char*						mtl_file_name,
	std::vector<MeshMaterial>&		material )
{
	std::ifstream mtl( mtl_file_name );	//	load mtl failed
	if( !mtl.is_open() ){
		#ifndef BE_QUIET
		std::cout << "readMaterialFromMtl, cannot open " << mtl_file_name << std::endl;
		#endif
		return 0;
	}

	std::string dir = utils::getDirectoryFromString( mtl_file_name );

	material.clear();

	MeshMaterial* new_mtl_ptr = NULL;

	int error_log_flag = 0;

	std::string line;
	int line_idx = 0;
	while( std::getline(mtl, line) ){
		line_idx ++;
		std::stringstream ss;
		std::string cmd;
		ss << line;
		ss >> cmd;

		if( cmd == "newmtl" ){		//	a new material
			if( new_mtl_ptr ){	//	a material already exist, record the old one
				material.push_back( *new_mtl_ptr );
				delete new_mtl_ptr;
				new_mtl_ptr = NULL;
			}

			new_mtl_ptr = new MeshMaterial;
			ss >> new_mtl_ptr->name;
		}
		else if( cmd == "Ka" ){
			double3 rgb;
			readMtlKadsFromString(rgb, line);
			new_mtl_ptr->Ka = rgb;
		}
		else if( cmd == "Kd" ){
			double3 rgb;
			readMtlKadsFromString(rgb, line);
			new_mtl_ptr->Kd = rgb;					
		}
		else if( cmd == "Ks" ){
			double3 rgb;
			readMtlKadsFromString(rgb, line);
			new_mtl_ptr->Ks = rgb;					
		}
		else if( cmd == "Tf" ){
			double3 rgb;
			readMtlKadsFromString(rgb, line);
			new_mtl_ptr->Tf = rgb;					
		}
		else if( cmd == "illum" ){
			ss >> new_mtl_ptr->illum;
		}
		else if( cmd == "d" ){
			ss >> cmd;
			if( cmd == "-halo" ){
				new_mtl_ptr->halo = 1;
				ss >> new_mtl_ptr->d_factor;
			}
			else{
				ss << line;
				ss >> cmd;
				ss >> new_mtl_ptr->d_factor;
			}
		}
		else if( cmd == "Ns" ){
			ss >> new_mtl_ptr->Ns;
		}
		else if( cmd == "Sharpness" ){
			ss >> new_mtl_ptr->sharpness;
		}
		else if( cmd == "Ni" ){
			ss >> new_mtl_ptr->Ni;
		}
		else if( cmd == "map_Ka" ){
			std::string file_name;
			readMtlMapFromString(file_name, line);
			new_mtl_ptr->map_Ka.file_dir = dir.c_str();
			new_mtl_ptr->map_Ka.file_name = file_name.c_str();
		}
		else if( cmd == "map_Kd" ){
			std::string file_name;
			readMtlMapFromString(file_name, line);
			new_mtl_ptr->map_Kd.file_dir = dir.c_str();
			new_mtl_ptr->map_Kd.file_name = file_name.c_str();
		}
		else if( cmd == "map_Ks" ){
			std::string file_name;
			readMtlMapFromString(file_name, line);
			new_mtl_ptr->map_Ks.file_dir = dir.c_str();
			new_mtl_ptr->map_Ks.file_name = file_name.c_str();
		}
	}

	//	add the last material into mtllib
	if( new_mtl_ptr ){
		material.push_back( *new_mtl_ptr );
		delete new_mtl_ptr;
		new_mtl_ptr = NULL;
	}

	mtl.close();

	//	read maps
	for(int i=0; i<material.size(); i++){
		material[i].map_Ka.ReadTexImage();
		material[i].map_Kd.ReadTexImage();
		material[i].map_Ks.ReadTexImage();
	}

	return 1;
}

/**
	Write material library to .mtl file
*/
inline int writeMaterialToMtl(
	const char*						mtl_file_name,
	std::vector<MeshMaterial>&		material )
{
	std::ofstream mtl( mtl_file_name );	//	load mtl failed
	if( !mtl.is_open() ){
		#ifndef BE_QUIET
		std::cout << "writeMaterialToMtl, cannot open " << mtl_file_name << std::endl;
		#endif
		return 0;
	}

	mtl << "# mtl file generated in yzLib" << std::endl << std::endl;

	for(int i=0; i<material.size(); i++){
		mtl << "newmtl " << material[i].name << std::endl;
		mtl << "Ka " << material[i].Ka[0] << ' ' << material[i].Ka[1] << ' ' << material[i].Ka[2] << std::endl;
		mtl << "Kd " << material[i].Kd[0] << ' ' << material[i].Kd[1] << ' ' << material[i].Kd[2] << std::endl;
		mtl << "Ks " << material[i].Ks[0] << ' ' << material[i].Ks[1] << ' ' << material[i].Ks[2] << std::endl;
		mtl << "Tf " << material[i].Tf[0] << ' ' << material[i].Tf[1] << ' ' << material[i].Tf[2] << std::endl;

		if( ! material[i].map_Ka.file_name.empty() ){
			std::string filename = utils::getFileNameWithoutDirectory( 
				material[i].map_Ka.file_name.c_str() );
			mtl << "map_Ka " << filename.c_str() << std::endl;
			material[i].map_Ka.WriteTexImage();
		}
		if( ! material[i].map_Kd.file_name.empty() ){
			std::string filename = utils::getFileNameWithoutDirectory( 
				material[i].map_Kd.file_name.c_str() );
			mtl << "map_Kd " << filename.c_str() << std::endl;
			material[i].map_Kd.WriteTexImage();
		}
		if( ! material[i].map_Ks.file_name.empty() ){
			std::string filename = utils::getFileNameWithoutDirectory( 
				material[i].map_Ks.file_name.c_str() );
			mtl << "map_Ks " << filename.c_str() << std::endl;
			material[i].map_Ks.WriteTexImage();
		}

		mtl << std::endl;
	}

	mtl.close();

	return 1;
}

/**
	Read Triangle Mesh with vector, convert non triangle face to triangle

	All component are pointers, if you don't want certain component, set it to NULL.

	Extension: added material compared with old version, currently allow one material

	\param	obj_file_name	file name
	\param	v				vertex list
	\param	f				face list
	\param	vn				vertex normal list
	\param	vn_idx			vertex normal index list
	\param	vt				vertex texture coordinate list
	\param	vt_idx			vertex texture index list
	\param	material_ptr	the material of this mesh
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromObjExt(
	const char*				obj_file_name, 
	std::vector<Vec3<T>>*	v		= NULL, 
	std::vector<int3>*		f		= NULL, 
	std::vector<Vec3<T>>*	vn		= NULL,
	std::vector<int3>*		vn_idx	= NULL,
	std::vector<Com2<T>>*	vt		= NULL,
	std::vector<int3>*		vt_idx	= NULL,
	MeshMaterial*			material_ptr = NULL,
	std::string*			mtl_file_name_ptr = NULL )
{
	std::ifstream obj(obj_file_name);	//	load obj failed, do not touch old data
	if (!obj.is_open()){
#ifndef BE_QUIET
		std::cout << "error: readTriMeshFromObj, cannot open " << obj_file_name << std::endl;
#endif
		return 0;
	}

	if (v)				v->clear();
	if (f)				f->clear();
	if (vn)				vn->clear();
	if (vn_idx)			vn_idx->clear();
	if (vt)				vt->clear();
	if (vt_idx)			vt_idx->clear();
	if (material_ptr)	material_ptr->name.clear();


	int error_log_flag = 0;

	std::vector<MeshMaterial>	mtllib;
	int face_format = 0;	//	0:unkonwn;  1:v;  2:v/vt;  3:v/vt/vn;  4:v//vn;  5:v//		5 is not correct file format, but seen in some files
	std::string line;
	int line_idx = 0;
	while (std::getline(obj, line)){
		line_idx++;
		std::stringstream ss;
		std::string cmd;
		ss << line;
		ss >> cmd;

		if (cmd == "mtllib"){			//	read a material library
			std::string dir = utils::getDirectoryFromString(obj_file_name);
			std::string mtllib_file_name;
			ss >> mtllib_file_name;
			if (mtl_file_name_ptr)
				*mtl_file_name_ptr = mtllib_file_name;
			mtllib_file_name = utils::getFileNameCombineDirentory( dir.c_str(), mtllib_file_name.c_str() );
			readMaterialFromMtl( mtllib_file_name.c_str(), mtllib );
		}
		else if (cmd == "usemtl"){		//	apply a material
			ss >> cmd;
			for (std::vector<MeshMaterial>::iterator iter = mtllib.begin(); iter != mtllib.end(); iter++){
				if (cmd == iter->name && material_ptr != NULL && material_ptr->name.empty()){
					*material_ptr = *iter;
					break;
				}
			}
		}
		else if (cmd == "v"){			//	got a vertex, insert into v
			if (v){
				T xyz[3] = { 0, 0, 0 };
				int i = 0;
				while (i < 3 && ss >> xyz[i])
					i++;
				if (i < 3){
					if (!error_log_flag){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert v anyway" << std::endl;
				}
				v->push_back(Vec3<T>(xyz));
			}
		}
		else if (cmd == "vn"){			//	got a vertex normal, insert into vn
			if (vn){
				T xyz[3] = { 0, 0, 0 };
				int i = 0;
				while (i < 3 && ss >> xyz[i])
					i++;
				if (i < 3){
					if (!error_log_flag){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vn anyway" << std::endl;
				}
				vn->push_back(Vec3<T>(xyz));
			}
		}
		else if (cmd == "vt"){			//	got a vertex texture, insert into vt
			if (vt){
				T uv[2] = { 0, 0 };
				int i = 0;
				while (i < 2 && ss >> uv[i])
					i++;
				if (i < 2){
					if (!error_log_flag){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vt anyway" << std::endl;
				}
				vt->push_back(Com2<T>(uv));
			}
		}
		else if (cmd == "f"){		//	got a face, check its format, then insert into different parts
			std::string data;
			if (f){
				int fv_idx[3] = { 0, 0, 0 }, fvt_idx[3] = { 0, 0, 0 }, fvn_idx[3] = { 0, 0, 0 };
				int i = 0, count = 0;
				while (ss >> data) {	//	we don't know how much vertices this face contain
					if (face_format == 0){		//	face format unknown yet
						face_format = 1;			//	v
						size_t found = data.find("//");
						if (found != std::string::npos){
							if (found == data.length() - 2)
								face_format = 5;	//	v//, an error case, but acceptable
							else
								face_format = 4;	//	v//vn
						}
						else{
							found = data.find("/");
							if (found != std::string::npos)
								face_format = 2;	//	v/vt
							found = data.find("/", found + 1);
							if (found != std::string::npos)
								face_format = 3;	//	v/vt/vn
						}
					}
					int scan_check = 0;
					if (face_format == 1 || face_format == 5)
						scan_check = sscanf(data.c_str(), "%d", &fv_idx[i]);
					else if (face_format == 2)
						scan_check = sscanf(data.c_str(), "%d/%d", &fv_idx[i], &fvt_idx[i]);
					else if (face_format == 3)
						scan_check = sscanf(data.c_str(), "%d/%d/%d", &fv_idx[i], &fvt_idx[i], &fvn_idx[i]);
					else if (face_format == 4)
						scan_check = sscanf(data.c_str(), "%d//%d", &fv_idx[i], &fvn_idx[i]);
					if (!scan_check){
						std::cout << "line " << line_idx << " : " << line << std::endl;
						std::cout << "this vertex of this face is not read correctly due format or other error" << std::endl;
						std::cout << "don't use this vertex" << std::endl;
						face_format = 0;
					}
					else{
						if (fv_idx[i] < 0)	//	negative index is relative index
							fv_idx[i] = v->size() + fv_idx[i];
						else
							fv_idx[i] --;	//	start from 1, so we have to minus 1 to make it start from 0

						if (fvt_idx[i] < 0)
							fvt_idx[i] = vt->size() + fvt_idx[i];
						else
							fvt_idx[i] --;

						if (fvn_idx[i] < 0)
							fvn_idx[i] = vn->size() + fvn_idx[i];
						else
							fvn_idx[i] --;

						if (i == 2){
							f->push_back(int3(fv_idx));
							if (vt_idx && (face_format == 2 || face_format == 3))
								vt_idx->push_back(int3(fvt_idx));
							if (vn_idx && (face_format == 3 || face_format == 4))
								vn_idx->push_back(int3(fvn_idx));
							fv_idx[1] = fv_idx[2];
							fvt_idx[1] = fvt_idx[2];
							fvn_idx[1] = fvn_idx[2];
							fv_idx[2] = 0;
							fvt_idx[2] = 0;
							fvn_idx[2] = 0;
						}
						else
							i++;
						count++;
					}
				}
				if (count < 3){
					if (!error_log_flag){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "too few vertices, don't insert f" << std::endl;
				}
			}
		}
	}

	obj.close();

	return 1;
}


/**
	Write Triangle Mesh using vector to .obj file. vn should be one-one projection to v

	\param	obj_file_name	file name
	\param	v				vertex_list vector ptr
	\param	f				face_list vector ptr
	\param	vn				vertex_normal list vector ptr
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTextureTriMeshToObj(
	const char*					obj_file_name, 
	const std::vector<Vec3<T>>*	v	= NULL, 
	const std::vector<int3>*	f	= NULL,
	const std::vector<Vec3<T>>*	vn	= NULL,
	const std::vector<Vec2<T>>*	vt	= NULL )
{

	std::ofstream obj( obj_file_name );
	if(!obj.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: writeTriMeshToObj, cannot open " << obj_file_name << std::endl;
		#endif
		return 0;
	}

	if( vn && v ){	//	check whether vn size equal v
		if( (*vn).size() != (*v).size() ){
			std::cout << "error: writeTriMeshToObj, vn size not equal v, ignore vn" << std::endl;
			vn = NULL;
		}
	}

	if( vt && v ){	//	check whether vt size equal v
		if( (*vt).size() != (*v).size() ){
			std::cout << "error: writeTriMeshToObj, vt size not equal v, ignore vt" << std::endl;
			vt = NULL;
		}
	}

	obj << "# obj file generated in yzLib" << std::endl;

	if( v ){	//	write v
		for( int i=0; i<(*v).size(); i++ )
			obj << "v " << (*v)[i].x << " " << (*v)[i].y << " " << (*v)[i].z << std::endl;
	}

	if( vn ){
		for( int i=0; i<(*vn).size(); i++ )
			obj << "vn " << (*vn)[i].x << " " << (*vn)[i].y << " " << (*vn)[i].z << std::endl;
	}

	if( vt ){
		for( int i=0; i<(*vt).size(); i++ )
			obj << "vt " << (*vt)[i].x << " " << (*vt)[i].y << std::endl;
	}

	if( f ){	//	write f
		if( v && vn && vt ){		//	v/vt/vn
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << '/' << (*f)[i].x+1 << '/' << (*f)[i].x+1 << ' ' 
					<< (*f)[i].y+1 << '/' << (*f)[i].y+1 << '/' << (*f)[i].y+1 << ' ' 
					<< (*f)[i].z+1 << '/' << (*f)[i].z+1 << '/' << (*f)[i].z+1 << std::endl;
			}
		}
		else if( v && vn ){		//	v//vn
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << "//" << (*f)[i].x+1 << ' ' 
					<< (*f)[i].y+1 << "//" << (*f)[i].y+1 << ' '
					<< (*f)[i].z+1 << "//" << (*f)[i].z+1 << std::endl;
			}
		}
		else if( v && vt ){		//	v/vt
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << '/' << (*f)[i].x+1 << ' '
					<< (*f)[i].y+1 << '/' << (*f)[i].y+1 << ' ' 
					<< (*f)[i].z+1 << '/' << (*f)[i].z+1 << std::endl;
			}
		}
		else{				//	v
			for( int i=0; i<(*f).size(); i++ )
				obj << "f " << (*f)[i].x+1 << ' ' << (*f)[i].y+1 << ' ' << (*f)[i].z+1 << std::endl;
		}
	}

	obj.close();

	return 1;
}

/**
	Write Triangle Mesh using vector to .obj file. vn should be one-one projection to v

	\param	obj_file_name	file name
	\param	v				vertex_list vector ptr
	\param	f				face_list vector ptr
	\param	vn				vertex_normal list vector ptr
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTextureTriMeshToObjExt(
	const char*					obj_file_name, 
	const std::vector<Vec3<T>>*	v	= NULL, 
	const std::vector<int3>*	f	= NULL,
	const std::vector<Vec3<T>>*	vn	= NULL,
	const std::vector<Vec2<T>>*	vt	= NULL,
	MeshMaterial*				material_ptr = NULL )
{

	std::ofstream obj( obj_file_name );
	if(!obj.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: writeTriMeshToObj, cannot open " << obj_file_name << std::endl;
		#endif
		return 0;
	}

	if( vn && v ){	//	check whether vn size equal v
		if( (*vn).size() != (*v).size() ){
			std::cout << "error: writeTriMeshToObj, vn size not equal v, ignore vn" << std::endl;
			vn = NULL;
		}
	}

	if( vt && v ){	//	check whether vt size equal v
		if( (*vt).size() != (*v).size() ){
			std::cout << "error: writeTriMeshToObj, vt size not equal v, ignore vt" << std::endl;
			vt = NULL;
		}
	}

	obj << "# obj file generated in yzLib" << std::endl;

	if( material_ptr ){
		std::vector<MeshMaterial> mtllib;
		mtllib.push_back( *material_ptr );

		std::string mtl_file_dir = utils::getDirectoryFromString(obj_file_name);
		std::string mtl_file_name = utils::getFileNameChangeExtension(obj_file_name, "mtl");
		std::string mtl_file_name_no_dir = utils::getFileNameWithoutDirectory(mtl_file_name.c_str());

		obj << "mtllib " << mtl_file_name_no_dir.c_str() << std::endl;


		if( mtllib.front().name.empty() ){
			mtllib.front().name = "material1";
		}
		mtllib.front().map_Ka.file_dir = mtl_file_dir;
		mtllib.front().map_Kd.file_dir = mtl_file_dir;
		mtllib.front().map_Ks.file_dir = mtl_file_dir;		
		
		writeMaterialToMtl(mtl_file_name.c_str(), mtllib);
	}

	if( v ){	//	write v
		for( int i=0; i<(*v).size(); i++ )
			obj << "v " << (*v)[i].x << " " << (*v)[i].y << " " << (*v)[i].z << std::endl;
	}

	if( vn ){
		for( int i=0; i<(*vn).size(); i++ )
			obj << "vn " << (*vn)[i].x << " " << (*vn)[i].y << " " << (*vn)[i].z << std::endl;
	}

	if( vt ){
		for( int i=0; i<(*vt).size(); i++ )
			obj << "vt " << (*vt)[i].x << " " << (*vt)[i].y << std::endl;
	}

	if( f ){	//	write f
		if( material_ptr ){
			obj << "usemtl " << material_ptr->name.c_str() << std::endl;
		}

		if( v && vn && vt ){		//	v/vt/vn
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << '/' << (*f)[i].x+1 << '/' << (*f)[i].x+1 << ' ' 
					<< (*f)[i].y+1 << '/' << (*f)[i].y+1 << '/' << (*f)[i].y+1 << ' ' 
					<< (*f)[i].z+1 << '/' << (*f)[i].z+1 << '/' << (*f)[i].z+1 << std::endl;
			}
		}
		else if( v && vn ){		//	v//vn
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << "//" << (*f)[i].x+1 << ' ' 
					<< (*f)[i].y+1 << "//" << (*f)[i].y+1 << ' '
					<< (*f)[i].z+1 << "//" << (*f)[i].z+1 << std::endl;
			}
		}
		else if( v && vt ){		//	v/vt
			for( int i=0; i<(*f).size(); i++ ){
				obj << "f " 
					<< (*f)[i].x+1 << '/' << (*f)[i].x+1 << ' '
					<< (*f)[i].y+1 << '/' << (*f)[i].y+1 << ' ' 
					<< (*f)[i].z+1 << '/' << (*f)[i].z+1 << std::endl;
			}
		}
		else{				//	v
			for( int i=0; i<(*f).size(); i++ )
				obj << "f " << (*f)[i].x+1 << ' ' << (*f)[i].y+1 << ' ' << (*f)[i].z+1 << std::endl;
		}
	}

	obj.close();

	return 1;
}




}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_TEXTURE_H__