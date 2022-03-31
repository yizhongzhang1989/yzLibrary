/***********************************************************/
/**	\file
	\brief		Read and Write Mesh From File
	\details	the main entry of is readTriMeshFromFile(). and 
				writeTriMeshToFile(). These function calls other 
				functions according to the file extension provided.

				Data in vector and data in memory pointer use different
				functions in this file to solve conflict.
	\author		Yizhong Zhang
	\date		6/2/2012
*/
/***********************************************************/
#ifndef __YZ_MESH_RW_H__
#define __YZ_MESH_RW_H__

#pragma warning (disable:4996)  /* sscanf warnings. */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_utils/yz_string_utils.h"

namespace yz{	namespace geometry{

//	========================================
///@{
/**	@name Read & Write Wavefront .obj File
*/
//	========================================

/**
	Scan .obj file, get the number of each component of a mesh

	If you want to get the number of certain component, pass the address
	to the pointer, or just set NULL to it.

	\param	obj_file_name	file name
	\param	vertex_number	pointer to vertex number
	\param	face_number		pointer to face number
	\param	vn_number		pointer to vn number
	\param	vt_number		pointer to vt number
	\param	face_format		pointer to face_format, 
							0:unkonwn;  001:v;  011:v/vt;  111:v/vt/vn;  101:v//vn;  1001:v//
	\return					1: succeed,		0: failed
*/
inline int scanTriMeshFromObj(	
	const char*	obj_file_name,
	int*		vertex_number	= NULL,
	int*		face_number		= NULL,
	int*		vn_number		= NULL,
	int*		vt_number		= NULL,
	int*		face_format		= NULL
) {
	std::ifstream obj( obj_file_name );
	if(!obj.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: scanTriMeshFromObj, cannot open " << obj_file_name << std::endl;
		#endif		
		return 0;
	}
							
	int v_num	= 0;
	int f_num	= 0;
	int vn_num	= 0;
	int vt_num	= 0;
	int format	= 0;	//	0:unkonwn;  001:v;  011:v/vt;  111:v/vt/vn;  101:v//vn;  1001:v//		1001 is not correct file format, but seen in some files

	std::string line;
	std::string format_checker;
	while( std::getline(obj, line) ){
		std::stringstream ss;
		std::string cmd;
		ss << line;
		ss >> cmd;

		if( cmd == "v" )
			v_num ++;
		else if( cmd == "vn" )
			vn_num ++;
		else if( cmd == "vt" )
			vt_num ++;
		else if( cmd == "f" ){	//	convert face to triangle face
			ss >> line;
			if( !format )
				format_checker = line;
			ss >> line;
			while( ss >> line )
				f_num ++;
		}
	}

	obj.close();

	//	check format
	if (!format_checker.empty()) {
		format |= 0x01;			//	v
		size_t found = format_checker.find("//");
		if (found != std::string::npos) {
			if (found == format_checker.length() - 2)
				format |= 0x08;	//	v//, an error case
			else
				format |= 0x04;	//	v//vn
		}
		else {
			found = format_checker.find("/");
			if (found != std::string::npos)
				format |= 0x02;	//	v/vt
			found = format_checker.find("/", found + 1);
			if (found != std::string::npos)
				format |= 0x04;	//	v/vt/vn
		}
	}

	if( vertex_number )	*vertex_number	= v_num;
	if( face_number )	*face_number	= f_num;
	if( vn_number )		*vn_number		= vn_num;
	if( vt_number )		*vt_number		= vt_num;
	if( face_format )	*face_format	= format;

	return 1;
}

/**
Read Mesh with vector. Face can be any polygon

All component are pointers, if you don't want certain component, set it to NULL.

\param	obj_file_name	file name
\param	v				vertex list
\param	f				vertex id of faces
\param	f_start			start of each face vertex, f_start[0] - f_start[1] represents the first face
\param	vn				vertex normal list
\param	vn_idx			vertex normal index list, correspond to f_start
\param	vt				vertex texture coordinate list
\param	vt_idx			vertex texture index list, correspond to f_start
\return					1: succeed,		0: failed
*/
template<typename T>
inline int readMeshFromObj(
	const char*				obj_file_name,
	std::vector<Vec3<T>>*	v		= NULL,
	std::vector<int>*		f		= NULL,
	std::vector<int>*		f_start = NULL,
	std::vector<Vec3<T>>*	vn		= NULL,
	std::vector<int>*		vn_idx	= NULL,
	std::vector<Vec2<T>>*	vt		= NULL,
	std::vector<int>*		vt_idx	= NULL
) {
	std::ifstream obj(obj_file_name);	//	load obj failed, do not touch old data
	if (!obj.is_open()) {
#ifndef BE_QUIET
		std::cout << "error: readMeshFromObj, cannot open " << obj_file_name << std::endl;
#endif
		return 0;
	}

	if (v)			v->clear();
	if (f)			f->clear();
	if (f_start)	f_start->clear();
	if (vn)			vn->clear();
	if (vn_idx)		vn_idx->clear();
	if (vt)			vt->clear();
	if (vt_idx)		vt_idx->clear();

	int error_log_flag = 0;

	int face_format = 0;	//	0:unkonwn;  001:v;  011:v/vt;  111:v/vt/vn;  101:v//vn;  1001:v//		1001 is not correct file format, but seen in some files
	std::string line;
	int line_idx = 0;
	while (std::getline(obj, line)) {
		line_idx++;
		std::stringstream ss;
		std::string cmd;
		ss << line;
		ss >> cmd;

		if (cmd == "v") {			//	got a vertex, insert into v
			if (v) {
				T xyz[3] = { 0, 0, 0 };
				int i = 0;
				while (i<3 && ss >> xyz[i])
					i++;
				if (i < 3) {
					if (!error_log_flag) {
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert v anyway" << std::endl;
				}
				v->push_back(Vec3<T>(xyz));
			}
		}
		if (cmd == "vn") {			//	got a vertex normal, insert into vn
			if (vn) {
				T xyz[3] = { 0, 0, 0 };
				int i = 0;
				while (i<3 && ss >> xyz[i])
					i++;
				if (i < 3) {
					if (!error_log_flag) {
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vn anyway" << std::endl;
				}
				vn->push_back(Vec3<T>(xyz));
			}
		}
		if (cmd == "vt") {			//	got a vertex texture, insert into vt
			if (vt) {
				T uv[2] = { 0, 0 };
				int i = 0;
				while (i<2 && ss >> uv[i])
					i++;
				if (i < 2) {
					if (!error_log_flag) {
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vt anyway" << std::endl;
				}
				vt->push_back(Com2<T>(uv));
			}
		}
		else if (cmd == "f") {		//	got a face, check its format, then insert into different parts
			std::string data;
			if (f && f_start) {
				int i = 0, count = 0;
				while (ss >> data) {	//	we don't know how much vertices this face contain
					if (face_format == 0) {			//	face format unknown yet
						face_format |= 0x01;		//	v
						size_t found = data.find("//");
						if (found != std::string::npos) {
							if (found == data.length() - 2)
								face_format |= 0x08;	//	v//, an error case, but acceptable
							else
								face_format |= 0x04;	//	v//vn
						}
						else {
							found = data.find("/");
							if (found != std::string::npos)
								face_format |= 0x02;	//	v/vt
							found = data.find("/", found + 1);
							if (found != std::string::npos)
								face_format |= 0x04;	//	v/vt/vn
						}
						if (f_start->empty())
							f_start->push_back(0);	//	init f_start
					}
					int scan_check = 0;
					int fv_idx = 0, fvt_idx = 0, fvn_idx = 0;
					if (face_format == 0x01 || face_format == 0x09)	//	just v
						scan_check = sscanf(data.c_str(), "%d", &fv_idx);
					else if (face_format == 0x03)					//	v/vt
						scan_check = sscanf(data.c_str(), "%d/%d", &fv_idx, &fvt_idx);
					else if (face_format == 0x07)					//	v/vt/vn
						scan_check = sscanf(data.c_str(), "%d/%d/%d", &fv_idx, &fvt_idx, &fvn_idx);
					else if (face_format == 0x05)					//	v//vn
						scan_check = sscanf(data.c_str(), "%d//%d", &fv_idx, &fvn_idx);
					if (!scan_check) {
						std::cout << "line " << line_idx << " : " << line << std::endl;
						std::cout << "this vertex of this face is not read correctly due format or other error" << std::endl;
						std::cout << "don't use this vertex" << std::endl;
						face_format = 0;
					}
					else {
						if (fv_idx < 0)	//	negative index is relative index
							fv_idx = v->size() + fv_idx;
						else
							fv_idx--;	//	start from 1, so we have to minus 1 to make it start from 0

						if (fvt_idx < 0)
							fvt_idx = vt->size() + fvt_idx;
						else
							fvt_idx--;

						if (fvn_idx < 0)
							fvn_idx = vn->size() + fvn_idx;
						else
							fvn_idx--;

						if (face_format & 0x01)
							f->push_back(fv_idx);
						if (vt_idx && (face_format & 0x02))
							vt_idx->push_back(fvt_idx);
						if (vn_idx && (face_format & 0x04))
							vn_idx->push_back(fvn_idx);

						count++;
					}
				}
				if (count < 3) {
					if (!error_log_flag) {
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "too few vertices, don't insert f" << std::endl;
					f->resize(f_start->back());
					if (vt_idx)
						vt_idx->resize(f_start->back());
					if (vn_idx)
						vn_idx->resize(f_start->back());
				}
				else
					f_start->push_back(f->size());
			}
		}
	}

	obj.close();

	return 1;
}

/**
	Read Triangle Mesh with vector, convert non triangle face to triangle

	All component are pointers, if you don't want certain component, set it to NULL.

	\param	obj_file_name	file name
	\param	v				vertex list
	\param	f				face list
	\param	vn				vertex normal list
	\param	vn_idx			vertex normal index list
	\param	vt				vertex texture coordinate list
	\param	vt_idx			vertex texture index list
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromObj(	
	const char*				obj_file_name,
	std::vector<Vec3<T>>*	v		= NULL,
	std::vector<int3>*		f		= NULL,
	std::vector<Vec3<T>>*	vn		= NULL,
	std::vector<int3>*		vn_idx	= NULL,
	std::vector<Vec2<T>>*	vt		= NULL,
	std::vector<int3>*		vt_idx	= NULL
) {
	std::ifstream obj( obj_file_name );	//	load obj failed, do not touch old data
	if(!obj.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: readTriMeshFromObj, cannot open " << obj_file_name << std::endl;
		#endif
		return 0;
	}

	if(v)		v->clear();
	if(f)		f->clear();
	if(vn)		vn->clear();
	if(vn_idx)	vn_idx->clear();
	if(vt)		vt->clear();
	if(vt_idx)	vt_idx->clear();

	int error_log_flag = 0;

	int face_format = 0;	//	0:unkonwn;  001:v;  011:v/vt;  111:v/vt/vn;  101:v//vn;  1001:v//		1001 is not correct file format, but seen in some files
	std::string line;
	int line_idx = 0;
	while( std::getline(obj, line) ){
		line_idx ++;
		std::stringstream ss;
		std::string cmd;
		ss << line;
		ss >> cmd;

		if( cmd == "v" ){			//	got a vertex, insert into v
			if( v ){
				T xyz[3] = {0, 0, 0};
				int i = 0;
				while(i<3 && ss>>xyz[i])
					i ++;
				if( i < 3 ){
					if( !error_log_flag ){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert v anyway" << std::endl;
				}
				v->push_back( Vec3<T>(xyz) );
			}
		}
		if( cmd == "vn" ){			//	got a vertex normal, insert into vn
			if( vn ){
				T xyz[3] = {0, 0, 0};
				int i = 0;
				while(i<3 && ss>>xyz[i])
					i ++;
				if( i < 3 ){
					if( !error_log_flag ){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vn anyway" << std::endl;
				}
				vn->push_back( Vec3<T>(xyz) );
			}
		}
		if( cmd == "vt" ){			//	got a vertex texture, insert into vt
			if( vt ){
				T uv[2] = {0, 0};
				int i = 0;
				while(i<2 && ss>>uv[i])
					i ++;
				if( i < 2 ){
					if( !error_log_flag ){
						std::cout << "read " << obj_file_name << " error:" << std::endl;
						error_log_flag = 1;
					}
					std::cout << "line " << line_idx << " : " << line << std::endl;
					std::cout << "insert vt anyway" << std::endl;
				}
				vt->push_back( Com2<T>(uv) );
			}
		}
		else if( cmd == "f" ){		//	got a face, check its format, then insert into different parts
			std::string data;
			if( f ){
				int fv_idx[3] = { 0, 0, 0 }, fvt_idx[3] = { 0, 0, 0 }, fvn_idx[3] = { 0, 0, 0 };
				int i = 0, count = 0;
				while(ss >> data) {	//	we don't know how much vertices this face contain
					if( face_format == 0 ){			//	face format unknown yet
						face_format |= 0x01;		//	v
						size_t found = data.find("//");
						if( found != std::string::npos ){
							if( found == data.length()-2 )
								face_format |= 0x08;	//	v//, an error case, but acceptable
							else
								face_format |= 0x04;	//	v//vn
						}
						else{
							found = data.find("/");
							if( found != std::string::npos )
								face_format |= 0x02;	//	v/vt
							found = data.find("/", found+1);	
							if( found != std::string::npos )
								face_format |= 0x04;	//	v/vt/vn
						}
					}
					if (i == 2){	//	more than 3 vertices, we need to clear the data
						fv_idx[i] = fvt_idx[i] = fvn_idx[i] = 0;
					}
					int scan_check = 0;
					if (face_format == 0x01 || face_format == 0x09)	//	just v
						scan_check = sscanf(data.c_str(), "%d", &fv_idx[i]);
					else if (face_format == 0x03)					//	v/vt
						scan_check = sscanf(data.c_str(), "%d/%d", &fv_idx[i], &fvt_idx[i]);
					else if (face_format == 0x07)					//	v/vt/vn
						scan_check = sscanf(data.c_str(), "%d/%d/%d", &fv_idx[i], &fvt_idx[i], &fvn_idx[i]);
					else if (face_format == 0x05)					//	v//vn
						scan_check = sscanf(data.c_str(), "%d//%d", &fv_idx[i], &fvn_idx[i]);
					if( !scan_check ){
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

						if( i==2 ){
							if (face_format & 0x01)		
								f->push_back(int3(fv_idx));
							if (vt_idx && (face_format & 0x02))
								vt_idx->push_back(int3(fvt_idx));
							if (vn_idx && (face_format & 0x04))
								vn_idx->push_back( int3(fvn_idx) );
							fv_idx[1]	= fv_idx[2];
							fvt_idx[1]	= fvt_idx[2];
							fvn_idx[1]	= fvn_idx[2];
							fv_idx[2] = 0;
							fvt_idx[2] = 0;
							fvn_idx[2] = 0;
						}
						else
							i++;
						count ++;
					}
				}
				if( count<3 ){
					if( !error_log_flag ){
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
	Read Triangle Mesh, convert non triangle face to triangle

	If you don't want certain component, set its pointer to NULL.

	it is assumed that the array is big enough to hold the data

	\param	obj_file_name	file name
	\param	v				vertex ptr
	\param	v_num			pointer to vertex number
	\param	f				face ptr
	\param	f_num			pointer to face number
	\param	vn				vertex normal ptr
	\param	vn_num			pointer to vertex normal number
	\param	vn_idx			vertex normal index ptr
	\param	vt				vertex texture coordinate ptr
	\param	vt_num			pointer to vertex texture coordinate number
	\param	vt_idx			vertex texture index ptr
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromObjPtr(
	const char*	obj_file_name,
	T*			v		= NULL,
	int*		v_num	= NULL,
	int*		f		= NULL,
	int*		f_num	= NULL,
	T*			vn		= NULL,
	int*		vn_num	= NULL,
	int*		vn_idx	= NULL,
	T*			vt		= NULL,
	int*		vt_num	= NULL,
	int*		vt_idx	= NULL
) {
	std::vector<Vec3<T>>	v_list;
	std::vector<int3>		f_list;
	std::vector<Vec3<T>>	vn_list;
	std::vector<int3>		vn_idx_list;
	std::vector<Com2<T>>	vt_list;
	std::vector<int3>		vt_idx_list;

	int ret_val = readTriMeshFromObj(obj_file_name, &v_list, &f_list, &vn_list, &vn_idx_list, &vt_list, &vt_idx_list);

	if( ret_val ){
		if( v ){
			memcpy(v, &v_list[0], sizeof(Vec3<T>)*v_list.size());
			if( v_num )	*v_num = v_list.size();
		}
		if( f ){
			memcpy(f, &f_list[0], sizeof(int3)*f_list.size());
			if( f_num )	*f_num = f_list.size();
		}
		if( vn ){
			memcpy(vn, &vn_list[0], sizeof(Vec3<T>)*vn_list.size());
			if( vn_num ) *vn_num = vn_list.size();
			if( vn_idx ) memcpy(vn_idx, &vn_idx_list[0], sizeof(int3)*vn_idx_list.size());
		}
		if( vt ){
			memcpy(vt, &vt_list[0], sizeof(Com2<T>)*vt_list.size());
			if( vt_num ) *vt_num = vt_list.size();
			if( vt_idx ) memcpy(vt_idx, &vt_idx_list[0], sizeof(int3)*vt_idx_list.size());
		}
	}

	return ret_val;
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
inline int writeTriMeshToObj(
	const char*					obj_file_name,
	const std::vector<Vec3<T>>*	v = NULL,
	const std::vector<int3>*	f = NULL,
	const std::vector<Vec3<T>>*	vn = NULL,
	const std::vector<Vec2<T>>* vt = NULL,
	const std::vector<int3>*	ft = NULL
) {
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

	if (f && ft) {	//	check whether f size equal ft
		if ((*f).size() != (*ft).size()) {
			std::cout << "error: writeTriMeshToObj, f size not equal ft, ignore ft" << std::endl;
			ft = NULL;
		}
	}
	if (!ft || !vt) {	//	ft, vt must exist at the same time
		ft = NULL;
		vt = NULL;
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

	if (vt) {
		for (int i = 0; i<(*vt).size(); i++)
			obj << "vt " << (*vt)[i].x << " " << (*vt)[i].y << std::endl;
	}

	if( f ){	//	write f
		if (v && vn && vt) {		//	v/vt/vn
			for (int i = 0; i<(*f).size(); i++) {
				obj << "f "
					<< (*f)[i].x + 1 << '/' << (*ft)[i].x + 1 << '/' << (*f)[i].x + 1 << ' '
					<< (*f)[i].y + 1 << '/' << (*ft)[i].y + 1 << '/' << (*f)[i].y + 1 << ' '
					<< (*f)[i].z + 1 << '/' << (*ft)[i].z + 1 << '/' << (*f)[i].z + 1 << std::endl;
			}
		}
		else if (v && vn) {		//	v//vn
			for (int i = 0; i<(*f).size(); i++) {
				obj << "f "
					<< (*f)[i].x + 1 << "//" << (*f)[i].x + 1 << ' '
					<< (*f)[i].y + 1 << "//" << (*f)[i].y + 1 << ' '
					<< (*f)[i].z + 1 << "//" << (*f)[i].z + 1 << std::endl;
			}
		}
		else if (v && vt) {		//	v/vt
			for (int i = 0; i<(*f).size(); i++) {
				obj << "f "
					<< (*f)[i].x + 1 << '/' << (*ft)[i].x + 1 << ' '
					<< (*f)[i].y + 1 << '/' << (*ft)[i].y + 1 << ' '
					<< (*f)[i].z + 1 << '/' << (*ft)[i].z + 1 << std::endl;
			}
		}
		else {				//	v
			for (int i = 0; i<(*f).size(); i++)
				obj << "f " << (*f)[i].x + 1 << ' ' << (*f)[i].y + 1 << ' ' << (*f)[i].z + 1 << std::endl;
		}
	}

	obj.close();

	return 1;
}

/**
	Write Triangle Mesh to .obj file. vn should be one-one projection to v

	\param	obj_file_name	file name
	\param	v				vertex list ptr
	\param	v_number		vertex number
	\param	f				face list ptr
	\param	f_number		face number
	\param	vn				vertex normal list ptr
	\param	vn_number		vertex normal number, if not same as v_number, it will be ignored
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTriMeshToObjPtr(
	const char*	obj_file_name,
	const T*	v			= NULL,
	int			v_number	= 0,
	const int*	f			= NULL,
	int			f_number	= 0,
	const T*	vn			= NULL,
	int			vn_number	= 0
) {
	std::ofstream obj(obj_file_name);
	if (!obj.is_open()) {
#ifndef BE_QUIET
		std::cout << "error: writeTriMeshToObjPtr, cannot open " << obj_file_name << std::endl;
#endif
		return 0;
	}

	if( vn && v ){	//	check whether vn size equal v
		if( vn_number != v_number ){
			std::cout << "error: writeTriMeshToObjPtr, vn size not equal v, ignore vn" << std::endl;
			vn = NULL;
			vn_number = 0;
		}
	}

	obj << "# obj file generated in yzLib" << std::endl;

	if( v ){	//	write v
		for( int i=0; i<v_number; i++ )
			obj << "v " << v[i*3] << " " << v[i*3+1] << " " << v[i*3+2] << std::endl;
	}

	if( vn ){
		for( int i=0; i<vn_number; i++ )
			obj << "vn " << vn[i*3] << " " << vn[i*3+1] << " " << vn[i*3+2] << std::endl;
	}

	if( f ){	//	write f
		if( v && vn ){		//	v//vn
			for( int i=0; i<f_number; i++ )
				obj << "f " 
					<< f[i*3]+1 << "//" << f[i*3]+1 << " " 
					<< f[i*3+1]+1 << "//" << f[i*3+1]+1 << " " 
					<< f[i*3+2]+1 << "//" << f[i*3+2]+1 << std::endl;
		}
		else{				//	v
			for( int i=0; i<f_number; i++ )
				obj << "f " << f[i*3]+1 << " " << f[i*3+1]+1 << " " << f[i*3+2]+1 << std::endl;
		}
	}

	obj.close();

	return 1;
}

/**
	Read Triangle Mesh 2D with vector, convert non triangle face to triangle

	All component are pointers, if you don't want certain component, set it to NULL.

	\param	obj_file_name		file name
	\param	v					vertex list
	\param	f					face list
	\param	ignore_dimension	the dimension to ignore, 'x', 'y', 'z'
	\return						whether read succeed
*/
template<typename T>
inline int readTriMeshFromObj2D(
	const char*				obj_file_name,
	std::vector<Vec2<T>>*	v = NULL,
	std::vector<int3>*		f = NULL,
	char					ignore_dimension = 'z'
) {
	std::vector<Vec3<T>> tmp_v;
	if (!readTriMeshFromObj(obj_file_name, &tmp_v, f))
		return 0;

	v->clear();
	for( int i=0; i<tmp_v.size(); i++ ){
		Vec2<T> v2d(tmp_v[i].x, tmp_v[i].y);
		if( ignore_dimension == 'x' || ignore_dimension == 'X' )
			v2d.x = tmp_v[i].y, v2d.y = tmp_v[i].z;
		else if( ignore_dimension == 'y' || ignore_dimension == 'Y' )
			v2d.y = tmp_v[i].z;
		v->push_back(v2d);
	}
	return 1;
}

/**
	Write Triangle Mesh 2D using vector to .obj file. Only v, f.

	\param	obj_file_name			file name
	\param	v						vertex_list vector ptr
	\param	f						face_list vector ptr
	\param	additional_dimension	the additional dimension to be added to 0, 'x', 'y', 'z'
	\return							1: succeed,		0: failed
*/
template<typename T>
inline int writeTriMeshToObj2D(
	const char*					obj_file_name,
	const std::vector<Vec2<T>>*	v					= NULL,
	const std::vector<int3>*	f					= NULL,
	char						additional_dimension = 'z'
) {
	std::vector<Vec3<T>> tmp_v;

	for (int i = 0; i < v->size(); i++) {
		Vec3<T> v3d((*v)[i].x, (*v)[i].y, 0);
		if( additional_dimension == 'x' || additional_dimension == 'X' )
			v3d.z = v3d.y, v3d.y = v3d.x, v3d.x = 0;
		else if( additional_dimension == 'y' || additional_dimension == 'Y' )
			v3d.z = v3d.y, v3d.y = 0;
		tmp_v.push_back(v3d);
	}

	return writeTriMeshToObj(obj_file_name, &tmp_v, f);
}

///@}


//	========================================
///@{
/**	@name Read & Write Wavefront .stl File
*/
//	========================================

/**
	Write Triangle Mesh using vector to .stl file

	\param	stl_file_name	file name
	\param	v				vertex_list vector ptr
	\param	f				face_list vector ptr
	\param	fn				face_normal list vector ptr
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTriMeshToBinaryStl(
	const char*					stl_file_name,
	const std::vector<Vec3<T>>*	v = NULL,
	const std::vector<int3>*	f = NULL,
	const std::vector<Vec3<T>>*	fn = NULL
) {
	std::ofstream stl( stl_file_name, std::ios::binary );
	if(!stl.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: writeTriMeshToBinaryStl, cannot open " << stl_file_name << std::endl;
		#endif
		return 0;
	}

	if( fn && f ){	//	check whether vn size equal v
		if( (*fn).size() != (*f).size() ){
			std::cout << "error: writeTriMeshToBinaryStl, fn size not equal f, ignore fn" << std::endl;
			fn = NULL;
		}
	}

	if( !f || !v ){
		std::cout << "error: writeTriMeshToBinaryStl, no vertex or face" << std::endl;
		return 0;
	}

	//	write header
	char header[80] = "\0\0\0\0\0";
	memset(header, 0, 80);
	stl.write( header, 80 );

	//	write face number
	int f_num = (*f).size();
	stl.write( (char*)&f_num, 4 );

	if( f_num <= 0 )
		return 1;

	//	create buffer for write
	const int buffer_max_capacity = 10000;
	int buf_cap = myMin(f_num, buffer_max_capacity);
	char* buffer = new char[buf_cap*50];
	if( ! buffer ){
		std::cout << "error: writeTriMeshToBinaryStl, create buffer failed" << std::endl;
		return 0;
	}

	//	write each face
	int buf_count = 0;
	for(int i=0; i<(*f).size(); i++){
		char* data = buffer + buf_count*50;
		//	normal
		if( fn )
			*(Vec3f*)&data[0] = (*fn)[i];
		else
			*(Vec3f*)&data[0] = yz::Vec3f(0,0,0);
		//	each vertex
		*(Vec3f*)&data[12] = (*v)[ (*f)[i].x ];
		*(Vec3f*)&data[24] = (*v)[ (*f)[i].y ];
		*(Vec3f*)&data[36] = (*v)[ (*f)[i].z ];

		//	the buffer is full, or we have reached the last face, write and clear buffer
		buf_count ++;
		if( buf_count == buf_cap || i+1==(*f).size() ){
			stl.write( buffer, 50*buf_count );
			buf_count = 0;
		}
	}

	delete[] buffer;

	stl.close();

	return 1;
}


/**
	Read triangle mesh from a binary stl file

	All component are pointers, if you don't want certain component, set it to NULL.

	\param	stl_file_name	file name
	\param	v				vertex list
	\param	f				face list
	\param	fn				face normal list
	\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromBinaryStl(
	const char*				stl_file_name,
	std::vector<Vec3<T>>*	v = NULL,
	std::vector<int3>*		f = NULL,
	std::vector<Vec3<T>>*	fn = NULL
) {
	std::ifstream stl( stl_file_name, std::ios::binary );	//	load stl failed, do not touch old data
	if(!stl.is_open()){
		#ifndef BE_QUIET
			std::cout << "error: readTriMeshFromBinaryStl, cannot open " << stl_file_name << std::endl;
		#endif
		return 0;
	}

	//	read header
	char header[80];
	stl.read( header, 80 );
	//if( header[0]=='s' && header[1]=='o' && header[2]=='l' && header[3]=='i' && header[4]=='d' ){
	//	std::cout << "error: readTriMeshFromBinaryStl, ascii stl file" << std::endl;
	//	stl.close();
	//	return 0;
	//}
		
	if(v)		v->clear();
	if(f)		f->clear();
	if(fn)		fn->clear();

	//	read face number
	int f_num;
	stl.read( (char*)&f_num, 4 );

	if( f_num <= 0 ){
		stl.close();
		return 1;
	}

	//	setup space
	if(v)		v->resize(f_num*3);
	if(f)		f->resize(f_num);
	if(fn)		fn->resize(f_num);

	//	create buffer for read
	const int buffer_max_capacity = 10000;
	int buf_cap = myMin(f_num, buffer_max_capacity);
	char* buffer = new char[buf_cap*50];
	if( ! buffer ){
		std::cout << "error: readTriMeshFromBinaryStl, create buffer failed" << std::endl;
		stl.close();
		return 0;
	}

	//	read each face
	int buf_i = buf_cap;
	int buf_end = buf_cap;
	for(int i=0; i<f_num; i++, buf_i++){
		//	read data to buffer
		if( buf_i == buf_cap ){
			buf_i = 0;
			int buf_count = myMin(f_num-i, buf_cap);
			stl.read( buffer, 50*buf_count );
		}

		char* data = buffer + buf_i*50;
		//	normal
		if( fn )
			(*fn)[i] = *(Vec3f*)&data[0];

		//	each vertex
		if( v ){
			(*v)[i*3  ] = *(Vec3f*)&data[12];
			(*v)[i*3+1] = *(Vec3f*)&data[24];
			(*v)[i*3+2] = *(Vec3f*)&data[36];
		}

		if( f ){
			(*f)[i].x = i*3;
			(*f)[i].y = i*3+1;
			(*f)[i].z = i*3+2;
		}
	}

	delete[] buffer;

	stl.close();

	return 1;
}

///@}

//	========================================
///@{
/**	@name Read & Write Mesh File of Any Format
*/
//	========================================

/**
Read Triangle Mesh From File of Any Format

Check the file format, then read the file according to the file format

Currently supported format:	.obj, .stl

\param	file_name		file name
\param	v				vertex list
\param	f				face list
\param	vn				vertex normal list
\param	vn_idx			vertex normal index list
\param	vt				vertex texture coordinate list
\param	vt_idx			vertex texture index list
\param	fn				face normal list
\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromFile(
	const char*				file_name,
	std::vector<Vec3<T>>*	v		= NULL,
	std::vector<int3>*		f		= NULL,
	std::vector<Vec3<T>>*	vn		= NULL,
	std::vector<int3>*		vn_idx	= NULL,
	std::vector<Vec2<T>>*	vt		= NULL,
	std::vector<int3>*		vt_idx	= NULL,
	std::vector<Vec3<T>>*	fn		= NULL
) {
	std::string extension = utils::getFileExtensionFromString(file_name);
	std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

	//	call read functions according to file extension
	if (extension == "obj")
		return readTriMeshFromObj(file_name, v, f, vn, vn_idx, vt, vt_idx);
	else if (extension == "stl")
		return readTriMeshFromBinaryStl(file_name, v, f, fn);
	else
		std::cout << "Failed to read mesh " << file_name << ", unsupported file format: ." << extension << std::endl;

	return 0;
}

/**
Read Triangle Mesh From File of Any Format

Check the file format, then read the file according to the file format.
Memory pointer are used to store data directly.

Currently supported format:	.obj

\param	file_name		file name
\param	v				vertex ptr
\param	v_num			pointer to vertex number
\param	f				face ptr
\param	f_num			pointer to face number
\param	vn				vertex normal ptr
\param	vn_num			pointer to vertex normal number
\param	vn_idx			vertex normal index ptr
\param	vt				vertex texture coordinate ptr
\param	vt_num			pointer to vertex texture coordinate number
\param	vt_idx			vertex texture index ptr
\return					1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromFilePtr(
	const char*	file_name,
	T*			v		= NULL,
	int*		v_num	= NULL,
	int*		f		= NULL,
	int*		f_num	= NULL,
	T*			vn		= NULL,
	int*		vn_num	= NULL,
	int*		vn_idx	= NULL,
	T*			vt		= NULL,
	int*		vt_num	= NULL,
	int*		vt_idx	= NULL
) {
	std::string extension = utils::getFileExtensionFromString(file_name);
	std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

	//	call read functions according to file extension
	if (extension == "obj")
		return readTriMeshFromObjPtr(file_name, v, v_num, f, f_num, vn, vn_num, vn_idx, vt, vt_num, vt_idx);
	else
		std::cout << "Failed to read mesh " << file_name << ", unsupported file format: ." << extension << std::endl;

	return 0;
}

/**
Write Triangle Mesh to File.

choose file format accroding to file extension

\param	file_name		file name
\param	v				vertex_list vector ptr
\param	f				face_list vector ptr
\param	vn				vertex_normal_list vector ptr
\param	fn				face_normal vector ptr, each vertex has only one normal
\param	vt				texture coordinate
\param	ft				face texture coordinate index
\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTriMeshToFile(
	const char*					file_name,
	const std::vector<Vec3<T>>*	v = NULL,
	const std::vector<int3>*	f = NULL,
	const std::vector<Vec3<T>>*	vn = NULL,
	const std::vector<Vec3<T>>*	fn = NULL,
	const std::vector<Vec2<T>>* vt = NULL,
	const std::vector<int3>*	ft = NULL
) {
	std::string extension = utils::getFileExtensionFromString(file_name);
	std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

	//	call write functions according to file extension
	if (extension == "obj")
		return writeTriMeshToObj(file_name, v, f, vn, vt, ft);
	else if (extension == "stl")
		return writeTriMeshToBinaryStl(file_name, v, f, fn);
	else
		std::cout << "Failed to write mesh " << file_name << ", unsupported file format: ." << extension << std::endl;

	return 0;
}

/**
Write Triangle Mesh to File.

choose file format accroding to file extension.
This function takes memory pointer as argument.

\param	file_name		file name
\param	v				vertex list ptr
\param	v_number		vertex number
\param	f				face list ptr
\param	f_number		face number
\return					1: succeed,		0: failed
*/
template<typename T>
inline int writeTriMeshToFilePtr(
	const char*	file_name,
	const T*	v			= NULL,
	int			v_number	= 0,
	const int*	f			= NULL,
	int			f_number	= 0
) {
	std::string extension = utils::getFileExtensionFromString(file_name);
	std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

	//	call write functions according to file extension
	if (extension == "obj")
		return writeTriMeshToObjPtr(file_name, v, v_number, f, f_number);
	else {
#ifndef BE_QUIET
		std::cout << "Failed to write mesh " << file_name << ", unsupported file format: ." << extension << std::endl;
#endif		
	}

	return 0;
}

/**
Read Triangle Mesh From File of Any Format

Check the file format, then read the file according to the file format

Currently supported format:	.obj

\param	file_name			file name
\param	v					vertex list
\param	f					face list
\param	ignore_dimension	the dimension to ignore, 'x', 'y', 'z'
\return						1: succeed,		0: failed
*/
template<typename T>
inline int readTriMeshFromFile2D(
	const char*				file_name,
	std::vector<Vec2<T>>*	v					= NULL,
	std::vector<int3>*		f					= NULL,
	char					ignore_dimension	= 'z'
) {
	std::string extension = utils::getFileExtensionFromString(file_name);
	std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

	//	call read functions according to file extension
	if (extension == "obj")
		return readTriMeshFromObj2D(file_name, v, f, ignore_dimension);
	else
		std::cout << "Failed to read mesh 2D " << file_name << ", unsupported file format: ." << extension << std::endl;

	return 0;
}
///@}


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_RW_H__