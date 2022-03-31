/***********************************************************/
/**	\file
	\brief		Mesh UV
	\details	create UV
	\author		Yizhong Zhang
	\date		10/26/2014
*/
/***********************************************************/
#ifndef __YZ_MESH_UV_H__
#define __YZ_MESH_UV_H__

#include <iostream>
#include <vector>
#include <queue>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_geometry/yz_mesh_utils.h"
#include "yzLib/yz_geometry/yz_polygon/yz_polygon_deform.h"
#include "yzLib/yz_geometry/yz_mesh_deform.h"

namespace yz{	namespace geometry{

/**
	Unwrap a mesh patch 

	create UV for single patch mesh (only one close loop boundary)
*/
template <class T>
class MeshPatchUVunwrapper :
	public BaseTriMesh<T>,
	public MeshTextureCoordinate<T>
{
public:
	using BaseTriMesh<T>::vertex;
	using BaseTriMesh<T>::face;
	using MeshTextureCoordinate<T>::tex_coord;

public:
	/**
		Unwrap the patch mesh

		\return		whether wrap succeed
	*/
	int Unwrap(int max_iterations = 100) {
		//	check whether mesh exist
		if (vertex.empty()) {
			std::cout << "error: PatchUVunwrapper::Unwrap(), mesh not read yet" << std::endl;
			return 0;
		}

		//	check whether the mesh is a patch (has a single loop boundary)
		if (!CreateBoundary())
			return 0;

		//	flatten the mesh into UV
		Flatten(max_iterations);

		//	refine the UV
		//Refine(max_iterations*10);

		return 1;
	}

	/**
		draw patch UV
	*/
	void Draw2D() {
#ifdef YZ_glut_h
		glColor3f(0, 0, 1);
		yz::opengl::drawMeshEdgeFromFace2D(tex_coord, face);
#else
		std::cout << "error:MeshPatchUVunwrapper, glut.h must be included to enable Draw2D()" << std::endl;
#endif
	}

protected:
	/**
		create boundary and deformation parameters

		\return		whether the patch has only 1 boundary curve
	*/
	int CreateBoundary() {
		//	get boundary curve
		std::vector<std::vector<int>>	curve_vertex;
		std::vector<int>				loop_flag;
		getBoundaryCurveVerticesFromMesh(curve_vertex, loop_flag, face);

		//	check whether the mesh is a patch
		if (curve_vertex.size() != 1) {	//	the mesh should have only one curve
			std::cout << "error: MeshPatchUVunwrapper::CreateBoundary(), the mesh has " << curve_vertex.size() << " boundary curves" << std::endl;
			return 0;
		}
		if (loop_flag[0] != 1) {	//	the single curve should be a loop
			std::cout << "error: MeshPatchUVunwrapper::CreateBoundary(), the boundary curve is not a loop" << std::endl;
			return 0;
		}

		boundary_curve_vertex.swap(curve_vertex[0]);

		//	label boundary vertices
		boundary_vertex_flag.clear();
		boundary_vertex_flag.resize(vertex.size(), 0);
		for (int i = 0; i < boundary_curve_vertex.size(); i++) {
			boundary_vertex_flag[boundary_curve_vertex[i]] = 1;
		}

		//	calculate target length
		target_length.resize(boundary_curve_vertex.size());
		for (int i = 0; i < target_length.size(); i++) {
			int vid0 = boundary_curve_vertex[i];
			int vid1 = boundary_curve_vertex[(i + 1) % boundary_curve_vertex.size()];
			target_length[i] = (vertex[vid0] - vertex[vid1]).Length();
		}

		//	create a lookup table
		std::unordered_map<int, int> boundary_loop_idx;
		boundary_loop_idx.reserve(boundary_curve_vertex.size());
		for (int i = 0; i < boundary_curve_vertex.size(); i++) {
			boundary_loop_idx.insert(std::pair<int, int>(boundary_curve_vertex[i], i));
		}

		//	calculate target angle
		target_angle.clear();
		target_angle.resize(boundary_curve_vertex.size(), 0);
		for (int f = 0; f < face.size(); f++) {
			for (int i = 0; i < 3; i++) {
				int vid0 = face[f][i];
				if (!boundary_vertex_flag[vid0])	//	this vertex is not boundary, skip
					continue;
				int vid1 = face[f][(i + 1) % 3];
				int vid2 = face[f][(i + 2) % 3];
				Vec3<T> r1 = vertex[vid1] - vertex[vid0];
				Vec3<T> r2 = vertex[vid2] - vertex[vid0];
				T angle = angleRadBetweenVectors(r1, r2);

				auto iter = boundary_loop_idx.find(vid0);
				target_angle[iter->second] += angle;
			}
		}

		return 1;
	}

	/**
		flatten boundary, then laplacian the inner
	*/
	void Flatten(int iterations) {
		//	flatten boundary
		polygon::PolygonDeformer<T> poly_deformer;
		poly_deformer.SetTarget(target_length, target_angle);
		poly_deformer.Deform();

		tex_coord.resize(vertex.size());
		for (int i = 0; i < boundary_curve_vertex.size(); i++) {
			int vid = boundary_curve_vertex[i];
			tex_coord[vid] = poly_deformer.v[i];
		}

		//	set inner coord, using laplacian
		std::unordered_multimap<int, int> vv_hash;
		createVVHashFromTriFace(vv_hash, face);	
		while (iterations--) {
			//	copy original tex_coord
			std::vector<Vec2<T>>	tmp_tex_coord;
			tmp_tex_coord.swap(tex_coord);
			tex_coord.resize(tmp_tex_coord.size());

			//	for each vertex, set its position as average
			for (int i = 0; i < tex_coord.size(); i++) {
				if (boundary_vertex_flag[i])	//	boundary vertex, fix it
					tex_coord[i] = tmp_tex_coord[i];
				else {							//	inner vertex, set as average
					Vec2<T> avg_coord;
					int count = 0;
					auto range = vv_hash.equal_range(i);
					for (auto iter = range.first; iter != range.second; iter++) {
						int nvid = iter->second;
						avg_coord += tmp_tex_coord[nvid];
						count++;
					}
					if (count)
						avg_coord /= count;
					tex_coord[i] = avg_coord;
				}
			}

			//	check whether laplacian should terminate
		}
	}

	void Refine(int iterations) {
		MeshPatchDeformer2D<T> deformer;
		deformer.vertex = tex_coord;
		deformer.face = face;

		for (int fid = 0; fid < face.size(); fid++) {
			for (int i = 0; i < 3; i++) {
				int2 edge(face[fid][i], face[fid][(i + 1) % 3]);
				T edge_length = (vertex[edge.x] - vertex[edge.y]).Length();
				deformer.SetTargetEdgeLength(edge, edge_length);
			}
		}

		deformer.Deform(iterations);

		tex_coord = deformer.vertex;
	}

protected:
	std::vector<int>		boundary_curve_vertex;	//	vertices on the boundary curve
	std::vector<int>		boundary_vertex_flag;	//	label whether each vertex is boundary

	//	polygon deformation
	std::vector<T>			target_length;
	std::vector<T>			target_angle;
};


/**
	UV unwrapper of triangle mesh

	create patch for the triangle mesh seperated by seams, then flatten each patch as UV

	How to use:	\n
	1, call ReadMeshFromFile(file_name) to read the mesh, existing UV will be read (call ClearUV to clear) \n
	2, call SetSeam() to set seam \n
	3, call Unwrap() to unwrap the mesh
	4, call TranslatePatchUV(), RotatePatchUV() to refine the layout of UV patches
	5, call WriteMeshToFile() to write the result

	display functions are provided for manual layout each patch
*/
template <class T>
class UVunwrapper : 
	public TriMesh<T>,
	public MeshTextureCoordinate<T>,
	public MeshTextureFace
{
public:
	using TriMesh<T>::vertex;
	using TriMesh<T>::face;
	using TriMesh<T>::vertex_normal;
	using TriMesh<T>::face_normal;
	using MeshTextureCoordinate<T>::tex_coord;

public:
	/**
	read the mesh
	*/
	int ReadMeshFromFile(const char* file_name) {
		Clear();

		int succ = readTriMeshFromFile(
			file_name,
			&vertex,
			&face,
			(std::vector<Vec3<T>>*)NULL,
			(std::vector<int3>*)NULL,
			&tex_coord,
			&tex_face
		);
		if (!succ)
			return 0;

		this->CalculateNormals();

		//	initialize topology and seam
		CreateTopology();
		SetUVFromMesh();

		//	texture read from the file, create label
		if (!tex_face.empty())
			CreateFaceLabel();

		return 1;
	}

	/**
	set multiple edges to be seam

	\param	seam_edge	list of seam edges
	\return				number of seam edges added
	*/
	int SetSeam(const std::vector<int2>& seam_edge) {
		seam_unordered_set.reserve(seam_unordered_set.size() + seam_edge.size());
		seam.reserve(seam.size() + seam_edge.size());

		int count = 0;
		for (unsigned int i = 0; i != seam_edge.size(); i++) {
			count += SetSeam(seam_edge[i]);
		}

		return count;
	}

	/**
	set a single edge to be seam

	\param	seam_edge	the edge to be labeled
	\return				whether this edge is labeled
	*/
	int SetSeam(int2 seam_edge) {
		if (seam_edge.x > seam_edge.y)
			mySwap(seam_edge.x, seam_edge.y);

		if (seam_unordered_set.find(seam_edge) != seam_unordered_set.end())	//	this edge is already labeled as seam
			return 0;
		if (ef_multimap.find(seam_edge) == ef_multimap.end())	//	this edge doesn't exist
			return 0;
		seam.push_back(seam_edge);
		seam_unordered_set.insert(seam_edge);
		return 1;
	}

	/**
	remove a seam edge, but don't update patch

	\param	seam_edge	the seam edge to unlabel
	\return				whether this seam edge is unlabeled
	*/
	int RemoveSeam(int2 seam_edge) {
		if (seam_edge.x > seam_edge.y)
			mySwap(seam_edge.x, seam_edge.y);

		auto hash_iter = seam_unordered_set.find(seam_edge);
		if (hash_iter == seam_unordered_set.end())	//	this edge is not a seam
			return 0;

		//	erase the hash
		seam_unordered_set.erase(hash_iter);

		//	erase the vector
		auto seam_iter = std::find(seam.begin(), seam.end(), seam_edge);
		if (seam_iter == seam.end()) {	//	this case should not happen
			std::cout << "error: UVunwrapper::RemoveSeam, seam in hash, but not in vector" << std::endl;
			return 0;
		}
		*seam_iter = seam.back();
		seam.pop_back();

		return 1;
	}

	/**
	Clear existing UV
	*/
	void ClearUV() {
		tex_coord.clear();
		tex_face.clear();

		face_label.clear();
		patch_face.clear();

		tex_v_mapping.clear();
		v_tex_multimap.clear();
	}

	/**
	Unwrap the mesh
	*/
	int Unwrap() {
		ClearUV();
		CreateFaceLabel();
		CreateTexFace();
		int succ = CreateTexCoord();

		if (!succ) {
			std::cout << "error: UVunwrapper::Unwrap, unwrap failed, seam need to be added" << std::endl;
			return 0;
		}

		return 1;
	}

	/**
	Write the mesh to file, with UV
	*/
	int WriteMeshToFile(const char* file_name) {
		return writeTriMeshToFile(
			file_name,
			&vertex,
			&face,
			(std::vector<Vec3<T>>*)NULL,
			(std::vector<Vec3<T>>*)NULL,
			tex_coord.empty() ? (std::vector<Vec2<T>>*)NULL : &tex_coord,
			tex_face.empty() ? (std::vector<int3>*)NULL : &tex_face
		);
	}

	/**
	offset a single patch uv
	*/
	void TranslatePatchUV(int label_id, Vec2<T> offset) {
		std::unordered_set<int> tex_hash;
		tex_hash.reserve(tex_coord.size());

		for (int i = 0; i < patch_face[label_id].size(); i++) {
			int fid = patch_face[label_id][i];
			for (int j = 0; j < 3; j++) {
				int tvid = tex_face[fid][j];
				if (tex_hash.find(tvid) != tex_hash.end())
					continue;
				tex_coord[tvid] += offset;
				tex_hash.insert(tvid);
			}
		}
	}

	/**
	rotate a single patch uv
	*/
	void RotatePatchUV(int label_id, T angle_deg) {
		//	calculate patch center
		Vec2<T> cen(0, 0);
		for (int i = 0; i < patch_face[label_id].size(); i++) {
			int fid = patch_face[label_id][i];
			for (int j = 0; j < 3; j++) {
				int tvid = tex_face[fid][j];
				cen += tex_coord[tvid];
			}
		}
		cen /= (patch_face[label_id].size() * 3);

		//	rotate
		std::unordered_set<int> tex_hash;
		tex_hash.reserve(tex_coord.size());

		for (int i = 0; i < patch_face[label_id].size(); i++) {
			int fid = patch_face[label_id][i];
			for (int j = 0; j < 3; j++) {
				int tvid = tex_face[fid][j];
				if (tex_hash.find(tvid) != tex_hash.end())
					continue;

				Vec2<T> r = tex_coord[tvid] - cen;
				r.SetRotateDeg(angle_deg);
				tex_coord[tvid] = cen + r;
				tex_hash.insert(tvid);
			}
		}
	}

	/**
	scale a single patch uv
	*/
	void ScalePatchUV(int label_id, T scale) {
		//	calculate patch center
		Vec2<T> cen(0, 0);
		for (int i = 0; i < patch_face[label_id].size(); i++) {
			int fid = patch_face[label_id][i];
			for (int j = 0; j < 3; j++) {
				int tvid = tex_face[fid][j];
				cen += tex_coord[tvid];
			}
		}
		cen /= (patch_face[label_id].size() * 3);

		//	scale respect to center
		std::unordered_set<int> tex_hash;
		tex_hash.reserve(tex_coord.size());

		for (int i = 0; i < patch_face[label_id].size(); i++) {
			int fid = patch_face[label_id][i];
			for (int j = 0; j < 3; j++) {
				int tvid = tex_face[fid][j];
				if (tex_hash.find(tvid) != tex_hash.end())
					continue;

				Vec2<T> r = tex_coord[tvid] - cen;
				r *= scale;
				tex_coord[tvid] = cen + r;
				tex_hash.insert(tvid);
			}
		}
	}

	/**
	number of created patches
	*/
	int PatchNumber() {
		return patch_face.size();
	}

public:		//	display related functions
	void Draw3D(double seam_radius = 0.01) {
#ifdef YZ_glut_h
		//	draw the mesh
		if (face_label.empty()) {	//	if the face is not labeled, just draw the mesh
			glColor3f(1, 1, 1);
			opengl::drawSmoothShadingTriMesh(vertex, face, vertex_normal);
			glColor3f(0, 0, 1);
			opengl::drawMeshEdgeFromFace(vertex, face, vertex_normal);
		}
		else {	//	if the face already labeled, draw mesh face with face label as color
			std::vector<yz::uchar3> face_color;
			face_color.resize(face_label.size());
			for (int i = 0; i < face_color.size(); i++) {
				opengl::getSequentialDisplayColor(&face_color[i].x, face_label[i]);
			}

			glColor3f(1, 1, 1);
			glDisable(GL_LIGHTING);
			opengl::drawFlatColorTriMesh(vertex, face, face_color);
			glEnable(GL_LIGHTING);

			//	draw mesh edge
			glColor3f(1, 1, 1);
			opengl::drawMeshEdgeFromFace(vertex, face, vertex_normal);
		}

		//	draw seam
		glColor3f(1, 0, 1);
		for (int i = 0; i < seam.size(); i++)
			yz::opengl::drawCylinder(vertex[seam[i].x], vertex[seam[i].y], seam_radius);
#else
		std::cout << "error:UVunwrapper, glut.h must be included to enable Draw3D()" << std::endl;
#endif
	}

	void Draw2D() {
#ifdef YZ_glut_h
		if (tex_coord.empty() || tex_face.empty())
			return;

		//	draw the mesh
		if (face_label.empty()) {	//	if the face is not labeled, just draw the mesh
			glColor3f(0, 0, 1);
			opengl::drawMeshEdgeFromFace2D(tex_coord, tex_face);
		}
		else {	//	if the face already labeled, draw mesh face with face label as color
			std::vector<yz::uchar3> face_color;
			face_color.resize(face_label.size());
			for (int i = 0; i < face_color.size(); i++) {
				opengl::getSequentialDisplayColor(&face_color[i].x, face_label[i]);
			}

			glColor3f(1, 1, 1);
			opengl::drawFlatColorTriMesh2D(tex_coord, tex_face, face_color);

			//	draw mesh edge
			glColor3f(0, 0, 0);
			opengl::drawMeshEdgeFromFace2D(tex_coord, tex_face);
		}
#else
		std::cout << "error:UVunwrapper, glut.h must be included to enable Draw2D()" << std::endl;
#endif
	}

	void PickingDraw2D() {
#ifdef YZ_glut_h
		if (tex_coord.empty() || tex_face.empty() || face_label.empty())
			return;

		glBegin(GL_TRIANGLES);
		for (int fid = 0; fid < tex_face.size(); fid++) {
			opengl::setPickingIndex(face_label[fid]);
			glVertex2d(tex_coord[tex_face[fid].x].x, tex_coord[tex_face[fid].x].y);
			glVertex2d(tex_coord[tex_face[fid].y].x, tex_coord[tex_face[fid].y].y);
			glVertex2d(tex_coord[tex_face[fid].z].x, tex_coord[tex_face[fid].z].y);
		}
		glEnd();
#else
		std::cout << "error:UVunwrapper, glut.h must be included to enable PickingDraw2D()" << std::endl;
#endif
	}

protected:
	/**
	Clear all data
	*/
	void Clear() {
		//	mesh
		vertex.clear();
		face.clear();
		vertex_normal.clear();
		face_normal.clear();
		tex_coord.clear();
		tex_face.clear();

		//	seam
		seam.clear();

		//	labeling
		face_label.clear();
		patch_face.clear();

		//	topology
		vf_multimap.clear();
		ef_multimap.clear();
		ff.clear();
		ff_start.clear();

		//	lookup table
		seam_unordered_set.clear();
		tex_v_mapping.clear();
		v_tex_multimap.clear();
	}

	/**
	create topology of the mesh
	*/
	void CreateTopology() {
		vf_multimap.clear();
		ef_multimap.clear();
		ff.clear();
		ff_start.clear();

		ff.reserve(face.size() * 3);
		ff_start.reserve(face.size() + 1);
		ff_start.push_back(0);

		//	first parse, create vf, ef multimap		
		for (unsigned int f = 0; f != face.size(); f++) {	//	for each face
			for (int i = 0; i < 3; i++) {						//	for each edge of the face
				yz::int2 edge(face[f][i], face[f][(i + 1) % 3]);
				if (edge.x > edge.y)
					yz::mySwap(edge.x, edge.y);

				//	insert vertex-face
				vf_multimap.insert(std::pair<int, int>(face[f][i], f));

				//	insert edge-face
				ef_multimap.insert(std::pair<yz::int2, int>(edge, f));
			}
		}

		//	second parse, extract face-face
		for (unsigned int f = 0; f != face.size(); f++) {	//	for each face
			for (int i = 0; i < 3; i++) {						//	for each edge of the face
				yz::int2 edge(face[f][i], face[f][(i + 1) % 3]);
				if (edge.x > edge.y)
					yz::mySwap(edge.x, edge.y);

				//	for each neighbor of this edge
				auto range = ef_multimap.equal_range(edge);
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

	/**
	set existing from the mesh directly

	non-manifold edges and edge seperated in UV are seams
	*/
	void SetUVFromMesh() {
		if (ef_multimap.empty())	//	must create topology first
			return;

		seam_unordered_set.clear();
		tex_v_mapping.clear();
		v_tex_multimap.clear();

		//	texture already created, setup v-tex mapping
		if (face.size() == tex_face.size()) {
			//	setup v - tex mapping
			tex_v_mapping.resize(tex_coord.size(), -1);
			for (unsigned int fid = 0; fid != face.size(); fid++) {
				for (int i = 0; i < 3; i++) {
					int tid = tex_face[fid][i];
					int vid = face[fid][i];

					//	record tex-v
					if (tex_v_mapping.size() <= tid)
						std::cout << "error: UVunwrapper::SetUVFromMesh, illegal texture index" << std::endl;
					else if (tex_v_mapping[tid] == -1)
						tex_v_mapping[tid] = vid;
					else if (tex_v_mapping[tid] != vid)
						std::cout << "error: UVunwrapper::SetUVFromMesh, tex_v multi projection" << std::endl;

					//	record v-tex
					v_tex_multimap.insert(std::pair<int, int>(vid, tid));
				}
			}

			//	
			int2 curr_edge(-1, -1);
			for (auto ef_iter = ef_multimap.begin(); ef_iter != ef_multimap.end(); ef_iter++) {
				//	this edge is already checked
				if (ef_iter->first == curr_edge)
					continue;

				//	a new edge
				curr_edge = ef_iter->first;
				if (ef_multimap.count(curr_edge) != 2) {	//	non manifold edge, must be seam
					SetSeam(curr_edge);
					continue;
				}

				//	an edge with 2 neighboring faces, check whether it is a single edge on UV
				auto ef_iter_next = ef_iter;
				ef_iter_next++;
				int fid0 = ef_iter->second;
				int fid1 = ef_iter_next->second;
				for (int i = 0; i < 2; i++) {
					int tid0 = -1, tid1 = -1;
					for (int j = 0; j < 3; j++) {
						if (face[fid0][j] == curr_edge[i]) {
							tid0 = tex_face[fid0][j];
							break;
						}
					}
					for (int j = 0; j < 3; j++) {
						if (face[fid1][j] == curr_edge[i]) {
							tid1 = tex_face[fid1][j];
							break;
						}
					}
					if (tid0 != tid1) {
						SetSeam(curr_edge);
						break;
					}
				}
			}

			return;
		}

		//	texture not created, just add non-manifold edges
		int2 curr_edge(-1, -1);
		for (auto ef_iter = ef_multimap.begin(); ef_iter != ef_multimap.end(); ef_iter++) {
			//	this edge is already checked
			if (ef_iter->first == curr_edge)
				continue;

			//	a new edge
			curr_edge = ef_iter->first;
			if (ef_multimap.count(curr_edge) != 2) {	//	non manifold edge, must be seam
				SetSeam(curr_edge);
				continue;
			}
		}
	}

	/**
	label each face, segmented by seam
	*/
	int CreateFaceLabel() {
		face_label.clear();
		face_label.resize(face.size(), -1);

		int label_id = 0;
		for (unsigned int f = 0; f != face.size(); f++) {
			if (face_label[f] != -1)	//	already labeled
				continue;

			//	set the seed
			std::queue<int> face_queue;
			face_queue.push(f);
			face_label[f] = label_id;

			//	floodfill
			while (!face_queue.empty()) {
				int fid = face_queue.front();
				face_queue.pop();

				//	check whether the label can extend through the edge
				for (int i = 0; i < 3; i++) {
					//	for each edge of the face
					yz::int2 edge(face[fid][i], face[fid][(i + 1) % 3]);
					if (edge.x > edge.y)
						yz::mySwap(edge.x, edge.y);

					//	if the edge is seem, cannot extend in this direction
					if (seam_unordered_set.find(edge) != seam_unordered_set.end())
						continue;

					//	for each neighbor face of the edge
					auto range = ef_multimap.equal_range(edge);
					for (auto it = range.first; it != range.second; ++it) {
						int nfid = it->second;
						if (face_label[nfid] == -1) {
							face_queue.push(nfid);
							face_label[nfid] = label_id;
						}
					}	//	for each neighbor face of the edge
				}	// for each edge
			}	//	while face_queue not empty

			label_id++;
		}	//	for each face

			//	seperate patch by label
		patch_face.clear();
		patch_face.resize(label_id);
		for (int fid = 0; fid < face_label.size(); fid++) {
			int label_id = face_label[fid];
			patch_face[label_id].push_back(fid);
		}

		return label_id;
	}

	/**
	Split the mesh along seam
	*/
	void CreateTexFace() {
		//	mapping split vertex to original
		tex_v_mapping.resize(vertex.size());
		for (int i = 0; i < tex_v_mapping.size(); i++)
			tex_v_mapping[i] = i;

		//	mapping original seam vertex to split
		v_tex_multimap.clear();

		tex_face = face;

		//	for each vertex, check whether it should be split
		for (unsigned int vid = 0; vid != vertex.size(); vid++) {
			//	get face groups around this vertex
			std::vector<int> fan_f, fan_f_start;
			int groups = GetFaceGroup(vid, fan_f, fan_f_start);
			if (groups < 2)	//	we don't split it if just one group (normal vertex)
				continue;

			//	if more than 1 group, this vertex should be split from the second group
			v_tex_multimap.insert(std::pair<int, int>(vid, vid));
			for (int i = 1; i < groups; i++) {	//	skip the first group
				for (int f = fan_f_start[i]; f < fan_f_start[i + 1]; f++) {
					int fid = fan_f[f];
					for (int j = 0; j < 3; j++) {
						if (tex_face[fid][j] == vid)
							tex_face[fid][j] = tex_v_mapping.size();
					}
				}
				v_tex_multimap.insert(std::pair<int, int>(vid, tex_v_mapping.size()));
				tex_v_mapping.push_back(vid);
			}
		}	//	end for each vid
	}

	/**
	unwrap UV
	*/
	int CreateTexCoord() {
		if (patch_face.empty())
			return 0;

		//	count tex_number
		int tex_coord_count = -1;
		for (int fid = 0; fid < tex_face.size(); fid++) {
			for (int i = 0; i < 3; i++)
				if (tex_coord_count < tex_face[fid][i])
					tex_coord_count = tex_face[fid][i];
		}
		tex_coord_count++;
		tex_coord.resize(tex_coord_count);

		std::vector<int> tex_patch_mapping;
		std::vector<int> patch_tex_mapping;
		double x_start = 0;

		//	for each patch, calculate local tex_coord, then merge into tex_coord
		for (int label_id = 0; label_id < patch_face.size(); label_id++) {
			//	create patch mesh face
			BaseTriMesh<T> patch_mesh;
			tex_patch_mapping.clear();
			tex_patch_mapping.resize(tex_coord_count, -1);
			for (int j = 0; j < patch_face[label_id].size(); j++) {
				int fid = patch_face[label_id][j];
				patch_mesh.face.push_back(tex_face[fid]);
				for (int i = 0; i < 3; i++) {
					tex_patch_mapping[tex_face[fid][i]] = 0;
				}
			}

			//	set tex patch mapping
			patch_tex_mapping.clear();
			int patch_v_count = 0;
			for (int i = 0; i < tex_patch_mapping.size(); i++) {
				if (tex_patch_mapping[i] == 0) {
					tex_patch_mapping[i] = patch_v_count;
					patch_tex_mapping.push_back(i);
					patch_v_count++;
				}
			}

			//	create patch mesh vertex
			patch_mesh.vertex.resize(patch_v_count);
			for (int i = 0; i < tex_patch_mapping.size(); i++) {
				if (tex_patch_mapping[i] == -1)
					continue;
				int pvid = tex_patch_mapping[i];
				int vid = tex_v_mapping[i];
				patch_mesh.vertex[pvid] = vertex[vid];
			}

			//	update patch face index
			for (int fid = 0; fid < patch_mesh.face.size(); fid++) {
				for (int i = 0; i < 3; i++) {
					patch_mesh.face[fid][i] = tex_patch_mapping[patch_mesh.face[fid][i]];
				}
			}

			//==================

			//	flating
			MeshPatchUVunwrapper<T> flatter;
			flatter.vertex.swap(patch_mesh.vertex);
			flatter.face.swap(patch_mesh.face);
			int wrap_succ = flatter.Unwrap();

			if (!wrap_succ) {
				tex_coord.clear();
				tex_v_mapping.clear();
				v_tex_multimap.clear();
				return 0;
			}

			//	calculate texture offset
			Vec2<T> offset;
			AABB2D<T> aabb;
			aabb.GetAABBCoef(flatter.tex_coord);
			offset.x = x_start - aabb.bb_min.x;
			offset.y = -aabb.bb_min.y;

			x_start += aabb.bb_max.x - aabb.bb_min.x;

			//	record texture
			for (int i = 0; i < patch_tex_mapping.size(); i++) {
				int tvid = patch_tex_mapping[i];
				tex_coord[tvid] = flatter.tex_coord[i] + offset;
			}
		}

		return 1;
	}

	/**
	get face neighboring a vertex, segmented by seams into groups

	\param	vid				the vertex
	\param	fan_f			all faces neighboring vid, arranged in group order
	\param	fan_f_start		start position in fan_f of each group
	\return					number of groups
	*/
	int GetFaceGroup(int vid, std::vector<int>& fan_f, std::vector<int>& fan_f_start) {
		fan_f.clear();
		fan_f_start.clear();
		fan_f_start.push_back(0);

		//	get all faces neighbor this vertex
		std::map<int, int> f_idx;	//	face - index
		auto range = vf_multimap.equal_range(vid);
		for (auto it = range.first; it != range.second; ++it) {
			f_idx.insert(std::pair<int, int>(it->second, -1));
		}
		if (f_idx.empty())
			return 0;

		//	perform floodfill on fan faces
		int fan_idx = 0;
		for (auto iter = f_idx.begin(); iter != f_idx.end(); iter++) {
			if (iter->second != -1)
				continue;

			//	a seed face
			std::queue<std::map<int, int>::iterator> fan_queue;
			iter->second = fan_idx;
			fan_queue.push(iter);

			//	floodfill on the fan
			while (!fan_queue.empty()) {
				auto f_iter = fan_queue.front();
				fan_queue.pop();

				int fid = f_iter->first;

				//	check whether the label can extend through the edge
				for (int i = 0; i < 3; i++) {
					//	get the edge from vid
					if (face[fid][i] == vid)
						continue;
					yz::int2 edge(vid, face[fid][i]);
					if (edge.x > edge.y)
						yz::mySwap(edge.x, edge.y);

					//	if the edge is seem, cannot extend in this direction
					if (seam_unordered_set.find(edge) != seam_unordered_set.end())
						continue;

					//	for each neighbor face of the edge
					auto range = ef_multimap.equal_range(edge);
					for (auto it = range.first; it != range.second; ++it) {
						int nfid = it->second;
						auto nei_iter = f_idx.find(nfid);
						if (nei_iter == f_idx.end()) {
							std::cout << "this should not happen" << std::endl;
							continue;
						}
						if (nei_iter->second == -1) {	//	this fan face is not labeled yet
							nei_iter->second = fan_idx;
							fan_queue.push(nei_iter);
						}
					}	//	for each neighbor face of the edge
				}	// for each edge
			}	//	queue empty
			fan_idx++;
		}

		//	group the fan
		std::vector<std::pair<int, int>> idx_f;
		for (auto iter = f_idx.begin(); iter != f_idx.end(); iter++) {
			idx_f.push_back(std::pair<int, int>(iter->second, iter->first));
		}
		std::sort(idx_f.begin(), idx_f.end());

		fan_f.reserve(idx_f.size());
		fan_f_start.reserve(fan_idx + 1);
		fan_f.push_back(idx_f[0].second);
		for (int i = 1; i < idx_f.size(); i++) {
			if (idx_f[i].first != idx_f[i - 1].first)	//	a new label
				fan_f_start.push_back(fan_f.size());
			fan_f.push_back(idx_f[i].second);
		}
		fan_f_start.push_back(fan_f.size());

		return fan_idx;
	}

protected:
	//	seam
	std::vector<int2>			seam;

	//	labeling, seperated by seam
	std::vector<int>				face_label;
	std::vector<std::vector<int>>	patch_face;

	//	mesh topology
	std::multimap<int, int>		vf_multimap;
	std::multimap<int2, int>	ef_multimap;	///< edge is represented as (x, y), x < y
	std::vector<int>			ff;
	std::vector<int>			ff_start;

	//	lookup table
	std::unordered_set<int2, utils::BitwiseHasher<int2>>	seam_unordered_set;
	std::vector<int>			tex_v_mapping;
	std::multimap<int, int>		v_tex_multimap;
};


}}	//	namespace yz::geometry

#endif	//	__YZ_MESH_UV_H__