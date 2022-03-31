#include <iostream>
#include <GL/glut.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <yzLib/yz_lib.h>


namespace yz {  namespace geometry {

/**
Deform mesh using key points
*/
template<class T>
class MeshDeformer {
public:
	int SetMesh(
		int vertex_number,
		T*	vertex_xyz_ptr,
		int face_number,
		int* face_xyz_ptr)
	{
		if (vertex_number < 1 || face_number < 1) {
			std::cout << "error: MeshDeformer::SetMesh, empty mesh" << std::endl;
			return 0;
		}

		//	setup mesh
		mesh.vertex.resize(vertex_number);
		for (int i = 0; i < vertex_number; i++) {
			mesh.vertex[i].x = vertex_xyz_ptr[i * 3];
			mesh.vertex[i].y = vertex_xyz_ptr[i * 3 + 1];
			mesh.vertex[i].z = vertex_xyz_ptr[i * 3 + 2];
		}

		mesh.face.resize(face_number);
		for (int i = 0; i < face_number; i++) {
			mesh.face[i].x = face_xyz_ptr[i * 3];
			mesh.face[i].y = face_xyz_ptr[i * 3 + 1];
			mesh.face[i].z = face_xyz_ptr[i * 3 + 2];
		}

		//	setup constraint
		constraint_point.clear();
		constraint_point.resize(mesh.vertex.size());
		constraint_weight.clear();
		constraint_weight.resize(mesh.vertex.size(), 0);

		//	create topoloty
		mesh.CalculateNormals();
		vf.CreateVF(mesh.vertex.size(), mesh.face);
		boundary.MarkBoundaryVertex(mesh.vertex.size(), mesh.face);

		//	create Laplacian matrix
		yz::geometry::createLaplacianMatrixForMesh(Laplacian, mesh.vertex, mesh.face);

		//	set matrix
		int dim = mesh.vertex.size();
		L.resize(dim, dim);
		L.reserve(Laplacian.NNZ());
		for (int i = 0; i < Laplacian.NNZ(); i++) {
			int row_id = Laplacian.row_id[i];
			int col_id = Laplacian.col_id[i];
			T value = Laplacian.value[i];
			L.insert(row_id, col_id) = value;
		}
		L.makeCompressed();
		LT = L.transpose();
		LTL = LT * L;

		//	create Laplacian normal
		Laplacian_normal.resize(mesh.vertex.size());
		yz::DenseVector<T>	X(mesh.vertex.size());
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < X.Dim(); j++) {
				X[j] = mesh.vertex[j][i];
			}
			X = Laplacian * X;
			for (int j = 0; j < X.Dim(); j++) {
				Laplacian_normal[j][i] = X[j];
			}
		}

		Laplacian_normal_length.resize(Laplacian_normal.size());
		for (int i = 0; i < Laplacian_normal.size(); i++) {
			Laplacian_normal_length[i] = Laplacian_normal[i].Length();
		}

		CalculateLaplacianNormalWeight();

		return 1;
	}

	int SetConstraint(
		int constraint_number,
		int* constraint_vertex_index_ptr,
		T* constraint_point_xyz_ptr,
		T* constraint_weight_ptr)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformer::SetConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < constraint_number; i++) {
			int vertex_idx = constraint_vertex_index_ptr[i];
			constraint_point[vertex_idx].x = constraint_point_xyz_ptr[i * 3 + 0];
			constraint_point[vertex_idx].y = constraint_point_xyz_ptr[i * 3 + 1];
			constraint_point[vertex_idx].z = constraint_point_xyz_ptr[i * 3 + 2];
			constraint_weight[vertex_idx] = constraint_weight_ptr[i];
		}

		//	set A and b
		SetA();
		Setb();

		return 1;
	}

	int SetConstraint(
		T* constraint_point_xyz_ptr,
		T* constraint_weight_ptr)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformer::SetConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < mesh.vertex.size(); i++) {
			constraint_point[i].x = constraint_point_xyz_ptr[i * 3 + 0];
			constraint_point[i].y = constraint_point_xyz_ptr[i * 3 + 1];
			constraint_point[i].z = constraint_point_xyz_ptr[i * 3 + 2];
			constraint_weight[i] = constraint_weight_ptr[i];
		}

		//	set A and b
		SetA();
		Setb();

		return 1;
	}

	int InitDeformer() {
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformer::InitDeformer, empty mesh" << std::endl;
			return 0;
		}

		//	set matrix A
		SetA();

		//	set right hand side b
		Setb();

		return 1;
	}

	int DeformStep() {
		if (!A.rows() || !b.rows() || A.rows() != b.rows()) {
			std::cout << "error: MeshDeformer::DeformStep, A & b are not set" << std::endl;
			return 0;
		}

		//	solve matrix
		X = solver.solve(b);

		//	refresh mesh
		for (int i = 0; i < X.rows(); i++) {
			mesh.vertex[i].x = X(i, 0);
			mesh.vertex[i].y = X(i, 1);
			mesh.vertex[i].z = X(i, 2);
		}

		mesh.CalculateNormals();

		UpdateLaplacianNormal();

		return 1;
	}

	int AddConstraint(int vertex_idx, T new_point_xyz[3], T weight) {
		//	weight will change in this function, so A need to be recalculated
		if (vertex_idx < 0 || vertex_idx >= mesh.vertex.size()) {
			std::cout << "error: MeshDeformer::AddConstraint, invalid vertex index" << std::endl;
			return 0;
		}

		constraint_point[vertex_idx].x = new_point_xyz[0];
		constraint_point[vertex_idx].y = new_point_xyz[1];
		constraint_point[vertex_idx].z = new_point_xyz[2];
		constraint_weight[vertex_idx] = weight;

		//	set A and b
		SetA();
		Setb();

		return 1;
	}

	int UpdateConstraint(int vertex_idx, T new_point_xyz[3]) {
		//	weight doesn't change in this function, so A doesn't change
		if (vertex_idx < 0 || vertex_idx >= mesh.vertex.size()) {
			std::cout << "error: MeshDeformer::UpdateConstraint, invalid vertex index" << std::endl;
			return 0;
		}

		//	it is assumed that this vertex is constraint originally, so it doesn't change the format of the matrix
		constraint_point[vertex_idx].x = new_point_xyz[0];
		constraint_point[vertex_idx].y = new_point_xyz[1];
		constraint_point[vertex_idx].z = new_point_xyz[2];

		//	update b
		Setb();

		return 1;
	}

protected:
	virtual int SetA() {
		if (mesh.vertex.empty())
			return 0;

		//	set matrix
		A = LTL;
		for (int i = 0; i < constraint_weight.size(); i++) {
			A.coeffRef(i, i) += constraint_weight[i];
		}

		//	prepare solver
		solver.compute(A);

		return 1;
	}

	virtual int Setb() {
		if (mesh.vertex.empty())
			return 0;

		int dim = mesh.vertex.size();

		//	set right hand side b
		b.resize(dim, 3);
		for (int i = 0; i < dim; i++) {
			b(i, 0) = Laplacian_normal[i].x;
			b(i, 1) = Laplacian_normal[i].y;
			b(i, 2) = Laplacian_normal[i].z;
		}

		b = LT * b;

		for (int i = 0; i < dim; i++) {
			b(i, 0) += constraint_weight[i] * constraint_point[i].x;
			b(i, 1) += constraint_weight[i] * constraint_point[i].y;
			b(i, 2) += constraint_weight[i] * constraint_point[i].z;
		}

		return 1;
	}

	int CalculateLaplacianNormalWeight() {
		if (mesh.vertex.empty())
			return 0;

		normal_weight.resize(vf.vf.size(), 1);

		for (int v = 0; v < mesh.vertex.size(); v++) {
			//	M * Y = N, solve with (MT * M + I) * Y = MT * N
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	M;
			Eigen::Matrix<T, Eigen::Dynamic, 1>					Y;
			Eigen::Matrix<T, 3, 1>								N;

			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	MTMPLUSI;
			Eigen::Matrix<T, Eigen::Dynamic, 1>					MTN;

			int nn = vf.vf_start[v + 1] - vf.vf_start[v];
			M.resize(3, nn);
			Y.resize(nn);
			MTMPLUSI.resize(nn, nn);
			MTN.resize(nn);

			//	set M
			for (int i = 0, j = vf.vf_start[v]; i < nn; i++, j++) {
				int f = vf.vf[j];
				yz::Vec3<T>	fn = mesh.face_normal[f];
				M(0, i) = fn.x;
				M(1, i) = fn.y;
				M(2, i) = fn.z;
			}

			//	set N
			yz::Vec3<T>	ln = Laplacian_normal[v].Normalize();
			N(0) = ln.x;
			N(1) = ln.y;
			N(2) = ln.z;

			//	set MT * M + I
			MTMPLUSI = M.transpose() * M;
			for (int i = 0; i < nn; i++)
				MTMPLUSI(i, i) += 0.1;

			//	set MT * N
			MTN = M.transpose() * N;

			//	calculate Y
			Y = MTMPLUSI.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(MTN);

			//	validate Y
			T abs_sum = 0;
			for (int i = 0; i < nn; i++) {
				abs_sum += fabs(Y(i));
			}
			if (abs_sum != abs_sum || abs_sum < 1e-3) {
				Y.setOnes();
			}

			//	write to weight array
			for (int i = 0, j = vf.vf_start[v]; i < nn; i++, j++) {
				normal_weight[j] = Y(i);

				//std::cout << Y(i) << ", ";
			}
			//std::cout << std::endl;
		}

		return 1;
	}

	int UpdateLaplacianNormal() {
		for (int i = 0; i < mesh.vertex.size(); i++) {
			yz::Vec3<T>	n(0, 0, 0);
			for (int j = vf.vf_start[i], k = 0; j < vf.vf_start[i + 1]; j++, k++) {
				int f = vf.vf[j];
				T w = normal_weight[j];
				n += mesh.face_normal[f] * w;
			}
			n.SetNormalize();
			Laplacian_normal[i] = n * Laplacian_normal_length[i];
		}

		Setb();

		return 1;
	}

public:
	yz::geometry::TriMesh<T>							mesh;
	std::vector<yz::Vec3<T>>							constraint_point;
	std::vector<T>										constraint_weight;

protected:
	yz::geometry::TriMeshVF								vf;
	yz::geometry::TriMeshBoundaryVertex					boundary;

	yz::SparseMatrixCoo<T>								Laplacian;					//	Laplacian matrix
	std::vector<yz::Vec3<T>>							Laplacian_normal;			//	= Laplacian * mesh_vertex
	std::vector<T>										Laplacian_normal_length;	//	= Laplacian * mesh_vertex
	std::vector<T>										normal_weight;				//	weight of laplacian_normal direction with one ring neighbor face normal

	Eigen::SparseMatrix<T>								L;		//	Laplacian matrix
	Eigen::SparseMatrix<T>								LT;		//	transpose(Laplacian matrix)
	Eigen::SparseMatrix<T>								LTL;	//	transpose(Laplacian matrix) * Laplacian matrix

																//	elements for A X = b
	Eigen::SparseMatrix<T>								A;		//	A = LTL + w * I
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	X;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	b;		//	b = LT * Laplacian_normal + w * constraint_point
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<T> >		solver;
};


/**
deforme mesh with point plane
*/
template<class T>
class MeshDeformerPointPlane : public MeshDeformer<T> {
public:
	int SetMesh(
		int vertex_number,
		T*	vertex_xyz_ptr,
		int face_number,
		int* face_xyz_ptr)
	{
		MeshDeformer<T>::SetMesh(vertex_number, vertex_xyz_ptr, face_number, face_xyz_ptr);

		//	setup plane constraint
		plane_constraint_point.clear();
		plane_constraint_point.resize(mesh.vertex.size());
		plane_constraint_normal.clear();
		plane_constraint_normal.resize(mesh.vertex.size());
		plane_constraint_weight.clear();
		plane_constraint_weight.resize(mesh.vertex.size(), 0);

		CalculateNNT();

		//	set matrix (L is originally set in SetMesh, now change it)
		int dim = mesh.vertex.size() * 3;
		L.resize(dim, dim);
		L.reserve(Laplacian.NNZ() * 3);
		for (int i = 0; i < Laplacian.NNZ(); i++) {
			int row_id = Laplacian.row_id[i] * 3;
			int col_id = Laplacian.col_id[i] * 3;
			T value = Laplacian.value[i];
			L.insert(row_id, col_id) = value;
			L.insert(row_id + 1, col_id + 1) = value;
			L.insert(row_id + 2, col_id + 2) = value;
		}
		L.makeCompressed();
		LT = L.transpose();
		LTL = LT * L;

		return 1;
	}

	int SetPlaneConstraint(
		int plane_constraint_number,
		int* plane_constraint_vertex_index_ptr,
		T* plane_constraint_point_xyz_ptr,
		T* plane_constraint_normal_xyz_ptr,
		T* plane_constraint_weight_ptr,
		int set_matrix_flag = 1)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformerPointPlane::SetPlaneConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < plane_constraint_number; i++) {
			int vertex_idx = plane_constraint_vertex_index_ptr[i];
			plane_constraint_point[vertex_idx].x = plane_constraint_point_xyz_ptr[i * 3 + 0];
			plane_constraint_point[vertex_idx].y = plane_constraint_point_xyz_ptr[i * 3 + 1];
			plane_constraint_point[vertex_idx].z = plane_constraint_point_xyz_ptr[i * 3 + 2];
			plane_constraint_normal[vertex_idx].x = plane_constraint_normal_xyz_ptr[i * 3 + 0];
			plane_constraint_normal[vertex_idx].y = plane_constraint_normal_xyz_ptr[i * 3 + 1];
			plane_constraint_normal[vertex_idx].z = plane_constraint_normal_xyz_ptr[i * 3 + 2];
			plane_constraint_weight[vertex_idx] = plane_constraint_weight_ptr[i];
		}

		if (set_matrix_flag) {
			//	set A and b
			CalculateNNT();
			SetA();
			Setb();
		}

		return 1;
	}

	int SetPlaneConstraint(
		T* plane_constraint_point_xyz_ptr,
		T* plane_constraint_normal_xyz_ptr,
		T* plane_constraint_weight_ptr,
		int set_matrix_flag = 1)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformerPointPlane::SetPlaneConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < mesh.vertex.size(); i++) {
			plane_constraint_point[i].x = plane_constraint_point_xyz_ptr[i * 3 + 0];
			plane_constraint_point[i].y = plane_constraint_point_xyz_ptr[i * 3 + 1];
			plane_constraint_point[i].z = plane_constraint_point_xyz_ptr[i * 3 + 2];
			plane_constraint_normal[i].x = plane_constraint_normal_xyz_ptr[i * 3 + 0];
			plane_constraint_normal[i].y = plane_constraint_normal_xyz_ptr[i * 3 + 1];
			plane_constraint_normal[i].z = plane_constraint_normal_xyz_ptr[i * 3 + 2];
			plane_constraint_weight[i] = plane_constraint_weight_ptr[i];
		}

		if (set_matrix_flag) {
			//	set A and b
			CalculateNNT();
			SetA();
			Setb();
		}

		return 1;
	}

	int SetConstraint(
		int constraint_number,
		int* constraint_vertex_index_ptr,
		T* constraint_point_xyz_ptr,
		T* constraint_weight_ptr,
		int set_matrix_flag = 1)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformer::SetConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < constraint_number; i++) {
			int vertex_idx = constraint_vertex_index_ptr[i];
			constraint_point[vertex_idx].x = constraint_point_xyz_ptr[i * 3 + 0];
			constraint_point[vertex_idx].y = constraint_point_xyz_ptr[i * 3 + 1];
			constraint_point[vertex_idx].z = constraint_point_xyz_ptr[i * 3 + 2];
			constraint_weight[vertex_idx] = constraint_weight_ptr[i];
		}

		if (set_matrix_flag) {
			//	set A and b
			SetA();
			Setb();
		}

		return 1;
	}

	int SetConstraint(
		T* constraint_point_xyz_ptr,
		T* constraint_weight_ptr,
		int set_matrix_flag = 1)
	{
		if (mesh.vertex.empty() || mesh.face.empty()) {
			std::cout << "error: MeshDeformer::SetConstraint, empty mesh" << std::endl;
			return 0;
		}

		for (int i = 0; i < mesh.vertex.size(); i++) {
			constraint_point[i].x = constraint_point_xyz_ptr[i * 3 + 0];
			constraint_point[i].y = constraint_point_xyz_ptr[i * 3 + 1];
			constraint_point[i].z = constraint_point_xyz_ptr[i * 3 + 2];
			constraint_weight[i] = constraint_weight_ptr[i];
		}

		if (set_matrix_flag) {
			//	set A and b
			SetA();
			Setb();
		}

		return 1;
	}

	int DeformStep() {
		if (!A.rows() || !b.rows() || A.rows() != b.rows()) {
			std::cout << "error: MeshDeformer::DeformStep, A & b are not set" << std::endl;
			return 0;
		}

		//	solve matrix
		X = solver.solve(b);

		//	refresh mesh
		for (int i = 0; i < mesh.vertex.size(); i++) {
			mesh.vertex[i].x = X(i * 3);
			mesh.vertex[i].y = X(i * 3 + 1);
			mesh.vertex[i].z = X(i * 3 + 2);
		}

		mesh.CalculateNormals();

		UpdateLaplacianNormal();

		return 1;
	}

	int AddPlaneConstraint(int vertex_idx, T new_plane_point_xyz[3], T new_plane_normal_xyz[3], T weight) {
		if (vertex_idx < 0 || vertex_idx >= mesh.vertex.size()) {
			std::cout << "error: MeshDeformerPointPlane::AddPlaneConstraint, invalid vertex index" << std::endl;
			return 0;
		}

		plane_constraint_point[vertex_idx].x = new_plane_point_xyz[0];
		plane_constraint_point[vertex_idx].y = new_plane_point_xyz[1];
		plane_constraint_point[vertex_idx].z = new_plane_point_xyz[2];
		plane_constraint_normal[vertex_idx].x = new_plane_normal_xyz[0];
		plane_constraint_normal[vertex_idx].y = new_plane_normal_xyz[1];
		plane_constraint_normal[vertex_idx].z = new_plane_normal_xyz[2];
		plane_constraint_weight[vertex_idx] = weight;

		//	set A and b
		CalculateNNT();
		SetA();
		Setb();

		return 1;
	}

	int SetMatrix() {
		//	set A and b
		CalculateNNT();
		SetA();
		Setb();

		return 1;
	}

protected:
	virtual int SetA() {
		if (mesh.vertex.empty())
			return 0;

		//	set matrix
		A = LTL + NNT;
		for (int i = 0; i < constraint_weight.size(); i++) {
			A.coeffRef(i * 3, i * 3) += constraint_weight[i];
			A.coeffRef(i * 3 + 1, i * 3 + 1) += constraint_weight[i];
			A.coeffRef(i * 3 + 2, i * 3 + 2) += constraint_weight[i];
		}

		//	prepare solver
		solver.compute(A);

		return 1;
	}

	virtual int Setb() {
		if (mesh.vertex.empty())
			return 0;

		//	set right hand side b
		b.resize(mesh.vertex.size() * 3, 1);
		for (int i = 0; i < mesh.vertex.size(); i++) {
			b(i * 3) = Laplacian_normal[i].x;
			b(i * 3 + 1) = Laplacian_normal[i].y;
			b(i * 3 + 2) = Laplacian_normal[i].z;
		}

		b = LT * b;

		for (int i = 0; i < mesh.vertex.size(); i++) {
			b(i * 3) += constraint_weight[i] * constraint_point[i].x;
			b(i * 3 + 1) += constraint_weight[i] * constraint_point[i].y;
			b(i * 3 + 2) += constraint_weight[i] * constraint_point[i].z;
		}


		for (int i = 0; i < mesh.vertex.size(); i++) {
			yz::Vec3<T> p = plane_constraint_point[i];
			yz::Vec3<T>	n = plane_constraint_normal[i];
			T w = plane_constraint_weight[i];
			b(i * 3) += w*(n[0] * n[0] * p[0] + n[0] * n[1] * p[1] + n[0] * n[2] * p[2]);
			b(i * 3 + 1) += w*(n[1] * n[0] * p[0] + n[1] * n[1] * p[1] + n[1] * n[2] * p[2]);
			b(i * 3 + 2) += w*(n[2] * n[0] * p[0] + n[2] * n[1] * p[1] + n[2] * n[2] * p[2]);
		}

		return 1;
	}

	int CalculateNNT() {
		if (mesh.vertex.empty())
			return 0;

		int dim = mesh.vertex.size() * 3;
		NNT.resize(dim, dim);
		NNT.reserve(dim * 3);
		for (int idx = 0; idx < plane_constraint_normal.size(); idx++) {
			yz::Vec3<T>	nor = plane_constraint_normal[idx];
			T w = plane_constraint_weight[idx];
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < 3; i++) {
					NNT.insert(idx * 3 + i, idx * 3 + j) = w * nor[j] * nor[i];
				}
			}
		}

		return 1;
	}

public:
	std::vector<yz::Vec3<T>>	plane_constraint_point;
	std::vector<yz::Vec3<T>>	plane_constraint_normal;
	std::vector<T>				plane_constraint_weight;

	Eigen::SparseMatrix<T>		NNT;
};


}}	//	namespace yz::geometry

int main() {

	return 0;
}