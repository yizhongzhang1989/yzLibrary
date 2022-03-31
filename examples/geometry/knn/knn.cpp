#include <iostream>
#include <GL/glut.h>
#include "yzLib/yz_lib.h"

template <typename T>
int createClusterOnPointCloud(
	std::vector<int>&					cluster_id,
	const std::vector<yz::Vec3<T>>&		vertex,
	const yz::geometry::AABBTree3D<T>&	vertex_aabb,
	T									dist_thre,
	int									min_point_thre = 50,
	int									search_radius = 10
) {
	//	preprocess cluster_id
	if (cluster_id.size() != vertex.size()) {
		cluster_id.clear();
		cluster_id.resize(vertex.size(), 0);
	}
	else {
		//	negative value means this point should be ignored
		for (int i = 0; i < cluster_id.size(); i++)
			cluster_id[i] = cluster_id[i] >= 0 ? 0 : -1;
	}

	std::vector<yz::Vec3<T>>	knp;
	std::vector<int>			knpid;
	std::vector<int>			patch_vid;
	std::queue<int>				patch_queue;
	int patch_id = 1;
	T squ_dist_thre = dist_thre * dist_thre;
	for (int vid = 0; vid < vertex.size(); vid++) {
		if (cluster_id[vid] != 0)
			continue;

		patch_vid.clear();
		patch_vid.push_back(vid);
		patch_queue.push(vid);
		cluster_id[vid] = patch_id;
		while (!patch_queue.empty()) {
			int curr_vid = patch_queue.front();
			T max_dist = yz::geometry::getKNearestPointsOnVertices(
				knp, knpid, 10, vertex[curr_vid], vertex, vertex_aabb);
			patch_queue.pop();

			for (int i = 0; i < knpid.size(); i++) {
				if (cluster_id[knpid[i]])
					continue;
				T squ_dist = (vertex[curr_vid] - vertex[knpid[i]]).SquareLength();
				if (squ_dist > squ_dist_thre)
					continue;

				patch_vid.push_back(knpid[i]);
				patch_queue.push(knpid[i]);
				cluster_id[knpid[i]] = patch_id;
			}
		}

		//	if patch is too small, ignore it
		if (patch_vid.size() < min_point_thre) {
			for (int i = 0; i < patch_vid.size(); i++)
				cluster_id[patch_vid[i]] = 0;
		}
		else {
			patch_id++;
		}
	}

	//	regularize
	for (int i = 0; i < cluster_id.size(); i++) {
		if (cluster_id[i] >= 0)
			cluster_id[i]--;
	}

	return patch_id - 1;
}



template <typename T>
int fitPlaneFromPointsRegionGrowing(
	yz::Vec3<T>&						plane_org,
	yz::Vec3<T>&						plane_x,
	yz::Vec3<T>&						plane_y,
	yz::Vec3<T>&						plane_z,
	yz::Vec3<T>&						eigen_xyz,
	std::vector<int>&					inlier_mask,
	const std::vector<yz::Vec3<T>>&		point,
	const yz::geometry::AABBTree3D<T>&	point_aabb,
	T									dist_thre,
	int									seed_point_idx
) {
	if (seed_point_idx < 0 || seed_point_idx >= point.size())
		return 0;
	if (point.size() < 20)
		return 0;
	
	inlier_mask.clear();
	inlier_mask.resize(point.size(), 0);

	//	create a dist from the initial 
	std::vector<yz::Vec3<T>>	knp;
	std::vector<int>			knpid;
	T max_dist = yz::geometry::getKNearestPointsOnVertices(
		knp, knpid, 20, point[seed_point_idx], point, point_aabb);
	if (max_dist < dist_thre * 3) {

	}

}


template <class T>
class ClusterSeparater {
public:
	void Separate() {

	}

	void LabelPlanes(T dist_thre) {
		yz::Vec3<T>			plane_org, plane_x, plane_y, plane_z, eigen_xyz;
		std::vector<int>	inlier_mask;
		int fit_point_num = yz::geometry::fitPlaneFromPointsRANSAC(
			plane_org, plane_x, plane_y, plane_z, eigen_xyz, inlier_mask, vertex, dist_thre);
	}

public:
	std::vector<yz::Vec3<T>>		vertex;

	yz::geometry::AABBTree3D<T>		vertex_aabb;

	std::vector<int>				label;
};


template <class T>
class PointCluster {
public:
	PointCluster(const std::vector<yz::Vec3<T>>& v) : vertex(v) {
		vertex_aabb.BuildVertexAABBTree(vertex);

		cluster_id.resize(vertex.size(), 0);
		cluster_rgb.resize(vertex.size(), yz::uchar3(255, 255, 0));
	}

	void ClusterByDistance(T dist_thre) {
		createClusterOnPointCloud(cluster_id, vertex, vertex_aabb, dist_thre);

		CreateClusterRGB();
	}

	void MergeClusterByAABB() {
		int cluster_num = 0;
		for (int i = 0; i < cluster_id.size(); i++) {
			if (cluster_id[i] >= cluster_num)
				cluster_num = cluster_id[i] + 1;
		}

		yz::geometry::AABB2D<double> aabb_init;
		aabb_init.bb_min = yz::Vec2d(1e6, 1e6);
		aabb_init.bb_max = yz::Vec2d(-1e6, -1e6);
		std::vector<yz::geometry::AABB2D<double>> aabb;
		aabb.resize(cluster_num, aabb_init);

		//	calculate aabb of each cluster
		for (int i = 0; i < cluster_id.size(); i++) {
			int id = cluster_id[i];
			if (id >= 0)
				aabb[id].Expand(yz::Vec2d(vertex[i].x, vertex[i].y));
		}

		//	merge aabb
		std::vector<int> reid;
		for (int i = 0; i < aabb.size(); i++)
			reid.push_back(i);

		int merge_flag = 0;
		do {
			merge_flag = 0;

			for (int i = 0; i < aabb.size(); i++) {
				if (reid[i] != i)
					continue;
				for (int idx = i + 1; idx < aabb.size(); idx++) {
					int j = reid[idx];
					while (reid[j] != j)
						j = reid[j];
					if (i >= j)
						continue;

					if (yz::geometry::isAABBsIntersect(aabb[i].bb_min, aabb[i].bb_max, aabb[j].bb_min, aabb[j].bb_max)) {
						aabb[i].Expand(aabb[j]);
						reid[idx] = i;
						merge_flag = 1;
					}
				}
			}
		} while (merge_flag);

		//	reid
		for (int i = 0; i < reid.size(); i++) {
			int j = reid[i];
			while (reid[j] != j)
				j = reid[j];
			reid[i] = j;
		}

		int count = 0;
		for (int i = 0; i < reid.size(); i++) {
			if (reid[i] == i) {
				reid[i] = count;
				count++;
			}
			else {
				reid[i] = reid[reid[i]];
			}
		}

		for (int i = 0; i < cluster_id.size(); i++) {
			cluster_id[i] = reid[cluster_id[i]];
		}

		CreateClusterRGB();
	}

	void CreateAABB() {
		int cluster_num = 0;
		for (int i = 0; i < cluster_id.size(); i++) {
			if (cluster_id[i] >= cluster_num)
				cluster_num = cluster_id[i] + 1;
		}

		yz::geometry::AABB3D<double> aabb_init;
		aabb_init.bb_min = yz::Vec3d(1e6, 1e6, 1e6);
		aabb_init.bb_max = yz::Vec3d(-1e6, -1e6, -1e6);
		aabb.resize(cluster_num, aabb_init);

		//	calculate aabb of each cluster
		for (int i = 0; i < cluster_id.size(); i++) {
			int id = cluster_id[i];
			if (id >= 0)
				aabb[id].Expand(vertex[i]);
		}
	}

	void Draw() {
		if (vertex.empty())
			return;

		static std::vector<int> disp_indice;
		for (size_t i = disp_indice.size(); i < vertex.size(); i++)
			disp_indice.push_back(i);

		glPointSize(3);
		yz::opengl::drawColorPointCloud(
			(T*)&vertex[0],
			(unsigned char*)&cluster_rgb[0],
			vertex.size(),
			&disp_indice[0]);

		for (int i = 0; i < aabb.size(); i++) {
			yz::opengl::setSequentialDisplayColor(i);

			yz::opengl::drawAABB(aabb[i].bb_min, aabb[i].bb_max);

			//glDisable(GL_LIGHTING);
			//yz::opengl::drawAABBWire(aabb[i].bb_min, aabb[i].bb_max);
		}		
	}

	void CreateClusterRGB() {
		for (int i = 0; i < cluster_id.size(); i++) {
			if (cluster_id[i] < 0) {
				cluster_rgb[i] = yz::uchar3(0, 0, 0);
			}
			else {
				yz::uchar3 color;
				yz::opengl::getSequentialDisplayColor(&color.x, cluster_id[i]);
				cluster_rgb[i] = color;
			}
		}
	}

public:
	const std::vector<yz::Vec3<T>>& vertex;
	yz::geometry::AABBTree3D<T>		vertex_aabb;

	std::vector<int>				cluster_id;
	std::vector<yz::uchar3>			cluster_rgb;

	std::vector<yz::geometry::AABB3D<double>> aabb;
};

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow3D		win3d;

yz::geometry::TriMesh<double>	mesh;
PointCluster<double>*			cluster_ptr = NULL;

void draw() {
	yz::opengl::drawXYZAxis();

	cluster_ptr->Draw();

	yz::geometry::AABB3D<double> ground_aabb = cluster_ptr->aabb[0];
	for (int i = 1; i < cluster_ptr->aabb.size(); i++)
		ground_aabb.Expand(cluster_ptr->aabb[i]);

	//glColor3f(0.5, 0.5, 0.5);
	//glBegin(GL_QUADS);
	//glNormal3f(0, 0, -1);
	//glVertex2f(ground_aabb.bb_min.x, ground_aabb.bb_min.y);
	//glVertex2f(ground_aabb.bb_max.x, ground_aabb.bb_min.y);
	//glVertex2f(ground_aabb.bb_max.x, ground_aabb.bb_max.y);
	//glVertex2f(ground_aabb.bb_min.x, ground_aabb.bb_max.y);
	//glEnd();
}

int main() {
	mesh.ReadMeshFromFile("non_ground.obj");

	cluster_ptr = new PointCluster<double>(mesh.vertex);
	cluster_ptr->ClusterByDistance(0.2);
	//cluster_ptr->MergeClusterByAABB();
	cluster_ptr->CreateAABB();

	win3d.back_ground_red = 0.5;
	win3d.back_ground_green = 0.5;
	win3d.back_ground_blue = 0.5;
	win3d.use_arcball_flag = 1;
	win3d.SetDraw(draw);
	win3d.CreateGLUTWindow();

	manager.EnterMainLoop();
}
