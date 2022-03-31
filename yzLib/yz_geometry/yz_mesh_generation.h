/***********************************************************/
/**	\file
	\brief		Mesh Generation
	\details	In this file, several functions are provided 
				to generate basic shape meshes
	\author		Yizhong Zhang
	\date		12/8/2016
*/
/***********************************************************/
#ifndef __YZ_MESH_GENERATION_H__
#define __YZ_MESH_GENERATION_H__

namespace yz {	namespace geometry {

/**
	Create a unit circle on XY plane, store the data from a given iterator

	\param	vertex		return the vertex list
	\param	iter_start	the iterator to store the circle data, must be in vertex
	\param	slices		how many slices of the ring
	\return				whether create succeed
*/
template<typename T>
int createUnitCircleXY(
	std::vector<Vec3<T>>&					vertex,
	typename std::vector<Vec3<T>>::iterator	iter_start,
	int										slices = 16
) {
	if (slices < 3) {
		return 0;
	}

	//	make sure there is enough space to store the generated data
	if (vertex.end() - iter_start < slices) {
		vertex.resize(slices - (vertex.end() - iter_start) + vertex.size());
	}

	//	create vertex
	const T theta_interval = 2 * YZ_PI / slices;
	T theta = 0;
	for (int i = 0; i < slices; i++, iter_start++) {
		(*iter_start) = Vec3<T>(cos(theta), sin(theta), 0);
		theta += theta_interval;
	}

	return 1;
}

/**
	Create a unit circle on XY plane

	\param	vertex		return the vertex list
	\param	slices		how many slices of the ring
	\return				whether create succeed
*/
template<typename T>
int createUnitCircleXY(
	std::vector<Vec3<T>>&	vertex,
	int						slices = 16
) {
	if (slices < 3) {
		vertex.clear();
		return 0;
	}

	//	create vertex
	vertex.resize(slices);
	return createUnitCircleXY(vertex, vertex.begin(), slices);
}

/**
	Create a unit tri-mesh cylinder,
	axis in +z direction with length 1, diameter start from +x, radius 1 (max radius)

	\param	vertex		return the vertex list
	\param	face		return the face list
	\param	slices		how many slices of the ring
	\return				whether create succeed
*/
template<typename T>
int createUnitCylinderTriMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face,
	int						slices = 16
) {
	if (slices < 3) {
		vertex.clear();
		face.clear();
		return 0;
	}

	//	create vertex
	vertex.reserve(slices * 2);
	createUnitCircleXY(vertex, slices);		//	create a planar circle, then duplicate
	vertex.insert(vertex.end(), vertex.begin(), vertex.end());
	for (int i = slices; i < vertex.size(); i++) {
		vertex[i].z = 1;
	}

	//	create face
	face.resize(slices * 2);
	for (int i = 0; i < slices - 1; i++) {
		face[i << 1] = int3(i, i + 1, i + 1 + slices);						//	i*2
		face[(i << 1) | 0x01] = int3(i, i + 1 + slices, i + slices);		//	i*2+1
	}
	face[(slices << 1) - 2] = int3(slices - 1, 0, slices);					//	loop the ring
	face[(slices << 1) - 1] = int3(slices - 1, slices, (slices << 1) - 1);

	return 1;
}

/**
	Create a TriMesh cylinder

	\param	vertex		return the vertex list
	\param	face		return the face list, corresponding to v (start from 0)
	\param	v1			the center of start circle of the cylinder
	\param	v2			the center of end circle of the cylinder
	\param	radius		the radius of the cylinder (max radius)
	\param	slices		how many slices to form a circle (mininum 3)
	\param	diameter_r	the vector where the circle begins, for a zero length vector, we generate a vector automatically
	\return				whether the creation succeed
*/
template<typename T>
int createCylinderTriMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face,
	Vec3<T>					v1,
	Vec3<T>					v2,
	T						radius,
	int						slices = 16,
	Vec3<T>					diameter_r = Vec3d(0, 0, 0)
) {
	if (slices < 3) {	//	cannot form a mesh
		vertex.clear();
		face.clear();
		return 0;
	}

	//	if the cylinder is too short, just return an empty mesh
	Vec3<T> rz = v2 - v1;
	if (rz.Length() < T(1e-5)) {
		vertex.clear();
		face.clear();
		return 1;
	}

	//	calculate transform matrix
	int need_random_diameter = 0;
	if (diameter_r.Length() < 1e-5)
		need_random_diameter = 1;
	if (!need_random_diameter) {
		T ang = angleDegBetweenVectors(rz, diameter_r);
		if (ang < T(0.1) || ang > T(179.9))
			need_random_diameter = 1;
	}
	if (need_random_diameter) {
		diameter_r = Vec3<T>(1, 0, 0);
		T ang = angleDegBetweenVectors(rz, diameter_r);
		if (ang < T(0.1) || ang > T(179.9))
			diameter_r = Vec3<T>(0, 1, 0);
	}
	Vec3<T> ry = cross(rz, diameter_r);
	ry.SetNormalize();
	ry *= radius;
	Vec3<T> rx = cross(ry, rz);
	rx.SetNormalize();
	rx *= radius;
	yz::Matrix3x3<T> R(rx, ry, rz);

	//	create the mesh
	createUnitCylinderTriMesh(vertex, face, slices);
	for (int i = 0; i < vertex.size(); i++) {
		vertex[i] = R * vertex[i] + v1;
	}

	return 1;
}

/**
	Create a tube from a curve

	We don't check whether the curve is legal. The curve should be smooth, or the mesh may contain self-intersection

	\param	vertex		return the vertex list
	\param	face		return the face list, corresponding to v (start from 0)
	\param	curve		vertex list of the curve, the input curve must be legal
	\param	radius		the radius of the cylinder at each curve point.
						if size of this array is smaller than curve point number, radius will be constant after the last point
	\param	slices		how many slices to form a circle (mininum 3)
	\param	diameter_r	the vector where the circle begins (start at the first point of the curve),
						for a zero length vector, we generate a vector automatically
	\return				whether the creation succeed
*/
template<typename T>
int createTubeTriMesh(
	std::vector<Vec3<T>>&		vertex,
	std::vector<int3>&			face,
	const std::vector<Vec3<T>>&	curve,
	const std::vector<T>&		radius,
	int							slices = 16,
	Vec3<T>						diameter_r = Vec3d(0, 0, 0)
) {
	if (slices < 3 || curve.size() < 2 || radius.empty())
		return 0;

	//	create unit circle for each curve node
	vertex.clear();
	vertex.reserve(curve.size()*slices);
	createUnitCircleXY(vertex, slices);
	for (int i = 1; i < curve.size(); i++) {
		vertex.insert(vertex.end(), vertex.begin(), vertex.begin() + slices);
	}

	//	find a legal diameter_r
	int need_random_diameter = 0;
	if (diameter_r.Length() < 1e-5)
		need_random_diameter = 1;
	if (!need_random_diameter) {
		T ang = angleDegBetweenVectors(curve[1] - curve[0], diameter_r);
		if (ang < T(0.1) || ang > T(179.9))
			need_random_diameter = 1;
	}
	if (need_random_diameter) {
		diameter_r = Vec3<T>(1, 0, 0);
		T ang = angleDegBetweenVectors(curve[1] - curve[0], diameter_r);
		if (ang < T(0.1) || ang > T(179.9))
			diameter_r = Vec3<T>(0, 1, 0);
	}

	//	update the circle for each joint of the curve
	for (int i = 0; i < curve.size(); i++) {
		Vec3<T> v_prev = (i == 0 ? curve[0] * 2 - curve[1] : curve[i - 1]);
		Vec3<T> v_curr = curve[i];
		Vec3<T> v_next = (i + 1 == curve.size() ? curve[i] * 2 - curve[i - 1] : curve[i + 1]);

		//	calculate plane normal (+z unit vector of local coordinate)
		Vec3<T> rz;
		Vec3<T> r1 = (v_prev - v_curr).Normalize();
		Vec3<T> r2 = (v_next - v_curr).Normalize();
		Vec3<T> r_plane = cross(r1, r2);
		if (r_plane.Length() < 1e-5) //	the two line segments are almost colinear
			rz = r2;
		else
			rz = cross(r_plane, r1 + r2).Normalize();

		//	calculate rx, ry perpendicular to r1, length is same as radius
		T curr_radius = (i >= radius.size() ? radius.back() : radius[i]);
		Vec3<T> rx = (diameter_r - r1 * dot(diameter_r, r1)).Normalize() * curr_radius;
		Vec3<T> ry = cross(rx, r1).Normalize() * curr_radius;

		//	project rx, ry to rz plane in r1 direction
		rx = rx - r1 * (dot(rx, rz) / dot(r1, rz));
		ry = ry - r1 * (dot(ry, rz) / dot(r1, rz));

		//	calculate the transform matrix
		yz::Matrix3x3<T> R(rx, ry, rz);

		//	transform the circle
		for (int j = 0; j < slices; j++) {
			vertex[i*slices + j] = R * vertex[i*slices + j] + v_curr;
		}

		//	update diameter_r
		diameter_r = rx;
	}

	//	create face list
	face.resize(slices * 2 * (curve.size() - 1));
	for (int j = 0; j < curve.size() - 1; j++) {
		for (int i = 0; i < slices - 1; i++) {
			face[(slices << 1) * j + (i << 1)] = int3(slices * j + i, slices * j + i + 1, slices * j + i + 1 + slices);					//	i*2
			face[(slices << 1) * j + ((i << 1) | 0x01)] = int3(slices * j + i, slices * j + i + 1 + slices, slices * j + i + slices);		//	i*2+1
		}
		face[(slices << 1) * j + (slices << 1) - 2] = int3(slices * j + slices - 1, slices * j, slices * j + slices);					//	loop the ring
		face[(slices << 1) * j + (slices << 1) - 1] = int3(slices * j + slices - 1, slices * j + slices, slices * j + (slices << 1) - 1);
	}

	return 1;
}

/**
	Create a tube from a curve

	We don't check whether the curve is legal. The curve should be smooth, or the mesh may contain self-intersection

	\param	vertex		return the vertex list
	\param	face		return the face list, corresponding to v (start from 0)
	\param	curve		vertex list of the curve, the input curve must be legal
	\param	radius		the radius of the cylinder (max radius)
	\param	slices		how many slices to form a circle (mininum 3)
	\param	diameter_r	the vector where the circle begins (start at the first point of the curve),
						for a zero length vector, we generate a vector automatically
	\return				whether the creation succeed
*/
template<typename T>
int createTubeTriMesh(
	std::vector<Vec3<T>>&		vertex,
	std::vector<int3>&			face,
	const std::vector<Vec3<T>>&	curve,
	T							radius,
	int							slices = 16,
	Vec3<T>						diameter_r = Vec3d(0, 0, 0)
) {
	std::vector<T> tmp_radius;
	tmp_radius.push_back(radius);

	std::vector<Vec3<T>> tmp_diameter_r;
	tmp_diameter_r.push_back(diameter_r);

	return createTubeTriMesh(vertex, face, curve, tmp_radius, slices, diameter_r);
}

/**
	Create a unit radius tri-mesh sphere

	\param	vertex		return the vertex list
	\param	face		return the face list
	\param	slices		how many slices of the sphere
	\param	stacks		how many stacks of the sphere
	\return				whether create succeed
*/
template<typename T>
int createUnitSphereTriMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face,
	int						slices = 10,
	int						stacks = 10
) {
	if (slices < 3 || stacks < 2) {
		vertex.clear();
		face.clear();
		return 0;
	}

	//	create vertex
	vertex.reserve(slices * (stacks - 1) + 2);
	createUnitCircleXY(vertex, slices);		//	the unit circle
	for (int j = 2; j < stacks; j++)		//	copy unit circle
		vertex.insert(vertex.end(), vertex.begin(), vertex.begin() + slices);
	vertex.push_back(Vec3<T>(0, 0, -1));	//	polar points of the sphere
	vertex.push_back(Vec3<T>(0, 0, 1));

	for (int j = 1; j < stacks; j++) {		//	update each ring
		T z = -cos(YZ_PI * j / stacks);
		T r = sin(YZ_PI * j / stacks);
		for (int i = 0; i < slices; i++) {
			vertex[slices * (j - 1) + i] *= r;
			vertex[slices * (j - 1) + i].z = z;
		}
	}

	//	create face
	face.resize(slices * (stacks - 1) * 2, 0);
	if (stacks > 2) {
		//	the base cylinder
		for (int i = 0; i < slices - 1; i++) {
			face[i << 1] = int3(i, i + 1, i + 1 + slices);						//	i*2
			face[(i << 1) | 0x01] = int3(i, i + 1 + slices, i + slices);		//	i*2+1
		}
		face[(slices << 1) - 2] = int3(slices - 1, 0, slices);					//	loop the ring
		face[(slices << 1) - 1] = int3(slices - 1, slices, (slices << 1) - 1);

		//	duplicate base cylinder
		for (int j = 1; j < stacks - 2; j++) {
			int idx_inc = slices * j;
			for (int i = 0; i < slices * 2; i++) {
				face[j*slices * 2 + i] = face[i];
				face[j*slices * 2 + i].x += idx_inc;
				face[j*slices * 2 + i].y += idx_inc;
				face[j*slices * 2 + i].z += idx_inc;
			}
		}
	}

	//	the bottom cone
	int face_offset = (stacks - 2)*slices * 2;
	for (int i = 0; i < slices - 1; i++)
		face[face_offset + i] = int3(i, vertex.size() - 2, i + 1);
	face[face_offset + slices - 1] = int3(slices - 1, vertex.size() - 2, 0);

	//	the top cone
	face_offset = (stacks - 2)*slices * 2 + slices;
	int vertex_offset = slices * (stacks - 2);
	for (int i = 0; i < slices - 1; i++)
		face[face_offset + i] = int3(vertex_offset + i, vertex_offset + i + 1, vertex.size() - 1);
	face[face_offset + slices - 1] = int3(vertex_offset + slices - 1, vertex_offset, vertex.size() - 1);

	return 1;
}

/**
	Create a radius tri-mesh sphere

	\param	vertex		return the vertex list
	\param	face		return the face list
	\param	center		center of the sphere
	\param	radius		radius of the sphere
	\param	slices		how many slices of the sphere
	\param	stacks		how many stacks of the sphere
	\return				whether create succeed
*/
template<typename T>
int createSphereTriMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face,
	Vec3<T>					center,
	T						radius,
	int						slices = 10,
	int						stacks = 10
) {
	if (!createUnitSphereTriMesh(vertex, face, slices, stacks))
		return 0;

	for (int i = 0; i < vertex.size(); i++) {
		vertex[i] *= radius;
		vertex[i] += center;
	}

	return 1;
}

/**
Create an axis aligned bonding box with triangle mesh

\param	vertex		return the vertex list
\param	face		return the face list, corresponding to v (start from 0)
\param	aabb_min	the min coordinate of aabb
\param	aabb_max	the max coordinate of aabb
\return				whether the creation succeed
*/
template<typename T>
int createAABBTriMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int3>&		face,
	Vec3<T>					aabb_min,
	Vec3<T>					aabb_max
) {
	vertex.resize(8);
	vertex[0] = aabb_min;
	vertex[1] = Vec3<T>(aabb_max.x, aabb_min.y, aabb_min.z);
	vertex[2] = Vec3<T>(aabb_max.x, aabb_max.y, aabb_min.z);
	vertex[3] = Vec3<T>(aabb_min.x, aabb_max.y, aabb_min.z);
	vertex[4] = Vec3<T>(aabb_min.x, aabb_min.y, aabb_max.z);
	vertex[5] = Vec3<T>(aabb_max.x, aabb_min.y, aabb_max.z);
	vertex[6] = aabb_max;
	vertex[7] = Vec3<T>(aabb_min.x, aabb_max.y, aabb_max.z);

	face.resize(12);
	face[0] = int3(1, 2, 6);
	face[1] = int3(1, 6, 5);
	face[2] = int3(0, 4, 7);
	face[3] = int3(0, 7, 3);
	face[4] = int3(2, 3, 7);
	face[5] = int3(2, 7, 6);
	face[6] = int3(0, 1, 5);
	face[7] = int3(0, 5, 4);
	face[8] = int3(4, 5, 6);
	face[9] = int3(4, 6, 7);
	face[10] = int3(0, 3, 2);
	face[11] = int3(0, 2, 1);

	return 1;
}

/**
Create an axis aligned bonding box with quad mesh

\param	vertex		return the vertex list
\param	face		return the face list, corresponding to v (start from 0)
\param	aabb_min	the min coordinate of aabb
\param	aabb_max	the max coordinate of aabb
\return				whether the creation succeed
*/
template<typename T>
int createAABBQuadMesh(
	std::vector<Vec3<T>>&	vertex,
	std::vector<int4>&		face,
	Vec3<T>					aabb_min,
	Vec3<T>					aabb_max
) {
	vertex.resize(8);
	vertex[0] = aabb_min;
	vertex[1] = Vec3<T>(aabb_max.x, aabb_min.y, aabb_min.z);
	vertex[2] = Vec3<T>(aabb_max.x, aabb_max.y, aabb_min.z);
	vertex[3] = Vec3<T>(aabb_min.x, aabb_max.y, aabb_min.z);
	vertex[4] = Vec3<T>(aabb_min.x, aabb_min.y, aabb_max.z);
	vertex[5] = Vec3<T>(aabb_max.x, aabb_min.y, aabb_max.z);
	vertex[6] = aabb_max;
	vertex[7] = Vec3<T>(aabb_min.x, aabb_max.y, aabb_max.z);

	face.resize(6);
	face[0] = int4(1, 2, 6, 5);
	face[1] = int4(0, 4, 7, 3);
	face[2] = int4(2, 3, 7, 6);
	face[3] = int4(0, 1, 5, 4);
	face[4] = int4(4, 5, 6, 7);
	face[5] = int4(0, 3, 2, 1);

	return 1;
}


}}	//	namespace yz::geometry


#endif	//	__YZ_MESH_GENERATION_H__