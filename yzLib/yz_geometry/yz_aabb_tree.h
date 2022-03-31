/***********************************************************/
/**	\file
	\brief		AABB Tree
	\author		Yizhong Zhang
	\date		9/3/2012
*/
/***********************************************************/
#ifndef __YZ_AABB_TREE_H__
#define __YZ_AABB_TREE_H__

#include <iostream>
#include <vector>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_geometry/yz_bonding_box.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"		//	include this file for display AABB Tree
#endif

namespace yz{	namespace geometry{

/**
AABB Tree

Template is required to specify vector representation, such as Vec3<double>, Vec2<float>

An AABB tree is full binary tree. If the tree has N leaves,
then there are 2N-1 nodes in total.

We store node in the following order:	\n
0 ~ N-1,	leaf nodes	\n
N,			root node	\n
N+1 ~ 2N-1,	other non-leaf nodes
*/
template<class VEC_NT>
class AABBTree {
public:
	AABBTree() {
		leaf_number = 0;
	}

	/**
	Build AABB Tree for faces of a mesh

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the mesh
	\param	face			face of the mesh
	\param	expand_value	how much should each AABB expand
	*/
	inline void BuildTriangleAABBTree(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int3>&	face,
		typename VEC_NT::Type		expand_value = 0
	) {
		node.clear();
		SetupLeavesForTriangles(vertex, face, expand_value);
		BuildAABBTreeFromLeafNodes();
	}

	/**
	Build AABB Tree for edges of a mesh

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the mesh
	\param	edge			edge of the mesh
	\param	expand_value	how much should each AABB expand
	*/
	inline void BuildEdgeAABBTree(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int2>&	edge,
		typename VEC_NT::Type		expand_value = 0
	) {
		node.clear();
		SetupLeavesForEdges(vertex, edge, expand_value);
		BuildAABBTreeFromLeafNodes();
	}

	/**
	Build AABB Tree for a point cloud

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the point cloud
	\param	expand_value	how much should each AABB expand
	*/
	inline void BuildVertexAABBTree(
		const std::vector<VEC_NT>&	vertex,
		typename VEC_NT::Type		expand_value = 0
	) {
		node.clear();
		SetupLeavesForVertices(vertex, expand_value);
		BuildAABBTreeFromLeafNodes();
	}

	/**
	Update AABB Tree for faces of a mesh

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the mesh
	\param	face			face of the mesh
	\param	expand_value	how much should each AABB expand
	*/
	inline void UpdateTriangleAABBTree(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int3>&	face,
		typename VEC_NT::Type		expand_value = 0
	) {
		if (face.size() != leaf_number) {
#ifndef BE_QUIET
			std::cout << "error: AABBTree::UpdateTriangleAABBTree, face number doesn't match aabb face number" << std::endl;
#endif
			return;
		}

		SetupLeavesForTriangles(vertex, face, expand_value);

		if (leaf_number <= 1)
			return;

		UpdateHierarchy(leaf_number);
	}

	/**
	Update AABB Tree for edges of a mesh

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the mesh
	\param	edge			edge of the mesh
	\param	expand_value	how much should each AABB expand
	*/
	inline void UpdateEdgeAABBTree(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int2>&	edge,
		typename VEC_NT::Type		expand_value = 0
	) {
		if (edge.size() != leaf_number) {
#ifndef BE_QUIET
			std::cout << "error: AABBTree::UpdateEdgeAABBTree, edge number doesn't match aabb edge number" << std::endl;
#endif
			return;
		}

		SetupLeavesForEdges(vertex, edge, expand_value);

		if (leaf_number <= 1)
			return;

		UpdateHierarchy(leaf_number);
	}

	/**
	Update AABB Tree for point cloud

	by specifying expand value, each bonding box can be a little larger
	than the minimal bonding box

	\param	vertex			vertex of the point cloud
	\param	expand_value	how much should each AABB expand
	*/
	inline void UpdateVertexAABBTree(
		const std::vector<VEC_NT>&	vertex,
		typename VEC_NT::Type		expand_value = 0
	) {
		if (vertex.size() != leaf_number) {
#ifndef BE_QUIET
			std::cout << "error: AABBTree::UpdateVertexAABBTree, vertex number doesn't match aabb leaf number" << std::endl;
#endif
			return;
		}

		SetupLeavesForVertices(vertex, expand_value);

		if (leaf_number <= 1)
			return;

		UpdateHierarchy(leaf_number);
	}

	/**
	Setup leaves by triangles of a triangle mesh

	\param	vertex			vertex list of the mesh
	\param	face			face vertex index list
	\param	expand_value	how much should each AABB expand
	*/
	inline void SetupLeavesForTriangles(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int3>&	face,
		typename VEC_NT::Type		expand_value = 0
	) {
		leaf_number = face.size();
		if (node.size() < leaf_number)	//	node must be big enough
			node.resize(leaf_number);

		std::vector<VEC_NT> tmp;
		tmp.resize(3);
		for (int i = 0; i<face.size(); i++) {
			tmp[0] = vertex[face[i].x];
			tmp[1] = vertex[face[i].y];
			tmp[2] = vertex[face[i].z];
			node[i].GetAABBCoef(tmp);
			for (int j = 0; j < VEC_NT::dims; j++) {
				node[i].bb_min[j] -= expand_value;
				node[i].bb_max[j] += expand_value;
			}
		}
	}

	/**
	Setup leaves by edge of a triangle mesh

	\param	vertex			vertex list of the mesh
	\param	edge			edge vertex index list
	\param	expand_value	how much should each AABB expand
	*/
	inline void SetupLeavesForEdges(
		const std::vector<VEC_NT>&	vertex,
		const std::vector<int2>&	edge,
		typename VEC_NT::Type		expand_value = 0
	) {
		leaf_number = edge.size();
		if (node.size() < leaf_number)	//	node must be big enough
			node.resize(leaf_number);

		for (int i = 0; i<edge.size(); i++) {
			node[i].bb_min = vertex[edge[i].x];
			node[i].bb_max = vertex[edge[i].y];
			for (int j = 0; j < VEC_NT::dims; j++) {
				if (node[i].bb_min[j] > node[i].bb_max[j])
					mySwap(node[i].bb_min[j], node[i].bb_max[j]);
			}
			for (int j = 0; j < VEC_NT::dims; j++) {
				node[i].bb_min[j] -= expand_value;
				node[i].bb_max[j] += expand_value;
			}
		}
	}

	/**
	Setup leaves of point cloud

	\param	vertex			vertex list of point cloud
	\param	expand_value	how much should each AABB expand
	*/
	inline void SetupLeavesForVertices(
		const std::vector<VEC_NT>&	vertex,
		typename VEC_NT::Type		expand_value = 0
	) {
		leaf_number = vertex.size();
		if (node.size() < leaf_number)	//	node must be big enough
			node.resize(leaf_number);

		for (int i = 0; i < vertex.size(); i++) {
			node[i].bb_min = vertex[i];
			node[i].bb_max = vertex[i];
			for (int j = 0; j < VEC_NT::dims; j++) {
				node[i].bb_min[j] -= expand_value;
				node[i].bb_max[j] += expand_value;
			}
		}
	}

	/**
	Build AABB Tree after bonding boxes of leaf nodes have been created
	*/
	inline int BuildAABBTreeFromLeafNodes() {
		if (leaf_number <= 1)	//	no need to build tree
			return 1;
		if (node.size() < leaf_number) {
#ifndef BE_QUIET
			std::cout << "leaf node of AABB tree not setup" << std::endl;
#endif
			return 0;
		}

		//	setup vector
		node.resize(leaf_number);		//	 just save leaf nodes

		std::vector<AABBTreeNode> tmp_node = node;
		for (int i = 0; i<tmp_node.size(); i++)
			tmp_node[i].left = i;	//	record the original position of this node

		node.reserve(leaf_number * 2 - 1);

		//	recursively split the bonding box
		BuildHierarchy(&tmp_node[0], tmp_node.size());

		return 1;
	}

	/**
	Display bonding boxes of a given depth with wire cube

	\param	depth	the depth of node to display, 0: root
	*/
	inline int Display(int depth = 0) {
#ifdef YZ_gl_h
		DisplayWithDepth(leaf_number, depth);
		return 1;
#else
#ifndef BE_QUIET
		std::cout << "gl.h has to be included in order to use Display() in AABBTree3D" << std::endl;
#endif
		return 0;
#endif
	}
public:

	/**
	A node in AABB Tree 3D
	*/
	class AABBTreeNode : public AABB<VEC_NT> {
	public:
		int left, right;	//	two children

		AABBTreeNode() :left(-1), right(-1) {}	//	-1 indicate no child
	};

	int leaf_number;
	std::vector<AABBTreeNode> node;

protected:
	/**
	Recursively build the hierarchy of given leaf node array

	This function is used to construct AABB tree after
	bonding boxes of the leaf nodes have been calculated.

	\param	node_array		array of nodes
	\param	node_number		nodes contained in the array
	\return					id of current root
	*/
	inline int BuildHierarchy(AABBTreeNode* node_array, int node_number) {
		if (node_number == 1) {		//	leaf node
			return node_array[0].left;
		}

		//	calculate bonding box of the node array and push back the father node
		AABBTreeNode father_node = node_array[0];
		for (int i = 1; i<node_number; i++)
			father_node.Expand(node_array[i]);
		int father_node_id = node.size();
		node.push_back(father_node);

		//	find the longest edge of the bonding box
		int axis_id = 0;	//	0: x, 1: y, 2: z
		typename VEC_NT::Type length = 0;
		for (int i = 0; i < VEC_NT::dims; i++) {
			if (father_node.bb_max[i] - father_node.bb_min[i] > length) {
				length = father_node.bb_max[i] - father_node.bb_min[i];
				axis_id = i;
			}
		}

		//	cut the bonding box by its longest edge
		//	split the array into two parts, then recursively call this function
		typename VEC_NT::Type split_plane = father_node.bb_min[axis_id] + length * 0.5;
		int f = 0, r = node_number - 1;

		while (f < r) {	//	scan the array from head and tail, swap them if needed
			while (f < r && (node_array[f].bb_min[axis_id] + node_array[f].bb_max[axis_id]) * 0.5 < split_plane)
				f++;
			while (f < r && (node_array[r].bb_min[axis_id] + node_array[r].bb_max[axis_id]) * 0.5 > split_plane)
				r--;
			if (f < r) {
				mySwap(node_array[f], node_array[r]);
				f++;
				r--;
			}
		}

		if (f == 0 || r == node_number - 1) {	//	split anyway
			f = node_number / 2 - 1;
			r = f + 1;
		}
		else {
			if (f == r) {
				if ((node_array[f].bb_min[axis_id] + node_array[f].bb_max[axis_id]) * 0.5 < split_plane)
					r++;
				else
					f--;
			}
			else if (f > r) {
				f--;
				r++;
			}
		}

		//	recursive call this function
		node[father_node_id].left = BuildHierarchy(node_array, f + 1);
		node[father_node_id].right = BuildHierarchy(node_array + f + 1, node_number - r);

		return father_node_id;
	}

	/**
	Recursively update the hierarchy of given leaf node array
	*/
	inline void UpdateHierarchy(int node_id) {
		if (node_id < leaf_number)
			return;

		int left_id = node[node_id].left;
		int right_id = node[node_id].right;
		UpdateHierarchy(left_id);
		UpdateHierarchy(right_id);
		for (int i = 0; i < VEC_NT::dims; i++) {
			node[node_id].bb_min[i] = myMin(node[left_id].bb_min[i], node[right_id].bb_min[i]);
			node[node_id].bb_max[i] = myMax(node[left_id].bb_max[i], node[right_id].bb_max[i]);
		}
	}

	/**
	resursive call display. This function is only valid if gl.h included

	\param	node_id		current node id to draw
	\param	depth		current depth to draw
	*/
	inline void DisplayWithDepth(int node_id, int depth) {
#ifdef YZ_gl_h
		if (node_id < 0 || node_id >= node.size() || depth < 0)		//	input must be legal
			return;

		if (depth == 0)		//	draw this bonding box
			opengl::drawAABBWire(node[node_id].bb_min, node[node_id].bb_max);
		else {					//	draw its children
			DisplayWithDepth(node[node_id].left, depth - 1);
			DisplayWithDepth(node[node_id].right, depth - 1);
		}
#endif
	}
};


/**
2D AABB Tree
*/
template<class T>
class AABBTree2D : public AABBTree<Vec2<T>> {
};


/**
3D AABB Tree
*/
template<class T>
class AABBTree3D : public AABBTree<Vec3<T>> {
};


//	------------------------------------
//	type define
//	------------------------------------
typedef AABBTree2D<float>	AABBTree2Df;
typedef AABBTree2D<double>	AABBTree2Dd;

typedef AABBTree3D<float>	AABBTree3Df;
typedef AABBTree3D<double>	AABBTree3Dd;


}}	//	namespace yz::geometry

#endif	//	__YZ_AABB_TREE_H__