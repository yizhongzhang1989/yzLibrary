/***********************************************************/
/**	\file
	\brief		Sphere Tree
	\author		Yizhong Zhang
	\date		10/23/2012
*/
/***********************************************************/
#ifndef __YZ_SPHERE_TREE_H__
#define __YZ_SPHERE_TREE_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "yzLib/yz_setting.h"
#include "yzLib/yz_math/yz_numerical_utils.h"
#include "yzLib/yz_geometry/yz_sphere.h"
#include "yzLib/yz_utils/yz_string_utils.h"

#ifdef YZ_gl_h
#	include "yzLib/yz_opengl/yz_vector_opengl_utils.h"	
#endif

namespace yz{	namespace geometry{

/**
*/
template<class T>
class SphereTree{
public:
	/**
		constructer, initial depth is -1 
	*/
	SphereTree() : tree_depth(-1), branching_factor(0) {};

	/**
		read sphere tree from .sph file

		.sph file is defined by Gareth Bradshaw.

		Refer to paper: \n
		Sphere-tree construction using dynamic medial axis approximation,
		Gareth Bradshaw, Carol O'Sullivan, SCA 2002
	*/
	int ReadSphereTree(const char* file_name){
		std::string extension = utils::getFileExtensionFromString(file_name);
		std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

		if( extension != "sph" ){
			#ifndef BE_QUIET
				std::cout << "error: SphereTree::ReadSphereTree, unsupported file format" << std::endl;
			#endif
			return 0;
		}

		std::ifstream sph( file_name );	//	load obj failed, do not touch old data
		if(!sph.is_open()){
			#ifndef BE_QUIET
				std::cout << "error: SphereTree::ReadSphereTree, cannot open " << file_name << std::endl;
			#endif
			return 0;
		}

		sph >> tree_depth >> branching_factor;
		tree_depth --;		//	root is treated to be depth 0, but read data is the depth count, so we need to minus 1
		int node_number = (power(branching_factor, tree_depth+1) - 1) / (branching_factor - 1);
		node.resize( node_number );

		for(int i=0; i<node_number; i++){
			T	importance;
			sph >> node[i].center[0] >> node[i].center[1] >> node[i].center[2] >> node[i].radius >> importance;
		}

		sph.close();

		return 1;
	}

	/**
		write to .sph file
	*/
	int WriteSphereTree(const char* file_name){
		std::string extension = utils::getFileExtensionFromString(file_name);
		std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

		if( extension != "sph" ){
			#ifndef BE_QUIET
				std::cout << "error: SphereTree::WriteSphereTree, unsupported file format" << std::endl;
			#endif
			return 0;
		}

		std::ofstream sph( file_name );	//	load obj failed, do not touch old data
		if(!sph.is_open()){
			#ifndef BE_QUIET
				std::cout << "error: SphereTree::WriteSphereTree, cannot open " << file_name << std::endl;
			#endif
			return 0;
		}

		sph << tree_depth+1 << ' ' << branching_factor << '\n';

		for(int i=0; i<node.size(); i++){
			sph << node[i].center[0] << ' '
				<< node[i].center[1] << ' '
				<< node[i].center[2] << ' '
				<< node[i].radius << ' '
				<< 1.0 << '\n';
		}

		sph.close();

		return 1;
	}

	/**
	*/
	inline void MarkIllegalSpheres(T threshold){
		for(int i=0; i<NodeNumber(); i++){
			if( node[i].radius < threshold )
				node[i].legal_flag = 0;
			else
				node[i].legal_flag = 1;
		}
	}

	/**
	*/
	inline int NodeNumber(){
		return node.size();
	}

	/**
	*/
	inline int FirstChildID(int father_id){
		assert(father_id>=0 && father_id<NodeNumber());
		return father_id * branching_factor + 1;
	}

	/**
	*/
	inline int FatherID(int child_id){
		assert(child_id>=0 && child_id<NodeNumber());
		return (child_id - 1) / branching_factor;
	}

	/**
		Display a certain depth of the sphere-tree

		\param	depth		depth id to display
		\param	shape_mode	0: solid shpere		\n
							1: wire sphere
		\param	color_mode	0: default color, color set before this function	\n
							1: same random color for nodes of same father		\n
							2: node of the same father have different color, but color pattern repeat	\n
							3: all random color
	*/
	int Display(int depth=0, int shape_mode=0, int color_mode=0){
		#ifdef YZ_gl_h
			if( depth<0 || depth>tree_depth )
				return 1;

			int start_id	= ( depth==0 ? 0 : (power(branching_factor, depth) - 1) / (branching_factor - 1) );
			int end_id		= (power(branching_factor, depth+1) - 1) / (branching_factor - 1);

			if( color_mode == 0 ){			//	default color
				for(int i=start_id; i<end_id; i++){
					if( shape_mode )	opengl::drawPointAsWireSphere(node[i].center, node[i].radius);
					else				opengl::drawPointAsSphere(node[i].center, node[i].radius);
				}
			}
			else if( color_mode == 1 ){		//	same color of same father
				opengl::srandColor(NodeNumber());
				int count = 0;
				for(int i=start_id; i<end_id; i++){
					if( count++ % branching_factor == 0 ){
						float rgb[3];
						opengl::randColor(rgb);
						glColor3fv(rgb);
					}
					if( shape_mode )	opengl::drawPointAsWireSphere(node[i].center, node[i].radius);
					else				opengl::drawPointAsSphere(node[i].center, node[i].radius);
				}
			}
			else if( color_mode == 2 ){		//	repeat color pattern
				opengl::srandColor(NodeNumber());
				//	create color array
				float* rgb = new float[3*branching_factor];
				for(int i=0; i<branching_factor; i++)
					opengl::randColor(rgb+i*3);

				for(int i=start_id; i<end_id; i++){
					glColor3fv(rgb+ i%branching_factor * 3 );
					if( shape_mode )	opengl::drawPointAsWireSphere(node[i].center, node[i].radius);
					else				opengl::drawPointAsSphere(node[i].center, node[i].radius);
				}
				delete[] rgb;
			}
			else if( color_mode == 3 ){		//	all random color
				opengl::srandColor(NodeNumber());
				for(int i=start_id; i<end_id; i++){
					float rgb[3];
					opengl::randColor(rgb);
					glColor3fv(rgb);
					if( shape_mode )	opengl::drawPointAsWireSphere(node[i].center, node[i].radius);
					else				opengl::drawPointAsSphere(node[i].center, node[i].radius);
				}
			}

			return 1;
		#else
			#ifndef BE_QUIET
				std::cout << "gl.h has to be included in order to use Display() in SphereTree" << std::endl;
			#endif
			return 0;
		#endif
	}

public:
	/**
		we add legal flag to each node, because the complete tree
		may contain illegal nodes
	*/
	class SphereTreeNode : public Sphere<T>{
	public:
		int	legal_flag;

		SphereTreeNode() : legal_flag(1) {};
	};

	int	tree_depth;
	int	branching_factor;
	std::vector<SphereTreeNode>	node;


};


}}	//	namespace yz::geometry

#endif	//	__YZ_SPHERE_TREE_H__