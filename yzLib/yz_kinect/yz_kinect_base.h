/***********************************************************/
/**	\file
	\brief		Base class for Kinect
	\details	There are several versions of kinect driver,
				such as OpenNI, freenect. All of them have
				same basic functions. This file is the basic 
				abstraction of all the kinect libraries.
	\author		Yizhong Zhang
	\date		11/7/2012
*/
/***********************************************************/
#ifndef __YZ_KINECT_BASE_H__
#define __YZ_KINECT_BASE_H__

namespace yz{	namespace kinect{

/**
	Kinect Base

	Control parameters are included in this class. However, 
	it is possible that the driver doesn't support certain
	parameters. In this case, that parameter is just ignored.

	\todo	the whole class is a mass
*/
template<class T>
class Kinect{
public:
	//	color  map
	unsigned char*	color_map;		//	24 bit image
	int	color_width, color_height;

	//	depth map
	unsigned short*	depth_map;		//	depth in millimeters
	int depth_width, depth_height;
	int	mapped_flag;

	//	point cloud in 3D space arranged in xyz
	T*	point_cloud;				//	dimension is the same as depth map
	unsigned char* point_cloud_color;

public:
	Kinect(){
		color_width		= 640;
		color_height	= 480;
		color_map		= new unsigned char[color_width*color_height*3];
		depth_width		= 640;
		depth_height	= 480;
		mapped_flag		= 1;
		depth_map		= new unsigned short[depth_width*depth_height];
		point_cloud		= new T[depth_width*depth_height*3];

		point_cloud_color	= new unsigned char[depth_width*depth_height*3];

		SetParametersToDefault();
	}

	/**
		Initialize kinect
		
		multi-kinect is allowed, in this case we initialize the kinect with given index

		\param	kinect_idx		index of kinect we want to initialize, start from 0 to 
								(number of kinects - 1)
		\return					whether initialization succeed, 1: succeed, 0: failed
	*/
	virtual int InitKinect(int kinect_idx = 0) = 0;

	/**
	*/
	virtual int GetPointCloud() = 0;

public:
	/**
		set parameters to default

		by default, we just use kinect as a RGBD camera
	*/
	void SetParametersToDefault(){
		color_flag			= 1;
		depth_flag			= 1;
		skeleton_flag		= 0;
		audio_flag			= 0;
		player_index_flag	= 0;
		hd_color_flag		= 0;
	}

protected:
	//	control parameters, we set them to be protected because
	//	more steps have to be taken to change the state of kinect
	int	color_flag;			///<	whether enable RGB image
	int depth_flag;			///<	whether enable depth map
	int skeleton_flag;		///<	whether enable skeleton
	int audio_flag;			///<	whether enable audio
	int player_index_flag;	///<	whether enable player segmentation on depth map
	int hd_color_flag;		///<	whether enable high quality RGB image, at the cost of drop frame rate


};


}}	//	namespace yz::kinect

#endif //	__YZ_KINECT_BASE_H__