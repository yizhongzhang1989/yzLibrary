/***********************************************************/
/**	\file
	\brief		Microsoft Kinect with Microsoft Kinect SDK
	\author		Yizhong Zhang
	\date		11/7/2012
*/
/***********************************************************/
#ifndef __YZ_MICROSOFT_KINECT_H__
#define __YZ_MICROSOFT_KINECT_H__

#ifndef YZ_NuiApi_h
#	error yz_microsoft_kinect.h must be included after NuiApi.h
#endif
#ifndef YZ_windows_h	//	windows is used since the sdk is developed by microsoft
#	error yz_microsoft_kinect.h must be included after windows.h
#endif

#pragma comment(lib, "Kinect10.lib")

#include <iostream>
#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_kinect/yz_kinect_base.h"


namespace yz{	namespace kinect{

/**
	Microsoft Kinect
*/
template<class T>
class MicrosoftKinect : public Kinect<T>{
public:
	MicrosoftKinect(){
		pNuiSensor	= NULL;

		pColorStreamHandle = INVALID_HANDLE_VALUE;
		pDepthStreamHandle = INVALID_HANDLE_VALUE;
		hNextColorFrameEvent = INVALID_HANDLE_VALUE;
		hNextDepthFrameEvent = INVALID_HANDLE_VALUE;
		hNextSkeletonEvent = INVALID_HANDLE_VALUE;

		rgbx		= NULL;
	}

	~MicrosoftKinect(){
		if (NULL != pNuiSensor) {
			pNuiSensor->NuiShutdown();
			pNuiSensor->Release();
		}

		CloseHandle(hNextDepthFrameEvent);
		CloseHandle(hNextColorFrameEvent);
		CloseHandle(hNextSkeletonEvent);
	}
	/**
		Initialize kinect
		
		multi-kinect is allowed, in this case we initialize the kinect with given index

		\param	kinect_idx		index of kinect we want to initialize, start from 0 to 
								(number of kinects - 1)
		\return					whether initialization succeed, 1: succeed, 0: failed
	*/
	int InitKinect(int kinect_idx = 0){
		skeleton_flag		= 1;
		player_index_flag	= 1;

		//	check the number of connected sensor
		int iSensorCount = 0;
		HRESULT hr = NuiGetSensorCount(&iSensorCount);
		if( IsFailed(hr, "InitKinect: unable to count connected") )
			return 0;

		//	check whether given index is legal
		if( kinect_idx < 0 || kinect_idx >= iSensorCount ){	
			#ifndef BE_QUIET
				std::cout << "error: MicrosoftKinect, InitKinect:  " 
					<< iSensorCount << " kinects connected, index : "
					<< kinect_idx << " is not legal" << std::endl;
			#endif
			kinect_idx = 0;
		}

		//	create kinect by given index
		hr = NuiCreateSensorByIndex(kinect_idx, &pNuiSensor);
		if( IsFailed(hr, "InitKinect: unable to connect to kinect by index") )
			return 0;

		// Get the status of the sensor, and if connected, then we can initialize it
		hr = pNuiSensor->NuiStatus();
		if (S_OK != hr){
			pNuiSensor->Release();
			#ifndef BE_QUIET
				std::cout << "error: MicrosoftKinect, InitKinect: "
					<< "connection invalid, error code: " << hr << std::endl;
			#endif
			return 0;
		}

		//	initialize parameters 
		int init_flag = GetInitializeFlag();
		hr = pNuiSensor->NuiInitialize(init_flag); 
		if( IsFailed(hr, "InitKinect: initialization failed") )
			return 0;

		//	create event and streams
		///<	\todo	this part is under construction, currently we just open color and depth stream
		hNextDepthFrameEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
		hr = pNuiSensor->NuiImageStreamOpen(
			NUI_IMAGE_TYPE_DEPTH,
			NUI_IMAGE_RESOLUTION_640x480,
			0,
			2,
			hNextDepthFrameEvent,
			&pDepthStreamHandle);
		if( IsFailed(hr, "InitKinect: open depth frame") )
			return 0;

		hNextColorFrameEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
		hr = pNuiSensor->NuiImageStreamOpen(
			NUI_IMAGE_TYPE_COLOR,
			NUI_IMAGE_RESOLUTION_640x480,
			0,
			2,
			hNextColorFrameEvent,
			&pColorStreamHandle );
		if( IsFailed(hr, "InitKinect: open color frame") )
			return 0;

		hNextSkeletonEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
		hr = pNuiSensor->NuiSkeletonTrackingEnable(
			hNextSkeletonEvent, 
			0); 
		if( IsFailed(hr, "InitKinect: open enalbe skeleton tracking") )
			return 0;


		return 1;
	}

	/**
		stop kinect and release every resources
	*/
	int StopKinect(){
		if (NULL != pNuiSensor) {
			pNuiSensor->NuiShutdown();
			pNuiSensor->Release();
		}
		CloseHandle(hNextDepthFrameEvent);
		CloseHandle(hNextColorFrameEvent);
		CloseHandle(hNextSkeletonEvent);

		pNuiSensor	= NULL;

		pColorStreamHandle = INVALID_HANDLE_VALUE;
		pDepthStreamHandle = INVALID_HANDLE_VALUE;
		hNextColorFrameEvent = INVALID_HANDLE_VALUE;
		hNextDepthFrameEvent = INVALID_HANDLE_VALUE;
		hNextSkeletonEvent = INVALID_HANDLE_VALUE;

		rgbx		= NULL;

		return 1;
	}
	/**
	*/
	int GetColorMap(){
		NUI_IMAGE_FRAME imageFrame;

		HRESULT hr = pNuiSensor->NuiImageStreamGetNextFrame(pColorStreamHandle, 0, &imageFrame);
		if( hr == E_NUI_FRAME_NO_DATA ){	//	new frame not arrived yet, this is not an error
			return -1;
		}
		if( IsFailed(hr, "GetColorMap: get next frame") )
			return 0;

		NUI_LOCKED_RECT LockedRect;
		hr = imageFrame.pFrameTexture->LockRect(0, &LockedRect, NULL, 0);
		if( IsFailed(hr, "GetColorMap: lock rect") )
			return 0;

		if(!rgbx)	rgbx = new unsigned char[color_width*color_height*4];
		memcpy(rgbx, LockedRect.pBits, LockedRect.size);

		//	copy to color_map
		for(int i=0; i<color_width*color_height; i++){
			color_map[i*3  ] = rgbx[i*4+2];
			color_map[i*3+1] = rgbx[i*4+1];
			color_map[i*3+2] = rgbx[i*4  ];
		}

		hr = imageFrame.pFrameTexture->UnlockRect(0);
		if( IsFailed(hr, "GetColorMap: unlock rect") )
			return 0;

		hr = pNuiSensor->NuiImageStreamReleaseFrame(pColorStreamHandle, &imageFrame);
		if( IsFailed(hr, "GetColorMap: release frame") )
			return 0;

		return 1;
	}

	/**
	*/
	int GetDepthMap(){
		NUI_IMAGE_FRAME imageFrame;

		HRESULT hr = pNuiSensor->NuiImageStreamGetNextFrame(pDepthStreamHandle, 0, &imageFrame);
		if( hr == E_NUI_FRAME_NO_DATA ){	//	new frame not arrived yet, this is not an error
			return -1;
		}
		if( IsFailed(hr, "GetDepthMap: get next frame") )
			return 0;

		NUI_LOCKED_RECT LockedRect;
		hr = imageFrame.pFrameTexture->LockRect(0, &LockedRect, NULL, 0);
		if( IsFailed(hr, "GetDepthMap: lock rect") )
			return 0;

		memcpy(depth_map, LockedRect.pBits, LockedRect.size);

		hr = imageFrame.pFrameTexture->UnlockRect(0);
		if( IsFailed(hr, "GetDepthMap: unlock rect") )
			return 0;

		hr = pNuiSensor->NuiImageStreamReleaseFrame(pDepthStreamHandle, &imageFrame);
		if( IsFailed(hr, "GetDepthMap: release frame") )
			return 0;

		if( mapped_flag )
			MapDepthToColor();

		return 1;
	}

	/**
	*/
	int GetSkeleton(){
		NUI_SKELETON_FRAME skeletonFrame = {0};

		HRESULT hr = pNuiSensor->NuiSkeletonGetNextFrame(0, &skeletonFrame);
		if( E_NUI_FRAME_NO_DATA == hr )
			return 0;
		if( IsFailed(hr, "GetSkeleton: read from frame") ){
			return 0;
		}

		// smooth out the skeleton data
		pNuiSensor->NuiTransformSmooth(&skeletonFrame, NULL);

		//	find skeleton
		skeleton_quality = 0;
		for (int i = 0 ; i < NUI_SKELETON_COUNT; ++i){
			NUI_SKELETON_TRACKING_STATE trackingState = skeletonFrame.SkeletonData[i].eTrackingState;

			if (NUI_SKELETON_TRACKED == trackingState){
				skeleton_quality = 1;

				joint_pos.resize(NUI_SKELETON_POSITION_COUNT);
				joint_tracking_state.resize(NUI_SKELETON_POSITION_COUNT);

				NUI_SKELETON_DATA & skel = skeletonFrame.SkeletonData[i];
				for (int j = 0; j < NUI_SKELETON_POSITION_COUNT; ++j){
					//joint_pos[j].x = skel.SkeletonPositions[j].x + skel.Position.x;
					//joint_pos[j].y = skel.SkeletonPositions[j].y + skel.Position.y;
					//joint_pos[j].z = skel.SkeletonPositions[j].z + skel.Position.z;
					joint_pos[j].x = -skel.SkeletonPositions[j].x;
					joint_pos[j].y = skel.SkeletonPositions[j].y;
					joint_pos[j].z = -skel.SkeletonPositions[j].z;
					joint_tracking_state[j] = skel.eSkeletonPositionTrackingState[j];
				}

				break;	//	 just one skeleton
			}
		}

		return skeleton_quality;
	}

	/**
	*/
	int CalculatePointCloudColor(){
		static int		coord_count = 0;
		static long*	coord_ptr	= NULL;
		if( coord_count != depth_width*depth_height ){
			if( coord_ptr ){
				delete coord_ptr;
				coord_ptr = NULL;
			}
			coord_count	= depth_width*depth_height;
			coord_ptr	= new long[coord_count*2];
		}

		pNuiSensor->NuiImageGetColorPixelCoordinateFrameFromDepthPixelFrameAtResolution(
			NUI_IMAGE_RESOLUTION_640x480,
			NUI_IMAGE_RESOLUTION_640x480,
			coord_count,
			depth_map,
			coord_count*2,
			coord_ptr );

		for(int j=0; j<depth_height; j++)
			for(int i=0; i<depth_width; i++){
				int depthIndex = j*depth_width+i;
				LONG colorInDepthX = coord_ptr[depthIndex * 2];
				LONG colorInDepthY = coord_ptr[depthIndex * 2 + 1];

				if( colorInDepthX >= 0 && colorInDepthX < color_width && colorInDepthY >= 0 && colorInDepthY < color_height ){
					// calculate index into color array
					LONG colorIndex = colorInDepthX + colorInDepthY * color_width;

					//	write depth to the pixel that match color
					point_cloud_color[depthIndex*3  ] = color_map[colorIndex*3  ];
					point_cloud_color[depthIndex*3+1] = color_map[colorIndex*3+1];
					point_cloud_color[depthIndex*3+2] = color_map[colorIndex*3+2];
				}
				else{
					point_cloud_color[depthIndex*3  ] = 0;
					point_cloud_color[depthIndex*3+1] = 0;
					point_cloud_color[depthIndex*3+2] = 0;
				}

			}

		return 1;
	}

	/**
	*/
	int GetPointCloud(){
		const T DegreesToRadians = 3.14159265359f / 180.0f;
		float fov = mapped_flag ? NUI_CAMERA_COLOR_NOMINAL_HORIZONTAL_FOV : NUI_CAMERA_DEPTH_NOMINAL_HORIZONTAL_FOV;
		//float fov = NUI_CAMERA_DEPTH_NOMINAL_HORIZONTAL_FOV;
		T	xyScale = tanf(fov * DegreesToRadians * 0.5f) / (depth_width * 0.5f);
		int	half_width	= depth_width / 2;
		int	half_height	= depth_height / 2;
		for(int j=0; j<depth_height; j++)
			for(int i=0; i<depth_width; i++){
				int idx = j*depth_width+i;
				unsigned short pixel_depth = depth_map[idx] >> 3;
				T	depth = - pixel_depth * 0.001;	//	unit in meters
				//point_cloud[idx*3  ] = (i + 0.5 - half_width) * xyScale * depth + (mapped_flag ? 0.025 : 0);
				point_cloud[idx*3  ] = (i + 0.5 - half_width) * xyScale * depth;
				point_cloud[idx*3+1] = (j + 0.5 - half_height) * xyScale * depth;
				point_cloud[idx*3+2] = depth;		//	in OpenGL coordinate
			}

		//imdebug("rgb b=32f w=%d h=%d %p", depth_width, depth_height, point_cloud);

		return 1;
	}

	void Prj(int& u, int& v, T x, T y, T z){
		const T DegreesToRadians = 3.14159265359f / 180.0f;
		float fov = mapped_flag ? NUI_CAMERA_COLOR_NOMINAL_HORIZONTAL_FOV : NUI_CAMERA_DEPTH_NOMINAL_HORIZONTAL_FOV;
		T	xyScale = tanf(fov * DegreesToRadians * 0.5f) / (depth_width * 0.5f);
		int	half_width	= depth_width / 2;
		int	half_height	= depth_height / 2;

		u = (x ) / (xyScale * z) + half_width;
		//u = (x - (mapped_flag ? 0.025 : 0)) / (xyScale * z) + half_width;
		v = y / (xyScale * z) + half_height;
	}

protected:
	//	Kinect SDK Interface
	INuiSensor* pNuiSensor;
	HANDLE		pColorStreamHandle;
	HANDLE		pDepthStreamHandle;
	HANDLE		hNextColorFrameEvent;
	HANDLE		hNextDepthFrameEvent;
	HANDLE		hNextSkeletonEvent;

	//	intermediate values
	unsigned char*			rgbx;

public:
	int						skeleton_quality;
	std::vector<Vec3<T>>	joint_pos;
	std::vector<int>		joint_tracking_state;

//protected:
	/**
	*/
	int MapDepthToColor(){
		//ExplicitMapDepthToColor(depth_map, depth_width, color_width);
		//mapped_flag = 1;
		//return 1;

		static int		coord_count = 0;
		static long*	coord_ptr	= NULL;
		static unsigned short*	mapped_depth_map = NULL;
		if( coord_count != depth_width*depth_height ){
			if( coord_ptr ){
				delete coord_ptr;
				coord_ptr = NULL;
			}
			if( mapped_depth_map ){
				delete mapped_depth_map;
				mapped_depth_map = NULL;
			}
			coord_count	= depth_width*depth_height;
			coord_ptr	= new long[coord_count*2];
			mapped_depth_map = new unsigned short[coord_count];
		}

		memset(mapped_depth_map, 0, sizeof(unsigned short)*coord_count);

		pNuiSensor->NuiImageGetColorPixelCoordinateFrameFromDepthPixelFrameAtResolution(
			NUI_IMAGE_RESOLUTION_640x480,
			NUI_IMAGE_RESOLUTION_640x480,
			coord_count,
			depth_map,
			coord_count*2,
			coord_ptr );

		for(int j=0; j<depth_height; j++)
			for(int i=0; i<depth_width; i++){
				int depthIndex = j*depth_width+i;
				LONG colorInDepthX = coord_ptr[depthIndex * 2];
				LONG colorInDepthY = coord_ptr[depthIndex * 2 + 1];

				if( colorInDepthX >= 0 && colorInDepthX < color_width && colorInDepthY >= 0 && colorInDepthY < color_height ){
					// calculate index into color array
					LONG colorIndex = colorInDepthX + colorInDepthY * color_width;

					//	write depth to the pixel that match color
					if( !mapped_depth_map[colorIndex] ||									//	the pixel is not referenced
						mapped_depth_map[colorIndex]>>3 > depth_map[depthIndex]>>3 ){		//	current value is bigger
							mapped_depth_map[colorIndex] = depth_map[depthIndex];
					}
				}

			}

		//std::vector<float> prj;
		//prj.resize(coord_count*3);
		//for(int j=0; j<depth_height; j++){
		//	for(int i=0; i<depth_width; i++){
		//		int idx = j*depth_width+i;
		//		prj[idx*3] = coord_ptr[idx*2] - i;
		//		prj[idx*3+1] = coord_ptr[idx*2+1] - j;
		//		prj[idx*3+2] = 0;
		//	}
		//}

		memcpy(depth_map, mapped_depth_map, sizeof(unsigned short)*coord_count);

		mapped_flag = 1;

		return 1;
	}
	/**
	*/
	void ExplicitMapDepthToColor(unsigned short* depth_image, const int depth_width, const int color_width){
		if( (depth_width!=640 && depth_width!=320 && depth_width!=80) || 
			(color_width!=1280 && color_width!=640 && color_width!=320 && color_width!=80) ){
				printf("error: DumpDataProcesser::MapDepthToColor, illegal image size\n");
				return;
		}
		const int	depth_height	= (depth_width==640? 480 : (depth_width==320? 240 : 60));
		const int	color_height	= (color_width==1280? 960 : (color_width==640? 480 : (color_width==320? 240 : 60)));
		const float depth_horizontal_fov	= 58.5f;	//	data from Kinect SDK
		const float color_horizontal_fov	= 62.0f;
		const float	offset					= 25.0f;	//	distance between RGB and Depth camera 25mm

		unsigned short* mapped_depth = new unsigned short[depth_width*depth_height];
		memset(mapped_depth, 0, depth_width*depth_height*sizeof(unsigned short));

		//	calculate parameters
		const int	half_depth_width	= depth_width / 2;
		const int	half_depth_height	= depth_height / 2;
		const int	half_color_width	= color_width / 2;
		const int	half_color_height	= color_height / 2;
		const float DegreesToRadians = 3.14159265359f / 180.0f;
		const float	depth_xyScale = tanf(depth_horizontal_fov * DegreesToRadians * 0.5f) / (depth_width * 0.5f);
		const float inv_color_xyScale = (color_width * 0.5f) / tanf(color_horizontal_fov * DegreesToRadians * 0.5f);

		//	mapping
		for(int j=0; j<depth_height; j++){
			for(int i=0; i<depth_width; i++){
				int idx = j * depth_width + i;
				unsigned short value = depth_image[idx];
				if( !value )	continue;	//	this pixel don't have depth value, we skip it

				//	calculate coordinate in depth camera coordinate
				float	z = value >> 3;
				float	x = (i + 0.5f - half_depth_width) * depth_xyScale * z;
				float	y = (j + 0.5f - half_depth_height) * depth_xyScale * z;
				//	transform to color camera coordinate
				x += offset;
				//	project onto color image
				int		u = (x / z) * inv_color_xyScale + half_color_width;
				int		v = (y / z) * inv_color_xyScale + half_color_height;
				//	write to mapped depth map
				if( u>=0 && u<depth_width && v>=0 && v<depth_height ){
					int new_idx = v * color_width + u;
					if( !mapped_depth[new_idx] || value>>3 < mapped_depth[new_idx]>>3 )
						mapped_depth[new_idx] = value;
				}
			}	
		}

		//	over write old data
		memcpy(depth_image, mapped_depth, depth_width * depth_height * sizeof(unsigned short));

		delete mapped_depth;
	}

	/**
	*/
	int MapSkeletonToColor(){
	}

	int GetInitializeFlag(){
		int flag = 0;
		//	color image
		if(color_flag)			flag |= NUI_INITIALIZE_FLAG_USES_COLOR;
		if(hd_color_flag)		flag |= NUI_INITIALIZE_FLAG_USES_HIGH_QUALITY_COLOR;
		//	depth
		if(depth_flag)			flag |= NUI_INITIALIZE_FLAG_USES_DEPTH;
		if(player_index_flag)	flag |= NUI_INITIALIZE_FLAG_USES_DEPTH_AND_PLAYER_INDEX;
		//	audio and skeleton
		if(audio_flag)			flag |= NUI_INITIALIZE_FLAG_USES_AUDIO;
		if(skeleton_flag)		flag |= NUI_INITIALIZE_FLAG_USES_SKELETON;

		return flag;
	}
	inline int IsFailed(HRESULT hr, const char* str){
		if (FAILED(hr) ) { 
			#ifndef BE_QUIET
				std::cout << "error: MicrosoftKinect, " 
					<< str << ", error code: " << std::hex << hr << std::endl;
			#endif
			return 1; 
		}
		return 0;
	}
};




}}	//	namespace yz::kinect





#endif	//	__YZ_MICROSOFT_KINECT_H__