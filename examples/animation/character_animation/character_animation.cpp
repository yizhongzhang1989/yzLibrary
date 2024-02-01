/***********************************************************/
/**	\file
	\brief		Example of Character animation
	\author		Yizhong Zhang
	\date		7/1/2012
*/
/***********************************************************/

//	This example shows how to use Character Animation

#ifdef ENABLE_OPENMP
#	include <omp.h>
#endif
#include <iostream>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <mkl_types.h>
//#include <mkl_dss.h>
#include <yzLib/yz_lib.h>
using namespace std;

#include "data_path.h"

yz::opengl::DemoWindowManager		manager;
yz::opengl::DemoWindow3D			win3d;

typedef double REAL;

float		fps;
yz::animation::SkeletonAnimation<REAL>		James_skeleton;
yz::geometry::SmoothShadingTriMesh<REAL>	James_skin;
yz::animation::RiggingLBS<REAL>				James_rigging;
int frame_id=0;
int moving_flag;
int step_flag;
int draw_skeleton_flag;

//	over-ridden functions
void calculateFps(){
	static yz::utils::FPSCalculator fps_calc;
	fps = fps_calc.GetFPS();
}

void win3d_showInfo(){
	glColor3f(1, 0, 1);
	const char* cmd1 = "press 's' to start";
	const char* cmd2 = "press 's' to stop";
	const char* cmd3 = "press 'r' to back to rest pose";
	yz::opengl::printInfo(0, 0, "fps: %f\n%s\n%s", fps, (moving_flag?cmd2:cmd1), cmd3);
}

void win3d_draw(){
	yz::opengl::drawXYZAxis();

	if( draw_skeleton_flag )
		James_skeleton.Display(1);

	glColor3f(1, 1, 1);
	James_skin.Display();
}

void win3d_keyboard(unsigned char key, int x, int y){
	switch( key ){
		case 27:
			exit(0);
		case 's':
			moving_flag = !moving_flag;
			break;
		case 'r':
			moving_flag = 0;
			James_skeleton.TransformSkeletonToRest();
			break;
		case 'd':
			draw_skeleton_flag = !draw_skeleton_flag;
			break;
		case ' ':
			moving_flag = 0;
			step_flag = 1;
			break;
	}
}

void idle(){
	if( moving_flag || step_flag ){
		James_skeleton.TransformSkeletonToRest();
		yz::animation::transformSkeletonAndMeshLBS(James_skeleton.skeleton,
			James_skin.vertex, James_skeleton.transformer_data.frame[frame_id], James_rigging);

		James_skin.CalculateVertexNormal();

		frame_id ++;
		if( frame_id >= James_skeleton.transformer_data.FrameNumber() )
			frame_id = 0;

		if( step_flag )
			step_flag = 0;
	}
}

int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	char bvh_filename[1024], obj_filename[1024];
	sprintf(bvh_filename, "%s/James_motion.bvh", data_path);
	sprintf(obj_filename, "%s/James_simplify.obj", data_path);

	if( !James_skeleton.ReadSkeletonAndFramesFromBVH(bvh_filename, 0.025) ){
		cout << "read James skeleton failed" << endl;
		return 0;
	}
	if( !James_skin.ReadMeshFromFile(obj_filename) ){
		cout << "read James skin failed" << endl;
		return 0;
	}
	yz::geometry::removeIsolatedVerticesFromMesh(James_skin.vertex, James_skin.face);
	James_skin.CalculateVertexNormal();

	James_skeleton.TransformSkeletonToRest();
	yz::animation::riggingBoneHeat(James_rigging, James_skeleton.skeleton, James_skeleton.transformer_data.frame[0],
		James_skin.vertex, James_skin.face, 0.001, 0.001);
	James_skeleton.ResetFrameInterval(0.1);

	win3d.SetDraw(win3d_draw);
	win3d.SetDrawAppend(win3d_showInfo);
	win3d.keyboardFunc = win3d_keyboard;
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(calculateFps);	//	add each idle function in sequence
	manager.AddIdleFunc(idle);
	manager.AddIdleFunc(win3d.idleFunc);

	manager.EnterMainLoop();			//	then enter main loop
}