#include <iostream>
// #include <GL/glew.h>
// #include <GL/glut.h>
#include <yzLib/yz_lib.h>


yz::opengl::DemoWindowManager	manager;
yz::opengl::GLUTWindow3D<0>		win3d;

yz::physics::MassSpring<double> mass_spring;
int simulate_flag = 0;


void draw() {
	yz::opengl::drawXYZAxis();
	mass_spring.Draw();
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	case ' ':
		simulate_flag = !simulate_flag;
		break;
	case 's':
		mass_spring.SimulateStep(0.01);
		break;
	case 'q':
		mass_spring.SimulateStep(-1);
		break;
	}
}

void idle() {
	if (simulate_flag) {
		//simulate_flag = 0;
		//mass_spring.SimulateStepExplicit(0.01);
		mass_spring.SimulateStep(1);
	}
}

int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	std::vector<yz::Vec3d> node;
	node.push_back(yz::Vec3d(0, 0, 0));
	node.push_back(yz::Vec3d(1, 0, 0));
	node.push_back(yz::Vec3d(0, 1, 0));
	node.push_back(yz::Vec3d(0, 0, 1));
	std::vector<yz::int2> spring;
	spring.push_back(yz::int2(0, 1));
	spring.push_back(yz::int2(0, 2));
	spring.push_back(yz::int2(0, 3));
	spring.push_back(yz::int2(1, 2));
	spring.push_back(yz::int2(1, 3));
	spring.push_back(yz::int2(2, 3));

	mass_spring.SetupMassSpring(node, spring);
	mass_spring.InitPhysicaParameters(1, 10, 1, 0);
	//mass_spring.AddAttachConstraint(0);
	mass_spring.rest_length[0] *= 10;
	//mass_spring.v[0] = yz::Vec3d(10, 0, 0);
	//mass_spring.v[1] = yz::Vec3d(-10, 0, 0);

	//std::vector<yz::Vec3d> df;
	//std::vector<yz::Vec3d> dx;
	//dx.push_back(yz::Vec3d(0.01, 0, 0));
	//dx.push_back(yz::Vec3d(0, 0.01, 0));
	//dx.push_back(yz::Vec3d(0, 0, 0.01));
	//dx.push_back(yz::Vec3d(0, 0, 0));
	//mass_spring.ComputeElasticDifferentialForce(df, dx);
	//for (int i = 0; i < df.size(); i++)
	//	std::cout << df[i] << std::endl;

	//mass_spring.ClearForce();
	//mass_spring.AddElasticForce();
	//df = mass_spring.f;
	//for (int i = 0; i < dx.size(); i++)
	//	mass_spring.x[i] += dx[i];
	//mass_spring.ClearForce();
	//mass_spring.AddElasticForce();
	//for (int i = 0; i < df.size(); i++)
	//	df[i] = mass_spring.f[i] - df[i];
	//for (int i = 0; i < df.size(); i++)
	//	std::cout << df[i] << std::endl;


	win3d.keyboardFunc = keyboard;
	win3d.SetDraw(draw);
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(win3d.idleFunc);
	manager.AddIdleFunc(idle);
	manager.EnterMainLoop();
}
