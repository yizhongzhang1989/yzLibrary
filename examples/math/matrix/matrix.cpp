/***********************************************************/
/**	
\file
\brief		Example of Matrix
\details	This file illustrates how to use matrix in yzLib
\author		Yizhong Zhang
\date		3/15/2020
*/
/***********************************************************/
#include <iostream>
#include <yzLib/yz_math.h>

int main() {
	{
		std::cout << "==============================" << std::endl;
		std::cout << "Test Matrix2x2:" << std::endl;
		yz::Matrix2x2<double>	M1, M2;
		M1.SetIdentity();
		M2.SetRotationDeg(90);
		std::cout << "M1: " << std::endl << M1;
		std::cout << "M2: " << std::endl << M2;
		std::cout << "M1 * M2: " << std::endl << M1 * M2;
		std::cout << "--------------------" << std::endl;
	}

	//std::cout << "==============================" << std::endl;
	//std::cout << "Test Matrix3x3:" << std::endl;
	//yz::Matrix3x3d	M3x3;
	//M3x3.SetRotationDeg(yz::Vec3d(1, 1, 0), 90);
	//std::cout << "Rotate around (1,1,0) 90deg: " << std::endl << M3x3 << std::endl;

	//std::cout << "==============================" << std::endl;
	//std::cout << "Test Matrix4x4:" << std::endl;
	//yz::Matrix4x4d	M4x4;
	//M4x4.SetRotationDeg(yz::Vec3d(1, 1, 0), 90);
	//std::cout << "Rotate around (1,1,0) 90deg: " << std::endl << M4x4 << std::endl;

	return 0;
}