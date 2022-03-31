/***********************************************************/
/**	\file
	\brief		Ascii displayer
	\author		Yizhong Zhang
	\date		4/24/2013
*/
/***********************************************************/
#ifndef __YZ_ASCII_DISPLAYER_H__
#define __YZ_ASCII_DISPLAYER_H__

#include "yzLib/yz_setting.h"
#include <iostream>

#include "yzLib/yz_math/yz_vector.h"
#include "yzlib/yz_utils/yz_ascii_table.h"

namespace yz{	namespace opengl{


/**
	Display ascii characters with quads
*/
class AsciiDisplayer{
public:
	/**
		create displayer content from ascii table pixel mask

		\param	ascii		pixel mask of each character
		\param	width		pixels in width of each character, should be 8 * x
		\param	height		pixels in height of each character
		\param	length		length of the ascii table, typically 128 or 256
	*/
	void Create(const unsigned char* ascii, int width, int height, int length){
		//	create vertex
		vertex.clear();
		vertex.reserve( (width+1)*(height+1) );
		//float x_scale_coef = float(width)/height;
		float x_scale_coef = 1.0f;
		for(int j=0; j<=height; j++)
			for(int i=0; i<=width; i++)
				vertex.push_back( yz::Vec2f(float(i)/width * x_scale_coef, 1.0f - float(j)/height) );

		//	create face
		quad.clear();
		quad_start.clear();
		quad_start.resize(length+1);
		int bytes_per_line = (width+7) / 8;
		int bytes_per_char = bytes_per_line * height;
		for(int c=0; c<length; c++){
			quad_start[c] = quad.size();
			for(int j=0; j<height; j++)
				for(int i=0; i<width; i++){
					int byte_id = i / 8;
					unsigned char line_seg = ascii[c*bytes_per_char + j*bytes_per_line + byte_id];
					line_seg <<= (i - byte_id * 8);
					if( line_seg & 0x80 ){		//	this pixel should be displayed
						int base_id = j * (width + 1) + i;
						quad.push_back( yz::int4(base_id, base_id+width+1, base_id+width+2, base_id+1) );
					}
				}
		}
		quad_start[length] = quad.size();
	}

	/**
		setup 8x16 ascii table
	*/
	void Setup8x16(){
		Create(utils::ascii8x16, 8, 16, 256);
	}

	/**
		draw the given character in xy plane as a (0,0)-(1,1) quad
	*/
	void Draw(unsigned char ascii){
		if( ascii > quad_start.size()-1 )
			return;

		glNormal3f(0.0f, 0.0f, 1.0f);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2, GL_FLOAT, 0, (float*)&vertex[0]);

		glDrawElements( GL_QUADS, (quad_start[ascii+1] - quad_start[ascii])*4, GL_UNSIGNED_INT, (int*)&quad[quad_start[ascii]] );

		glDisableClientState(GL_VERTEX_ARRAY);
	}

	/**
		draw the given character in xy plane as given quad
	*/
	void Draw(unsigned char ascii, float x_min, float y_min, float x_max, float y_max){
		if( ascii > quad_start.size()-1 )
			return;

		glPushMatrix();
		glTranslatef(x_min, y_min, 0.0f);
		glScalef(x_max-x_min, y_max-y_min, 1.0f);

		Draw(ascii);

		glPopMatrix();
	}

public:
	std::vector<yz::Vec2f>	vertex;			//	normalize each character into (0,0) - (1,1)
	std::vector<yz::int4>	quad;			//	quad of each pixel
	std::vector<int>		quad_start;		//	start quad of each character

};


}}	//	end namespace yz::opengl


#endif	//	__YZ_ASCII_DISPLAYER_H__