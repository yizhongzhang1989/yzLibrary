/***********************************************************/
/**	\file
	\brief		OpenGL Utilities
	\details	util functions to make using OpenGL easier. 
				must include after <glut.h>.
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_OPENGL_UTILS_H__
#define __YZ_OPENGL_UTILS_H__

#pragma  warning(disable:4996)

#include "yzLib/yz_setting.h"

#if !(defined(YZ_glut_h) || defined(YZ_freeglut_h) )
#	error yz_opengl_utils.h must be included after glut.h or freeglut.h
#endif

#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <sstream>
#include <math.h>
#include <stdarg.h>
#include "yzLib/yz_math/yz_lookup_table.h"
#include "yzLib/yz_math/yz_numerical_utils.h"


namespace yz{	namespace opengl{
//	========================================
///@{
/**	@name OpenGL Check and Utils
*/
//	========================================
/**
	If gl error exist, print the error on console

	\return		whether error exist
*/
inline int printGLError(const char* info = NULL){
	if( GLenum err = glGetError() ){
		std::cout << info << ", " << "glError " << err << " : " << gluErrorString(err) << std::endl;
		return 1;
	}
	return 0;
}

/**
	Check whether the matrix stack is full

	\param	matrix		the matrix to check, GL_MODELVIEW / GL_PROJECTION / GL_TEXTURE / GL_COLOR(glew)
	\return		1: is full, 0: not full, -1: invalid input
*/
inline int isMatrixStackFull(GLint matrix = GL_MODELVIEW){
	int curr_dep, max_dep;
	switch(matrix){
		case GL_MODELVIEW:
			glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &curr_dep);
			glGetIntegerv(GL_MAX_MODELVIEW_STACK_DEPTH, &max_dep);
			return curr_dep == max_dep;
		case GL_PROJECTION:
			glGetIntegerv(GL_PROJECTION_STACK_DEPTH, &curr_dep);
			glGetIntegerv(GL_MAX_PROJECTION_STACK_DEPTH, &max_dep);
			return curr_dep == max_dep;
		case GL_TEXTURE:
			glGetIntegerv(GL_TEXTURE_STACK_DEPTH, &curr_dep);
			glGetIntegerv(GL_MAX_TEXTURE_STACK_DEPTH, &max_dep);
			return curr_dep == max_dep;
#ifdef YZ_glew_h
		case GL_COLOR:
			glGetIntegerv(GL_COLOR_MATRIX_STACK_DEPTH, &curr_dep);
			glGetIntegerv(GL_MAX_COLOR_MATRIX_STACK_DEPTH, &max_dep);
			return curr_dep == max_dep;
#endif
		default:
			return -1;
	}
}

/**
	Check whether the attribute stack is full

	\return		1: is full, 0: not full
*/
inline int isAttributeStackFull(){
	int curr_dep, max_dep;
	glGetIntegerv(GL_ATTRIB_STACK_DEPTH, &curr_dep);
	glGetIntegerv(GL_MAX_ATTRIB_STACK_DEPTH, &max_dep);
	return curr_dep == max_dep;
}

/**
	Check whether the client attribute stack is full

	\return		1: is full, 0: not full
*/
inline int isClientAttributeStackFull(){
	int curr_dep, max_dep;
	glGetIntegerv(GL_CLIENT_ATTRIB_STACK_DEPTH, &curr_dep);
	glGetIntegerv(GL_MAX_CLIENT_ATTRIB_STACK_DEPTH, &max_dep);
	return curr_dep == max_dep;
}

/**
	push all attributes and matrices

	must use in pairs with popAllAttributesAndMatrices

	in Debug mode, stack depth is checked to avoid overflow
*/
inline void pushAllAttributesAndMatrices(){
	assert( !isMatrixStackFull(GL_PROJECTION) && !isMatrixStackFull(GL_MODELVIEW) && !isMatrixStackFull(GL_TEXTURE) );
#ifdef YZ_glew_h
	assert( !isMatrixStackFull(GL_COLOR) );
#endif
	assert( !isAttributeStackFull() && !isClientAttributeStackFull() );

	//	push all attributes
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

	//	push all matrices
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();		//	push projection matrix
	glMatrixMode(GL_TEXTURE);
	glPushMatrix();		//	push texture matrix
	glMatrixMode(GL_COLOR);
	glPushMatrix();		//	push color matrix
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();		//	push modelview matrix
}

/**
	pop all attributes and matrices

	must use in pairs with pushAllAttributesAndMatrices
*/
inline void popAllAttributesAndMatrices(){
	//	pop all matrices
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();		//	pop modelview matrix
	glMatrixMode(GL_COLOR);
	glPopMatrix();		//	pop color matrix
	glMatrixMode(GL_TEXTURE);
	glPopMatrix();		//	pop texture matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();		//	pop projection matrix

	//	pop all attributes
	glPopClientAttrib();
	glPopAttrib();
}

/**
	get projection matrix of current OpenGL context
*/
inline Matrix4x4d getProjectionMatrixRowMajor(){
	Matrix4x4d m;
	glGetDoublev(GL_PROJECTION_MATRIX, m.data[0]);
	m.SetTranspose();
	return m;
}

/**
	get modelview matrix of current OpenGL context
*/
inline Matrix4x4d getModelviewMatrixRowMajor(){
	Matrix4x4d m;
	glGetDoublev(GL_MODELVIEW_MATRIX, m.data[0]);
	m.SetTranspose();
	return m;
}

/**
	A helper function of glDrawElements
*/
template<typename TV, typename TN, typename TC, typename TI>
void DrawElements(
	GLenum			elem_mode, 
	const TV*		vertex,
	const TN*		normal,
	const TC*		color,
	const TI*		indices,
	unsigned int	indices_number,
	unsigned int	color_size = 3,
	unsigned int	vertex_size = 3
	)
{
	if (!indices_number)
		return;

	if (vertex) {
		glEnableClientState(GL_VERTEX_ARRAY);
		if (Is_float<TV>::check_type)
			glVertexPointer(vertex_size, GL_FLOAT, 0, vertex);
		else if (Is_double<TV>::check_type)
			glVertexPointer(vertex_size, GL_DOUBLE, 0, vertex);
		else
			std::cout << "error: DrawElements, invalid vertex data type" << std::endl;
	}
	if (normal) {
		glEnableClientState(GL_NORMAL_ARRAY);
		if (Is_float<TN>::check_type)
			glNormalPointer(GL_FLOAT, 0, normal);
		else if (Is_double<TN>::check_type)
			glNormalPointer(GL_DOUBLE, 0, normal);
		else
			std::cout << "error: DrawElements, invalid normal data type" << std::endl;
	}
	if (color) {
		glEnableClientState(GL_COLOR_ARRAY);
		if (Is_float<TC>::check_type)
			glColorPointer(color_size, GL_FLOAT, 0, color);
		else if (Is_unsigned_char<TC>::check_type || Is_char<TC>::check_type)
			glColorPointer(color_size, GL_UNSIGNED_BYTE, 0, color);
		else
			std::cout << "error: DrawElements, invalid color data type" << std::endl;
	}

	glDrawElements(elem_mode, indices_number, GL_UNSIGNED_INT, indices);

	if (vertex)
		glDisableClientState(GL_VERTEX_ARRAY);
	if (normal)
		glDisableClientState(GL_NORMAL_ARRAY);
	if (color)
		glDisableClientState(GL_COLOR_ARRAY);
}

///@}

//	========================================
///@{
/**	@name Set Random Color
*/
//	========================================
/**
	set seed for rand() of stdlib
*/
inline void srandColor(unsigned int seed){
	srand(seed);
}
/**
	get color randomly
*/
inline void randColor(float& r, float& g, float& b){
	r = randNumber(0.3f, 1);
	g = randNumber(0.3f, 1);
	b = randNumber(0.3f, 1);
}

/**
	get color randomly
*/
inline void randColor(float rgb[3]){
	randColor(rgb[0], rgb[1], rgb[2]);
}

/**
	Get Sequential Color

	given a number, return a uniqie color to represent this number

	The color representing the number is specifically designed,
	it is unique for each number, so it can be used for color picking
*/
inline void getSequentialColor(unsigned char rgba[4], int index){
	unsigned char rgba_idx[3] = {0, 0, 0};
	rgba[3] = (~index & 0xFF000000) >> 24;
	for( int i=0; i<8; i++ ){
		rgba_idx[0] |= (index & 0x01<<(i*3)  ) >> (i*2);
		rgba_idx[1] |= (index & 0x01<<(i*3+1)) >> (i*2+1);
		rgba_idx[2] |= (index & 0x01<<(i*3+2)) >> (i*2+2);
	}
	for( int i=0; i<3; i++ )
		rgba[i] = ShuffledSequenceTableUC[ rgba_idx[i] ];
}

/**
	Get Index of Sequential Color

	For each RGBA color, get the unique index
*/
inline int getIndexOfSequentialColor(unsigned char rgba[4]){
	unsigned char rgba_idx[3] = {
		InverseShuffledSequenceTableUC[rgba[0]],
		InverseShuffledSequenceTableUC[rgba[1]],
		InverseShuffledSequenceTableUC[rgba[2]] };

	int index = int(~rgba[3]) << 24;
	for( int i=0; i<8; i++ ){
		index |= ((rgba_idx[0] & 0x01<<i) << i*2) | 
			((rgba_idx[1] & 0x01<<i) << (i*2+1)) | 
			((rgba_idx[2] & 0x01<<i) << (i*2+2));
	}

	return index;
}

/**
	Get 24 bits Sequential Color, alpha is not included

	given a number, return a uniqie color to represent this number

	The color representing the number is specifically designed,
	it is unique for each number, so it can be used for color picking
*/
inline void getSequentialColor24Bits(unsigned char rgba[3], int index){
	unsigned char rgba_idx[3] = {0, 0, 0};
	for( int i=0; i<8; i++ ){
		rgba_idx[0] |= (index & 0x01<<(i*3)  ) >> (i*2);
		rgba_idx[1] |= (index & 0x01<<(i*3+1)) >> (i*2+1);
		rgba_idx[2] |= (index & 0x01<<(i*3+2)) >> (i*2+2);
	}
	for( int i=0; i<3; i++ )
		rgba[i] = ShuffledSequenceTableUC[ rgba_idx[i] ];
}

/**
	Get Index of 24 bits Sequential Color, alpha is not included

	For each RGBA color, get the unique index
*/
inline int getIndexOfSequentialColor24Bits(unsigned char rgba[3]){
	unsigned char rgba_idx[3] = {
		InverseShuffledSequenceTableUC[rgba[0]],
		InverseShuffledSequenceTableUC[rgba[1]],
		InverseShuffledSequenceTableUC[rgba[2]] };

	int index = 0;
	for( int i=0; i<8; i++ ){
		index |= ((rgba_idx[0] & 0x01<<i) << i*2) | 
			((rgba_idx[1] & 0x01<<i) << (i*2+1)) | 
			((rgba_idx[2] & 0x01<<i) << (i*2+2));
	}

	return index;
}

/**
	Get Sequential Color for display

	the color table contain 256 colors, so if index is bigger than 255, 
	it will start from beginning
*/
inline void getSequentialDisplayColor(unsigned char rgb[3], int index){
	index = index & 0xFF;	//	loop every 256
	int color = SequentialDisplayColor[index];
	rgb[0] = (color & 0xFF0000) >> 16;
	rgb[1] = (color & 0x00FF00) >> 8;
	rgb[2] = color & 0x0000FF;
}

/**
	Get Sequential Color for display

	the color table contain 256 colors, so if index is bigger than 255, 
	it will start from beginning
*/
inline void getSequentialDisplayColor(float rgb[3], int index){
	unsigned char ucrgb[3];
	getSequentialDisplayColor(ucrgb, index);
	rgb[0] = ucrgb[0] / 255.0f;
	rgb[1] = ucrgb[1] / 255.0f;
	rgb[2] = ucrgb[2] / 255.0f;
}

/**
	Set Sequential Color for display
*/
inline void setSequentialDisplayColor(int index){
	unsigned char rgb[3];
	getSequentialDisplayColor(rgb, index);
	glColor3ubv(rgb);
}

///@}

//	========================================
///@{
/**	@name Picking Related Functions
*/
//	========================================
/**
	set picking index by set color

	set to sequential color by index

	since we use -1 for background color, so -1 cannot be used as index
*/
inline void setPickingIndex(int index){
	unsigned char rgba[4];
	getSequentialColor(rgba, index);
	glColor4ubv(rgba);
}

/**
	get picking index from rgba color
*/
inline int getPickingIndex(unsigned char rgba[4]){
	return getIndexOfSequentialColor(rgba);
}

///@}

//	========================================
///@{
/**	@name Calculate Coordinate
*/
//	========================================

/**
	given 3d coordinate, calculate window coordinate

	\param	pixel_x		x on screen
	\param	pixel_y		y on screen
	\param	x			x in space
	\param	y			y in space
	\param	z			z in space
	\return				GL Error
*/
inline GLuint getWindowCoordinate(int &pixel_x, int &pixel_y, GLdouble x, GLdouble y, GLdouble z){
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLdouble win_x, win_y, win_z;

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	GLuint ret = gluProject( x, y, z, modelview, projection, viewport, &win_x, &win_y, &win_z);

	pixel_x = win_x;
	pixel_y = win_y;
	return ret;
}

/**
	given 3d coordinate, calculate window coordinate and depth

	\param	pixel_x		x on screen
	\param	pixel_y		y on screen
	\param	depth		depth of this point
	\param	x			x in space
	\param	y			y in space
	\param	z			z in space
	\return				GL Error
*/
inline GLuint getWindowCoordinate(double &pixel_x, double &pixel_y, double &depth, GLdouble x, GLdouble y, GLdouble z){
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLdouble win_x, win_y, win_z;

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	GLuint ret = gluProject( x, y, z, modelview, projection, viewport, &win_x, &win_y, &win_z);

	pixel_x = win_x;
	pixel_y = win_y;
	depth	= win_z;
	return ret;
}

/**
	given window coordinate, calculate world coordinate

	\param	x			x in space
	\param	y			y in space
	\param	z			z in space
	\param	pixel_x		x on screen
	\param	pixel_y		y on screen
	\param	depth		depth value of (x, y)
	\return				GL Error
*/
template <typename T>
inline GLuint getWorldCoordinate(T& x, T& y, T& z, GLdouble pixel_x, GLdouble pixel_y, GLdouble depth){
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLdouble world_x, world_y, world_z;

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	GLuint ret = gluUnProject(pixel_x, viewport[3]-pixel_y, depth, modelview, projection, viewport, &world_x, &world_y, &world_z);
	x = world_x;
	y = world_y;
	z = world_z;

	return ret;
}

/**
	transform window coordinate from old to new

	consider two windows, one with top-left being the origin,
	the other with bottom-right being the origin. We want to
	get the corresponding coordinate of the other window, then
	we can call this function.

	\param	old_width		width of old window
	\param	old_height		height of old window
	\param	old_coor_x		coordinate x in old window
	\param	old_coor_y		coordinate y in old window
	\param	new_width		width of new window
	\param	new_height		height of new window
	\param	new_coor_x		return coordinate in new window
	\param	new_coor_y		return coordinate in new window
	\param	vertical_flip	whether flip the transform vertically
	\param	horizontal_flip	whether flip the transform horizontally
*/
inline void getTransformedWindowCoordinate(int& new_coor_x, int& new_coor_y, int new_width, int new_height,
										   int old_coor_x, int old_coor_y, int old_width, int old_height, 
										   int vertical_flip = 0, int horizontal_flip = 0){
	float x = (old_coor_x + 0.5f) / old_width;
	float y = (old_coor_y + 0.5f) / old_height;
	if( vertical_flip )		y = 1 - y;
	if( horizontal_flip )	x = 1 - x;

	new_coor_x = floor( x*new_width );
	new_coor_y = floor( y*new_height );
}

/**
	transform texture coordinate from old to new

	different from window coordinate, (0, 0) is on the bondary
	of the texture, while (0, 0) is in the middle of the first
	pixel. 

	\param	old_width		width of old window
	\param	old_height		height of old window
	\param	old_coor_x		coordinate x in old window
	\param	old_coor_y		coordinate y in old window
	\param	new_width		width of new window
	\param	new_height		height of new window
	\param	new_coor_x		return coordinate in new window
	\param	new_coor_y		return coordinate in new window
	\param	vertical_flip	whether flip the transform vertically
	\param	horizontal_flip	whether flip the transform horizontally
*/
inline void getTransformedTextureCoordinate(int& new_coor_x, int& new_coor_y, int new_width, int new_height,
										   int old_coor_x, int old_coor_y, int old_width, int old_height, 
										   int vertical_flip = 0, int horizontal_flip = 0){
	float x = float(old_coor_x) / (old_width-1);
	float y = float(old_coor_y) / (old_height-1);
	if( vertical_flip )		y = 1 - y;
	if( horizontal_flip )	x = 1 - x;

	new_coor_x = roundToClosestInteger( x*(new_width-1) );
	new_coor_y = roundToClosestInteger( y*(new_height-1) );
}

///@}

//	========================================
///@{
/**	@name Show String On Screen
*/
//	========================================

/**
	Font size used in printString

	Three fonts are provided to select, small, medium, large
*/
enum{
	PRINT_STRING_FONT_SIZE_SMALL,
	PRINT_STRING_FONT_SIZE_MEDIUM,
	PRINT_STRING_FONT_SIZE_LARGE
};

/**
	Print string on window, left up corner (0, 0)

	called: PushAttrib once, PushMatrix(PROJECTION) once, PushMatrix(MODEL_VIEW) once

	We use window size to be coordinate size. So if you are printing the string
	to FBO that is not the same size as window, the position of string maynot be 
	the same as wanted

	\param	x			start point of x on window
	\param	y			start point of y on window
	\param	width		how many characters is it able to hold each line
	\param	font_size	font size used to display, should be one of 
						PRINT_STRING_FONT_SIZE_[SMALL][MEDIUM][LARGE]
	\param	vertical_inc	incremental of y to display each line
	\param	format		string format
*/
#ifndef _WINDOWS
	inline int _vscprintf (const char * format, va_list pargs) { 
      int retval; 
      va_list argcopy; 
      va_copy(argcopy, pargs); 
      retval = vsnprintf(NULL, 0, format, argcopy); 
      va_end(argcopy); 
      return retval; 
	}
#endif

inline void printString(float x, float y, int width, int font_size, int vertical_inc, const char* format, ...){	
	//	given window coordinate, print information on the screen just like normal printf
	void* font;
	if		(font_size == PRINT_STRING_FONT_SIZE_SMALL)		font = GLUT_BITMAP_TIMES_ROMAN_10;
	else if	(font_size == PRINT_STRING_FONT_SIZE_MEDIUM)	font = GLUT_BITMAP_9_BY_15;
	else if (font_size == PRINT_STRING_FONT_SIZE_LARGE)		font = GLUT_BITMAP_TIMES_ROMAN_24;
	else													font = GLUT_BITMAP_8_BY_13;

	//	setup planar show
	assert( !isMatrixStackFull(GL_PROJECTION) && !isMatrixStackFull(GL_MODELVIEW) );
	assert( !isAttributeStackFull() );

	glPushAttrib(GL_TRANSFORM_BIT | GL_LIGHTING_BIT | GL_DEPTH_BITS);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();		//	be careful pushing projection matrix, only 2 entries. 
	glLoadIdentity();
	gluOrtho2D(0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), 0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	//	get string
	va_list args;
	va_start( args, format );
	int len = _vscprintf( format, args ) + 1;
	char* buffer = new char[len];
	vsprintf( buffer, format, args );
	va_end( args );

	//	show
	y += vertical_inc;
	glRasterPos2i(x, y);
	int lines = 0, line_c_count = 0;
	for( int i=0; i<len; i++ ){
		char c = buffer[i];
		if( c == '\n' ){
			glRasterPos2i(x, y+vertical_inc*(++lines));
			line_c_count = 0;
		}
		else{
			glutBitmapCharacter( font, c );
			line_c_count ++;
			if( line_c_count == width ){
				glRasterPos2i(x, y+vertical_inc*(++lines));
				line_c_count = 0;
			}
		}
	}

	delete[] buffer;

	//	roll back OpenGL
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();
}

/**
	Print string on window with large font, left up corner (0, 0)

	We use window size to be coordinate size. So if you are printing the string
	to FBO that is not the same size as window, the position of string maynot be 
	the same as wanted

	\param	x		x on window
	\param	y		y on window
	\param	format	string format
*/
inline void printInfo(float x, float y, const char* format, ...){	
	//	get string
	va_list args;
	va_start( args, format );
	int len = _vscprintf( format, args ) + 1;
	char* buffer = new char[len];
	vsprintf( buffer, format, args );
	va_end( args );

	//	print string
	printString(x, y, 1000, PRINT_STRING_FONT_SIZE_LARGE, 30, buffer);

	delete[] buffer;
}

/**
	Print string list on the screen, simulating the output of console

	\param	str_list	list of string
	\param	x			start point of x on window
	\param	y			start point of y on window
	\param	font_size	font size used to display, should be one of 
						PRINT_STRING_FONT_SIZE_[SMALL][MEDIUM][LARGE]
	\param	width		how many characters each line can hold
	\param	height		how many lines can the screen hold
*/
inline void printStringList(const std::list<std::string>& str_list, 
							int x=0, 
							int y=0, 
							int font_size = PRINT_STRING_FONT_SIZE_MEDIUM, 
							int width=80, 
							int height=20 ){
	//	create a combined string, insert '\n' between each string
	std::string display_str;
	for(std::list<std::string>::const_iterator iter = str_list.begin(); iter!=str_list.end(); iter++ ){
		if( !display_str.empty() )	display_str += "\n";
		display_str += *iter;
	}

	//	calculate start position of each line
	std::vector<int> start_idx;
	int line_c_count = 0;
	for( int i=0; i<display_str.size(); i++ ){
		if( line_c_count == 0 )	// first character of a line
			start_idx.push_back(i);

		if( display_str[i] == '\n' ){
			line_c_count = 0;
			continue;
		}
		else{
			line_c_count ++;
		}
		if( line_c_count == width )
			line_c_count = 0;
	}

	//	the list has exceeded height limit
	if( start_idx.size() > height ){
		display_str.erase(0, start_idx[start_idx.size()-height]);
	}

	//	print the string
	int vertical_inc[4] = {14, 19, 29, 17};
	if( font_size < 0 )	font_size = 0;
	if( font_size > 3 ) font_size = 3;
	printString(x, y, width, font_size, vertical_inc[font_size], display_str.c_str());
}

/**
	Draw string on window by space coordinate

	called: PushAttrib once

	\param	str		string to display
	\param	x		x in space
	\param	y		y in space
	\param	z		z in space
*/
inline void drawString(std::string str, float x=0, float y=0, float z=0){
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);

#ifdef YZ_glew_h	//	able to start new line with glew
	int win_x, win_y;
	getWindowCoordinate(win_x, win_y, x, y, z);
	glWindowPos2i( win_x, win_y );
	int lines = 0;
	for( unsigned int i=0; i<str.length(); i++ ){
		char c = str[i];
		if( c == '\n' )
			glWindowPos2i( win_x, win_y-12*(++lines) );
		else
			glutBitmapCharacter( GLUT_BITMAP_9_BY_15, c );
	}
#else	//	unable to start a new line, we just print '/' character
	glRasterPos3f(x, y, z);
	for( unsigned int i=0; i<str.length(); i++ ){
		char c = str[i];
		glutBitmapCharacter( GLUT_BITMAP_9_BY_15, (c=='\n'? '/' : c) );
	}
#endif

	glPopAttrib();
}

/**
	Draw a number on window by space coordinate.
	"cout << " should be redefined on that number

	\param	value	the number to display
	\param	x		x in space
	\param	y		y in space
	\param	z		z in space
*/
template<typename T>
inline void drawNumber(T value, float x=0, float y=0, float z=0){
	glRasterPos3f(x, y, z);
	std::ostringstream tmp_str;		//	print value to string
	tmp_str << value;
	drawString(tmp_str.str(), x, y, z);
}
///@}

//	========================================
///@{
/**	@name Draw Texture
*/
//	========================================

/**
	Draw the whole texture (0,0)-(1,1) on 2D canvas

	We assume that the texture is already bind

	\param	x_min		min x on canvas
	\param	y_min		min y on canvas
	\param	x_max		max x on canvas
	\param	y_max		max y on canvas
	\param	flip_flag	whether vertical flip the texture, 1: flip, 0: don't flip
*/
inline void drawWholeTexture(float x_min, float y_min, float x_max, float y_max, int flip_flag=0){
	if( flip_flag ){
		glBegin( GL_QUADS );
			glTexCoord2f(0, 0);	glVertex2f(x_min, y_max);
			glTexCoord2f(1, 0);	glVertex2f(x_max, y_max);
			glTexCoord2f(1, 1);	glVertex2f(x_max, y_min);
			glTexCoord2f(0, 1);	glVertex2f(x_min, y_min);
		glEnd();
	}
	else{
		glBegin( GL_QUADS );
			glTexCoord2f(0, 0);	glVertex2f(x_min, y_min);
			glTexCoord2f(1, 0);	glVertex2f(x_max, y_min);
			glTexCoord2f(1, 1);	glVertex2f(x_max, y_max);
			glTexCoord2f(0, 1);	glVertex2f(x_min, y_max);
		glEnd();
	}
}

///@}

//	========================================
///@{
/**	@name Draw Certain Shape
*/
//	========================================

/**
	Draw XYZ axis of coordinate system from original point 
	using red, green, blue color for each axis

	called: PushAttrib once

	\param	length		the length of the axis to draw
*/
inline void drawXYZAxis(float length = 1){
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable( GL_LIGHTING );
	glBegin( GL_LINES );
		//	x
		glColor3f( 1, 0, 0 );
		glVertex3f( 0, 0, 0 );
		glVertex3f( length, 0, 0 );
		//	y
		glColor3f( 0, 1, 0 );
		glVertex3f( 0, 0, 0 );
		glVertex3f( 0, length, 0 );
		//	z
		glColor3f( 0, 0, 1 );
		glVertex3f( 0, 0, 0 );
		glVertex3f( 0, 0, length );
	glEnd();
	glPopAttrib();
}

/**
	Draw rotation ring using red, green, blue for each ring

	called: PushAttrib once

	\param	radius		the radius of the ring
	\param	slices		slices of the ring
*/
inline void drawRotationRings(float radius = 1, int slices = 32){
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable( GL_LIGHTING );

	//	X ring
	glColor3f(1, 0, 0);
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(0, x, y);
		}
	glEnd();

	//	Y ring
	glColor3f(0, 1, 0);
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(x, 0, y);
		}
	glEnd();

	//	Z ring
	glColor3f(0, 0, 1);
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(x, y, 0);
		}
	glEnd();

	glPopAttrib();
}

/**
	Draw rings at the origin, perpendicular to X, Y, Z axis \n
	Similar to drawRotationRings but don't have color, don't change
	lighting condition

	\param	radius		the radius of the ring
	\param	slices		slices of the ring
*/
inline void drawXYZRings(float radius = 1, int slices = 32){
	//	X ring
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(0, x, y);
		}
	glEnd();

	//	Y ring
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(x, 0, y);
		}
	glEnd();

	//	Z ring
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			float angle = i*2*YZ_PI/slices;
			float x = radius * cos(angle);
			float y = radius * sin(angle);
			glVertex3f(x, y, 0);
		}
	glEnd();
}

/**
	Draw a camera in camera coordinate system (x right, y down, z front)

	called: PushAttrib once

	\param	z_near		near plane of camera
	\param	z_far		far plane of camera
	\param	fov_x		field of view of x direction of camera
	\param	fov_y		field of view of y direction of camera
*/
inline void drawCamera(float z_near=0, float z_far=1, float fov_x=64, float fov_y=48){

	float x_near = fabs( tanf( fov_x/2 * 3.1415926/180 ) * z_near );
	float y_near = fabs( tanf( fov_y/2 * 3.1415926/180 ) * z_near );
	float x_far = fabs( tanf( fov_x/2 * 3.1415926/180 ) * z_far );
	float y_far	= fabs( tanf( fov_y/2 * 3.1415926/180 ) * z_far );

	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
		//	near
		glVertex3f( x_near, y_near, z_near );
		glVertex3f( -x_near, y_near, z_near );

		glVertex3f( -x_near, y_near, z_near );
		glVertex3f( -x_near, -y_near, z_near );

		glVertex3f( -x_near, -y_near, z_near );
		glVertex3f( x_near, -y_near, z_near );

		glVertex3f( x_near, -y_near, z_near );
		glVertex3f( x_near, y_near, z_near );

		//	near to far
		glVertex3f( x_near, y_near, z_near );
		glVertex3f( x_far, y_far, z_far );

		glVertex3f( -x_near, y_near, z_near );
		glVertex3f( -x_far, y_far, z_far );

		glVertex3f( -x_near, -y_near, z_near );
		glVertex3f( -x_far, -y_far, z_far );

		glVertex3f( x_near, -y_near, z_near );
		glVertex3f( x_far, -y_far, z_far );

		//	far
		glVertex3f( x_far, y_far, z_far );
		glVertex3f( -x_far, y_far, z_far );

		glVertex3f( -x_far, y_far, z_far );
		glVertex3f( -x_far, -y_far, z_far );

		glVertex3f( -x_far, -y_far, z_far );
		glVertex3f( x_far, -y_far, z_far );

		glVertex3f( x_far, -y_far, z_far );
		glVertex3f( x_far, y_far, z_far );

		//	up
		glVertex3f( x_far*0.5f, -y_far, z_far );
		glVertex3f( 0, -y_far*1.5f, z_far );

		glVertex3f( -x_far*0.5f, -y_far, z_far );
		glVertex3f( 0, -y_far*1.5f, z_far );
	glEnd();
	glPopAttrib();
}

/**
	Draw a flat solid bone of length 1 from original point to +y direction, front +z
*/
inline void drawFlatSolidBone(){
	float vertex[] = {
		0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.2f, 0.3f, 0.0f,
		0.0f, 0.3f, 0.2f,
		-0.2f, 0.3f, 0.0f
	};
	int face[] = {
		0, 2, 3,
		2, 1, 3,
		0, 3, 4,
		3, 1, 4,
		0, 4, 2,
		1, 2, 4
	};

	float face_normal[] = {
		0.639602f,	-0.426401f,	0.639602f,
		0.693103f,	0.19803f,	0.693103f,
		-0.639602f,	-0.426401f,	0.639602f,
		-0.693103f,	0.19803f,	0.693103f,
		0.0f,		0.0f,		-1.0f,
		0.0f,		-0.0f,		-1.0f
	};

	glBegin( GL_TRIANGLES );
		for( int i=0; i<6; i++ ){
			glNormal3f( face_normal[i*3], face_normal[i*3+1], face_normal[i*3+2] );
			glVertex3f( vertex[face[i*3  ]*3], vertex[face[i*3  ]*3+1], vertex[face[i*3  ]*3+2] );
			glVertex3f( vertex[face[i*3+1]*3], vertex[face[i*3+1]*3+1], vertex[face[i*3+1]*3+2] );
			glVertex3f( vertex[face[i*3+2]*3], vertex[face[i*3+2]*3+1], vertex[face[i*3+2]*3+2] );
		}
	glEnd();
}

/**
	Draw a bone as a cube of length 1 from original point to +y direction

	color of front: blue + green	\n
	color of back: blue	\n
	color of side: red	\n
	color of ends: green

	\param	width	width of other two edges, default 0.3
*/
inline void drawColorCubeBone(float width = 0.3f){
	float hw = width * 0.5f;
	float vertex[] = {
		hw,		0.0f,	hw,
		hw,		0.0f,	-hw,
		-hw,	0.0f,	-hw,
		-hw,	0.0f,	hw,
		hw,		1.0f,	hw,
		hw,		1.0f,	-hw,
		-hw,	1.0f,	-hw,
		-hw,	1.0f,	hw
	};
	int face[] = {
		0, 4, 7, 3,
		1, 5, 4, 0,
		2, 6, 5, 1,
		3, 7, 6, 2,
		4, 5, 6, 7,
		0, 3, 2, 1
	};
	float face_normal[] = {
		0.0f,	0.0f,	1.0f,
		1.0f,	0.0f,	0.0f,
		0.0f,	0.0f,	-1.0f,
		-1.0f,	0.0f,	0.0f,
		0.0f,	1.0f,	0.0f,
		0.0f,	-1.0f,	0.0f
	};
	unsigned char face_color[] = {
		0,		255,	255,
		255,	0,		0,
		0,		0,		255,
		255,	0,		0,
		0,		255,	0,
		0,		255,	0
	};

	
	glBegin( GL_QUADS );
	for( int i=0; i<6; i++ ){
		glNormal3fv( face_normal+i*3 );
		glColor3ubv( face_color+i*3 );
		glVertex3f( vertex[face[i*4  ]*3], vertex[face[i*4  ]*3+1], vertex[face[i*4  ]*3+2] );
		glVertex3f( vertex[face[i*4+1]*3], vertex[face[i*4+1]*3+1], vertex[face[i*4+1]*3+2] );
		glVertex3f( vertex[face[i*4+2]*3], vertex[face[i*4+2]*3+1], vertex[face[i*4+2]*3+2] );
		glVertex3f( vertex[face[i*4+3]*3], vertex[face[i*4+3]*3+1], vertex[face[i*4+3]*3+2] );
	}
	glEnd();
}

/**
	Draw a solid bone of length 1 from original point to +y direction, front +z

	This function assumes the current matrix mode is model view

	called: PushAttrib once, PushMatrix(MODEL_VIEW) twice
*/
inline void drawSolidBone(){
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glutSolidSphere(0.05, 8, 5);
	glPushMatrix();
	glTranslatef(0, 1, 0);
	glutSolidSphere(0.05, 8, 5);
	glPopMatrix();

	glPushMatrix();
	glRotatef(-90, 1, 0, 0);
	glRotatef(45, 0, 0, 1);

	assert( !isMatrixStackFull(GL_MODELVIEW) );	//	push twice

	glPushMatrix();
	glTranslatef(0, 0, 0.1);
	glutSolidCone(0.141421, 0.9, 4, 1);
	glRotatef(180, 1, 0, 0);
	glutSolidCone(0.141421, 0.1, 4, 1);
	glPopMatrix();

	glPopMatrix();
}

/**
	Draw a wire bone of length 1 from original point to +y direction, front +z

	This function assumes the current matrix mode is model view

	called: PushAttrib once, PushMatrix(MODEL_VIEW) twice
*/
inline void drawWireBone(){
	assert( !isMatrixStackFull(GL_MODELVIEW) );
	assert( !isAttributeStackFull() );

	glPushAttrib(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);

	drawXYZRings(0.05);
	glPushMatrix();
	glTranslatef(0, 1, 0);
	drawXYZRings(0.05);
	glPopMatrix();

	glPushMatrix();
	glRotatef(-90, 1, 0, 0);
	glRotatef(45, 0, 0, 1);

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();
	glTranslatef(0, 0, 0.1);
	glutWireCone(0.141421, 0.9, 4, 1);
	glRotatef(180, 1, 0, 0);
	glutWireCone(0.141421, 0.1, 4, 1);
	glPopMatrix();

	glPopMatrix();

	glPopAttrib();
}
///@}

//	========================================
///@{
/**	@name Draw Curve
*/
//	========================================

/**
	Draw a long curve on 2D canvas

	\param	curve			the curve to be displayed, xy_xy_
	\param	point_number	number of points in the curve
*/
template<typename T>
inline void drawCurve2D(const T* curve, int point_number){
	glBegin(GL_LINE_STRIP);
		for( int i=0; i<point_number; i++ ){
			glVertex2f(curve[i*2], curve[i*2+1]);
		}
	glEnd();
}
/**
	Draw a long curve in 3D space

	\param	curve			the curve to be displayed, xyz_xyz_
	\param	point_number	number of points in the curve
*/
template<typename T>
inline void drawCurve3D(const T* curve, int point_number){
	glBegin(GL_LINE_STRIP);
		for( int i=0; i<point_number; i++ ){
			glVertex3f(curve[i*3], curve[i*3+1], curve[i*3+2]);
		}
	glEnd();
}

/**
	Draw a polygon on 2D canvas

	\param	point			the polygon points to be displayed, xy_xy_
	\param	point_number	number of points in the polygon
*/
template<typename T>
inline void drawPolygon2D(const T* point, int point_number){
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<point_number; i++ ){
			glVertex2f(point[i*2], point[i*2+1]);
		}
	glEnd();
}
/**
	Draw a polygon in 3D space

	\param	point			the polygon points to be displayed, xyz_xyz_
	\param	point_number	number of points in the curve
*/
template<typename T>
inline void drawPolygon3D(const T* point, int point_number){
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<point_number; i++ ){
			glVertex3f(point[i*3], point[i*3+1], point[i*3+2]);
		}
	glEnd();
}

///@}

//	========================================
///@{
/**	@name Draw Point Cloud
*/
//	========================================
/**
	Draw point cloud of single color

	\param	point_cloud		coordinate of each point in point cloud, xyz_xyz
	\param	point_number	number of points
	\param	indices			indices of each point to display
*/
template<typename T, typename TI>
void drawPointCloud(
	const T*		point_cloud, 
	unsigned int	point_number, 
	const TI*		indices)
{
	if (!point_number)
		return;

	if( !Is_float<T>::check_type && !Is_double<T>::check_type ){
		std::cout << "error: drawPointCloud, only float and double is allowed" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawPointCloud, only int or unsigned int is allowed for indices" << std::endl;
		return;
	}

	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);

	DrawElements(
		GL_POINTS,
		point_cloud,
		(float*)NULL,
		(float*)NULL,
		indices,
		point_number
	);

	//	setup array
	//glEnableClientState(GL_VERTEX_ARRAY);
	//if( Is_float<T>::check_type )
	//	glVertexPointer(3, GL_FLOAT, 0, point_cloud);
	//else if( Is_double<T>::check_type )
	//	glVertexPointer(3, GL_DOUBLE, 0, point_cloud);

	//glDrawElements(GL_POINTS, point_number, GL_UNSIGNED_INT, indices);

	//glDisableClientState(GL_VERTEX_ARRAY);

	glPopAttrib();
}

/**
	Draw color point cloud

	\param	point_cloud		coordinate of each point in point cloud, xyz_xyz_, type should be float or double
	\param	color			color of each point, rgb_rgb_, type should be float or unsigned char
	\param	point_number	number of points
	\param	indices			indices of each point to display
*/
template<typename TV, typename TC, typename TI>
inline void drawColorPointCloud(
	const TV*		point_cloud, 
	const TC*		color, 
	unsigned int	point_number, 
	const TI*		indices)
{
	if (!point_number)
		return;

	if( !Is_float<TV>::check_type && !Is_double<TV>::check_type ){
		std::cout << "error: drawColorPointCloud, only float or double is allowed for point cloud" << std::endl;
		return;
	}
	if( !Is_float<TC>::check_type && !Is_char<TC>::check_type && !Is_unsigned_char<TC>::check_type ){
		std::cout << "error: drawColorPointCloud, only float, char or unsigned char is allowed for color" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawColorPointCloud, only int or unsigned int is allowed for indices" << std::endl;
		return;
	}

	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);

	DrawElements(
		GL_POINTS,
		point_cloud,
		(float*)NULL,
		color,
		indices,
		point_number
	);

	//	setup array
	//glEnableClientState(GL_VERTEX_ARRAY);
	//if( Is_float<TV>::check_type )
	//	glVertexPointer(3, GL_FLOAT, 0, point_cloud);
	//else if( Is_double<TV>::check_type )
	//	glVertexPointer(3, GL_DOUBLE, 0, point_cloud);

	//glEnableClientState(GL_COLOR_ARRAY);
	//if( Is_float<TC>::check_type )
	//	glColorPointer(3, GL_FLOAT, 0, color);
	//else if( Is_unsigned_char<TC>::check_type )
	//	glColorPointer(3, GL_UNSIGNED_BYTE, 0, color);

	//glDrawElements(GL_POINTS, point_number, GL_UNSIGNED_INT, indices);

	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_COLOR_ARRAY);

	glPopAttrib();
}

/**
	Draw shading point cloud

	\param	point_cloud		coordinate of each point in point cloud, xyz_xyz_
	\param	normal			normal of each point, xyz_xyz_
	\param	point_number	number of points
*/
template<typename T>
inline void drawShadingPointCloud(
	const T*		point_cloud, 
	const T*		normal, 
	unsigned int	point_number)
{
	if (!point_number)
		return;

	static std::vector<unsigned int> indices;
	if (indices.size() < point_number) {
		indices.reserve(point_number);
		unsigned int count = point_number - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	drawShadingPointCloud(point_cloud, normal, point_number, &indices[0]);
}

/**
	Draw shading point cloud

	\param	point_cloud		coordinate of each point in point cloud, xyz_xyz_
	\param	normal			normal of each point, xyz_xyz_
	\param	point_number	number of points
	\param	indices			indices of each point to display
*/
template<typename T, typename TI>
inline void drawShadingPointCloud(
	const T*		point_cloud, 
	const T*		normal, 
	unsigned int	point_number, 
	const TI*		indices)
{
	if (!point_number)
		return;

	if( !Is_float<T>::check_type && !Is_double<T>::check_type ){
		std::cout << "error: drawShadingPointCloud, only float or double is allowed for point cloud" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawShadingPointCloud, only int or unsigned int is allowed for indices" << std::endl;
		return;
	}

	DrawElements(
		GL_POINTS,
		point_cloud,
		normal,
		(float*)NULL,
		indices,
		point_number
	);

	//glEnableClientState(GL_VERTEX_ARRAY);
	//if( Is_float<T>::check_type )
	//	glVertexPointer(3, GL_FLOAT, 0, point_cloud);
	//else if( Is_double<T>::check_type )
	//	glVertexPointer(3, GL_DOUBLE, 0, point_cloud);

	//glEnableClientState(GL_NORMAL_ARRAY);
	//if( Is_float<T>::check_type )
	//	glNormalPointer(GL_FLOAT, 0, normal);
	//else if( Is_double<T>::check_type )
	//	glNormalPointer(GL_DOUBLE, 0, normal);

	//glDrawElements(GL_POINTS, point_number, GL_UNSIGNED_INT, indices);

	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);
}

/**
	Draw shading point cloud

	\param	point_cloud		coordinate of each point in point cloud, xyz_xyz_
	\param	normal			normal of each point, xyz_xyz_
	\param	color			color of each point, rgb_rgb_
	\param	point_number	number of points
	\param	indices			indices of each point to display
*/
template<typename TV, typename TC, typename TI>
inline void drawShadingColorPointCloud(
	const TV*		point_cloud, 
	const TV*		normal, 
	const TC*		color, 
	unsigned int	point_number, 
	const TI*		indices)
{
	if( !Is_float<TV>::check_type && !Is_double<TV>::check_type ){
		std::cout << "error: drawShadingColorPointCloud, only float or double is allowed for point cloud" << std::endl;
		return;
	}
	if( !Is_float<TC>::check_type && !Is_char<TC>::check_type && !Is_unsigned_char<TC>::check_type ){
		std::cout << "error: drawShadingColorPointCloud, only float, char or unsigned char is allowed for color" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawShadingColorPointCloud, only int or unsigned int is allowed for indices" << std::endl;
		return;
	}

	DrawElements(
		GL_POINTS,
		point_cloud,
		normal,
		color,
		indices,
		point_number
	);

	//glEnableClientState(GL_VERTEX_ARRAY);
	//if( Is_float<TV>::check_type )
	//	glVertexPointer(3, GL_FLOAT, 0, point_cloud);
	//else if( Is_double<TV>::check_type )
	//	glVertexPointer(3, GL_DOUBLE, 0, point_cloud);

	//glEnableClientState(GL_NORMAL_ARRAY);
	//if( Is_float<TV>::check_type )
	//	glNormalPointer(GL_FLOAT, 0, normal);
	//else if( Is_double<TV>::check_type )
	//	glNormalPointer(GL_DOUBLE, 0, normal);

	//glEnableClientState(GL_COLOR_ARRAY);
	//if( Is_float<TC>::check_type )
	//	glColorPointer(3, GL_FLOAT, 0, color);
	//else if( Is_unsigned_char<TC>::check_type )
	//	glColorPointer(3, GL_UNSIGNED_BYTE, 0, color);

	//glDrawElements(GL_POINTS, point_number, GL_UNSIGNED_INT, indices);

	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);
	//glDisableClientState(GL_COLOR_ARRAY);
}

///@}

//	========================================
///@{
/**	@name Draw Mesh 3D
*/
//	========================================
/**
	Draw Triangle Mesh with Flat Shading

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	face_normal		face normal list, xyz_xyz_
*/
template<typename T, typename TI>
void drawFlatShadingTriMesh(
	const T*		vertex, 
	unsigned int	vertex_number, 
	const TI*		face, 
	unsigned int	face_number, 
	const T*		face_normal)
{
	if (!vertex_number || !face_number)
		return;

	if (!Is_float<T>::check_type && !Is_double<T>::check_type) {
		std::cout << "error: drawFlatShadingTriMesh, only float or double is allowed for vertex and normal" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawFlatShadingTriMesh, only int or unsigned int is allowed for face" << std::endl;
		return;
	}

	//	create split vertex array
	static std::vector<T> split_vertex;
	if (split_vertex.size() < face_number * 9) {	//	3 vertices each face * 3 digits (xyz) each vertex
		split_vertex.resize(face_number * 9);
	}
	static std::vector<T> vertex_normal;
	if (vertex_normal.size() < face_number * 9) {
		vertex_normal.resize(face_number * 9);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 3) {
		indices.reserve(face_number * 3);
		unsigned int count = face_number * 3 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 3;
			unsigned int acc2 = acc;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			vertex_normal[acc2++] = face_normal[idx_num];
			vertex_normal[acc2++] = face_normal[idx_num + 1];
			vertex_normal[acc2++] = face_normal[idx_num + 2];
		}
	}

	DrawElements(
		GL_TRIANGLES,
		&split_vertex[0],
		&vertex_normal[0],
		(float*)NULL,
		&indices[0],
		face_number * 3
	);
}

/**
	Draw Triangle Mesh with Smooth Shading

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices, this parameter is actually not used in this function
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	vertex_normal	vertex normal list, xyz_xyz_
*/
template<typename T, typename TI>
void drawSmoothShadingTriMesh(
	const T*		vertex, 
	unsigned int	vertex_number, 
	const TI*		face,
	unsigned int	face_number,
	const T*		vertex_normal)
{
	if( !Is_float<T>::check_type && !Is_double<T>::check_type ){
		std::cout << "error: drawSmoothShadingTriMesh, only float or double is allowed for vertices" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawSmoothShadingTriMesh, only int or unsigned int is allowed for face" << std::endl;
		return;
	}

	DrawElements(
		GL_TRIANGLES,
		vertex,
		vertex_normal,
		(float*)NULL,
		face,
		face_number*3
	);
}

/**
	Draw Triangle Mesh with Given Face Color

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	face_color		face color list, rgb_rgb_
	\param	color_size		3 rgb, 4 rgba
*/
template<typename TV, typename TI, typename TC>
void drawFlatColorTriMesh(
	const TV*		vertex, 
	unsigned int	vertex_number, 
	const TI*		face, 
	unsigned int	face_number, 
	const TC*		face_color,
	unsigned int	color_size = 3 )
{
	if (!vertex_number || !face_number)
		return;

	if (!Is_float<TV>::check_type && !Is_double<TV>::check_type) {
		std::cout << "error: drawFlatColorTriMesh, only float or double is allowed for vertex and normal" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawFlatColorTriMesh, only int or unsigned int is allowed for face" << std::endl;
		return;
	}
	if (!Is_float<TC>::check_type && !Is_char<TC>::check_type && !Is_unsigned_char<TC>::check_type) {
		std::cout << "error: drawFlatColorTriMesh, only float, char or unsigned char is allowed for color" << std::endl;
		return;
	}

	//	create split vertex array
	static std::vector<TV> split_vertex;
	if (split_vertex.size() < face_number * 9) {	//	3 vertices each face * 3 digits (xyz) each vertex
		split_vertex.resize(face_number * 9);
	}
	static std::vector<TC> vertex_color;
	if (vertex_color.size() < face_number * 3 * color_size) {
		vertex_color.resize(face_number * 3 * color_size);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 3) {
		indices.reserve(face_number * 3);
		unsigned int count = face_number * 3 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	unsigned int color_acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		unsigned int color_idx_num = f * color_size;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 3;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			for (unsigned int j = 0; j != color_size; j++) {
				vertex_color[color_acc++] = face_color[color_idx_num + j];
			}
		}
	}

	DrawElements(
		GL_TRIANGLES,
		&split_vertex[0],
		(TV*)NULL,
		&vertex_color[0],
		&indices[0],
		face_number * 3,
		color_size
	);

}

/**
	Draw Triangle Mesh with Given Vertex Color

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	vertex_color	vertex color list, rgb_rgb_
	\param	color_size		3 rgb, 4 rgba
*/
template<typename TV, typename TI, typename TC>
void drawSmoothColorTriMesh(
	const TV*		vertex, 
	unsigned int	vertex_number, 
	const TI*		face, 
	unsigned int	face_number, 
	const TC*		vertex_color,
	unsigned int	color_size = 3)
{
	if (!Is_float<TV>::check_type && !Is_double<TV>::check_type) {
		std::cout << "error: drawSmoothColorTriMesh, only float or double is allowed for vertex and normal" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawSmoothColorTriMesh, only int or unsigned int is allowed for face" << std::endl;
		return;
	}
	if (!Is_float<TC>::check_type && !Is_char<TC>::check_type && !Is_unsigned_char<TC>::check_type) {
		std::cout << "error: drawSmoothColorTriMesh, only float, char or unsigned char is allowed for color" << std::endl;
		return;
	}

	DrawElements(
		GL_TRIANGLES,
		vertex,
		(TV*)NULL,
		vertex_color,
		face,
		face_number * 3,
		color_size
	);
}

/**
	Draw Triangle Mesh with Sequential Color

	This function can draw mesh with different color of each face. 
	It can also be used for drawing picking image

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	start_index		index of start sequential color
*/
template<typename T>
void drawSequentialColorTriMesh(
	const T*		vertex, 
	unsigned int	vertex_number, 
	const int*		face, 
	unsigned int	face_number,
	int				start_index = 0)
{
	if (!vertex_number || !face_number)
		return;

	static std::vector<uchar4> face_color;
	if (face_color.size() < face_number)
		face_color.resize(face_number);
	for (unsigned int i = 0; i != face_number; i++)
		getSequentialColor(&face_color[i].x, start_index + i);

	drawFlatColorTriMesh(
		vertex,
		vertex_number,
		face,
		face_number,
		&face_color[0].x,
		4
	);
}

/**
	Draw Mesh vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	scale			scale of each normal
*/
template<typename T>
void drawVertexNormal(
	const T*		vertex, 
	unsigned int	vertex_number, 
	const T*		vertex_normal, 
	float			scale = 1.0f)
{
	if (!vertex || !vertex_number)
		return;

	//	setup array
	static std::vector<T> line_vertex;
	if (line_vertex.size() < vertex_number * 6) {
		line_vertex.resize(vertex_number * 6);
	}

	static std::vector<unsigned int> indices;
	if (indices.size() < vertex_number * 2) {
		indices.reserve(vertex_number * 2);
		unsigned int count = vertex_number * 2 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	create vertices
	unsigned int count1 = 0, count2 = 0;
	for (unsigned int i = 0; i < vertex_number; i++) {
		line_vertex[count1++] = vertex[count2++];
		line_vertex[count1++] = vertex[count2++];
		line_vertex[count1++] = vertex[count2++];
		count2 -= 3;
		line_vertex[count1++] = vertex[count2++];
		line_vertex[count1++] = vertex[count2++];
		line_vertex[count1++] = vertex[count2++];
		count1 -= 3;
		count2 -= 3;
		line_vertex[count1++] += vertex_normal[count2++] * scale;
		line_vertex[count1++] += vertex_normal[count2++] * scale;
		line_vertex[count1++] += vertex_normal[count2++] * scale;
	}

	GLboolean lighting_flag;
	glGetBooleanv(GL_LIGHTING, &lighting_flag);
	glDisable(GL_LIGHTING);

	DrawElements(
		GL_LINES,
		&line_vertex[0],
		(float*)NULL,
		(float*)NULL,
		&indices[0],
		vertex_number * 2
	);

	if (lighting_flag)
		glEnable(GL_LIGHTING);
}

/**
	Draw face normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	face_normal		normal of each face
	\param	scale			scale of each normal
*/
template<typename T>
void drawFaceNormal(
	const T*		vertex,
	const int*		face,
	unsigned int	face_number,
	const T*		face_normal = NULL,
	float			scale = 1.0f
) {
	if (!vertex || !face || !face_number)
		return;

	const Vec3<T>*	v_ptr = (const Vec3<T>*)vertex;
	const int3*		f_ptr = (const int3*)face;
	const Vec3<T>*	n_ptr = (const Vec3<T>*)face_normal;

	//	setup array
	static std::vector<Vec3<T>> line_vertex;
	if (line_vertex.size() < face_number * 2) {
		line_vertex.resize(face_number * 2);
	}

	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 2) {
		indices.reserve(face_number * 2);
		unsigned int count = face_number * 2 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	create vertices
	for (unsigned int i = 0; i != face_number; i++) {
		line_vertex[i * 2] = (v_ptr[f_ptr[i].x] + v_ptr[f_ptr[i].y] + v_ptr[f_ptr[i].z]) * 0.333333;
		if (n_ptr)
			line_vertex[i * 2 + 1] = line_vertex[i * 2] + n_ptr[i] * scale;
		else {
			Vec3<T> n = cross(v_ptr[f_ptr[i].y] - v_ptr[f_ptr[i].x], v_ptr[f_ptr[i].z] - v_ptr[f_ptr[i].x]);
			line_vertex[i * 2 + 1] = line_vertex[i * 2] + n.Normalize() * scale;
		}
	}

	GLboolean lighting_flag;
	glGetBooleanv(GL_LIGHTING, &lighting_flag);
	glDisable(GL_LIGHTING);

	DrawElements(
		GL_LINES,
		&line_vertex[0].x,
		(T*)NULL,
		(T*)NULL,
		&indices[0],
		face_number * 2
	);

	if (lighting_flag)
		glEnable(GL_LIGHTING);

}

/**
	Draw Mesh Edge, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
*/
template<typename T>
void drawMeshEdge(
	const T* vertex, int vertex_number, 
	const int* edge, int edge_number)
{
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshEdge(
	const T* vertex, int vertex_number, 
	const int* edge, int edge_number, 
	const T* vertex_normal, float offset = 0.001f)
{
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			float nor0[3] = {vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]};
			float nor1[3] = {vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]};
			float v0[3] = {vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]};
			float v1[3] = {vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]};
			for( int j=0; j<3; j++ ){
				v0[j] += nor0[j] * offset;
				v1[j] += nor1[j] * offset;
			}
			glNormal3fv( nor0 );
			glVertex3fv( v0 );
			glNormal3fv( nor1 );
			glVertex3fv( v1 );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given vertex color, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(
	const T* vertex, int vertex_number,
	const int* edge, int edge_number,
	const float* vertex_color)
{
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3f( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glColor3f( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given vertex color, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const T* vertex, int vertex_number,
							 const int* edge, int edge_number,
							 const unsigned char* vertex_color){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3ub( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glColor3ub( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given vertex color, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const T* vertex, int vertex_number, 
							 const int* edge, int edge_number, 
							 const float* vertex_color,
							 const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			float nor0[3] = {vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]};
			float nor1[3] = {vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]};
			float v0[3] = {vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]};
			float v1[3] = {vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]};
			for( int j=0; j<3; j++ ){
				v0[j] += nor0[j] * offset;
				v1[j] += nor1[j] * offset;
			}
			glColor3f( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glNormal3fv( nor0 );
			glVertex3fv( v0 );
			glColor3f( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glNormal3fv( nor1 );
			glVertex3fv( v1 );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given vertex color, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const T* vertex, int vertex_number, 
							 const int* edge, int edge_number, 
							 const unsigned char* vertex_color,
							 const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			float nor0[3] = {vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]};
			float nor1[3] = {vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]};
			float v0[3] = {vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]};
			float v1[3] = {vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]};
			for( int j=0; j<3; j++ ){
				v0[j] += nor0[j] * offset;
				v1[j] += nor1[j] * offset;
			}
			glColor3ub( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glNormal3fv( nor0 );
			glVertex3fv( v0 );
			glColor3ub( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glNormal3fv( nor1 );
			glVertex3fv( v1 );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given edge color, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_color		color of each edge
*/
template<typename T>
void drawFlatColorMeshEdge(const T* vertex, int vertex_number,
						   const int* edge, int edge_number,
						   const float* edge_color){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3fv( edge_color + i*3 );
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given edge color, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_color		color of each edge
*/
template<typename T>
void drawFlatColorMeshEdge(const T* vertex, int vertex_number,
						   const int* edge, int edge_number,
						   const unsigned char* edge_color){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3ubv( edge_color + i*3 );
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given edge color, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_color		color of each edge
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawFlatColorMeshEdge(const T* vertex, int vertex_number, 
						   const int* edge, int edge_number, 
						   const float* edge_color,
						   const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			Vec3f nor0(vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]);
			Vec3f nor1(vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]);
			Vec3f v0(vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]);
			Vec3f v1(vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]);
			v0 += nor0 * offset;
			v1 += nor1 * offset;
			Vec3f nor = (nor0 + nor1).Normalize();
			glColor3fv( edge_color + i*3 );
			glNormal3f( nor[0], nor[1], nor[2] );
			glVertex3f( v0[0], v0[1], v0[2] );
			glVertex3f( v1[0], v1[1], v1[2] );
		}
	glEnd();
}

/**
	Draw Mesh Edge with given edge color, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_color		color of each edge
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawFlatColorMeshEdge(const T* vertex, int vertex_number, 
						   const int* edge, int edge_number, 
						   const unsigned char* edge_color,
						   const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			Vec3f nor0(vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]);
			Vec3f nor1(vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]);
			Vec3f v0(vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]);
			Vec3f v1(vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]);
			v0 += nor0 * offset;
			v1 += nor1 * offset;
			Vec3f nor = (nor0 + nor1).Normalize();
			glColor3ubv( edge_color + i*3 );
			glNormal3f( nor[0], nor[1], nor[2] );
			glVertex3f( v0[0], v0[1], v0[2] );
			glVertex3f( v1[0], v1[1], v1[2] );
		}
	glEnd();
}

/**
	Draw Mesh Boundary from edge-face connectivity, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	ef				edge - face connectivity
*/
template<typename T>
void drawMeshBoundaryFromEF(const T* vertex, int vertex_number, 
							const int* edge, int edge_number,
							const int* ef){
	glBegin( GL_LINES );
	for( int i=0; i<edge_number; i++ ){
		if( ef[i*2] == -1 || ef[i*2+1] == -1 ){
			glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
			glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
		}
	}
	glEnd();
}
/**
	Draw Mesh Boundary from edge-face connectivity, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	ef				edge - face connectivity
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshBoundaryFromEF(const T* vertex, int vertex_number, 
							const int* edge, int edge_number, 
							const int* ef,
							const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			if( ef[i*2] == -1 || ef[i*2+1] == -1 ){
				float nor0[3] = {vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]};
				float nor1[3] = {vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]};
				float v0[3] = {vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]};
				float v1[3] = {vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]};
				for( int j=0; j<3; j++ ){
					v0[j] += nor0[j] * offset;
					v1[j] += nor1[j] * offset;
				}
				glNormal3fv( nor0 );
				glVertex3fv( v0 );
				glNormal3fv( nor1 );
				glVertex3fv( v1 );
			}
		}
	glEnd();
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_mark		boundary mark of each edge, 1: is boundary; 0: not boundary
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark(const T* vertex, int vertex_number, 
								  const int* edge, int edge_number,
								  const int* edge_mark){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			if( edge_mark[i] ){
				glVertex3f( vertex[edge[i*2  ]*3], vertex[edge[i*2  ]*3+1], vertex[edge[i*2  ]*3+2] );
				glVertex3f( vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2] );
			}
		}
	glEnd();
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_mark		boundary mark of each edge, 1: is boundary; 0: not boundary
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark(const T* vertex, int vertex_number, 
								  const int* edge, int edge_number, 
								  const int* edge_mark,
								  const T* vertex_normal, float offset = 0.001f){
	glBegin( GL_LINES );
	for( int i=0; i<edge_number; i++ ){
	  if( edge_mark[i] ){
		  float nor0[3] = {vertex_normal[edge[i*2]*3], vertex_normal[edge[i*2]*3+1], vertex_normal[edge[i*2]*3+2]};
		  float nor1[3] = {vertex_normal[edge[i*2+1]*3], vertex_normal[edge[i*2+1]*3+1], vertex_normal[edge[i*2+1]*3+2]};
		  float v0[3] = {vertex[edge[i*2]*3], vertex[edge[i*2]*3+1], vertex[edge[i*2]*3+2]};
		  float v1[3] = {vertex[edge[i*2+1]*3], vertex[edge[i*2+1]*3+1], vertex[edge[i*2+1]*3+2]};
		  for( int j=0; j<3; j++ ){
			  v0[j] += nor0[j] * offset;
			  v1[j] += nor1[j] * offset;
		  }
		  glNormal3fv( nor0 );
		  glVertex3fv( v0 );
		  glNormal3fv( nor1 );
		  glVertex3fv( v1 );
	  }
	}
	glEnd();
}


/**
	Draw Mesh Edge from face, no offset from vertices

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
*/
template<typename T>
void drawMeshEdgeFromFace(
	const T*		vertex,
	unsigned int	vertex_number,
	const int*		face,
	unsigned int	face_number
) {
	//	create split vertex array
	static std::vector<T> split_vertex;
	if (split_vertex.size() < face_number * 18) {	//	3 edges each face * 2 vertices each edge * 3 digits (xyz) each vertex
		split_vertex.resize(face_number * 18);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 6) {
		indices.reserve(face_number * 6);
		unsigned int count = face_number * 6 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 3;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			v = face[idx_num + (i + 1) % 3] * 3;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
		}
	}

	DrawElements(
		GL_LINES,
		&split_vertex[0],
		(T*)NULL,
		(T*)NULL,
		&indices[0],
		face_number * 6
	);
}

/**
	Draw Mesh Edge from face, offset vertices according to vertex normal

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	vertex_normal	vertex normal list, xyz_xyz_
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshEdgeFromFace(
	const T*		vertex,
	unsigned int	vertex_number,
	const int*		face,
	unsigned int	face_number,
	const T*		vertex_normal,
	float			offset = 0.001f
) {
	//	create split vertex array
	static std::vector<T> split_vertex;
	if (split_vertex.size() < face_number * 18) {	//	3 edges each face * 2 vertices each edge * 3 digits (xyz) each vertex
		split_vertex.resize(face_number * 18);
	}
	static std::vector<T> split_vertex_normal;
	if (split_vertex_normal.size() < face_number * 18) {	//	3 edges each face * 2 vertices each edge * 3 digits (xyz) each vertex
		split_vertex_normal.resize(face_number * 18);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 6) {
		indices.reserve(face_number * 6);
		unsigned int count = face_number * 6 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	unsigned int nor_acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 3;
			split_vertex[acc++] = vertex[v] + vertex_normal[v] * offset;
			split_vertex[acc++] = vertex[v + 1] + vertex_normal[v + 1] * offset;
			split_vertex[acc++] = vertex[v + 2] + vertex_normal[v + 2] * offset;
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
			v = face[idx_num + (i + 1) % 3] * 3;
			split_vertex[acc++] = vertex[v] + vertex_normal[v] * offset;
			split_vertex[acc++] = vertex[v + 1] + vertex_normal[v + 1] * offset;
			split_vertex[acc++] = vertex[v + 2] + vertex_normal[v + 2] * offset;
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
			split_vertex_normal[nor_acc++] = vertex_normal[v++];
		}
	}

	DrawElements(
		GL_LINES,
		&split_vertex[0],
		&split_vertex_normal[0],
		(T*)NULL,
		&indices[0],
		face_number * 6
	);
}

/**
	Draw mesh edge from face, with flat shading and color

	\param	vertex			vertex coordinate list, xyz_xyz_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	face_color		face color list, rgb_rgb_
	\param	face_normal		face normal list, xyz_xyz_
	\param	color_size		3 rgb, 4 rgba
*/
template<typename TV, typename TI, typename TC>
void drawFlatColorTriMeshEdgeFromFace(
	const TV*		vertex,
	unsigned int	vertex_number,
	const TI*		face,
	unsigned int	face_number,
	const TC*		face_color,
	const TV*		face_normal,
	unsigned int	color_size = 3
) {
	if (!vertex_number || !face_number)
		return;

	if (!Is_float<TV>::check_type && !Is_double<TV>::check_type) {
		std::cout << "error: drawFlatColorTriMeshEdgeFromFace, only float or double is allowed for vertex and normal" << std::endl;
		return;
	}
	if (!Is_int<TI>::check_type && !Is_unsigned_int<TI>::check_type) {
		std::cout << "error: drawFlatColorTriMeshEdgeFromFace, only int or unsigned int is allowed for face" << std::endl;
		return;
	}
	if (!Is_float<TC>::check_type && !Is_char<TC>::check_type && !Is_unsigned_char<TC>::check_type) {
		std::cout << "error: drawFlatColorTriMeshEdgeFromFace, only float, char or unsigned char is allowed for color" << std::endl;
		return;
	}

	//	create split vertex array
	static std::vector<TV> split_vertex;
	if (split_vertex.size() < face_number * 18) {	//	3 edges each face * 2 vertices each edge * 3 digits (xyz) each vertex
		split_vertex.resize(face_number * 18);
	}
	static std::vector<TV> vertex_normal;
	if (vertex_normal.size() < face_number * 18) {	//	3 edges each face * 2 vertices each edge * 3 digits (xyz) each vertex
		vertex_normal.resize(face_number * 18);
	}
	static std::vector<TC> vertex_color;
	if (vertex_color.size() < face_number * 6 * color_size) {	//	3 edges each face * 2 vertices each edge * color_size
		vertex_color.resize(face_number * 6 * color_size);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 6) {
		indices.reserve(face_number * 6);
		unsigned int count = face_number * 6 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		for (unsigned int i = 0; i != 3; i++) {
			//	create split vertex
			unsigned int v = face[idx_num + i] * 3;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			v = face[idx_num + (i + 1) % 3] * 3;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];

			//	create split vertex normal
			acc -= 6;
			vertex_normal[acc++] = face_normal[idx_num];
			vertex_normal[acc++] = face_normal[idx_num + 1];
			vertex_normal[acc++] = face_normal[idx_num + 2];
			vertex_normal[acc++] = face_normal[idx_num];
			vertex_normal[acc++] = face_normal[idx_num + 1];
			vertex_normal[acc++] = face_normal[idx_num + 2];

			//	create split vertex color
			acc -= 6;
			vertex_color[acc++] = face_color[idx_num];
			vertex_color[acc++] = face_color[idx_num + 1];
			vertex_color[acc++] = face_color[idx_num + 2];
			vertex_color[acc++] = face_color[idx_num];
			vertex_color[acc++] = face_color[idx_num + 1];
			vertex_color[acc++] = face_color[idx_num + 2];
		}
	}

	DrawElements(
		GL_LINES,
		&split_vertex[0],
		&vertex_normal[0],
		&vertex_color[0],
		&indices[0],
		face_number * 6,
		color_size
	);
}

/**
	draw index of each vertex of a mesh

	\param	vertex			vertex position of the mesh
	\param	vertex_number	number of vertices
*/
template<typename T>
void drawMeshVertexIndex(const T*	vertex,
						 int		vertex_number){
	for( int i=0; i<vertex_number; i++ )
		drawNumber(i, vertex[i*3], vertex[i*3+1], vertex[i*3+2]);
}

/**
	draw number on each vertex of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	vertex_number	number of vertices
*/
template<typename TN, typename TV>
void drawNumberOnMeshVertices(const TN* number,
							  const TV* vertex,
							  int vertex_number){
	for( int i=0; i<vertex_number; i++ )
		drawNumber(number[i], vertex[i*3], vertex[i*3+1], vertex[i*3+2]);
}

/**
	draw number on each edge of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	edge			edges list
	\param	edge_number		number of edges
*/
template<typename TN, typename TV>
void drawNumberOnMeshEdges(const TN* number,
						   const TV* vertex,
						   const int* edge,
						   int edge_number){
	for( int i=0; i<edge_number; i++ ){
		int v1 = edge[i*2];
		int v2 = edge[i*2+1];
		drawNumber(number[i],
			(vertex[v1*3  ]+vertex[v2*3  ]) * 0.5,
			(vertex[v1*3+1]+vertex[v2*3+1]) * 0.5,
			(vertex[v1*3+2]+vertex[v2*3+2]) * 0.5 );
	}
}

///@}

//	========================================
///@{
/**	@name Draw Mesh 2D
*/
//	========================================

/**
	Draw Triangle Mesh on 2D plane

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
*/
template<typename T>
void drawTriMesh2D(const T* vertex, int vertex_number, 
				   const int* face, int face_number ){
	glBegin( GL_TRIANGLES );
		for( int i=0; i<face_number; i++ ){
			glVertex2f( vertex[face[i*3  ]*2], vertex[face[i*3  ]*2+1] );
			glVertex2f( vertex[face[i*3+1]*2], vertex[face[i*3+1]*2+1] );
			glVertex2f( vertex[face[i*3+2]*2], vertex[face[i*3+2]*2+1] );
		}
	glEnd();
}

/**
	Draw Triangle Mesh with Given Face Color on 2D

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	face_color		face color list, rgb_rgb_
*/
template<typename TV, typename TI, typename TC>
void drawFlatColorTriMesh2D(
	const TV*		vertex,
	int				vertex_number,
	const TI*		face,
	int				face_number,
	const TC*		face_color,
	unsigned int	color_size = 3)
{
	//	create split vertex array
	static std::vector<TV> split_vertex;
	if (split_vertex.size() < face_number * 6) {	//	3 vertices each face * 2 digits (xyz) each vertex
		split_vertex.resize(face_number * 6);
	}
	static std::vector<TC> vertex_color;
	if (vertex_color.size() < face_number * 3 * color_size) {
		vertex_color.resize(face_number * 3 * color_size);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 3) {
		indices.reserve(face_number * 3);
		unsigned int count = face_number * 3 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	unsigned int color_acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		unsigned int color_idx_num = f * color_size;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 2;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			for (unsigned int j = 0; j != color_size; j++) {
				vertex_color[color_acc++] = face_color[color_idx_num + j];
			}
		}
	}

	DrawElements(
		GL_TRIANGLES,
		&split_vertex[0],
		(TV*)NULL,
		&vertex_color[0],
		&indices[0],
		face_number * 3,
		color_size,
		2
	);
}

/**
	Draw Triangle Mesh with Given Vertex Color on 2D

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	vertex_color	vertex color list, rgb_rgb_
*/
template<typename T>
void drawSmoothColorTriMesh2D(const T* vertex, int vertex_number, 
							  const int* face, int face_number, 
							  const float* vertex_color){
	glBegin( GL_TRIANGLES );
		for( int i=0; i<face_number; i++ ){
			glColor3f( vertex_color[face[i*3  ]*3], vertex_color[face[i*3  ]*3+1], vertex_color[face[i*3  ]*3+2] );
			glVertex2f( vertex[face[i*3  ]*2], vertex[face[i*3  ]*2+1] );
			glColor3f( vertex_color[face[i*3+1]*3], vertex_color[face[i*3+1]*3+1], vertex_color[face[i*3+1]*3+2] );
			glVertex2f( vertex[face[i*3+1]*2], vertex[face[i*3+1]*2+1] );
			glColor3f( vertex_color[face[i*3+2]*3], vertex_color[face[i*3+2]*3+1], vertex_color[face[i*3+2]*3+2] );
			glVertex2f( vertex[face[i*3+2]*2], vertex[face[i*3+2]*2+1] );
		}
	glEnd();
}

/**
	Draw Triangle Mesh with Given Vertex Color on 2D

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
	\param	vertex_color	vertex color list, rgb_rgb_
*/
template<typename T>
void drawSmoothColorTriMesh2D(const T* vertex, int vertex_number, 
							  const int* face, int face_number, 
							  const unsigned char* vertex_color){
	glBegin( GL_TRIANGLES );
		for( int i=0; i<face_number; i++ ){
			glColor3ub( vertex_color[face[i*3  ]*3], vertex_color[face[i*3  ]*3+1], vertex_color[face[i*3  ]*3+2] );
			glVertex2f( vertex[face[i*3  ]*2], vertex[face[i*3  ]*2+1] );
			glColor3ub( vertex_color[face[i*3+1]*3], vertex_color[face[i*3+1]*3+1], vertex_color[face[i*3+1]*3+2] );
			glVertex2f( vertex[face[i*3+1]*2], vertex[face[i*3+1]*2+1] );
			glColor3ub( vertex_color[face[i*3+2]*3], vertex_color[face[i*3+2]*3+1], vertex_color[face[i*3+2]*3+2] );
			glVertex2f( vertex[face[i*3+2]*2], vertex[face[i*3+2]*2+1] );
		}
	glEnd();
}

/**
	Draw Mesh Edge 2D on plane

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
*/
template<typename T>
void drawMeshEdge2D(const T* vertex, int vertex_number, 
					const int* edge, int edge_number){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glVertex2f( vertex[edge[i*2  ]*2], vertex[edge[i*2  ]*2+1] );
			glVertex2f( vertex[edge[i*2+1]*2], vertex[edge[i*2+1]*2+1] );
		}
	glEnd();
}

/**
	Draw 2D Mesh Edge with given vertex color

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge2D(const T* vertex, int vertex_number,
							   const int* edge, int edge_number,
							   const float* vertex_color){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3f( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glVertex2f( vertex[edge[i*2  ]*2], vertex[edge[i*2  ]*2+1] );
			glColor3f( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glVertex2f( vertex[edge[i*2+1]*2], vertex[edge[i*2+1]*2+1] );
		}
	glEnd();
}

/**
	Draw 2D Mesh Edge with given vertex color

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge2D(const T* vertex, int vertex_number,
							   const int* edge, int edge_number,
							   const unsigned char* vertex_color){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			glColor3ub( vertex_color[edge[i*2  ]*3], vertex_color[edge[i*2  ]*3+1], vertex_color[edge[i*2  ]*3+2] );
			glVertex2f( vertex[edge[i*2  ]*2], vertex[edge[i*2  ]*2+1] );
			glColor3ub( vertex_color[edge[i*2+1]*3], vertex_color[edge[i*2+1]*3+1], vertex_color[edge[i*2+1]*3+2] );
			glVertex2f( vertex[edge[i*2+1]*2], vertex[edge[i*2+1]*2+1] );
		}
	glEnd();
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	ef				edge - face connectivity
*/
template<typename T>
void drawMeshBoundaryFromEF2D(const T* vertex, int vertex_number, 
							  const int* edge, int edge_number,
							  const int* ef){
	glBegin( GL_LINES );
	for( int i=0; i<edge_number; i++ ){
		if( ef[i*2] == -1 || ef[i*2+1] == -1 ){
			glVertex2f( vertex[edge[i*2  ]*2], vertex[edge[i*2  ]*2+1] );
			glVertex2f( vertex[edge[i*2+1]*2], vertex[edge[i*2+1]*2+1] );
		}
	}
	glEnd();
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_number		number of edge
	\param	edge_mark		boundary mark of each edge, 1: is boundary; 0: not boundary
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark2D(const T* vertex, int vertex_number, 
									const int* edge, int edge_number,
									const int* edge_mark){
	glBegin( GL_LINES );
		for( int i=0; i<edge_number; i++ ){
			if( edge_mark[i] ){
				glVertex2f( vertex[edge[i*2  ]*2], vertex[edge[i*2  ]*2+1] );
				glVertex2f( vertex[edge[i*2+1]*2], vertex[edge[i*2+1]*2+1] );
			}
		}
	glEnd();
}

/**
	Draw 2D Mesh Edge from face

	\param	vertex			vertex coordinate list, xy_xy_
	\param	vertex_number	number of vertices
	\param	face			face index list, v0v1v2_v0v1v2_
	\param	face_number		number of face
*/
template<typename T>
void drawMeshEdgeFromFace2D(
	const T*	vertex,
	int			vertex_number,
	const int*	face,
	int			face_number)
{
	//	create split vertex array
	static std::vector<T> split_vertex;
	if (split_vertex.size() < face_number * 12) {	//	3 edges each face * 2 vertices each edge * 2 digits (xy) each vertex
		split_vertex.resize(face_number * 12);
	}
	static std::vector<unsigned int> indices;
	if (indices.size() < face_number * 6) {
		indices.reserve(face_number * 6);
		unsigned int count = face_number * 6 - indices.size();
		do {
			indices.push_back(indices.size());
		} while (--count);
	}

	//	split each face
	unsigned int acc = 0;
	for (unsigned int f = 0; f != face_number; f++) {
		unsigned int idx_num = f * 3;
		for (unsigned int i = 0; i != 3; i++) {
			unsigned int v = face[idx_num + i] * 2;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
			v = face[idx_num + (i + 1) % 3] * 2;
			split_vertex[acc++] = vertex[v++];
			split_vertex[acc++] = vertex[v++];
		}
	}

	DrawElements(
		GL_LINES,
		&split_vertex[0],
		(T*)NULL,
		(T*)NULL,
		&indices[0],
		face_number * 6,
		3,
		2
	);
}

/**
	draw index of each vertex of a mesh

	\param	vertex			vertex position of the mesh
	\param	vertex_number	number of vertices
*/
template<typename T>
void drawMeshVertexIndex2D(const T*	vertex,
						   int		vertex_number){
	for( int i=0; i<vertex_number; i++ )
		drawNumber(i, vertex[i*2], vertex[i*2+1]);
}

/**
	draw number on each vertex of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	vertex_number	number of vertices
*/
template<typename TN, typename TV>
void drawNumberOnMeshVertices2D(const TN* number,
								const TV* vertex,
								int vertex_number){
	for( int i=0; i<vertex_number; i++ )
		drawNumber(number[i], vertex[i*2], vertex[i*2+1]);
}

/**
	draw number on each edge of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	edge			edges list
	\param	edge_number		number of edges
*/
template<typename TN, typename TV>
void drawNumberOnMeshEdges2D(const TN* number,
							 const TV* vertex,
							 const int* edge,
							 int edge_number){
	for( int i=0; i<edge_number; i++ ){
		int v1 = edge[i*2];
		int v2 = edge[i*2+1];
		drawNumber(number[i],
			(vertex[v1*2  ]+vertex[v2*2  ]) * 0.5,
			(vertex[v1*2+1]+vertex[v2*2+1]) * 0.5 );
	}
}

///@}

}}	//	end namespace yz::opengl


#endif	//	__YZ_OPENGL_UTILS_H__