/***********************************************************/
/**	\file
	\brief		FBO with Texture and Depth Buffer
	\details	must be included after <glew.h>
	\author		Yizhong Zhang
	\date		5/25/2012
*/
/***********************************************************/
#ifndef __YZ_FBO_H__
#define __YZ_FBO_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_glew_h
#	error yz_fbo.h must be included after glew.h
#endif

#include <iostream>

namespace yz{	namespace opengl{

/**
	FBO of OpenGL

	Frame Buffer Object is used to redirect display of 
	OpenGL from screen to FBO. To use FBO, you should:

	1,	create an FBO instance and call InitFBO to set dimension

	2,	call BeginRender() before you draw something and EndRender()
		after, then your display will be directed to this FBO
	
	When this FBO will never be used, call DestroyFBO to release resource
*/
class FBO{
public:
	int tex_width, tex_height;

	GLuint tex_id;
	GLuint fbo_id;
	GLuint rbo_id;	//	depth buffer

public:
	//	constructor and destructor
	FBO():tex_width(0), tex_height(0), tex_id(0), fbo_id(0), rbo_id(0), status(0){}
	~FBO(){
		Reset();
	}

	/**
		Initialize FBO according to given dimension

		This function must be called right after glutCreateWindow()
		of the window you want to render to FBO.

		\return		1: succeed;		0: failed
	*/
	inline int InitFBO(int width, int height){
		if( status & FBO_STATUS_BIT_INITIALIZED ){
			std::cout << "fbo : " << fbo_id << " already initialized" << std::endl;
			return 0;
		}

		tex_width	= width;
		tex_height	= height;

		//	get IDs
		glGenTextures( 1, &tex_id );
		glGenFramebuffers( 1, &fbo_id );
		glGenRenderbuffers( 1, &rbo_id );

		//	create and bind FBO
			//	fbo
		glBindFramebuffer( GL_FRAMEBUFFER, fbo_id );

			//	rbo
		glBindRenderbuffer( GL_RENDERBUFFER, rbo_id );
		glRenderbufferStorage( GL_RENDERBUFFER, GL_DEPTH_COMPONENT, tex_width, tex_height );
		glBindRenderbuffer( GL_RENDERBUFFER, 0 );

			//	texture
		glBindTexture( GL_TEXTURE_2D, tex_id );
		glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA32F, tex_width, tex_height, 0, GL_RGBA, GL_FLOAT, NULL );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture( GL_TEXTURE_2D, 0 );

			//	bind
		glFramebufferTexture2D( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_id, 0 );
		glFramebufferRenderbuffer( GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo_id);

			//	check status
		int fbo_code = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		if( fbo_code != GL_FRAMEBUFFER_COMPLETE ){
			std::cout << "FBO init incomplete, err code: " << std::hex << fbo_code << std::endl;
		}

		glBindFramebuffer( GL_FRAMEBUFFER, 0 );

		if( fbo_code != GL_FRAMEBUFFER_COMPLETE ){	//	initialize fbo failed
			Reset();
			return 0;
		}
		else{
			status |= FBO_STATUS_BIT_INITIALIZED;
			return 1;
		}
	}

	/**
		Begin Render to FBO, call before display

		You must be careful with context of OpenGL when rendering to FBO

		\snippet	glut_window.cpp	Render to FBO

		\return		1: succeed;		0: failed
	*/	
	inline int BeginRender(){
		if( ! (status & FBO_STATUS_BIT_INITIALIZED) ){
			std::cout << "fbo : " << fbo_id << " not initialized" << std::endl;
			return 0;
		}
		if( status & FBO_STATUS_BIT_BIND_FBO_ON ){
			std::cout << "fbo : " << fbo_id << " has started render, cannot start again" << std::endl;
			return 0;
		}

		glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

		status |= FBO_STATUS_BIT_BIND_FBO_ON;
		return 1;
	}

	/**
		End Render to FBO, call after display

		\return		1: succeed;		0: failed
	*/	
	inline int EndRender(){
		if( ! (status & FBO_STATUS_BIT_INITIALIZED) ){
			std::cout << "fbo : " << fbo_id << " not initialized" << std::endl;
			return 0;
		}
		if( ! (status & FBO_STATUS_BIT_BIND_FBO_ON) ){
			std::cout << "fbo : " << fbo_id << " is not start, cannot end render" << std::endl;
			return 0;
		}

		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		status &= ~FBO_STATUS_BIT_BIND_FBO_ON;
		return 1;
	}

	/**
		Destroy FBO, call when this FBO will never be used

		\return		1: succeed;		0: failed
	*/	
	inline int DestroyFBO(){
		if( status & FBO_STATUS_BIT_BIND_FBO_ON ){
			std::cout << "fbo : " << fbo_id << " is binded to buffer, cannot destroy" << std::endl;
			return 0;
		}

		Reset();

		return 1;
	}

	/**
		Print Help Information
	*/
	inline void Help(){
		std::cout
			<< "FBO usage: (OpenGL context must be ready)\n"
			<< "\t1, InitFBO to given dimension\n"
			<< "\t2, before and after display function, call BeginRender and EndRender\n"
			<< "\t   the display will be redirected to texture: tex_id\n"
			<< std::endl;

	}

protected:
	int	status;		//	show the current status of FBO
	enum{			//	control fbo status by bit
		FBO_STATUS_BIT_INITIALIZED	= 0x01,
		FBO_STATUS_BIT_BIND_FBO_ON	= 0x02
	};

protected:
	/**
		Reset the class
	*/
	inline void Reset(){
		if( tex_id )	glDeleteTextures( 1, &tex_id );
		if( fbo_id )	glDeleteFramebuffers( 1, &fbo_id );
		if( rbo_id )	glDeleteRenderbuffers( 1, &rbo_id );

		tex_width	= 0;
		tex_height	= 0;
		tex_id		= 0;
		fbo_id		= 0;
		rbo_id		= 0;
		status		= 0;
	}

};


}}	//	end namespace yz::opengl


#endif	//	__YZ_FBO_H__