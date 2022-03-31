/***********************************************************/
/**	\file
	\brief		VBO 
	\details	must be included after <glew.h>
	\author		Yizhong Zhang
	\date		5/4/2013
*/
/***********************************************************/
#ifndef __YZ_VBO_H__
#define __YZ_VBO_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_glew_h
#	error yz_vbo.h must be included after glew.h
#endif

namespace yz{	namespace opengl{

/**
	VBO that is used for rendering

	This class holds several vbos for a rendering scene, including
	vertex_buffer, normal_buffer, color_buffer, index_buffer

	default parameters are used for rendering smooth shading triangle mesh

	To use VBODisplayer, \n
	1,	set buffers, at least vertex buffern and index buffer should be set		\n
	2,	call Display, then object will be able to draw
*/
class VBODisplayer{
public:
	VBODisplayer(){
		Reset();
	}

	~VBODisplayer(){
		DeleteBuffers();
	}

	/**
		set vertex buffer

		\param	ptr				pointer to the vertex array
		\param	size_in_bytes	size of the array in bytes
		\param	coordinates_per_vertex		how many components in each vertex, typically 2, 3, 4, first param of glVertexPointer
		\param	type			type of the vertex array, GL_SHORT, GL_INT, GL_FLOAT, GL_DOUBLE, second param of glVertexPointer
		\param	stride			byte offset between consecutive vertices, third param of glVertexPointer
		\param	usage			type of the buffer, the last parameter of glBufferData
	*/
	void SetVertexBuffer(	void*	ptr, 
							int		size_in_bytes, 
							GLint  	coordinates_per_vertex	= 3,
							GLenum	type					= GL_FLOAT,
							GLsizei stride					= 0,
							GLenum	usage					= GL_STATIC_DRAW ){
		if( !vertex_vbo_id ){	//	buffer not created yet
			glGenBuffers(1, &vertex_vbo_id);
		}

		glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo_id);
		glBufferData(GL_ARRAY_BUFFER, size_in_bytes, ptr, usage);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		vertex_size		= coordinates_per_vertex;
		vertex_type		= type;
		vertex_stride	= stride;
	}

	/**
		set vertex normal buffer

		\param	ptr				pointer to the vertex normal array
		\param	size_in_bytes	size of the array in bytes
		\param	type			type of the vertex array, GL_SHORT, GL_INT, GL_FLOAT, GL_DOUBLE, first param of glNormalPointer
		\param	stride			byte offset between consecutive vertices, second param of glNormalPointer
		\param	usage			type of the buffer, the last parameter of glBufferData
	*/
	void SetVertexNormalBuffer(	void*	ptr, 
								int		size_in_bytes, 
								GLenum	type					= GL_FLOAT,
								GLsizei stride					= 0,
								GLenum	usage					= GL_STATIC_DRAW ){
		if( !normal_vbo_id ){	//	buffer not created yet
			glGenBuffers(1, &normal_vbo_id);
		}

		glBindBuffer(GL_ARRAY_BUFFER, normal_vbo_id);
		glBufferData(GL_ARRAY_BUFFER, size_in_bytes, ptr, usage);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		normal_type		= type;
		normal_stride	= stride;
	}

	/**
		set vertex color buffer

		\param	ptr				pointer to the vertex color array
		\param	size_in_bytes	size of the array in bytes
		\param	coordinates_per_vertex		how many components in each vertex, typically 2, 3, 4, first param of glColorPointer
		\param	type			type of the vertex array, GL_SHORT, GL_INT, GL_FLOAT, GL_DOUBLE, second param of glColorPointer
		\param	stride			byte offset between consecutive vertices, third param of glColorPointer
		\param	usage			type of the buffer, the last parameter of glBufferData
	*/
	void SetVertexColorBuffer(	void*	ptr, 
								int		size_in_bytes, 
								GLint  	coordinates_per_vertex	= 3,
								GLenum	type					= GL_FLOAT,
								GLsizei stride					= 0,
								GLenum	usage					= GL_STATIC_DRAW ){
		if( !color_vbo_id ){	//	buffer not created yet
			glGenBuffers(1, &color_vbo_id);
		}

		glBindBuffer(GL_ARRAY_BUFFER, color_vbo_id);
		glBufferData(GL_ARRAY_BUFFER, size_in_bytes, ptr, usage);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		color_size		= coordinates_per_vertex;
		color_type		= type;
		color_stride	= stride;
	}

	/**
		set element index buffer

		\param	ptr				pointer to the face array
		\param	size_in_bytes	size of the array in bytes
		\param	count			number of elements, count*sizeof(type) == size_in_bytes
		\param	type			type of the index array, GL_UNSIGNED_BYTE, GL_UNSIGNED_SHORT, or GL_UNSIGNED_INT, third param of glDrawElements
		\param	usage			type of the buffer, the last parameter of glBufferData
	*/
	void SetIndexBuffer(void*	ptr, 
						int		size_in_bytes, 
						GLsizei	count, 
						GLenum	type	= GL_UNSIGNED_INT,
						GLenum	usage	= GL_STATIC_DRAW ){
		if( !index_vbo_id ){	//	buffer not created yet
			glGenBuffers(1, &index_vbo_id);
		}

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo_id);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size_in_bytes, ptr, usage);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		index_count	= count;
		index_type	= type;
	}

	/**
		Draw elements with previously set vbos
	*/
	void Display(GLenum  draw_mode = GL_TRIANGLES){
		//	at least vertex and face should exist
		if( !vertex_vbo_id || !index_vbo_id )
			return;

		//	bind vertex
		glBindBuffer(GL_ARRAY_BUFFER, vertex_vbo_id);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(vertex_size, vertex_type, vertex_stride, NULL);

		//	bind vertex normal
		if( normal_vbo_id ){
			glBindBuffer(GL_ARRAY_BUFFER, normal_vbo_id);
			glEnableClientState(GL_NORMAL_ARRAY);
			glNormalPointer(normal_type, normal_stride, NULL); 
		}

		//	bind vertex color
		if( color_vbo_id ){
			glBindBuffer(GL_ARRAY_BUFFER, color_vbo_id);
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(color_size, color_type, color_stride, NULL); 
		}

		//	draw the mesh
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo_id);
		glDrawElements(draw_mode, index_count, index_type, 0);

		//	disable all
		if( color_vbo_id )
			glDisableClientState(GL_COLOR_ARRAY);
		if( normal_vbo_id )
			glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		//	unbind
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}


public:
	//	vertex
	GLuint	vertex_vbo_id;
	GLint  	vertex_size;
	GLenum	vertex_type;
	GLsizei vertex_stride;

	//	vertex normal
	GLuint	normal_vbo_id;
	GLenum	normal_type;
	GLsizei normal_stride;

	//	vertex color
	GLuint	color_vbo_id;
	GLint  	color_size;
	GLenum	color_type;
	GLsizei color_stride;

	//	element index
	GLuint	index_vbo_id;
	GLsizei	index_count;
	GLenum  index_type;

private:
	void Reset(){
		vertex_vbo_id	= 0;
		vertex_size		= 0;
		vertex_stride	= 0;

		normal_vbo_id	= 0;
		normal_stride	= 0;

		color_vbo_id	= 0;
		color_size		= 0;
		color_stride	= 0;

		index_vbo_id	= 0;
		index_count		= 0;
	}

	void DeleteBuffers(){
		if(vertex_vbo_id)		glDeleteBuffers(1, &vertex_vbo_id);
		if(normal_vbo_id)		glDeleteBuffers(1, &normal_vbo_id);
		if(color_vbo_id)		glDeleteBuffers(1, &color_vbo_id);
		if(index_vbo_id)		glDeleteBuffers(1, &index_vbo_id);
	}
};




}}	//	end namespace yz::opengl


#endif	//	__YZ_VBO_H__