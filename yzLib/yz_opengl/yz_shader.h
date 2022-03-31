/***********************************************************/
/**	\file
	\brief		Shader
	\details	must be included after <glew.h>
	\author		Yizhong Zhang
	\date		5/4/2013
*/
/***********************************************************/
#ifndef __YZ_SHADER_H__
#define __YZ_SHADER_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_glew_h
#	error yz_shader.h must be included after glew.h
#endif

namespace yz{	namespace opengl{

/**
	Shader manager

	To use shader, \n
	1,	call Setup or SetupFromFile to setup the shader	\n
	2,	call UseMe() before draw and UseDefault() after draw
*/
class Shader{
public:
	Shader(){
		Reset();
	}
	/**
		Setup shader by transfer string directly
	*/
	void Setup(const char* vert_str = NULL, const char* frag_str = NULL){
		if (glewIsSupported("GL_VERSION_3_3"))
			printf("Ready for OpenGL 3.3\n");
		else {
			printf("OpenGL 3.0 not supported\n");
			exit(1);
		}

		if( !vert_str && !frag_str )
			return;

		if( vert_str ){
			vert_shader = glCreateShader(GL_VERTEX_SHADER);
			glShaderSource(vert_shader, 1, &vert_str, NULL);
			glCompileShader(vert_shader);
			PrintShaderInfoLog(vert_shader);
		}

		if( frag_str ){
			frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
			glShaderSource(frag_shader, 1, &frag_str, NULL);
			glCompileShader(frag_shader);
			PrintShaderInfoLog(frag_shader);
		}

		prog = glCreateProgram();

		if( vert_str )
			glAttachShader(prog, vert_shader);
		if( frag_str )
			glAttachShader(prog, frag_shader);

		glBindFragDataLocation(prog, 0, "outputF");
		glLinkProgram(prog);
		PrintProgramInfoLog(prog);
	}

	/**
		Setup shader by read read shader string from file
	*/
	void SetupFromFile(const char* vert_file_name = NULL, const char* frag_file_name = NULL){
		if( !vert_file_name && !frag_file_name )
			return;

		std::string vert_str, frag_str;

		//	read vert to string
		std::ifstream vert_file(vert_file_name);
		if( vert_file.is_open() ){
			vert_file.seekg(0, std::ios::end);   
			vert_str.reserve(vert_file.tellg());
			vert_file.seekg(0, std::ios::beg);

			vert_str.assign((std::istreambuf_iterator<char>(vert_file)), std::istreambuf_iterator<char>());
		}

		//	read frag to string
		std::ifstream frag_file(frag_file_name);
		if( frag_file.is_open() ){
			frag_file.seekg(0, std::ios::end);   
			frag_str.reserve(frag_file.tellg());
			frag_file.seekg(0, std::ios::beg);

			frag_str.assign((std::istreambuf_iterator<char>(frag_file)), std::istreambuf_iterator<char>());
		}

		//	setup by string
		if( vert_str.size() || frag_str.size() ){
			Setup(vert_str.c_str(), frag_str.c_str());
		}

	}

	/**
		Use this shader
	*/
	void UseMe(){
		glUseProgram(prog);
	}

	/**
		Use default shader
	*/
	void UseDefault(){
		glUseProgram(0);
	}

	/**
		Get location of uniform variables in shader
	*/
	GLint GetUniformLocation(const char* uniform_name){
		GLint loc = glGetUniformLocation( prog, uniform_name ); 
		if( loc == -1 ){
			printf("cannot get location for uniform: %s", uniform_name);
		}

		return loc;
	}

	void Reset(){
		prog		= 0;
		vert_shader	= 0;
		frag_shader	= 0;
	}
public:
	GLuint prog, vert_shader, frag_shader;

private:
	/**
		Print error information for shader
	*/
	void PrintShaderInfoLog(GLuint obj){
		int infologLength = 0;
		int charsWritten  = 0;
		char *infoLog;

		glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

		if (infologLength > 0){
			infoLog = (char *)malloc(infologLength);
			glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
			printf("%s\n",infoLog);
			free(infoLog);
		}
	}

	/**
		Print error information for program
	*/
	void PrintProgramInfoLog(GLuint obj){
		int infologLength = 0;
		int charsWritten  = 0;
		char *infoLog;

		glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

		if (infologLength > 0){
			infoLog = (char *)malloc(infologLength);
			glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
			printf("%s\n",infoLog);
			free(infoLog);
		}
	}

};


}}	//	end namespace yz::opengl


#endif	//	__YZ_SHADER_H__
