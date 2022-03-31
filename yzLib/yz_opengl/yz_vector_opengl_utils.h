/***********************************************************/
/**	\file
	\brief		OpenGL Utilities of Vector
	\details	util functions taking vector as parameter to make 
				using OpenGL easier. 

				Must include after <glut.h>.
				yz_vector is included in this file, so if yz_opengl
				is included, yz_vector is automatically included.
				All the functions allow only yz::Vec, so you must 
				be careful if you have other vector implementation, 
				such as openCV
	\author		Yizhong Zhang
	\date		5/23/2012
*/
/***********************************************************/
#ifndef __YZ_VECTOR_OPENGL_UTILS_H__
#define __YZ_VECTOR_OPENGL_UTILS_H__

#include "yzLib/yz_setting.h"

#if !(	defined(YZ_glut_h) || defined(YZ_freeglut_h) )
#	error yz_vector_opengl_utils.h must be included after glut.h or freeglut.h
#endif

#include <vector>
#include <math.h>
#include "yzLib/yz_math/yz_vector.h"
#include "yzLib/yz_math/yz_interpolation.h"
#include "yzLib/yz_opengl/yz_opengl_utils.h"

namespace yz{	namespace opengl{

template<typename T1, typename T2>
inline void drawLineSegment(Com2<T1> v0, Com2<T2> v1);
template<typename T1, typename T2>
inline void drawLineSegment(Com3<T1> v0, Com3<T2> v1);
template<typename T1, typename T2>
inline void drawCylinder(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 4);

//	========================================
///@{
/**	@name Draw Point
*/
//	========================================
/**
	Draw a point on 2D canvas
*/
template<typename T>
inline void drawPoint(Com2<T> p){
	glBegin(GL_POINTS);
	glVertex2f(p.x, p.y);
	glEnd();
}

/**
	Draw a point in 3D space
*/
template<typename T>
inline void drawPoint(Com3<T> p){
	glBegin(GL_POINTS);
	glVertex3f(p.x, p.y, p.z);
	glEnd();
}

/**
	Draw a point with square in 2D space

	This function assumes the current matrix is model view

	\param	p		coordinate of point
	\param	size	length of the square edge
*/
template<typename T>
inline void drawPointAsSquare(Com2<T> p, float size){
	float half_size = size * 0.5;
	glBegin(GL_QUADS);
		glVertex2f(p.x - half_size, p.y - half_size);
		glVertex2f(p.x + half_size, p.y - half_size);
		glVertex2f(p.x + half_size, p.y + half_size);
		glVertex2f(p.x - half_size, p.y + half_size);
	glEnd();
}

/**
	Draw a point with circle in 2D space, not filled

	This function assumes the current matrix is model view

	\param	p			coordinate of point
	\param	radius		radius of the circle
	\param	slices		slices is the circle divided
*/
template<typename T>
inline void drawPointAsCircle(Com2<T> p, float radius, int slices = 16){
	drawCircle(p, radius, slices);
}

/**
	Draw a point with circle and a cross in 2D space, not filled

	This function assumes the current matrix is model view

	\param	p			coordinate of point
	\param	radius		radius of the circle
	\param	slices		slices is the circle divided
*/
template<typename T>
inline void drawPointAsCrossedCircle(Com2<T> p, float radius, int slices = 16){
	drawCircle(p, radius, slices);
	drawLineSegment(float2(p.x-radius, p.y), float2(p.x+radius, p.y));
	drawLineSegment(float2(p.x, p.y-radius), float2(p.x, p.y+radius));
}

/**
	Draw a point with ball in 2D space

	This function assumes the current matrix is model view

	\param	p			coordinate of point
	\param	radius		radius of the circle
	\param	slices		slices is the circle divided
*/
template<typename T>
inline void drawPointAsBall(Com2<T> p, float radius, int slices = 16){
	glBegin(GL_TRIANGLE_FAN);
		glVertex2f(p.x, p.y);
		for(int i=0; i<=slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex2f(p.x + r.x, p.y + r.y);
		}
	glEnd();
}

/**
	Draw a point with cube in 3D space

	This function assumes the current matrix is model view

	\param	p		coordinate of point
	\param	size	length of cube edge
*/
template<typename T>
inline void drawPointAsCube(Com3<T> p, float size){
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glPushMatrix();
	glTranslatef(p.x, p.y, p.z);
	glutSolidCube(size);
	glPopMatrix();
}

/**
	Draw a point with wire cube in 3D space

	This function assumes the current matrix is model view

	\param	p		coordinate of point
	\param	size	length of cube edge
*/
template<typename T>
inline void drawPointAsWireCube(Com3<T> p, float size){
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glPushMatrix();
	glTranslatef(p.x, p.y, p.z);
	glutWireCube(size);
	glPopMatrix();
}

/**
	Draw a point with sphere in 3D space

	This function assumes the current matrix mode is model view

	\param	p		coordinate of point
	\param	size	sphere diameter
	\param	slices	slices of the sphere
	\param	stacks	stacks of the sphere
*/
template<typename T>
inline void drawPointAsSphere(Com3<T> p, float size, int slices=10, int stacks=10){
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glPushMatrix();
	glTranslatef(p.x, p.y, p.z);
	glutSolidSphere(size, slices, stacks);
	glPopMatrix();
}

/**
	Draw a point with wire sphere in 3D space

	This function assumes the current matrix mode is model view

	\param	p		coordinate of point
	\param	size	sphere diameter
	\param	slices	slices of the sphere
	\param	stacks	stacks of the sphere
*/
template<typename T>
inline void drawPointAsWireSphere(Com3<T> p, float size, int slices=10, int stacks=10){
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glPushMatrix();
	glTranslatef(p.x, p.y, p.z);
	glutWireSphere(size, slices, stacks);
	glPopMatrix();
}

/**
	Draw number at a given point in 3D space

	\param	value	the number to draw
	\param	xyz		position of the point
*/
template<typename T>
inline void drawNumber(T value, float3 xyz){
	drawNumber(value, xyz.x, xyz.y, xyz.z);
}

/**
	Draw number at a given point on 2D plane

	\param	value	the number to draw
	\param	xy		position of the point
*/
template<typename T>
inline void drawNumber(T value, float2 xy){
	drawNumber(value, xy.x, xy.y);
}

///@}

//	========================================
///@{
/**	@name Draw Vector
*/
//	========================================
/**
	Draw vector on 2D canvas

	\param	start	start point of the vector
	\param	r		the vector
*/
template<typename T1, typename T2>
inline void drawVector(Com2<T1> start, Com2<T2> r){
	float2 start_f	= start;
	float2 end_f	= Vec2f(start) + Vec2f(r);
	drawLineSegment(start_f, end_f);
}

/**
	Draw vector in 3D space

	\param	start	start point of the vector
	\param	r		the vector
*/
template<typename T1, typename T2>
inline void drawVector(Com3<T1> start, Com3<T2> r){
	float3 start_f	= start;
	float3 end_f	= Vec3f(start) + Vec3f(r);
	drawLineSegment(start_f, end_f);
}

///@}

//	========================================
///@{
/**	@name Draw Curve
*/
//	========================================

/**
	Draw line Segment on 2D canvas

	\param	v0		coordinate of end point v0
	\param	v1		coordinate of end point v1
*/
template<typename T1, typename T2>
inline void drawLineSegment(Com2<T1> v0, Com2<T2> v1){
	float2 start_f	= v0;
	float2 end_f	= v1;
	glBegin(GL_LINES);
		glVertex2f( start_f.x, start_f.y );
		glVertex2f( end_f.x, end_f.y );
	glEnd();
}

/**
	Draw line segment on in 3D space

	\param	v0		coordinate of end point v0
	\param	v1		coordinate of end point v1
*/
template<typename T1, typename T2>
inline void drawLineSegment(Com3<T1> v0, Com3<T2> v1){
	float3 start_f	= v0;
	float3 end_f	= v1;
	glBegin(GL_LINES);
		glVertex3f( start_f.x, start_f.y, start_f.z );
		glVertex3f( end_f.x, end_f.y, end_f.z );
	glEnd();
}

/**
	Draw ray on 2D canvas, given origin and direction

	\param	origin		coordinate of origin point
	\param	next		direction = next - origin
	\param	scale		how many times longer is the ray displayed than direction length
*/
template<typename T1, typename T2>
inline void drawRay(Com2<T1> origin, Com2<T2> next, float scale=100){
	Vec2f direction = Vec2f(next) - Vec2f(origin);
	Vec2f end		= Vec2f(origin) + direction * scale;
	glBegin(GL_LINES);
		glVertex2f( origin.x, origin.y );
		glVertex2f( end.x, end.y );
	glEnd();
}

/**
	Draw ray in 3D space, given origin and direction

	\param	origin		coordinate of origin point
	\param	next		direction = next - origin
	\param	scale		how many times longer is the ray displayed than direction length
*/
template<typename T1, typename T2>
inline void drawRay(Com3<T1> origin, Com3<T2> next, float scale=100){
	Vec3f direction = Vec3f(next) - Vec3f(origin);
	Vec3f end		= Vec3f(origin) + direction * scale;
	glBegin(GL_LINES);
		glVertex3f( origin.x, origin.y, origin.z );
		glVertex3f( end.x, end.y, end.z );
	glEnd();
}

/**
	Draw line on 2D canvas, given two points on the line

	\param	v0		coordinate of point v0
	\param	v1		coordinate of point v1
	\param	scale	how many times longer is the line displayed than v0v1
*/
template<typename T1, typename T2>
inline void drawLine(Com2<T1> v0, Com2<T2> v1, float scale=100){
	Vec2f r = Vec2f(v1) - Vec2f(v0);
	Vec2f start = Vec2f(v0) - r * scale;
	Vec2f end	= Vec2f(v1) + r * scale;
	glBegin(GL_LINES);
		glVertex2f( start.x, start.y );
		glVertex2f( end.x, end.y );
	glEnd();
}

/**
	Draw line in 3D space, given two points on the line

	\param	v0		coordinate of point v0
	\param	v1		coordinate of point v1
	\param	scale	how many times longer is the line displayed than v0v1
*/
template<typename T1, typename T2>
inline void drawLine(Com3<T1> v0, Com3<T2> v1, float scale=100){
	Vec3f r = Vec3f(v1) - Vec3f(v0);
	Vec3f start = Vec3f(v0) - r * scale;
	Vec3f end	= Vec3f(v1) + r * scale;
	glBegin(GL_LINES);
		glVertex3f( start.x, start.y, start.z );
		glVertex3f( end.x, end.y, end.z );
	glEnd();
}

/**
	Draw ray on 2D canvas, given origin and direction

	\param	origin			coordinate of origin point of the arrow
	\param	end				the arrow is from origin to end
	\param	length_coef		length of the arrow part with respect of the whole arrow
	\param	width_coef		width of the arrow part with respect of the whole arrow
*/
template<typename T1, typename T2>
inline void drawArrow(Com2<T1> origin, Com2<T2> end, float length_coef = 0.25f, float width_coef = 0.1f){
	Vec2f dir = Vec2f(end) - Vec2f(origin);
	Vec2f p = Vec2f(origin) + dir * (1.0f - length_coef);
	dir.SetRotateRad(YZ_PI*0.5f);
	dir *= width_coef;

	glBegin(GL_LINES);
		glVertex2f( origin.x, origin.y );
		glVertex2f( end.x, end.y );
		glVertex2f( end.x, end.y );
		glVertex2f( p.x+dir.x, p.y+dir.y );
		glVertex2f( end.x, end.y );
		glVertex2f( p.x-dir.x, p.y-dir.y );
	glEnd();
}

/**
	draw 3D arrow

	\param	origin			coordinate of origin point of the arrow
	\param	end				the arrow is from origin to end
	\param	radius			radius of this arrow to draw
*/
template<typename T1, typename T2>
void drawArrow(Com3<T1> origin, Com3<T2> end, float radius){
	Vec3f z = Vec3f(end) - Vec3f(origin);

	drawCylinder(origin, Vec3f(origin) + z * 0.80f, radius);

	int max_dim = 0;
	float max_abs_val = fabs( z[max_dim] );
	for(int i=1; i<3; i++){
		if( fabs(z[i]) > max_abs_val ){
			max_dim = i;
			max_abs_val = fabs(z[i]);
		}
	}

	Vec3f x = yz::cross( z, (max_dim==0 ? Vec3f(0,1,0) : Vec3f(1,0,0) ) ).Normalize();
	Vec3f y = yz::cross( z, x ).Normalize();
	Matrix4x4f mat( Vec4f(x), Vec4f(y), Vec4f(z)*0.33f, Vec4f(Vec3f(origin) + z * 0.66f, 1) );
	mat.SetTranspose();

	glPushMatrix();
	glMultMatrixf(mat.data[0]);
	glutSolidCone(radius*2, 1, 16, 2);
	//glutSolidCone(radius*2, 1, 4, 2);
	glPopMatrix();
}


/**
	Draw a long curve on 2D canvas

	\param	curve		the curve to be displayed
*/
template<typename T>
inline void drawCurve(const std::vector<Vec2<T>>& curve){
	if( curve.empty() )
		return;
	drawCurve2D((T*)&curve[0], curve.size());
}

/**
	Draw a long curve in 3D space

	\param	curve		the curve to be displayed
*/
template<typename T>
inline void drawCurve(const std::vector<Vec3<T>>& curve){
	if( curve.empty() )
		return;
	drawCurve3D((T*)&curve[0], curve.size());
}

/**
	Draw cubic bezier curve on 2D canvas

	\param	v0		position of point 0
	\param	v1		position of point 1
	\param	v2		position of point 2
	\param	v3		position of point 3
	\param	slices	how many slices is the curve to be discretized
*/
template<typename T>
inline void drawCubicBezier(Vec2<T> v0, Vec2<T> v1, Vec2<T> v2, Vec2<T> v3, int slices = 16){
	std::vector<Vec2<T>> point;
	point.resize(slices+1);
	point[0] = v0;
	point[slices] = v3;
	for(int i=1; i<slices; i++){	
		point[i] = interpCubicBezier(v0, v1, v2, v3, double(i)/slices);
	}
	drawCurve(point);
}

/**
	Draw cubic bezier curve in 3D space

	\param	v0		position of point 0
	\param	v1		position of point 1
	\param	v2		position of point 2
	\param	v3		position of point 3
	\param	slices	how many slices is the curve to be discretized
*/
template<typename T>
inline void drawCubicBezier(Vec3<T> v0, Vec3<T> v1, Vec3<T> v2, Vec3<T> v3, int slices = 16){
	std::vector<Vec3<T>> point;
	point.resize(slices+1);
	point[0] = v0;
	point[slices] = v3;
	for(int i=1; i<slices; i++){	
		point[i] = interpCubicBezier(v0, v1, v2, v3, double(i)/slices);
	}
	drawCurve(point);
}
/**
	Draw a polygon on 2D canvas

	\param	point		the polygon points to be displayed
*/
template<typename T>
inline void drawPolygon(const std::vector<Vec2<T>>& point){
	if( point.empty() )
		return;
	drawPolygon2D((T*)&point[0], point.size());
}

/**
	Draw a polygon in 3D space

	\param	point		the polygon points to be displayed
*/
template<typename T>
inline void drawPolygon(const std::vector<Vec3<T>>& point){
	if( point.empty() )
		return;
	drawPolygon3D((T*)&point[0], point.size());
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

	\param	min_coord	min coordinate on canvas
	\param	max_coord	max coordinate on canvas
	\param	flip_flag	whether vertical flip the texture, 1: flip, 0: don't flip
*/
template<typename T1, typename T2>
inline void drawWholeTexture(Com2<T1> min_coord, Com2<T1> max_coord, int flip_flag=0){
	if( flip_flag ){
		glBegin( GL_QUADS );
			glTexCoord2f(0, 0);	glVertex2f(min_coord.x, max_coord.y);
			glTexCoord2f(1, 0);	glVertex2f(max_coord.x, max_coord.y);
			glTexCoord2f(1, 1);	glVertex2f(max_coord.x, min_coord.y);
			glTexCoord2f(0, 1);	glVertex2f(min_coord.x, min_coord.y);
		glEnd();
	}
	else{
		glBegin( GL_QUADS );
			glTexCoord2f(0, 0);	glVertex2f(min_coord.x, min_coord.y);
			glTexCoord2f(1, 0);	glVertex2f(max_coord.x, min_coord.y);
			glTexCoord2f(1, 1);	glVertex2f(max_coord.x, max_coord.y);
			glTexCoord2f(0, 1);	glVertex2f(min_coord.x, max_coord.y);
		glEnd();
	}
}

//	========================================
///@{
/**	@name Draw Certain Shape
*/
//	========================================
/**
	Draw a circle on 2D canvas

	\param	center		center coordinate of circle
	\param	radius		radius of the circle
	\param	slices		slices is the circle divided
*/
template<typename T>
inline void drawCircle(Com2<T> center, float radius, int slices = 32){
	glBegin(GL_LINE_LOOP);
		for( int i=0; i<slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex2f(center.x+r.x, center.y+r.y);
		}
	glEnd();
}

/**
	Draw a filled circle on 2D canvas

	\param	center		center coordinate of circle
	\param	radius		radius of the circle
	\param	slices		slices is the circle divided
*/
template<typename T>
inline void drawCircleFilled(Com2<T> center, float radius, int slices = 32){
	glBegin(GL_POLYGON);
		for( int i=0; i<slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex2f(center.x+r.x, center.y+r.y);
		}
	glEnd();
}
/**
	we draw a plane as a quad, orientation is not controlled

	\param	center		center of the quad
	\param	normal		normal of the quad
	\param	size		edge length of the quad
*/
template<typename T>
inline void drawPlane(Vec3<T> center, Vec3<T> normal, float size){
	normal.SetNormalize();

	Vec3f perp_dir = (normal.x>0.8f || normal.x<-0.8f) ? cross(normal, Vec3<T>(0, 1, 0)) : cross(normal, Vec3<T>(1, 0, 0));
	perp_dir.SetNormalize();
	perp_dir *= size * 0.707107f;
	Vec3f perp_dir2 = perp_dir.RotateDeg(normal, 90);

	Vec3f v[4] = {center+perp_dir, center+perp_dir2, center-perp_dir, center-perp_dir2};

	glBegin(GL_TRIANGLES);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(v[0].x, v[0].y, v[0].z);
		glVertex3f(v[1].x, v[1].y, v[1].z);
		glVertex3f(v[2].x, v[2].y, v[2].z);
		glVertex3f(v[2].x, v[2].y, v[2].z);
		glVertex3f(v[3].x, v[3].y, v[3].z);
		glVertex3f(v[0].x, v[0].y, v[0].z);
	glEnd();
}

/**
	we draw a plane as a quad, three axis are given, length of each edge are given by axis

	\param	center		center of the quad
	\param	normal		normal of the quad
	\param	axis_x		direction of axis x, length of vector is half the length of the edge
	\param	axis_y		direction of axis y, length of vector is half the length of the edge
*/
template<typename T>
inline void drawPlane(Vec3<T> center, Vec3<T> normal, Vec3<T> axis_x, Vec3<T> axis_y){
	normal.SetNormalize();

	if( dot( cross(axis_x, axis_y), normal ) < 0 )
		axis_y *= -1;

	Vec3f v[4] = {center+axis_x+axis_y, center-axis_x+axis_y, center-axis_x-axis_y, center+axis_x-axis_y};

	glBegin(GL_TRIANGLES);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(v[0].x, v[0].y, v[0].z);
		glVertex3f(v[1].x, v[1].y, v[1].z);
		glVertex3f(v[2].x, v[2].y, v[2].z);
		glVertex3f(v[2].x, v[2].y, v[2].z);
		glVertex3f(v[3].x, v[3].y, v[3].z);
		glVertex3f(v[0].x, v[0].y, v[0].z);
	glEnd();
}

/**
	we draw a quad as wire frames in x-y plane

	\param	size		edge length of the quad
	\param	slices		how many slices of the wire frame
*/
inline void drawPlaneXYWire(float size, int slices = 10){
	float offset = size * 0.5f;
	glBegin(GL_LINES);
		glNormal3f(0, 0, 1);
		for(int i=0; i<=slices; i++){
			glVertex3f(-offset, -offset+size*i/slices, 0);
			glVertex3f(offset, -offset+size*i/slices, 0);
		}
		for(int i=0; i<=slices; i++){
			glVertex3f(-offset+size*i/slices, -offset, 0);
			glVertex3f(-offset+size*i/slices, offset, 0);
		}
	glEnd();
}

/**
	we draw a plane as wire frames, orientation is not controlled

	\param	center		center of the quad
	\param	normal		normal of the quad
	\param	size		edge length of the quad
	\param	slices		how many slices do we want to draw
*/
template<typename T>
inline void drawPlaneWire(Vec3<T> center, Vec3<T> normal, float size, int slices = 10){
	normal.SetNormalize();

	Vec3f perp_dir = (normal.x>0.8f || normal.x<-0.8f) ? cross(normal, Vec3<T>(0, 1, 0)) : cross(normal, Vec3<T>(1, 0, 0));
	perp_dir.SetNormalize();
	Vec3f perp_dir2 = perp_dir.RotateDeg(normal, 90);

	//	calculate the transform
	float trans[16] = {
		perp_dir.x, perp_dir.y, perp_dir.z, 0, 
		perp_dir2.x, perp_dir2.y, perp_dir2.z, 0,
		normal.x, normal.y, normal.z, 0,
		center.x, center.y, center.z, 1
	};

	//	draw the wire plane
	glPushMatrix();
	glMultMatrixf(trans);
	drawPlaneXYWire(size, slices);
	glPopMatrix();
}


/**
	we draw a plane as wire frames, three axis are given, length of each edge are given by axis

	\param	center		center of the quad
	\param	normal		normal of the quad
	\param	axis_x		direction of axis x, length of vector is half the length of the edge
	\param	axis_y		direction of axis y, length of vector is half the length of the edge
	\param	slices		how many slices do we want to draw
*/
template<typename T>
inline void drawPlaneWire(Vec3<T> center, Vec3<T> normal, Vec3<T> axis_x, Vec3<T> axis_y, int slices = 10){
	normal.SetNormalize();

	if( dot( cross(axis_x, axis_y), normal ) < 0 )
		axis_y *= -1;

	//	calculate the transform
	float trans[16] = {
		axis_x.x*2, axis_x.y*2, axis_x.z*2, 0, 
		axis_y.x*2, axis_y.y*2, axis_y.z*2, 0,
		normal.x, normal.y, normal.z, 0,
		center.x, center.y, center.z, 1
	};

	//	draw the wire plane
	glPushMatrix();
	glMultMatrixf(trans);
	drawPlaneXYWire(1, slices);
	glPopMatrix();
}


/**
	Draw a cylinder in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinder(Com3<T1> v0, Com3<T2> v1, float radius, int slices){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	draw the faces
	glBegin(GL_QUADS);
		float h = ry_length;
		Vec2f r_old(radius, 0.0f);
		Vec2f nor_old = r_old.Normalize();
		for( int i=1; i<=slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor_old.x, 0, nor_old.y);
			glVertex3f(r_old.x, 0, r_old.y);
			glVertex3f(r_old.x, h, r_old.y);
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
			glVertex3f(r.x, 0, r.y);
			r_old	= r;
			nor_old	= nor;
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw a cylinder represented by wire in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinderWire(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 16){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	draw generatrix
	glBegin(GL_LINES);
		float h = ry_length;
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	//	draw two rings
	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw a closed cylinder (with two ends) in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinderClosed(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 16){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	the cylinder
	glBegin(GL_QUADS);
	float h = ry_length;
	Vec2f r_old(radius, 0.0f);
	Vec2f nor_old = r_old.Normalize();
	for( int i=1; i<=slices; i++ ){
		Vec2f r(radius, 0);
		r.SetRotateRad(i*2*YZ_PI/slices);
		Vec2f nor = r.Normalize();
		glNormal3f(nor_old.x, 0, nor_old.y);
		glVertex3f(r_old.x, 0, r_old.y);
		glVertex3f(r_old.x, h, r_old.y);
		glNormal3f(nor.x, 0, nor.y);
		glVertex3f(r.x, h, r.y);
		glVertex3f(r.x, 0, r.y);
		r_old	= r;
		nor_old	= nor;
	}
	glEnd();

	//	bottom circle
	glBegin(GL_TRIANGLE_FAN);
		glNormal3f(0, -1, 0);
		glVertex3f(0, 0, 0);
		for( int i=0; i<=slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	//	top circle
	glBegin(GL_TRIANGLE_FAN);
		glNormal3f(0, 1, 0);
		glVertex3f(0, h, 0);
		for( int i=0; i<=slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(-i*2*YZ_PI/slices);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw a closed cylinder (with two ends) represented by wire in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinderClosedWire(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 16){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	draw generatrix
	glBegin(GL_LINES);
		float h = ry_length;
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	//	draw two rings
	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	//	draw two ends
	glBegin(GL_LINES);
		glNormal3f(0, -1, 0);
		for( int i=0; i<slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex3f(0, 0, 0);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	glBegin(GL_LINES);
		glNormal3f(0, 1, 0);
		for( int i=0; i<slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex3f(0, h, 0);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw a half closed cylinder (with two ends) in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1, the closed end
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinderHalfClosed(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 16){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	the cylinder
	glBegin(GL_QUADS);
	float h = ry_length;
	Vec2f r_old(radius, 0.0f);
	Vec2f nor_old = r_old.Normalize();
	for( int i=1; i<=slices; i++ ){
		Vec2f r(radius, 0);
		r.SetRotateRad(i*2*YZ_PI/slices);
		Vec2f nor = r.Normalize();
		glNormal3f(nor_old.x, 0, nor_old.y);
		glVertex3f(r_old.x, 0, r_old.y);
		glVertex3f(r_old.x, h, r_old.y);
		glNormal3f(nor.x, 0, nor.y);
		glVertex3f(r.x, h, r.y);
		glVertex3f(r.x, 0, r.y);
		r_old	= r;
		nor_old	= nor;
	}
	glEnd();

	//	top circle
	glBegin(GL_TRIANGLE_FAN);
		glNormal3f(0, 1, 0);
		glVertex3f(0, h, 0);
		for( int i=0; i<=slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(-i*2*YZ_PI/slices);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw a half closed cylinder (with two ends) represented by wire in 3D space

	This function assumes the current matrix mode is model view

	\param	v0			coordinate of end point v0
	\param	v1			coordinate of end point v1, the closed end
	\param	radius		radius of the cylinder
	\param	slices		slices is the cylinder divided
*/
template<typename T1, typename T2>
inline void drawCylinderHalfClosedWire(Com3<T1> v0, Com3<T2> v1, float radius, int slices = 16){
	//	if the cylinder is too short or too thin, draw nothing
	Vec3f ry(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);	
	float ry_length = ry.Length();
	if( ry_length < 1e-5f || radius < 1e-5f )
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );
	glPushMatrix();

	//	find a axis that is farthest to r
	int axis_id = 0;
	if( fabs(ry[1]) < fabs(ry[axis_id]) )	axis_id = 1;
	if( fabs(ry[2]) < fabs(ry[axis_id]) )	axis_id = 2;

	//	find other two axis perpendicular to ry
	Vec3f tmp_r(0.0f, 0.0f, 0.0f);
	tmp_r[axis_id] = 1.0f;
	Vec3f rx = cross(ry, tmp_r).Normalize();
	Vec3f rz = cross(rx, ry).Normalize();

	//	calculate transform matrix, after = M * before
	Matrix4x4f before(Vec4f(v0, 1.0f), Vec4f(Vec3f(v0)+rx, 1.0f), Vec4f(v1, 1.0f), Vec4f(Vec3f(v0)+rz, 1.0f));
	Matrix4x4f after(Vec4f(0.0f, 0.0f, 0.0f, 1.0f), 
		Vec4f(1.0f, 0.0f, 0.0f, 1.0f), Vec4f(0.0f, ry_length, 0.0f, 1.0f), Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
	after = (before * after.Inverse()).Transpose();
	glMultMatrixf(after.data[0]);

	//	draw generatrix
	glBegin(GL_LINES);
		float h = ry_length;
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	//	draw two rings
	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, 0, r.y);
		}
	glEnd();

	glBegin(GL_LINE_LOOP);
		for(int i=0; i<slices; i++){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			Vec2f nor = r.Normalize();
			glNormal3f(nor.x, 0, nor.y);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	//	draw the closed end
	glBegin(GL_LINES);
		glNormal3f(0, 1, 0);
		for( int i=0; i<slices; i++ ){
			Vec2f r(radius, 0);
			r.SetRotateRad(i*2*YZ_PI/slices);
			glVertex3f(0, h, 0);
			glVertex3f(r.x, h, r.y);
		}
	glEnd();

	glPopMatrix();
}

/**
	Draw Axis Aligned Bonding Box with wires 2D

	given min and max coordinate of a bonding box, draw it
*/
template<typename T>
inline void drawAABBWire(Com2<T> bb_min, Com2<T> bb_max){
	Com2<T>	v1(bb_max.x, bb_min.y);
	Com2<T> v2(bb_min.x, bb_max.y);

	drawLineSegment(bb_min, v1);
	drawLineSegment(v1, bb_max);
	drawLineSegment(bb_max, v2);
	drawLineSegment(v2, bb_min);
}

/**
	Draw Axis Aligned Bonding Box with wires 3D

	given min and max coordinate of a bonding box, draw it
*/
template<typename T>
inline void drawAABBWire(Com3<T> bb_min, Com3<T> bb_max){
	Com3<T> v[8] = {
		bb_min,
		Com3<T>(bb_max.x, bb_min.y, bb_min.z),
		Com3<T>(bb_max.x, bb_max.y, bb_min.z),
		Com3<T>(bb_min.x, bb_max.y, bb_min.z),
		Com3<T>(bb_min.x, bb_min.y, bb_max.z),
		Com3<T>(bb_max.x, bb_min.y, bb_max.z),
		bb_max,
		Com3<T>(bb_min.x, bb_max.y, bb_max.z) };

	drawLineSegment(v[0], v[1]);
	drawLineSegment(v[4], v[5]);
	drawLineSegment(v[7], v[6]);
	drawLineSegment(v[3], v[2]);

	drawLineSegment(v[0], v[3]);
	drawLineSegment(v[1], v[2]);
	drawLineSegment(v[5], v[6]);
	drawLineSegment(v[4], v[7]);

	drawLineSegment(v[0], v[4]);
	drawLineSegment(v[1], v[5]);
	drawLineSegment(v[2], v[6]);
	drawLineSegment(v[3], v[7]);
}

/**
Draw Axis Aligned Bonding Box

given min and max coordinate of a bonding box, draw it
*/
template<typename T>
inline void drawAABB(Com3<T> bb_min, Com3<T> bb_max) {
	Com3<T> v[8] = {
		bb_min,
		Com3<T>(bb_max.x, bb_min.y, bb_min.z),
		Com3<T>(bb_max.x, bb_max.y, bb_min.z),
		Com3<T>(bb_min.x, bb_max.y, bb_min.z),
		Com3<T>(bb_min.x, bb_min.y, bb_max.z),
		Com3<T>(bb_max.x, bb_min.y, bb_max.z),
		bb_max,
		Com3<T>(bb_min.x, bb_max.y, bb_max.z) };

	glBegin(GL_QUADS);
	glNormal3f(1, 0, 0);
	glVertex3f(v[1].x, v[1].y, v[1].z);
	glVertex3f(v[2].x, v[2].y, v[2].z);
	glVertex3f(v[6].x, v[6].y, v[6].z);
	glVertex3f(v[5].x, v[5].y, v[5].z);
	glNormal3f(-1, 0, 0);
	glVertex3f(v[0].x, v[0].y, v[0].z);
	glVertex3f(v[4].x, v[4].y, v[4].z);
	glVertex3f(v[7].x, v[7].y, v[7].z);
	glVertex3f(v[3].x, v[3].y, v[3].z);
	glNormal3f(0, 1, 0);
	glVertex3f(v[2].x, v[2].y, v[2].z);
	glVertex3f(v[3].x, v[3].y, v[3].z);
	glVertex3f(v[7].x, v[7].y, v[7].z);
	glVertex3f(v[6].x, v[6].y, v[6].z);
	glNormal3f(0, -1, 0);
	glVertex3f(v[0].x, v[0].y, v[0].z);
	glVertex3f(v[1].x, v[1].y, v[1].z);
	glVertex3f(v[5].x, v[5].y, v[5].z);
	glVertex3f(v[4].x, v[4].y, v[4].z);
	glNormal3f(0, 0, 1);
	glVertex3f(v[4].x, v[4].y, v[4].z);
	glVertex3f(v[5].x, v[5].y, v[5].z);
	glVertex3f(v[6].x, v[6].y, v[6].z);
	glVertex3f(v[7].x, v[7].y, v[7].z);
	glNormal3f(0, 0, -1);
	glVertex3f(v[0].x, v[0].y, v[0].z);
	glVertex3f(v[3].x, v[3].y, v[3].z);
	glVertex3f(v[2].x, v[2].y, v[2].z);
	glVertex3f(v[1].x, v[1].y, v[1].z);
	glEnd();
}

/**
	Draw camera at with parametes

	This function assumes the current matrix mode is model view

	\param	eye		position of camera
	\param	center	look center of camera
	\param	up		up direction of camera
	\param	z_near		near plane of camera
	\param	z_far		far plane of camera
	\param	fov_x		field of view of x direction of camera
	\param	fov_y		field of view of y direction of camera
*/
template<typename T1, typename T2, typename T3>
inline void drawCamera(Com3<T1> eye, 
					   Com3<T2> center, 
					   Com3<T3> up, 
					   float z_near=0, 
					   float z_far=1, 
					   float fov_x=64, 
					   float fov_y=48){
	//	same as OpenGL gluLookAt parameter
	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glPushMatrix();

	Matrix4x4f mat;
	mat.SetLookAt(eye, center, up);
	mat.SetInverse();
	mat.SetTranspose();
	glMultMatrixf(mat.data[0]);

	drawCamera(-fabs(z_near), -fabs(z_far), fov_x, fov_y);

	glPopMatrix();
}

/**
	draw a flat solid bone with given parameter

	This function assumes the current matrix mode is model view

	\param	bone_origin		origin of the bone
	\param	bone_vector		bone direction
	\param	bone_front		front direction of the bone
*/
template<typename T1, typename T2, typename T3>
inline void drawFlatSolidBone(Vec3<T1> bone_origin, Vec3<T2> bone_vector, Vec3<T3> bone_front){
	float length = bone_vector.Length();
	if( length < 1e-6 )	//	an zero length bone
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glEnable(GL_NORMALIZE);
	glPushMatrix();

	glTranslatef(bone_origin.x, bone_origin.y, bone_origin.z);

	Vec3f		bone_vec = bone_vector.Normalize();
	Vec3f		bone_fro = bone_front.Normalize();
	Vec3f		normal	= cross(bone_vec, bone_fro);
	Matrix3x3f	trans_mat(normal, bone_vec, bone_fro);
	Matrix4x4f	gl_trans_mat = trans_mat.SetTranspose();

	glMultMatrixf(gl_trans_mat.data[0]);
	glScalef(length, length, length);

	drawFlatSolidBone();

	glPopMatrix();
}

/**
	draw a color cube bone with given parameter

	This function assumes the current matrix mode is model view

	\param	bone_origin		origin of the bone
	\param	bone_vector		bone direction
	\param	bone_front		front direction of the bone
	\param	width_scale		width of bone, with respect to bone of length 1, default value is 0.3
*/
template<typename T1, typename T2, typename T3>
inline void drawColorCubeBone(Vec3<T1> bone_origin, Vec3<T2> bone_vector, Vec3<T3> bone_front, float width_scale = 0.3f){
	float length = bone_vector.Length();
	if( length < 1e-6 )	//	an zero length bone
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glEnable(GL_NORMALIZE);
	glPushMatrix();

	glTranslatef(bone_origin.x, bone_origin.y, bone_origin.z);

	Vec3f		bone_vec = bone_vector.Normalize();
	Vec3f		bone_fro = bone_front.Normalize();
	Vec3f		normal	= cross(bone_vec, bone_fro);
	Matrix3x3f	trans_mat(normal, bone_vec, bone_fro);
	Matrix4x4f	gl_trans_mat = trans_mat.SetTranspose();

	glMultMatrixf(gl_trans_mat.data[0]);
	glScalef(length, length, length);

	drawColorCubeBone(width_scale);

	glPopMatrix();
}

/**
	draw a solid bone with given parameter

	This function assumes the current matrix mode is model view

	\param	bone_origin		origin of the bone
	\param	bone_vector		bone direction
	\param	bone_front		front direction of the bone
*/
template<typename T1, typename T2, typename T3>
inline void drawSolidBone(Vec3<T1> bone_origin, Vec3<T2> bone_vector, Vec3<T3> bone_front){
	float length = bone_vector.Length();
	if( length < 1e-6 )	//	an zero length bone
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glEnable(GL_NORMALIZE);
	glPushMatrix();

	glTranslatef(bone_origin.x, bone_origin.y, bone_origin.z);

	Vec3f		bone_vec = bone_vector.Normalize();
	Vec3f		bone_fro = bone_front.Normalize();
	Vec3f		normal	= cross(bone_vec, bone_fro);
	Matrix3x3f	trans_mat(normal, bone_vec, bone_fro);
	Matrix4x4f	gl_trans_mat = trans_mat.SetTranspose();

	glMultMatrixf(gl_trans_mat.data[0]);
	glScalef(length, length, length);

	drawSolidBone();

	glPopMatrix();
}

/**
	draw a wire bone with given parameter

	This function assumes the current matrix mode is model view

	\param	bone_origin		origin of the bone
	\param	bone_vector		bone direction
	\param	bone_front		front direction of the bone
*/
template<typename T1, typename T2, typename T3>
inline void drawWireBone(Vec3<T1> bone_origin, Vec3<T2> bone_vector, Vec3<T3> bone_front){
	float length = bone_vector.Length();
	if( length < 1e-6 )	//	an zero length bone
		return;

	assert( !isMatrixStackFull(GL_MODELVIEW) );

	glEnable(GL_NORMALIZE);
	glPushMatrix();

	glTranslatef(bone_origin.x, bone_origin.y, bone_origin.z);

	Vec3f		bone_vec = bone_vector.Normalize();
	Vec3f		bone_fro = bone_front.Normalize();
	Vec3f		normal	= cross(bone_vec, bone_fro);
	Matrix3x3f	trans_mat(normal, bone_vec, bone_fro);
	Matrix4x4f	gl_trans_mat = trans_mat.SetTranspose();

	glMultMatrixf(gl_trans_mat.data[0]);
	glScalef(length, length, length);

	drawWireBone();

	glPopMatrix();
}

///@}

//	========================================
///@{
/**	@name Draw Mesh with Vector as container
*/
//	========================================

/**
	Draw a single triangle
*/
template<typename T>
void drawTriangle(Vec3<T> v0, Vec3<T> v1, Vec3<T> v2){
	Vec3<T> nor = cross(v1-v0, v2-v2).Normalize();
	glBegin(GL_TRIANGLES);
		glNormal3f(nor.x, nor.y, nor.z);
		glVertex3f(v0.x, v0.y, v0.z);
		glVertex3f(v1.x, v1.y, v1.z);
		glVertex3f(v2.x, v2.y, v2.z);
	glEnd();
}

/**
	Draw Triangle Mesh with Flat Shading, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	face_normal		face normal container
*/
template<typename T>
void drawFlatShadingTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<Vec3<T>>& face_normal)
{
	if( vertex.empty() || face.empty() || face_normal.empty() )
		return;
	drawFlatShadingTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (T*)&face_normal[0]);
}

/**
	Draw Triangle Mesh with Smooth Shading, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_normal	vertex normal container
*/
template<typename T>
void drawSmoothShadingTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<Vec3<T>>&	vertex_normal)
{
	if( vertex.empty() || face.empty() || vertex_normal.empty()
		|| vertex.size() > vertex_normal.size())
		return;
	drawSmoothShadingTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (T*)&vertex_normal[0]);
}
/**
	Draw Triangle Mesh with Given Face Color, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
*/
template<typename T>
void drawFlatColorTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<float3>&	face_color)
{
	if( vertex.empty() || face.empty() || face.size() > face_color.size() )
		return;
	drawFlatColorTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (float*)&face_color[0]);
}

/**
	Draw Triangle Mesh with Given Face Color, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
*/
template<typename T>
void drawFlatColorTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<uchar3>&	face_color)
{
	if( vertex.empty() || face.empty() || face.size() > face_color.size() )
		return;
	drawFlatColorTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (unsigned char*)&face_color[0]);
}

/**
	Draw Triangle Mesh with Given Vertex Color, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_color	vertex color container
*/
template<typename T>
void drawSmoothColorTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<float3>&	vertex_color)
{
	if( vertex.empty() || face.empty() || vertex.size() > vertex_color.size() )
		return;
	drawSmoothColorTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (float*)&vertex_color[0]);
}

/**
	Draw Triangle Mesh with Given Vertex Color, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_color	vertex color container
*/
template<typename T>
void drawSmoothColorTriMesh(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<uchar3>&	vertex_color)
{
	if( vertex.empty() || face.empty() || vertex.size() > vertex_color.size() )
		return;
	drawSmoothColorTriMesh((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (unsigned char*)&vertex_color[0]);
}

/**
	Draw Vertex Normal

	\param	vertex			vertex container
	\param	vertex_normal	vertex normal container
	\param	scale			scale of the normal
*/
template<typename T>
void drawVertexNormal(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<Vec3<T>>&	vertex_normal,
	float						scale = 1.0f)
{
	if( vertex.empty() || vertex.size() > vertex_normal.size() )
		return;
	drawVertexNormal((T*)&vertex[0], vertex.size(), (T*)&vertex_normal[0], scale);
}

/**
	Draw Face Normal, face normals are calculated

	\param	vertex			vertex container
	\param	face			face container
	\param	scale			scale of the normal
*/
template<typename T>
void drawFaceNormal(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	float						scale = 1.0f
) {
	if (vertex.empty() || face.empty())
		return;

	drawFaceNormal(
		&vertex[0].x,
		&face[0].x,
		face.size(),
		(T*)NULL,
		scale
	);
}

/**
	Draw Face Normal

	\param	vertex			vertex container
	\param	face			face container
	\param	face_normal		face normal container	
	\param	scale			scale of the normal
*/
template<typename T>
void drawFaceNormal(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<Vec3<T>>&	face_normal,
	float						scale = 1.0f
) {
	if (vertex.empty() || face.empty())
		return;
	if (face.size() != face_normal.size())
		return;

	drawFaceNormal(
		&vertex[0].x,
		&face[0].x,
		face.size(),
		&face_normal[0].x,
		scale
	);
}

/**
	Draw Edge of Mesh, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
*/
template<typename T>
void drawMeshEdge(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int2>&	edge)
{
	if( vertex.empty() || edge.empty() )
		return;
	drawMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size());
}

/**
	Draw Edge of Mesh with vertex offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshEdge(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int2>&	edge,
	const std::vector<Vec3<T>>&	vertex_normal,
	float						offset = 0.001f)
{
	if( vertex.empty() || edge.empty() || vertex.size() > vertex_normal.size() )
		return;
	drawMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (T*)&vertex_normal[0], offset);
}

/**
	Draw Edge of Mesh with vertex color, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
							 const std::vector<int2>&		edge,
							 const std::vector<float3>&		vertex_color){
	if( vertex.empty() || edge.empty() || vertex.size() > vertex_color.size() )
		return;
	
	drawSmoothColorMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (float*)&vertex_color[0]);
}

/**
	Draw Edge of Mesh with vertex color, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
							 const std::vector<int2>&		edge,
							 const std::vector<uchar3>&		vertex_color){
	if( vertex.empty() || edge.empty() || vertex.size() > vertex_color.size() )
		return;
	
	drawSmoothColorMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (unsigned char*)&vertex_color[0]);
}

/**
	Draw Edge of Mesh with vertex color and offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	vertex_color	color of each vertex
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
							 const std::vector<int2>&		edge, 
							 const std::vector<float3>&		vertex_color,
							 const std::vector<Vec3<T>>&	vertex_normal, 
							 float							offset = 0.001f){
	if( vertex.empty() || edge.empty() || vertex.size() > vertex_color.size() 
		|| vertex.size() > vertex_normal.size() )
		return;
	drawSmoothColorMeshEdge((T*)&vertex[0], vertex.size(), 
		(int*)&edge[0], edge.size(), (float*)&vertex_color[0], (T*)&vertex_normal[0], offset);
}

/**
	Draw Edge of Mesh with vertex color and offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	vertex_color	color of each vertex
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
							 const std::vector<int2>&		edge, 
							 const std::vector<uchar3>&		vertex_color,
							 const std::vector<Vec3<T>>&	vertex_normal, 
							 float							offset = 0.001f){
	if( vertex.empty() || edge.empty() || vertex.size() > vertex_color.size() 
		|| vertex.size() > vertex_normal.size() )
		return;
	drawSmoothColorMeshEdge((T*)&vertex[0], vertex.size(), 
		(int*)&edge[0], edge.size(), (unsigned char*)&vertex_color[0], (T*)&vertex_normal[0], offset);
}

/**
	Draw Edge of Mesh with edge color, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	edge_color		color of each edge
*/
template<typename T>
void drawFlatColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
						   const std::vector<int2>&		edge,
						   const std::vector<float3>&	edge_color){
	if( vertex.empty() || edge.empty() || vertex.size() > edge_color.size() )
		return;
	
	drawFlatColorMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (float*)&edge_color[0]);
}

/**
	Draw Edge of Mesh with edge color, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	edge_color		color of each edge
*/
template<typename T>
void drawFlatColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
						   const std::vector<int2>&		edge,
						   const std::vector<uchar3>&	edge_color){
	if( vertex.empty() || edge.empty() || vertex.size() > edge_color.size() )
		return;
	
	drawFlatColorMeshEdge((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (unsigned char*)&edge_color[0]);
}

/**
	Draw Edge of Mesh with edge color and offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	edge_color		color of each edge
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawFlatColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
						   const std::vector<int2>&		edge, 
						   const std::vector<float3>&	edge_color,
						   const std::vector<Vec3<T>>&	vertex_normal, 
						   float						offset = 0.001f){
	if( vertex.empty() || edge.empty() || vertex.size() > edge_color.size() 
		|| vertex.size() > vertex_normal.size() )
		return;
	drawFlatColorMeshEdge((T*)&vertex[0], vertex.size(), 
		(int*)&edge[0], edge.size(), (float*)&edge_color[0], (T*)&vertex_normal[0], offset);
}

/**
	Draw Edge of Mesh with edge color and offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	edge_color		color of each edge
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawFlatColorMeshEdge(const std::vector<Vec3<T>>&	vertex, 
						   const std::vector<int2>&		edge, 
						   const std::vector<uchar3>&	edge_color,
						   const std::vector<Vec3<T>>&	vertex_normal, 
						   float						offset = 0.001f){
	if( vertex.empty() || edge.empty() || vertex.size() > edge_color.size() 
		|| vertex.size() > vertex_normal.size() )
		return;
	drawFlatColorMeshEdge((T*)&vertex[0], vertex.size(), 
		(int*)&edge[0], edge.size(), (unsigned char*)&edge_color[0], (T*)&vertex_normal[0], offset);
}

/**
	Draw Boundary of Mesh from edge-face, with vector as container

	\param	vertex		vertex container
	\param	edge		edge container
	\param	ef			edge - face connectivity
*/
template<typename T>
void drawMeshBoundaryFromEF(const std::vector<Vec3<T>>&	vertex, 
							const std::vector<int2>&	edge,
							const std::vector<int2>&	ef){
	if( vertex.empty() || edge.empty() || ef.empty() )
		return;
	drawMeshBoundaryFromEF((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (int*)&ef[0]);
}

/**
	Draw Boundary of Mesh from edge-face with vertex offset, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
	\param	ef			edge - face connectivity
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshBoundaryFromEF(const std::vector<Vec3<T>>&	vertex, 
							const std::vector<int2>&	edge, 
							const std::vector<int2>&	ef,
							const std::vector<Vec3<T>>&	vertex_normal, 
							float						offset = 0.001f){
	if( vertex.empty() || edge.empty() || ef.empty() || vertex_normal.empty() )
		return;
	drawMeshBoundaryFromEF((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (int*)&ef[0], (T*)&vertex_normal[0], offset);
}

/**
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark(const std::vector<Vec3<T>>&	vertex, 
								  const std::vector<int2>&		edge,
								  const std::vector<int>&		edge_mark){
	if( vertex.empty() || edge.empty() || edge_mark.empty() )
		return;
	drawMeshBoundaryFromEdgeMark((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (int*)&edge_mark[0]);
}

/**
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark(const std::vector<Vec3<T>>&	vertex, 
								  const std::vector<int2>&		edge, 
								  const std::vector<int>&		edge_mark,
								  const std::vector<Vec3<T>>&	vertex_normal, 
								  float							offset = 0.001f){
	if( vertex.empty() || edge.empty() || edge_mark.empty() || vertex_normal.empty() )
		return;
	drawMeshBoundaryFromEdgeMark((T*)&vertex[0], vertex.size(), 
		(int*)&edge[0], edge.size(), (int*)&edge_mark[0], (T*)&vertex_normal[0], offset);
}
/**
	Draw Edge of Mesh from face, with vector as container

	\param	vertex			vertex container
	\param	face			face container
*/
template<typename T>
void drawMeshEdgeFromFace(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face
) {
	if( vertex.empty() || face.empty() )
		return;
	drawMeshEdgeFromFace((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}

/**
	Draw Edge of Mesh from face with vertex offset, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_normal	vertex normal container
	\param	offset			offset length of each vertex
*/
template<typename T>
void drawMeshEdgeFromFace(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<Vec3<T>>&	vertex_normal,
	float						offset = 0.001f
) {
	if( vertex.empty() || face.empty() || vertex_normal.empty() )
		return;
	drawMeshEdgeFromFace((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (T*)&vertex_normal[0], offset);
}

/**
	Draw mesh edge from face, with flat shading and color

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
	\param	face_normal		face normal container
*/
template<typename T>
void drawFlatColorTriMeshEdgeFromFace(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<uchar3>&	face_color,
	const std::vector<Vec3<T>>& face_normal
) {
	if (vertex.empty() || face.empty() || face_color.empty() || face_normal.empty())
		return;

	drawFlatColorTriMeshEdgeFromFace(
		(T*)&vertex[0],
		vertex.size(),
		(int*)&face[0],
		face.size(),
		(unsigned char*)&face_color[0],
		(T*)&face_normal[0]
	);
}

/**
	Draw mesh edge from face, with flat shading and color

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
	\param	face_normal		face normal container
*/
template<typename T>
void drawFlatColorTriMeshEdgeFromFace(
	const std::vector<Vec3<T>>&	vertex,
	const std::vector<int3>&	face,
	const std::vector<float3>&	face_color,
	const std::vector<Vec3<T>>& face_normal
) {
	if (vertex.empty() || face.empty() || face_color.empty() || face_normal.empty())
		return;

	drawFlatColorTriMeshEdgeFromFace(
		(T*)&vertex[0],
		vertex.size(),
		(int*)&face[0],
		face.size(),
		(float*)&face_color[0],
		(T*)&face_normal[0]
	);
}

/**
	draw index of each vertex of a mesh

	\param	vertex			vertex position of the mesh
*/
template<typename T>
void drawMeshVertexIndex(const std::vector<Vec3<T>>& vertex){
	drawMeshVertexIndex((T*)&vertex[0], vertex.size());
}

/**
	draw number on each vertex of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
*/
template<typename TN, typename TV>
void drawNumberOnMeshVertices(const std::vector<TN>&		number,
							  const std::vector<Vec3<TV>>&	vertex ){
	if( number.empty() || vertex.empty() )
		return;
	drawNumberOnMeshVertices(&number[0], (TV*)&vertex[0], vertex.size());
}

/**
	draw number on each edge of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	edge			edge list
*/
template<typename TN, typename TV>
void drawNumberOnMeshEdges(const std::vector<TN>&		number,
						   const std::vector<Vec3<TV>>&	vertex,
						   const std::vector<int2>&		edge ){
	if( number.empty() || vertex.empty() || edge.empty() )
		return;
	drawNumberOnMeshEdges(&number[0], (TV*)&vertex[0], (int*)&edge[0], edge.size());
}

///@}

//	========================================
///@{
/**	@name Draw Mesh 2D
*/
//	========================================

/**
	Draw a single triangle
*/
template<typename T>
void drawTriangle(Vec2<T> v0, Vec2<T> v1, Vec2<T> v2){
	glBegin(GL_TRIANGLES);
		glVertex2f(v0.x, v0.y);
		glVertex2f(v1.x, v1.y);
		glVertex2f(v2.x, v2.y);
	glEnd();
}

/**
	Draw Triangle Mesh on 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
*/
template<typename T>
void drawTriMesh2D(const std::vector<Vec2<T>>&	vertex, 
				   const std::vector<int3>&		face) {
	if( vertex.empty() || face.empty() )
		return;
	drawTriMesh2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}

/**
	Draw Triangle Mesh with Given Face Color 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
*/
template<typename T>
void drawFlatColorTriMesh2D(const std::vector<Vec2<T>>&	vertex, 
							const std::vector<int3>&	face, 
							const std::vector<float3>&	face_color){
	if( vertex.empty() || face.empty() || face.size() > face_color.size() )
		return;
	drawFlatColorTriMesh2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (float*)&face_color[0]);
}

/**
	Draw Triangle Mesh with Given Face Color 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	face_color		face color container
*/
template<typename T>
void drawFlatColorTriMesh2D(const std::vector<Vec2<T>>&	vertex, 
							const std::vector<int3>&	face, 
							const std::vector<uchar3>&	face_color){
	if( vertex.empty() || face.empty() || face.size() > face_color.size() )
		return;
	drawFlatColorTriMesh2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (unsigned char*)&face_color[0]);
}

/**
	Draw Triangle Mesh with Given Vertex Color 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_color	vertex color container
*/
template<typename T>
void drawSmoothColorTriMesh2D(const std::vector<Vec2<T>>&	vertex, 
							  const std::vector<int3>&		face, 
							  const std::vector<float3>&	vertex_color){
	if( vertex.empty() || face.empty() || vertex.size() > vertex_color.size() )
		return;
	drawSmoothColorTriMesh2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (float*)&vertex_color[0]);
}

/**
	Draw Triangle Mesh with Given Vertex Color 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
	\param	vertex_color	vertex color container
*/
template<typename T>
void drawSmoothColorTriMesh2D(const std::vector<Vec2<T>>&	vertex, 
							  const std::vector<int3>&		face, 
							  const std::vector<uchar3>&	vertex_color){
	if( vertex.empty() || face.empty() || vertex.size() > vertex_color.size() )
		return;
	drawSmoothColorTriMesh2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size(), (unsigned char*)&vertex_color[0]);
}

/**
	Draw Edge of Mesh 2D, with vector as container

	\param	vertex			vertex container
	\param	edge			edge container
*/
template<typename T>
void drawMeshEdge2D(const std::vector<Vec2<T>>& vertex, 
					const std::vector<int2>&	edge){
	if( vertex.empty() || edge.empty() )
		return;
	drawMeshEdge2D((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size());
}

/**
	Draw 2D Mesh Edge with given vertex color, with vector as container

	\param	vertex			vertex coordinate list, xy_xy_
	\param	edge			edge index list, v0v1_v0v1_
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge2D(const std::vector<Vec2<T>>&	vertex, 
							   const std::vector<int2>&		edge,
							   const std::vector<float3>&	vertex_color){
	if( vertex.empty() || edge.empty() || vertex_color.empty() )
		return;
	drawSmoothColorMeshEdge2D((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (float*)&vertex_color[0]);
}

/**
	Draw 2D Mesh Edge with given vertex color, with vector as container

	\param	vertex			vertex coordinate list, xy_xy_
	\param	edge			edge index list, v0v1_v0v1_
	\param	vertex_color	color of each vertex
*/
template<typename T>
void drawSmoothColorMeshEdge2D(const std::vector<Vec2<T>>&	vertex, 
							   const std::vector<int2>&		edge,
							   const std::vector<uchar3>&	vertex_color){
	if( vertex.empty() || edge.empty() || vertex_color.empty() )
		return;
	drawSmoothColorMeshEdge2D((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (unsigned char*)&vertex_color[0]);
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity, with vector as container

	\param	vertex			vertex coordinate list, xy_xy_
	\param	edge			edge index list, v0v1_v0v1_
	\param	ef				edge - face connectivity
*/
template<typename T>
void drawMeshBoundaryFromEF2D(const std::vector<Vec2<T>>&	vertex, 
							  const std::vector<int2>&		edge,
							  const std::vector<int2>&		ef ){
	if( vertex.empty() || edge.empty() || ef.empty() )
		return;
	drawMeshBoundaryFromEF2D((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), (int*)&ef[0]);
}

/**
	Draw 2D Mesh Boundary from edge-face connectivity, with vector as container

	\param	vertex			vertex coordinate list, xy_xy_
	\param	edge			edge index list, v0v1_v0v1_
	\param	edge_mark		boundary mark of each edge, 1: is boundary; 0: not boundary
*/
template<typename T>
void drawMeshBoundaryFromEdgeMark2D(const std::vector<Vec2<T>>&	vertex, 
									const std::vector<int2>&	edge,
									const std::vector<int>&		edge_mark){
	if( vertex.empty() || edge.empty() || edge_mark.empty() )
		return;
	drawMeshBoundaryFromEdgeMark2D((T*)&vertex[0], vertex.size(), (int*)&edge[0], edge.size(), &edge_mark[0]);
}

/**
	Draw Edge of Mesh from face 2D, with vector as container

	\param	vertex			vertex container
	\param	face			face container
*/
template<typename T>
void drawMeshEdgeFromFace2D(const std::vector<Vec2<T>>& vertex, 
							const std::vector<int3>&	face){
	if( vertex.empty() || face.empty() )
		return;
	drawMeshEdgeFromFace2D((T*)&vertex[0], vertex.size(), (int*)&face[0], face.size());
}


/**
	draw index of each vertex of a mesh

	\param	vertex			vertex position of the mesh
*/
template<typename T>
void drawMeshVertexIndex2D(const std::vector<Vec2<T>>& vertex){
	if( vertex.empty() )
		return;
	drawMeshVertexIndex2D((T*)&vertex[0], vertex.size());
}

/**
	draw number on each vertex of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
*/
template<typename TN, typename TV>
void drawNumberOnMeshVertices2D(const std::vector<TN>&			number,
								const std::vector<Vec2<TV>>&	vertex ){
	if( number.empty() || vertex.empty() )
		return;
	drawNumberOnMeshVertices2D(&number[0], (TV*)&vertex[0], vertex.size());
}

/**
	draw number on each edge of a mesh

	\param	number			number to draw
	\param	vertex			vertex position of the mesh
	\param	edge			edge list
*/
template<typename TN, typename TV>
void drawNumberOnMeshEdges2D(const std::vector<TN>&			number,
							 const std::vector<Vec2<TV>>&	vertex,
							 const std::vector<int2>&		edge ){
	if( number.empty() || vertex.empty() || edge.empty() )
		return;
	drawNumberOnMeshEdges2D(&number[0], (TV*)&vertex[0], (int*)&edge[0], edge.size());
}

///@}

}}	//	end namespace yz::opengl


#endif	//	__YZ_VECTOR_OPENGL_UTILS_H__