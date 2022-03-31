/***********************************************************/
/**	\file
	\brief		Convert value to color
	\details	A lot of color conversion methods are implemented
				in this file. The following functions are provided:

				1,	convertToColor[X]		convert single value to color

				2,	convertToColorArray[X]	convert value array to color array

				[X] represent method of color conversion, method come
				from Matlab, the following can be chosen:
				Simple, Jet, HSV, Hot, Cool, Gray
	\author		Yizhong Zhang
	\date		5/12/2012
*/
/***********************************************************/
/***************************************************
convert value to color using different method, 
float is good enough to represent color

This header provide following functions
	1,	convertToColor[X]		convert single value to color
	2,	convertToColorArray[X]	convert value array to color array
[X] represent method of color conversion, method come
from Matlab, the following can be chosen:
Simple, Jet, HSV, Hot, Cool, Gray

resolution: how many discrete zones (1/res) is the bar divided
			0 means no discretization

Yizhong Zhang, 5/12/2012
***************************************************/
#ifndef __YZ_COLOR_BAR_H__
#define __YZ_COLOR_BAR_H__

namespace yz{	namespace utils{

//	========================================
///@{
/**	@name Convert Single Value to Single Color
*/
//	========================================

/**
	Simple Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorSimple(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	float color[3];
	convertToColorSimple(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorSimple<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float lamda = (value - min_value) / (max_value - min_value); 
	if( resolution > 0.001 )
		lamda = int(lamda / resolution) * resolution;

	float coef = lamda * 4.0f;

	if( coef <= 1.0f )
		rgb[0] = 0.0f, rgb[1] = coef, rgb[2] = 1.0f;
	else if( coef <= 2.0f )
		rgb[0] = 0.0f, rgb[1] = 1.0f, rgb[2] = 2.0f-coef;
	else if( coef <= 3.0f )
		rgb[0] = coef-2.0f, rgb[1] = 1.0f, rgb[2] = 0.0f;
	else
		rgb[0] = 1.0f, rgb[1] = 4.0f-coef, rgb[2] = 0.0f;
}

template<>
inline void convertToColorSimple<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorSimple(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}

/**
	Jet Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorJet(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){
	float color[3];
	convertToColorJet(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorJet<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float lamda = (value - min_value) / (max_value - min_value); 
	if( resolution > 0.001 )
		lamda = int(lamda / resolution) * resolution;

	float coef = lamda * 8.0f;

	if( coef <= 1.0f )
		rgb[0] = 0.0f, rgb[1] = 0.0f, rgb[2] = 0.5f*coef + 0.5f;
	else if( coef <= 3.0f )
		rgb[0] = 0.0f, rgb[1] = 0.5f*(coef-1.0f), rgb[2] = 1.0f;
	else if( coef <= 5.0f )
		rgb[0] = 0.5f*(coef-3.0f), rgb[1] = 1.0f, rgb[2] = 1.0f-0.5f*(coef-3.0f);
	else if( coef <= 7.0f )
		rgb[0] = 1.0f, rgb[1] = 1.0f-0.5f*(coef-5.0f), rgb[2] = 0.0f;
	else
		rgb[0] = 1.0f-0.5f*(coef-7.0f), rgb[1] = 0.0f, rgb[2] = 0.0f;
}

template<>
inline void convertToColorJet<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorJet(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}

/**
	HSV Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorHSV(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){
	float color[3];
	convertToColorHSV(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorHSV<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float lamda = (value - min_value) / (max_value - min_value); 
	if( resolution > 0.001 )
		lamda = int(lamda / resolution) * resolution;

	float coef = lamda * 6.0f;

	if( coef <= 1.0f )
		rgb[0] = 1.0f, rgb[1] = coef, rgb[2] = 0.0f;
	else if( coef <= 2.0f )
		rgb[0] = 1.0f-(coef-1.0f), rgb[1] = 1.0f, rgb[2] = 0.0f;
	else if( coef <= 3.0f )
		rgb[0] = 0.0f, rgb[1] = 1.0f, rgb[2] = coef-2.0f;
	else if( coef <= 4.0f )
		rgb[0] = 0.0f, rgb[1] = 1.0f-(coef-3.0f), rgb[2] = 1.0f;
	else if( coef <= 5.0f )
		rgb[0] = coef-4.0f, rgb[1] = 0.0f, rgb[2] = 1.0f;
	else
		rgb[0] = 1.0f, rgb[2] = 0.0f, rgb[2] = 1.0f-(coef-5.0f);
}

template<>
inline void convertToColorHSV<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorHSV(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}

/**
	Hot Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorHot(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){
	float color[3];
	convertToColorHot(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorHot<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float lamda = (value - min_value) / (max_value - min_value); 
	if( resolution > 0.001 )
		lamda = int(lamda / resolution) * resolution;

	float coef = lamda * 3.0f;

	if( coef <= 1.0f )
		rgb[0] = coef, rgb[1] = 0.0f, rgb[2] = 0.0f;
	else if( coef <= 2.0f )
		rgb[0] = 1.0f, rgb[1] = coef-1.0f, rgb[2] = 0.0f;
	else
		rgb[0] = 1.0f, rgb[1] = 1.0f, rgb[2] = coef-2.0f;
}

template<>
inline void convertToColorHot<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorHot(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}

/**
	Cool Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorCool(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){
	float color[3];
	convertToColorCool(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorCool<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float coef = (value - min_value) / (max_value - min_value);
	if( resolution > 0.001 )
		coef = int(coef / resolution) * resolution;

	rgb[0] = coef, rgb[1] = 1.0f-coef, rgb[2] = 1.0f;
}

template<>
inline void convertToColorCool<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorCool(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}

/**
	Gray Conversion method

	\param	rgb			color ptr
	\param	value		value to be converted
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename T>
inline void convertToColorGray(T rgb[3], float value, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){
	float color[3];
	convertToColorGray(color, value, min_value, max_value, resolution);
	rgb[0] = color[0] * 255;
	rgb[1] = color[1] * 255;
	rgb[2] = color[2] * 255;
}

template<>
inline void convertToColorGray<float>(float rgb[3], float value, float min_value, float max_value, float resolution){
	if( value < min_value )
		value = min_value;
	else if( value > max_value )
		value = max_value;

	float coef = (value - min_value) / (max_value - min_value);
	if( resolution > 0.001 )
		coef = int(coef / resolution) * resolution;

	rgb[0] = coef, rgb[1] = coef, rgb[2] = coef;
}

template<>
inline void convertToColorGray<double>(double rgb[3], float value, float min_value, float max_value, float resolution){
	float rgb_f[3];
	convertToColorGray(rgb_f, value, min_value, max_value, resolution);
	rgb[0] = rgb_f[0];
	rgb[1] = rgb_f[1];
	rgb[2] = rgb_f[2];
}
///@}

//	========================================
///@{
/**	@name Convert Value Array to Color Array
*/
//	========================================

/**
	Simple Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArraySimple(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorSimple(rgb+i*3, value[i], min_value, max_value, resolution);
}
/**
	Jet Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArrayJet(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorJet(rgb+i*3, value[i], min_value, max_value, resolution);
}
/**
	HSV Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArrayHSV(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorHSV(rgb+i*3, value[i], min_value, max_value, resolution);
}
/**
	Hot Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArrayHot(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorHot(rgb+i*3, value[i], min_value, max_value, resolution);
}
/**
	Cool Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArrayCool(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorCool(rgb+i*3, value[i], min_value, max_value, resolution);
}
/**
	Gray Conversion method

	\param	rgb			color array ptr
	\param	value		value array ptr
	\param	size		value number
	\param	min_value	value of minimal color
	\param	max_value	value of maximal color
	\param	resolution	how many discrete zones (1/res) is the bar divided. 
						0: no discretization.
*/
template<typename TC, typename TV>
inline void convertToColorArrayGray(TC* rgb, const TV* value, int size, float min_value=0.0f, float max_value=1.0f, float resolution=0.0f){	
	for( int i=0; i<size; i++ )
		convertToColorGray(rgb+i*3, value[i], min_value, max_value, resolution);
}
///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_COLOR_BAR_H__