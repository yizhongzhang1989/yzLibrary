/***********************************************************/
/**	\file
	\brief		Raster of Image
	\author		Yizhong Zhang
	\date		8/10/2019
*/
/***********************************************************/
#ifndef __YZ_IMAGE_RASTER_H__
#define __YZ_IMAGE_RASTER_H__

#include <queue>
#include "yzLib/yz_geometry/yz_clipping.h"
#include "yzlib/yz_utils/yz_ascii_table.h"


namespace yz {	namespace image {

/**
rasterize basic primitives on image

To use this class, image reference must be set at constructor.
Two constructors are provided, depending on whether stride exist.
The image array can be any type, but T_per_pixel must be correctly set (Tpp in constructor)

How to use:	\n
1, Create a Raster and setup image with constructor
2, call SetPenWidth() and SetPenColor()
3, call DrawXXX() to draw primitives on the image
*/
template <class T>
class Raster {
public:
	//	image reference
	const int		width, height, T_per_pixel;
	int				stride;			///< image stride in bytes (number of bytes each row)
	T*	const		image;			///< T can either be basic type, or user defined type, but must match  T_per_pixel

	//	pen control
	int				pen_width;		///< width of pen when displaying, legal value: 1 2 3
	T				pen_color[16];	///< foreground display color

public:
	Raster(int w, int h, int Tpp, T* img) :
		width(w), height(h), T_per_pixel(Tpp), image(img)
	{
		if (Tpp > 16) {
			std::cout << "error: Raster::Raster(), max T_per_pixel is 16, curr value: " << Tpp << std::endl;
			return;
		}

		stride = width*T_per_pixel * sizeof(T);

		pen_width = 1;
		std::fill(&pen_color[0], &pen_color[16], 0);
		//memset(pen_color, 0, sizeof(T) * 16);
	}

	Raster(int w, int h, int s, int Tpp, T* img) :
		width(w), height(h), stride(s), T_per_pixel(Tpp), image(img)
	{
		if (Tpp > 16) {
			std::cout << "error: Raster::Raster(), max T_per_pixel is 16, curr value: " << Tpp << std::endl;
			return;
		}

		if (!s)		//	illegal stride, treate as no stride
			stride = width*T_per_pixel * sizeof(T);

		pen_width = 1;
		std::fill(&pen_color[0], &pen_color[16], 0);
		//memset(pen_color, 0, sizeof(T) * 16);
	}

	inline void SetPenWidth(int pen_w) {
		pen_width = pen_w;
		if (pen_width < 1)
			pen_width = 1;
		else if (pen_width > 3)
			pen_width = 3;
	}

	inline void SetPenColor(T* color) {
		memcpy(pen_color, color, sizeof(T)*T_per_pixel);
	}

	inline void DrawPoint(int x, int y) {
		switch (pen_width){
		case 1:
			PlotPixel(x, y);
			break;
		case 2:
			PlotPixel(x, y);
			PlotPixel(x + 1, y);
			PlotPixel(x, y + 1);
			PlotPixel(x + 1, y + 1);
			break;
		case 3:
			PlotPixel(x - 1, y - 1);
			PlotPixel(x, y - 1);
			PlotPixel(x + 1, y - 1);
			PlotPixel(x - 1, y);
			PlotPixel(x, y);
			PlotPixel(x + 1, y);
			PlotPixel(x - 1, y + 1);
			PlotPixel(x, y + 1);
			PlotPixel(x + 1, y + 1);
			break;
		default:
			PlotPixel(x, y);
			break;
		}
	}

	inline void DrawLine(int x1, int y1, int x2, int y2) {
		int2 p1(x1, y1), p2(x2, y2);
		if (!geometry::clipCohenSutherland(p1, p2, int2(0, 0), int2(width - 1, height - 1)))
			return;

		x1 = p1.x;
		y1 = p1.y;
		x2 = p2.x;
		y2 = p2.y;

		int steep = 0;
		int sx = ((x2 - x1) > 0) ? 1 : -1;
		int sy = ((y2 - y1) > 0) ? 1 : -1;
		int dx = abs(x2 - x1);
		int dy = abs(y2 - y1);

		if (dy > dx) {
			mySwap(x1, y1);
			mySwap(dx, dy);
			mySwap(sx, sy);
			steep = 1;
		}

		int e = 2 * dy - dx;

		for (int i = 0; i < dx; ++i) {
			if (steep)
				DrawPoint(y1, x1);
			else
				DrawPoint(x1, y1);

			while (e >= 0) {
				y1 += sy;
				e -= (dx << 1);
			}

			x1 += sx;
			e += (dy << 1);
		}

		DrawPoint(x2, y2);
	}

	inline void DrawCircle(int center_x, int center_y, int radius) {
		int x = 0;
		int d = (1 - radius) << 1;

		while (radius >= 0) {
			DrawPoint(center_x + x, center_y + radius);
			DrawPoint(center_x + x, center_y - radius);
			DrawPoint(center_x - x, center_y + radius);
			DrawPoint(center_x - x, center_y - radius);

			if ((d + radius) > 0)
				d -= ((--radius) << 1) - 1;
			if (x > d)
				d += ((++x) << 1) + 1;
		}
	}

	inline void DrawCharacter8x16(unsigned char c, int x, int y){
		const unsigned char* mask = &utils::ascii8x16[c * 16];
		for (int j = 0; j < 16; j++) {
			unsigned char test_bit = 0x80;
			for (int i = 0; i < 8; i++) {
				if (mask[j] & test_bit)
					PlotPixel(x + i, y + j);
				test_bit >>= 1;
			}
		}
	}

	inline void DrawString8x16(const char* str, int x, int y) {
		int len = strlen(str);
		for (int i = 0; i < len; i++) {
			DrawCharacter8x16(str[i], x + i * 8, y);
		}
	}

	inline int FloodFill(int x, int y, int nei_count = 4) {
		//	seed pixel color
		T color[16], nei_color[16];
		if (!GetColor(color, x, y))
			return 0;
		if (ColorEqual(color, pen_color))
			return 0;

		int count = 1;
		int nei_offset_x[8] = { 0, -1, 1, 0, -1, 1, -1, 1 };
		int nei_offset_y[8] = { -1, 0, 0, 1, -1, -1, 1, 1 };

		//	setup seed
		std::queue<int2> queue_pixel;
		queue_pixel.push(int2(x, y));
		PlotPixel(x, y);

		//	floodfill
		while (!queue_pixel.empty()) {
			int2 curr_pixel = queue_pixel.front();
			queue_pixel.pop();

			for (int i = 0; i < nei_count; i++) {
				int2 nei_pixel(curr_pixel.x + nei_offset_x[i], curr_pixel.y + nei_offset_y[i]);
				if (!GetColor(nei_color, nei_pixel.x, nei_pixel.y))
					continue;
				if (ColorEqual(nei_color, color)) {
					queue_pixel.push(nei_pixel);
					PlotPixel(nei_pixel.x, nei_pixel.y);
					count++;
				}
			}
		}

		return count;
	}

protected:
	inline void PlotPixel(int x, int y) {
		if ((unsigned int)x >= width || (unsigned int)y >= height)
			return;
		int offset = y * stride + x * T_per_pixel * sizeof(T);
		memcpy((unsigned char*)image + offset, pen_color, sizeof(T)*T_per_pixel);
	}

	inline bool GetColor(T* color, int x, int y) {
		if ((unsigned int)x >= width || (unsigned int)y >= height)
			return false;
		int offset = y * stride + x * T_per_pixel * sizeof(T);
		memcpy(color, (unsigned char*)image + offset, sizeof(T)*T_per_pixel);
		return true;
	}

	inline bool ColorEqual(T* color1, T* color2) {
		for (int i = 0; i < T_per_pixel; i++)
			if (color1[i] != color2[i])
				return false;
		return true;
	}
};

}}


#endif	//	__YZ_IMAGE_RASTER_H__