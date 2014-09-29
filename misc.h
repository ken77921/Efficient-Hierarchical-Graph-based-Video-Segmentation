/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011, 2012 
  Chenliang Xu, Jason Corso.
Modifications Copyright (C) 2014
  Haw-Shiuan Chang
  
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

/* random stuff */

#ifndef MISC_H
#define MISC_H

#include <vector>

/*!
 * \brief Randomly permute the values 0 to n-1 and store the permuted list in a vector.
 *
 * \param n the upper bound on the range of values to generate and permute
 * \param v the vector to store the permuted value list in
 * \param pointer to a random seed, necessary for thread-safe random number generation
 */

// for color space transfer
const float ref_X = 95.047;
const float ref_Y = 100.000;
const float ref_Z = 108.883;

typedef unsigned char uchar;

typedef struct { uchar r, g, b; } rgb;
typedef struct { float L, a, b; } Lab;

inline bool operator==(const rgb &a, const rgb &b) {
  return ((a.r == b.r) && (a.g == b.g) && (a.b == b.b));
}

template <class T>
inline T abs(const T &x) { return (x > 0 ? x : -x); };

template <class T>
inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };

template <class T>
inline T square(const T &x) { return x*x; };

void my_RandomPermuteRange(unsigned int *v,const int Large_const,const int range)
{
	for(int i = 0; i < range; i++)
	{
		int index=int(rand()/(double)RAND_MAX*(Large_const-1));
		v[i]=v[index];
		v[index]=i;
	}
}

//void juRandomPermuteRange(int n, vector<int>& v, unsigned int* seed)
//{
//    v.clear();
//    v.resize(n);
//    int i, j;
//    v[0] = 0;
//	srand (*seed);
//    for (i = 1; i < n; i++)
//    {
//        j = rand() % (i + 1);
//        v[i] = v[j];
//        v[j] = i;
//    }
//}

bool isWhite(rgb color){
	return (color.r == 255 && color.g == 255 && color.b == 255);
}

// convert single int into unique rgb values
void setColor(rgb& color, int n){
	color.r = ((n>>16) % 256);
	color.g = ((n>>8) % 256);
	color.b = (n % 256);
}

// color space transfer
Lab RGBtoLab(rgb &color_in) {
	float var_R = color_in.r / 255.0;
	float var_G = color_in.g / 255.0;
	float var_B = color_in.b / 255.0;

	/*float var_R = p_node[a].color.r / 255.0;
	float var_G = p_node[a].color.g / 255.0;
	float var_B = p_node[a].color.b / 255.0;*/

	if (var_R > 0.04045)
		var_R = pow(((var_R + 0.055) / 1.055), 2.4);
	else
		var_R = var_R / 12.92;

	if (var_G > 0.04045)
		var_G = pow(((var_G + 0.055) / 1.055), 2.4);
	else
		var_G = var_G / 12.92;

	if (var_B > 0.04045)
		var_B = pow(((var_B + 0.055) / 1.055), 2.4);
	else
		var_B = var_B / 12.92;

	var_R = var_R * 100;
	var_G = var_G * 100;
	var_B = var_B * 100;

	//Observer. = 2¢X, Illuminant = D65
	float var_X = (var_R * 0.412453f + var_G * 0.357580f + var_B * 0.180423f)/ref_X;
	float var_Y = (var_R * 0.212671f + var_G * 0.715160f + var_B * 0.072169f)/ref_Y;
	float var_Z = (var_R * 0.019334f + var_G * 0.119193f + var_B * 0.950227f)/ref_Z;

	if (var_X > 0.008856)
		var_X = pow((double) var_X, 0.333333);
	else
		var_X = (7.787 * var_X) + (16 / 116);
	if (var_Y > 0.008856)
		var_Y = pow((double) var_Y, 0.333333);
	else
		var_Y = (7.787 * var_Y) + (16 / 116);
	if (var_Z > 0.008856)
		var_Z = pow((double) var_Z, 0.333333);
	else
		var_Z = (7.787 * var_Z) + (16 / 116);

	Lab color_out;
	color_out.L = (116 * var_Y) - 16; 
	color_out.a = 500 * (var_X - var_Y);
	color_out.b = 200 * (var_Y - var_Z);
	return color_out;
}

#endif
