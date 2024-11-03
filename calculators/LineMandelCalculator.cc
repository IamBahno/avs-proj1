/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xbahou00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 11.1.2024
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <malloc.h>
#include <immintrin.h>

#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(_mm_malloc(height * width * sizeof(int),64));
	real_arr = (float *)(_mm_malloc(width * sizeof(float),64));
	img_arr = (float *)(_mm_malloc(width * sizeof(float),64));
	x_axis = (float *)(_mm_malloc(width * sizeof(float),64));
	y_axis = (float *)(_mm_malloc((height/2) * sizeof(float),64));
}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	_mm_free(real_arr);
	_mm_free(img_arr);
	_mm_free(x_axis);
	_mm_free(y_axis);

}


int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;

	// precompute the tick on the axis
	for(int i = 0; i < height/2;i++){
		y_axis[i] = y_start + i * dy;
	}
	for(int i = 0; i < width;i++){
		x_axis[i] = x_start + i * dx;
	}

	//computint only for half of the picture
	for (int i = 0; i < height/2; i++)
	{
		// float y = y_start + i * dy; // current imaginary value
		float y = y_axis[i]; // current imaginary value

		//index to line in pdata
		int line_in_data = i*width;

		for (int k = 0; k < limit; ++k)
		{
			// variable that allows me to break when all the point in row, breaken out of the < 4 boundary
			int is_running = 0;
			#pragma omp simd reduction(+:is_running) simdlen(32)
			for (int j = 0; j < width; j++)
			{
				float x = x_axis[j]; // current real value

				float zReal = (k == 0) ? x : real_arr[j];
				float zImag = (k == 0) ? y : img_arr[j];

				float r2 = zReal * zReal;
				float i2 = zImag * zImag;

				if (r2 + i2 < 4.0f)
				{
					pdata[line_in_data+j] += 1;
					img_arr[j] = 2.0f * zReal * zImag + y;
					real_arr[j] = r2 - i2 + x;
					is_running += 1;
				}
			}
			
			if(is_running == 0){
				break;
			}
		}
	}

	//duplicate the second half of the matrix
	for (int i = 0; i < height / 2; i++) {
		for (int j = 0; j < width; j++) {
			pdata[(height - 1 - i) * width + j] = pdata[i * width + j];
		}
	}

	return pdata;
}
