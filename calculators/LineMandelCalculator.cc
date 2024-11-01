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


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	data = (int *)(calloc(height * width,sizeof(int)));
	real_arr = (float *)(calloc(width * 2, sizeof(float)));
	img_arr = (float *)(calloc(width * 2, sizeof(float)));

}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	free(real_arr);
	free(img_arr);

}


int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *real_array = real_arr;
	float *img_array = img_arr;

	//computint only for half of the picture
	for (int i = 0; i < height/2; i++)
	{
		float y = y_start + i * dy; // current imaginary value

		// variable that allows me to break when all the point in row, breaken out of the < 4 boundary
		int is_running = 0;
		for (int k = 0; k < limit; ++k)
		{
			#pragma omp simd reduction(+:is_running)
			for (int j = 0; j < width; j++)
			{
				float x = x_start + j * dx; // current real value

				float zReal = (k == 0) ? x : real_array[j];
				float zImag = (k == 0) ? y : img_array[j];

				float r2 = zReal * zReal;
				float i2 = zImag * zImag;

				if (r2 + i2 < 4.0f)
				{
					pdata[i*width+j] += 1;
					img_array[j] = 2.0f * zReal * zImag + y;
					real_array[j] = r2 - i2 + x;
					is_running += 1;
				}
			}
			if(is_running != width){
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
