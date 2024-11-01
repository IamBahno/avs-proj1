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
	// @TODO find out why 2 *
	data = (int *)(calloc(height * width,sizeof(int)));
	real_arr = (float *)(calloc(width, sizeof(float)));
	img_arr = (float *)(calloc(width, sizeof(float)));

}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	free(real_arr);
	free(img_arr);

}


int * LineMandelCalculator::calculateMandelbrot () {
	// @TODO implement the calculator & return array of integers
	int *pdata = data;
	float *real_array = real_arr;
	float *img_array = img_arr;

	for (int i = 0; i < height; i++)
	{
		float y = y_start + i * dy; // current imaginary value
		//counting iterations
		bool allOver;
		for (int k = 0; k < limit; ++k)
		{
			allOver = true;
			#pragma omp simd
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
					allOver = false;
				}
			}
			if(allOver){
				break;
			}
		}
	}
	return pdata;
}
