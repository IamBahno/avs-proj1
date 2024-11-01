/**
 * @file LineMandelCalculator.h
 * @author FULL NAME <xbahou0000@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 1.11.2024
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    int *data;
    float *real_arr;
	float *img_arr; 
};