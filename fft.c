/*
@file: fft.c
@author: ZZH
@time: 2022-11-05 17:10:13
@info: fft实现
*/
#include"fft.h"

const double PI_2 = 2 * M_PI;

void fft(double* output, double* input, size_t len)
{
}

void ifft(double* output, double* input, size_t len)
{
}

void dft(double* output, double* input, size_t len)
{
    for (size_t i = 0;i < len;i += 2)
    {
        double real = 0, image = 0;
        for (size_t j = 0;j < len;j++)
        {
            double arg = PI_2 * ((double) (i * j) / (double) len);
            real += input[j] * sin(arg);
            image += input[j] * cos(arg);
        }
        output[i] = real;
        output[i + 1] = image;
    }
}
