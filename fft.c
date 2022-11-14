/*
@file: fft.c
@author: ZZH
@time: 2022-11-05 17:10:13
@info: fft实现
*/
#include "fft.h"

const double PI_2 = 2 * M_PI;

void fft(Complex_t* output, number_t* input, size_t len)
{
}

void ifft(number_t* output, Complex_t* input, size_t len)
{
    
}

void dft(Complex_t* output, number_t* input, size_t len)
{
    number_t tmp = PI_2 / len;
    for (size_t i = 0;i < len / 2;i++)
    {
        number_t real = 0, image = 0;
        number_t w = tmp * i;
        for (size_t j = 0;j < len;j++)
        {
            number_t arg = j * w;
            real += input[j] * cos(arg);
            image -= input[j] * sin(arg);
        }
        output[i].real = real;
        output[i].image = image;
    }
}

void length(number_t* output, Complex_t* input, size_t len)
{
    for (size_t i = 0;i < len;i++)
        output[i] = sqrt(input[i].real * input[i].real + input[i].image * input[i].image);
}
