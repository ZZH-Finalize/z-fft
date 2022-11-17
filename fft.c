/*
@file: fft.c
@author: ZZH
@date: 2022-11-05
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
    for (size_t i = 0;i < len;i++)
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

void idft(number_t* output, Complex_t* input, size_t len)
{
    number_t tmp = PI_2 / len;
    for (size_t i = 0;i < len;i++)
    {
        number_t res = 0;
        number_t w = tmp * i;
        for (size_t j = 0;j < len;j++)
        {
            number_t arg = j * w;
            res += input[j].real * cos(arg) - input[j].image * sin(arg);
        }
        output[i] = res;
    }
}

void length(number_t* output, Complex_t* input, size_t len)
{
    for (size_t i = 0;i < len;i++)
        output[i] = sqrt(input[i].real * input[i].real + input[i].image * input[i].image);
}
