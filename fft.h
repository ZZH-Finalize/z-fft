/*
@file: fft.h
@author: ZZH
@time: 2022-11-05 17:09:49
@info: fft
*/
#pragma once
#include <stdint.h>
#include <math.h>

#ifndef BASIC_TYPE
#define BASIC_TYPE float
#endif

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    typedef BASIC_TYPE number_t;

    typedef struct
    {
        number_t real;
        number_t image;
    }Complex_t, * pComplex_t;

    void fft(Complex_t* output, number_t* input, size_t len);
    void ifft(number_t* output, Complex_t* input, size_t len);
    void dft(Complex_t* output, number_t* input, size_t len);
    void idft(number_t* output, Complex_t* input, size_t len);

    void length(number_t* output, Complex_t* input, size_t len);

#ifdef __cplusplus
}
#endif // __cplusplus
