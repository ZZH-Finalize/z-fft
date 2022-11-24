/*
@file: fft.h
@author: ZZH
@date: 2022-11-05
@info: fft
*/
#pragma once
#include <stdint.h>
#include <math.h>

//可以通过构建系统传入的宏定义替换此处的类型
#ifndef BASIC_TYPE
#define BASIC_TYPE float
#endif

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    typedef BASIC_TYPE Number_t;

    typedef struct
    {
        Number_t real;
        Number_t image;
    }Complex_t, * pComplex_t;

    void fft(Complex_t* output, Number_t* input, size_t len);
    void ifft(Number_t* output, Complex_t* input, size_t len);
    void dft(Complex_t* output, Number_t* input, size_t len);
    void idft(Number_t* output, Complex_t* input, size_t len);

    void length(Number_t* output, Complex_t* input, size_t len);

#ifdef __cplusplus
}
#endif // __cplusplus
