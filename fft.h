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

    //复数变换使用output的内存空间进行蝶形运算
    void fft(Complex_t* output, Complex_t* input, size_t len);
    void ifft(Complex_t* output, Complex_t* input, size_t len);

    //实数变换的正变换使用output的内存空间, 反变换使用input的内存空间. 也就是说实数反变换会改变输入的复数值
    void rfft(Complex_t* output, Number_t* input, size_t len);
    void irfft(Number_t* output, Complex_t* input, size_t len);

    void dft(Complex_t* output, Complex_t* input, size_t len);
    void idft(Complex_t* output, Complex_t* input, size_t len);

    void rdft(Complex_t* output, Number_t* input, size_t len);
    void irdft(Number_t* output, Complex_t* input, size_t len);

    void length(Number_t* output, Complex_t* input, size_t len);

#ifdef __cplusplus
}
#endif // __cplusplus
