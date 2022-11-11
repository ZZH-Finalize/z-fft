/*
@file: fft.h
@author: ZZH
@time: 2022-11-05 17:09:49
@info: fft
*/
#pragma once
#include <stdint.h>
#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


    void fft(void* output, void* input, size_t len);
    void ifft(void* output, void* input, size_t len);

#ifdef __cplusplus
}
#endif // __cplusplus
