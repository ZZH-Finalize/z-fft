/*
@file: benchmark.c
@author: ZZH
@date: 2022-11-24
@info: fft和dft性能差异对比, 此程序可能会运行长达几十秒, 取决于N的大小, 因为dft是非常低效的算法
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "fft.h"

#ifndef N
#define N 16384//此N值在我的电脑上, dft需要22秒多才能算出结果, 但同等条件fft仅需3毫秒
#endif // !N
#define FS 5000

int main(const int argc, const char** argv)
{
    //构建一个含有多频率成分的信号
    Number_t input[N];
    for (int i = 0;i < N;i++)
    {
        double arg = 2 * M_PI * i / FS;
        input[i] = sin(70 * arg) + sin(800 * arg) + sin(2200 * arg);
    }

    Complex_t out1[N], out2[N];

    clock_t dftStart = clock();
    rdft(out1, input, N);
    clock_t dftEnd = clock();

    clock_t fftStart = clock();
    rfft(out2, input, N);
    clock_t fftEnd = clock();

    clock_t dftCost = dftEnd - dftStart, fftCost = fftEnd - fftStart;

    printf("dft cost: %ld, fft cost: %ld\r\n", dftCost, fftCost);

    return dftCost > fftCost ? 0 : -1;
}
