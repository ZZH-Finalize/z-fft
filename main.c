/*
@file: main.c
@author: ZZH
@time: 2022-11-05 17:10:22
@info: fft测试案例
*/
#include <stdio.h>
#include <string.h>
#include "fft.h"

#define N 128
#define FS 5000

int main(const int argc, const char** argv)
{
    number_t input[N];
    for (int i = 0;i < N;i++)
    {
        double arg = 2 * M_PI * i / FS;
        input[i] = sin(70 * arg) + sin(800 * arg) + sin(2200 * arg);
    }

    Complex_t out[N];

    dft(out, input, N);

    number_t out2[N];

    length(out2, out, N);

    printf("length:\r\n");

    for (int i = 0;i < N / 2;i++)
    {
        printf("[%llf]: %llf\r\n", i * (double) FS / N, out2[i]);
    }

    // memset(&out[N / 2], 0, sizeof(out) / 2);

    idft(out2, out, N);

    for (int i = 0;i < N;i++)
    {
        printf("orig[%d]: %llf, idft[%d]: %llf\r\n", i, input[i], i, out2[i]);
    }

    return 0;
}
