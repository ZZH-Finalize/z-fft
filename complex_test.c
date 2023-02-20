/*
@file: complex_test.c
@author: ZZH
@date: 2023-01-06
@info: 复数变换测试案例
*/
#include <stdio.h>
#include <string.h>
#include "fft.h"

#define N 128
#define FS 5000

int main(const int argc, const char** argv)
{
    //构建一个含有多频率成分的复数信号
    Complex_t input[N];
    for (int i = 0;i < N;i++)
    {
        double arg = 2 * M_PI * i / FS;
        input[i].real = cos(70 * arg) + cos(800 * arg) + cos(2200 * arg);//I路数据
        input[i].image = sin(70 * arg) + sin(800 * arg) + sin(2200 * arg);//Q路数据
    }

    Complex_t out[N];

    //傅立叶正变换, 把input的信号变为复数频谱
    fft(out, input, N);

    Number_t out2[N];

    //对复数频谱求模得到实数, 这些实数就是每个频率下的信号强度
    length(out2, out, N);

    //打印幅度谱, 方框内将下标解算为实际频率, 后跟信号强度
    for (int i = 0;i < N / 2;i++)
    {
        printf("[%llf]: %llf\r\n", i * (double) FS / N, out2[i]);
    }

    // memset(&out[N / 2], 0, sizeof(out) / 2);

    //傅立叶反变换, 输入刚才得到的复数频谱, 输出复数信号, 这里的out2理论上应该和最初的input是一样的信号, 但是幅值被放大了N倍
    Complex_t out3[N];
    ifft(out3, out, N);

    for (int i = 0;i < N;i++)
    {
        out3[i].real /= N;
        out3[i].image /= N;
    }

    uint8_t res = 0;
    //输出原始信号和频谱反变换回来的信号, 这两个信号应该是完全一样的才正确
    for (int i = 0;i < N;i++)
    {
        printf("orig[%d]: %llf+%llfi, ifft[%d]: %llf+%llfi\r\n",
            i, input[i].real, input[i].image,
            i, out3[i].real, out3[i].image
        );

        if (-1 != res &&
                fabs(input[i].real - out3[i].real) > 0.001 ||
                fabs(input[i].image - out3[i].image) > 0.001
            )//误差不超过0.001
            res = -1;//否则认为是还原失败
    }

    return res;
}

