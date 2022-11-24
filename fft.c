/*
@file: fft.c
@author: ZZH
@date: 2022-11-05
@info: fft实现
*/
#include "fft.h"
#include <stdio.h>
#include <stdlib.h>

const double PI_2 = 2 * M_PI;

//二分法找到第一个为1的二进制位, 根据size_t的大小可做到32位64位自适应
static __attribute__((pure)) uint8_t ff1(size_t num)
{
    if (0 == num)
        return 0;

    uint8_t bitsOfVar = sizeof(num) * 8;//通过sizeof实现自适应
    uint8_t index = 1;

    //这里执行的时候变量的位数已经固定了, 所以实际上是O(1), 32位固定循环5次, 64位固定循环6次
    while (bitsOfVar)
    {
        bitsOfVar /= 2;
        size_t tmp = num >> bitsOfVar;
        if (tmp)
        {
            index += bitsOfVar;
            num = tmp;
        }
    }

    return index;
}

//位倒序
static __attribute__((pure)) size_t reverseBits(size_t num, uint8_t bitsOfVar)
{
    if (0 == num || bitsOfVar > (sizeof(num) * 8))
        return 0;

    const size_t topBit = 1 << (bitsOfVar - 1);
    bitsOfVar >>= 1;

    for (uint8_t i = 0;i < bitsOfVar;i++)
    {
        size_t lowMask = 1 << i;
        size_t highMask = topBit >> i;
        uint8_t lowBit = (num & lowMask) > 0;
        uint8_t highBit = (num & highMask) > 0;

        num &= ~(highMask | lowMask);
        num |= (highMask * lowBit) | (lowMask * highBit);
    }

    return num;
}

void fft(Complex_t* output, Number_t* input, size_t len)
{
    if (len < 2)
        return;

    uint8_t binLen = ff1(len) - 1;//首位二进制1的下标相当于求某个数的二进制长度, 例如1024的二进制长度是10
    if (len & (len - 1))//计算长度如果不是2的N次幂
        len = 1 << binLen;//求下一级的变换, 例如1023长度这里会当做512点数去计算, 1025会当做1024点数

    //重排序
    for (size_t i = 0;i < len;i++)
    {
        size_t rIndex = reverseBits(i, binLen);
        //交换顺序的同时输出到output, 后续在output上原地进行蝶形运算
        output[rIndex].real = input[i];
        output[rIndex].image = 0;
        output[i].real = input[rIndex];
        output[i].image = 0;
    }

    //进行蝶形运算
    //W(p, N) = e^-i*PI_2*p/N = cos(PI_2*p/N) - i*sin(PI_2*p/N) (i为虚数单位, N是变换点数, 相当于len变量)
    //蝶形运算为 A=A+B*W(p, N) B=A-B*W(p, N) (A,B 均为复数,且等号左侧的A,B为下次迭代的值, 等号右侧的为本次的值)
    for (size_t i = 0;i < binLen;i++)
    {
        size_t step = 1 << i;//数据间隔, 每层蝶形运算取数的间隔不同

        //相同旋转因子的数量和数据间隔是一致的, 比如第一层全是W(0, N), 因此相同旋转因子只有1种, 而第一层取数的间隔也是1 
        for (size_t j = 0;j < step;j++)//循环计算所有相同旋转因子
        {
            size_t k = 1 << binLen - i - 1;//旋转因子增量, 比如 W(0, N) 到 W(2, N) 增量为2
            size_t p = j * k;//第j个旋转指数

            Number_t arg = PI_2 * p / len;
            Number_t sv = sin(arg);
            Number_t cv = cos(arg);

            for (size_t z = 0;z < k;z++)//
            {
                size_t wStep = j + 2 * step * z;
                //参与蝶形运算的实际就三个数, A B W, 其中A和B都是上一次迭代的输出值, W是旋转因子, 也就是那个展开是sin cos的公式
                pComplex_t pA = &output[wStep], pB = pA + step;

                //这个就相当于是B*W(p, N), tr结果的实数部分, ti是结果的虚数部分
                Number_t tr = pB->real * cv + pB->image * sv;
                Number_t ti = pB->image * cv - pB->real * sv;

                //蝶形运算是A=A+B*WpN B=A-B*WpN,所以这里还得把A加进来
                pB->real = pA->real - tr;
                pB->image = pA->image - ti;
                pA->real = pA->real + tr;
                pA->image = pA->image + ti;
            }
        }
    }
}

void ifft(Number_t* output, Complex_t* input, size_t len)
{
}

void dft(Complex_t* output, Number_t* input, size_t len)
{
    Number_t tmp = PI_2 / len;
    for (size_t i = 0;i < len;i++)
    {
        Number_t real = 0, image = 0;
        Number_t w = tmp * i;
        for (size_t j = 0;j < len;j++)
        {
            Number_t arg = j * w;
            real += input[j] * cos(arg);
            image -= input[j] * sin(arg);
        }
        output[i].real = real;
        output[i].image = image;
    }
}

void idft(Number_t* output, Complex_t* input, size_t len)
{
    Number_t tmp = PI_2 / len;
    for (size_t i = 0;i < len;i++)
    {
        Number_t res = 0;
        Number_t w = tmp * i;
        for (size_t j = 0;j < len;j++)
        {
            Number_t arg = j * w;
            res += input[j].real * cos(arg) - input[j].image * sin(arg);
        }
        output[i] = res;
    }
}

void length(Number_t* output, Complex_t* input, size_t len)
{
    for (size_t i = 0;i < len;i++)
        output[i] = sqrt(input[i].real * input[i].real + input[i].image * input[i].image);
}
