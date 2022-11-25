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
    //蝶形运算所用到的变量有:
    //binLen - 总共需要进行几层运算(也就是蝶形图上有多少个列)
    //level - 当前进行的是第几层运算
    //group - 本层运算分为几组(也就是蝶形图中比较靠近的交点)
    //cg - current group, 也就是当前是第几组
    //inGroup - 每组内有多少个蝶形运算
    //cig - 当前是组内的第几个运算
    //step - 每个蝶形运算里, 参与运算的A和B之间的下标差
    //numOfFactors - 本层运算内包括多少个旋转因子, 也就是多少个不同的W(p ,N), W和N都是固定的, 其实就是多少个不同的p
    //cf - current factor, 当前是第几个旋转因子
    //p - W(p, N)里的p, N就是len, W是一个表达式, 带入p和N之后可以得到一个复数
    //numOfSameFactors - 相同的旋转因子出现了多少次, 比如第一层实际上只有一个旋转因子, 它就出现了group次
    for (size_t level = 1;level <= binLen;level++)//迭代计算所有的层
    {
        //以下规律通过观察蝶形图总结得到
        size_t group = 1 << (binLen - level);//每层内的组数是2^(总层数-当前层数)
        size_t inGroup = 1 << (level - 1);//每组内的运算数量是2^当前层数-1
        size_t step = inGroup;//取数的间隔和组内运算数量一致
        size_t numOfFactors = inGroup;//不同种类的旋转因子的数量和组内运算数量一致
        size_t numOfSameFactors = group;//同种类的旋转因子的数量和组的数量一致

        for (size_t cf = 0;cf < numOfFactors;cf++)//遍历使用所有不同的旋转因子参与计算
        {
            //这里通过观察蝶形图可以发现, p的变化规律就是每次增长group
            //前面也说过, 不同的旋转因子也就是不同的p
            size_t p = cf * group;
            Number_t omega = PI_2 * p / len;//这里就是2pi*p/N, 也就是准备传递给三角函数的参数omega

            //W(p, N)是指数表达式的简写, 也就是e^-j*omega
            //指数表达式可以用欧拉公式变为一般形式, 也就是cos(omega) - j*sin(omega), j不是数字而是虚数
            //所以所谓的W(p, N)也就是一个实部为cos(omega), 虚部为-sin(omega)的复数而已
            Complex_t w = { cos(omega), -sin(omega) };//计算W(p, N)

            //将当前旋转因子在本层内所需要参与的运算全部算出来, 这里就是用同一个w去取不同的A和B取计算
            for (size_t csf = 0; csf < numOfSameFactors; csf++)
            {
                //计算本次运算取的第一个数的下标
                size_t indexA = step * 2 * csf + cf;
                pComplex_t pA = &output[indexA], pB = &output[indexA + step];

                //蝶形运算是三个数参与算出两个数, 表达式为:
                //A(n) = A(n-1) + B(n-1) * W(p, N)
                //B(n) = A(n-1) + B(n-1) * W(p, N)
                //最终要计算的就是A(n)和B(n), W(p, N)上面已经算好了, 而A(n-1)和B(n-1)是上次算出来的值, 也是已经算好的
                //这里先算B(n-1) * W(p, N), 这是一个复数乘法, 所以要算实部tr和虚部ti
                //复数乘法的计算为B*W = (B.real*W.real - B.image*W.image) + (B.real*W.image + B.image*W.real)i
                Number_t tr = pB->real * w.real - pB->image * w.image;
                Number_t ti = pB->real * w.image + pB->image * w.real;

                //上面算的是B*W(p, N), 所以这里还得把A加进来
                pB->real = pA->real - tr;
                pB->image = pA->image - ti;//先B后A是因为A不能被覆盖, B公式里面的A代表的是本次的A而不是下一次的A
                pA->real = pA->real + tr;
                pA->image = pA->image + ti;
            }
        }
    }
}

void ifft(Number_t* output, Complex_t* input, size_t len)
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
        Complex_t tmp = input[i];
        input[i] = input[rIndex];
        input[rIndex] = tmp;
    }

    //进行蝶形运算
    //蝶形运算所用到的变量有:
    //binLen - 总共需要进行几层运算(也就是蝶形图上有多少个列)
    //level - 当前进行的是第几层运算
    //group - 本层运算分为几组(也就是蝶形图中比较靠近的交点)
    //cg - current group, 也就是当前是第几组
    //inGroup - 每组内有多少个蝶形运算
    //cig - 当前是组内的第几个运算
    //step - 每个蝶形运算里, 参与运算的A和B之间的下标差
    //numOfFactors - 本层运算内包括多少个旋转因子, 也就是多少个不同的W(p ,N), W和N都是固定的, 其实就是多少个不同的p
    //cf - current factor, 当前是第几个旋转因子
    //p - W(p, N)里的p, N就是len, W是一个表达式, 带入p和N之后可以得到一个复数
    //numOfSameFactors - 相同的旋转因子出现了多少次, 比如第一层实际上只有一个旋转因子, 它就出现了group次
    for (size_t level = 1;level <= binLen;level++)//迭代计算所有的层
    {
        //以下规律通过观察蝶形图总结得到
        size_t group = 1 << (binLen - level);//每层内的组数是2^(总层数-当前层数)
        size_t inGroup = 1 << (level - 1);//每组内的运算数量是2^当前层数-1
        size_t step = inGroup;//取数的间隔和组内运算数量一致
        size_t numOfFactors = inGroup;//不同种类的旋转因子的数量和组内运算数量一致
        size_t numOfSameFactors = group;//同种类的旋转因子的数量和组的数量一致

        for (size_t cf = 0;cf < numOfFactors;cf++)//遍历使用所有不同的旋转因子参与计算
        {
            //这里通过观察蝶形图可以发现, p的变化规律就是每次增长group
            //前面也说过, 不同的旋转因子也就是不同的p
            size_t p = cf * group;
            Number_t omega = PI_2 * p / len;//这里就是2pi*p/N, 也就是准备传递给三角函数的参数omega

            //W(-p, N)是指数表达式的简写, 也就是e^j*omega
            //指数表达式可以用欧拉公式变为一般形式, 也就是cos(omega) + j*sin(omega), j不是数字而是虚数
            //所以所谓的W(p, N)也就是一个实部为cos(omega), 虚部为sin(omega)的复数而已
            Complex_t w = { cos(omega), sin(omega) };//计算W(p, N)

            //将当前旋转因子在本层内所需要参与的运算全部算出来, 这里就是用同一个w去取不同的A和B取计算
            for (size_t csf = 0; csf < numOfSameFactors; csf++)
            {
                //计算本次运算取的第一个数的下标
                size_t indexA = step * 2 * csf + cf;
                pComplex_t pA = &input[indexA], pB = &input[indexA + step];

                //蝶形运算是三个数参与算出两个数, 表达式为:
                //A(n) = A(n-1) + B(n-1) * W(p, N)
                //B(n) = A(n-1) + B(n-1) * W(p, N)
                //最终要计算的就是A(n)和B(n), W(p, N)上面已经算好了, 而A(n-1)和B(n-1)是上次算出来的值, 也是已经算好的
                //这里先算B(n-1) * W(p, N), 这是一个复数乘法, 所以要算实部tr和虚部ti
                //复数乘法的计算为B*W = (B.real*W.real - B.image*W.image) + (B.real*W.image + B.image*W.real)i
                Number_t tr = pB->real * w.real - pB->image * w.image;
                Number_t ti = pB->real * w.image + pB->image * w.real;

                //上面算的是B*W(p, N), 所以这里还得把A加进来
                pB->real = pA->real - tr;
                pB->image = pA->image - ti;//先B后A是因为A不能被覆盖, B公式里面的A代表的是本次的A而不是下一次的A
                pA->real = pA->real + tr;
                pA->image = pA->image + ti;
            }
        }
    }
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
