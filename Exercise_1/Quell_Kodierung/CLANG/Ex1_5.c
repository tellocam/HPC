#include <stdio.h>
#include <math.h>
#include <time.h>

typedef unsigned int uint;
uint gray(uint num)
{
    return num ^ (num >> 1); 
}

uint yarg(uint num)
{
    uint mask = num;
    while (mask)
    {
        mask >>= 1;
        num ^= mask;
    }
    return num;
}

uint yarg32(uint num)
{
    num ^= num >> 16;
    num ^= num >> 8;
    num ^= num >> 4;
    num ^= num >> 2;
    num ^= num >> 1;
    return num;
}

int main()
{

    int d = 16;
    clock_t start = clock();
    for (int j = 1; j < pow(2, d); j++)
    {
        yarg(j);
    }
    clock_t res = clock() - start;
    printf("yarg(): %f s\n", (double)res / CLOCKS_PER_SEC);

    start = clock();
    for (int j = 1; j < pow(2, d); j++)
    {
        yarg32(j);
    }
    res = clock() - start;
    printf("yarg32(): %f s\n", (double)res / CLOCKS_PER_SEC);
    return 0;
}