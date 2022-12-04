#include <stdio.h>
#include <math.h>
#include <time.h>

typedef unsigned int uint;

// Gray code according to script p. 24
uint gray(uint i)
{
    return i ^ (i >> 1);
}

// inverse Gray code according to Algorithm 2, script p. 24
uint yarg(uint h)
{
    uint i = h;
    h = h >> 1;
    while (h != 0)
    {
        i = i ^ h;
        h = h >> 1;
    }
    return i;
}

uint yargFast(uint i)
{
    i ^= (i >> 16);
    i ^= (i >> 8);
    i ^= (i >> 4);
    i ^= (i >> 2);
    i ^= (i >> 1);
    return i;
}

uint naiveSucc(uint x)
{
    return gray(yargFast(x) + 1);
}

uint naivePred(uint x)
{
    return gray(yargFast(x) - 1);
}

/* https://www.johndcook.com/blog/2020/02/21/popcount/ */
int pop1(unsigned x)
{
    x -= ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0F0F0F0F;
    x += (x >> 8);
    x += (x >> 16);
    return x & 0x0000003F;
}

int parity(uint x)
{
    // int temp = __builtin_popcount(x);
    return pop1(x) % 2;
}

uint succ(uint x)
{
    if (parity(x) == 0)
    {
        return x ^ 1;
    }
    else
    {
        if (x == 1 << (32 - 1))
        {
            return 0;
        }
        // find most right one
        // by looking at James' picture
        // AND x with negative x
        uint y = x & (~x + 1);
        // flip bit left to the most right one
        return x ^ (y << 1);
    }
}

uint pred(uint x)
{

    if (parity(x) == 1)
    {
        return x ^ 1;
    }
    else
    {
        if (x == 0)
        {
            // return 1<<(32-1);
            return 0x80000000;
        }
        uint y = x & (~x + 1);
        return x ^ (y << 1);
    }
}

int main()
{

    uint input = 0;
    printf("\nRESULT: %u\n", yarg(succ(gray(input))));
    printf("\nRESULT: %u\n", yarg(pred(gray(input))));

    // dimension
    int d = 22;

    clock_t start = clock();
    for (int i = 0; i < pow(2, d); i++)
    {
        volatile uint temp = gray(i);
        for (int k = 0; k < 1e4; k++)
        {
            naiveSucc(temp);
            naivePred(temp);
        }
    }
    clock_t stop = clock() - start;

    printf("NAIVE: %f s\n", (double)stop / (CLOCKS_PER_SEC));

    start = clock();
    for (int i = 0; i < pow(2, d); i++)
    {
        volatile uint temp = gray(i);
        for (int k = 0; k < 1e4; k++)
        {
            succ(temp);
            pred(temp);
        }
    }
    stop = clock() - start;

    printf("FAST: %f s\n", (double)stop / (CLOCKS_PER_SEC));

    return 0;
}