#include <stdio.h>
#include <math.h>

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

    uint old = gray(0);

    for (int j = 1; j < pow(2, d); j++)
    {
        uint ans = gray(j);
        uint diff = old ^ ans;

        int check = 0;

        for (int k = 0; k < d; k++)
        {
            if (diff == (1 << k))
            {
                check++;
            }
        }

        if (check != 1)
        {
            printf("gray() error");
            break;
        }

        if (yarg(ans) != j)
        {
            printf("yarg() error");
            break;
        }
        old = ans;
    
    }
    printf("yarg and gray did not throw errors - Ex1.4 done :-)!\n");
    return 0;
}