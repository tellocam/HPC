#include <stdio.h>
#include <math.h>

typedef unsigned int uint;

// This function converts an unsigned binary number to reflected binary Gray code.
uint gray(uint num)
{
    return num ^ (num >> 1); // The operator >> is shift right. The operator ^ is exclusive or.
}

// This function converts a reflected binary Gray code number to a binary number.
uint yarg(uint num)
{
    uint mask = num;
    while (mask)
    { // Each Gray code bit is exclusive-ored with all more significant bits.
        mask >>= 1;
        num ^= mask;
    }
    return num;
}

// A more efficient version for Gray codes 32 bits or fewer through the use of SWAR (SIMD within a register) techniques.
// It implements a parallel prefix XOR function. The assignment statements can be in any order.
//
// This function can be adapted for longer Gray codes by adding steps.

uint yarg32(uint num)
{
    num ^= num >> 16;
    num ^= num >> 8;
    num ^= num >> 4;
    num ^= num >> 2;
    num ^= num >> 1;
    return num;
}
// A Four-bit-at-once variant changes a binary number (abcd)2 to (abcd)2 ^ (00ab)2, then to (abcd)2 ^ (00ab)2 ^ (0abc)2 ^ (000a)2.

int main()
{
    // printf("\n%u\n", gray(3));
    // printf("\n%u\n", yarg(0b111));

    int d = 20;

    uint old = gray(0);

    for (int j = 1; j < pow(2, d); j++)
    {
        uint ans = gray(j);
        // printf("gray(): %u\n", ans);
        uint diff = old ^ ans; // xor between old and ans (1 where old and ans are different - 0 otherwise)
        // printf("diff: %u\n", diff);
        int check = 0; // boolean

        for (int k = 0; k < d; k++)
        {
            // shift 1 by k
            if (diff == (1 << k))
            {
                check++;
            }
        }

        if (check != 1)
        {
            printf("ERROR gray()");
            break;
        }

        if (yarg(ans) != j)
        {
            printf("ERROR yarg()");
            break;
        }

        old = ans;
    }

    return 0;
}