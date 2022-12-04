#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

typedef unsigned int uint;
typedef struct
{
    uint edgeCoord;
    uint cycleCoord;
} node;

int customMod(int a, int b)
{
    return (a % b + b) % b;
}

int print_uint(uint value, int d)
{
    for (int i = d - 1; i >= 0; --i)
    {
        printf("%u", (value & (1 << i)) != 0);
    }
    printf("\n");
}

void printNode(node n, int d)
{
    printf("c: %2u, e: ", n.cycleCoord);
    print_uint(n.edgeCoord, d);
}

void path(node start, node end, uint d)
{
    // compute number of edges we have to traverse in cube (excluding cycles)
    uint xor = start.edgeCoord ^ end.edgeCoord;

    int r = customMod(start.cycleCoord - end.cycleCoord, d);
    int l = customMod(end.cycleCoord - start.cycleCoord, d);

    int count = 0;

    node curr = {start.edgeCoord, start.cycleCoord};

    printNode(curr, d);

    // should we flip the bit at start.cycleCoord?
    if ((xor&(1 << curr.cycleCoord)) != 0)
    {
        // bit flip if 1 is at this position
        curr.edgeCoord ^= (1 << curr.cycleCoord);
        ++count;
        printNode(curr, d);
    }

    int dir = l < r ? -1 : 1;

    for (uint i = 0; i < d - 1; ++i)
    {
        curr.cycleCoord = customMod(curr.cycleCoord + dir, d);
        ++count;
        printNode(curr, d);

        if ((xor&(1 << curr.cycleCoord)) != 0)
        {
            // bit flip if 1 is at this position
            curr.edgeCoord ^= (1 << curr.cycleCoord);
            ++count;
            printNode(curr, d);
        }
    }

    int stepsToEnd = l > r ? r - 1 : l - 1;

    if (stepsToEnd < 0) {
        dir = -dir;
    }

    for (uint i = 0; i < abs(stepsToEnd); ++i)
    {
        curr.cycleCoord = customMod(curr.cycleCoord - dir, d);
        ++count;
        printNode(curr, d);
    }
    printf("%d\n", count);
}

int pathLength(node start, node end, uint d)
{
    // compute number of edges we have to traverse in cube (excluding cycles)
    uint xor = start.edgeCoord ^ end.edgeCoord;

    int r = customMod(start.cycleCoord - end.cycleCoord, d);
    int l = customMod(end.cycleCoord - start.cycleCoord, d);

    int count = 0;

    // counts number of bits in the first d bits
    for (uint i = 0; i < d; ++i)
    {
        if ((xor&(1 << i)) != 0)
        {
            ++count;
        }
    }

    int stepsToEnd = l > r ? r : l;

    return d - 1 + abs(stepsToEnd - 1) + count;
}

#define d 10
#define arrSize 2 * d + d / 2

int main()
{
    // const uint d = 15;
    //  const arrSize = 2 * d + d / 2 - 2;
    //  node start = {0, 0}, end = {0xffffffff, 8};

    // path(start, end, d);

    // asymptotic complexity O(d) -> linear in d
    double sum = 0;

    // array to count path lengths
    int arr[arrSize] = {0};

    for (uint i = 0; i < pow(2, d); ++i)
    {
        if (i % d == 0)
        {
            printf("PROGRESS %.2f \n", 100 * i / (pow(2, d)));
        }
        for (uint k = 0; k < d; ++k)
        {
            for (uint j = 0; j < pow(2, d); ++j)
            {
                for (uint l = 0; l < d; ++l)
                {
                    node start = {i, k};
                    node end = {j, l};
                    int temp = pathLength(start, end, d);
                    sum += temp;
                    arr[temp]++;
                }
            }
        }
    }

    sum /= pow(d * pow(2, d), 2);
    printf("sum: %f\n", sum);

    uint arrSum = 0;

    for (int i = 0; i < arrSize; ++i)
    {
        printf("path length: %d, %d\n", i, arr[i]);
        arrSum += arr[i];
    }

    printf("array sum: %u\n", arrSum);

    node a = {0, 0};
    node b = {0b1111111, 0};
    path(a, b, d);
}