/* (C) Jesper Larsson Traff, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <mpi.h>

// #include <iostream>

/* Algorithm selection done at compile time */
// #ifndef ALLTOALL
// #define ALLTOALL Alltoall_fully
// // #define ALLTOALL Alltoall_telephone
// // #define ALLTOALL Alltoall_factor
// // #define ALLTOALL MPI_Alltoall
// #endif

// #define ALLTOALLTAG 7777

// Benchmarking parameters
#define WARMUP 8
#define REPEAT 40
#define MICRO 1000000.0

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

#define blockSize 4

int get_parent(int id)
{
    return id-1;
}

int get_childL(int id)
{
    return id+1;
}


int MY_Allreduce_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node)
{
    // workbuf is what will be worked on (Y from algorithm)
    // Start: workbuf = sendbuf
    tuwtype_t *workbuf;
    workbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    memcpy(&workbuf[0], sendbuf, count * sizeof(tuwtype_t));  

    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));

    MPI_Status status;
    int recvcount;

    int d = node; // level of node
    int d_max = size-1;
    int b = ceil((double)count/blockSize);
    int sendSize;

    for ( int j = 0; j < b+d_max; j++ )
    {
        if (node < size-1) // node != leaf
        {
            sendSize = count-j*blockSize >= count%blockSize ? blockSize : count%blockSize;
            sendSize = (j-(d+1)<0 || j-(d+1)>=b) ? 0 : sendSize;
            
            MPI_Sendrecv(workbuf + (j-(d+1))*blockSize, sendSize, MPI_DOUBLE, get_childL(node), j, tmpbuf, blockSize, MPI_DOUBLE, get_childL(node), j, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            for (int i = 0; i < recvcount; i++)
            {
                workbuf[j*blockSize+i] = tmpbuf[i] > workbuf[j*blockSize+i] ? tmpbuf[i] : workbuf[j*blockSize+i];
            }
            printf("tochild - rank: %d, j: %d, recv: %d, recvct: %d, send: %d, sendct: %d\n", node, j, j*blockSize, recvcount, (j-(d+1))*blockSize, sendSize);

        }
        if (node != 0) // node != root
        {
            sendSize = count-j*blockSize >= blockSize ? blockSize : count%blockSize;
            sendSize = (j<0 || j>=b) ? 0 : sendSize;
            MPI_Sendrecv(workbuf + j*blockSize, sendSize, MPI_DOUBLE, get_parent(node), j, tmpbuf, blockSize, MPI_DOUBLE, get_parent(node), j, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            if (j-d >= 0 && j-d < b && recvcount > 0)
            {
                memcpy(&workbuf[(j-d)*blockSize], tmpbuf, sizeof(tuwtype_t)*recvcount);
            }
            printf("toparnt - rank: %d, j: %d, recv: %d, recvct: %d, send: %d, sendct: %d\n", node, j, (j-d)*blockSize, recvcount, j*blockSize, sendSize);
        }
    }

    // End:  recvbuf = workbuf
    memcpy(&recvbuf[0], workbuf, count * sizeof(tuwtype_t));  
    free(workbuf);
    return MPI_SUCCESS;
}

int main(int argc, char *argv[])
{
    int rank, size;
    int count;
    tuwtype_t *sendbuf, *recvbuf, *testbuf;

    double start, stop;
    double runtime[REPEAT];

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) printf("----- \n Number of processors must be of size 2^n-1. \n----- \n");

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    count = 9;
   
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);

    for (int i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)i*(rank+1);
    for (int i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (int i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function
    // MY_Reduce_T(sendbuf, testbuf, count, size, rank);
    // MY_Bcast_T(testbuf, testbuf, count, size, rank);
    MY_Allreduce_T(sendbuf, testbuf, count, size, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("rank: %d -> %f, %f\n", rank, recvbuf[i], testbuf[i]); // for debugging
    }
 
    // TODO: add benchmarking here
 
    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
