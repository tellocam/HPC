/* (C) Jesper Larsson Traff, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <mpi.h>


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
    return floor((id-1)/2);
}

int get_childL(int id)
{
    return 2*id+1;
}

int get_childR(int id)
{
    return 2*id+2;
}




int MY_Allreduce_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node)
{
    // int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    // tuwtype_t *zerobuf = (tuwtype_t *)malloc(0); // for sending zero in case no valid information is available
    // MPI_getcount or MPI_getelements -> use this to understand wehen the proc is done 
    // tuwtype_t *tmpbufL_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    // tuwtype_t *tmpbufR_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    // tuwtype_t tmp;

    MPI_Status status;
    int recvcount;

    int d = floor(log2(node+1)); // level of node
    int d_max = floor(log2(size+1));
    int b = ceil(count/blockSize);
    int sendSize;

    for ( int j = 0; j < b+d_max; j++ )
    {
        if (node < (size-1)/2) // node != leaf
        {
            sendSize = (j-(d+1)<0 || j-(d+1)>=b) ? 0 : blockSize;        
            MPI_Sendrecv(recvbuf + (j-(d+1))*blockSize, sendSize, MPI_DOUBLE, get_childL(node), 0, tmpbuf, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            for (int i = 0; i < recvcount; i++)
            {
                recvbuf[j*blockSize+i] = tmpbuf[i] > recvbuf[(j-(d+1))*blockSize+i] ? tmpbuf[i] : recvbuf[(j-(d+1))*blockSize+i];
            }
            MPI_Sendrecv(recvbuf + (j-(d+1))*blockSize, sendSize, MPI_DOUBLE, get_childR(node), 0, tmpbuf, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            for (int i = 0; i < recvcount; i++)
            {
                recvbuf[j*blockSize+i] = tmpbuf[i] > recvbuf[(j-(d+1))*blockSize+i] ? tmpbuf[i] : recvbuf[(j-(d+1))*blockSize+i];
            }     
        }
        if (node != 0) // node != root
        {
            sendSize = (j<0 || j>=b) ? 0 : blockSize;
            MPI_Sendrecv(recvbuf + j*blockSize, sendSize, MPI_DOUBLE, get_parent(node), 0, tmpbuf, blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            if (j-d >= 0 && j-d < b && recvcount > 0)
            {
                memcpy(&recvbuf[(j-d)*blockSize], tmpbuf, sizeof(tuwtype_t)*recvcount);
            }
        }
        // printf("j = %d\n", j);
    }
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

    count = 6;
   
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);

    for (int i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)i;
    for (int i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (int i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function
    // MY_Reduce_T(sendbuf, testbuf, count, size, rank);
    // MY_Bcast_T(testbuf, testbuf, count, size, rank);
    MY_Allreduce_T(sendbuf, sendbuf, count, size, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < count; i++) {
        // assert(recvbuf[i] == sendbuf[i]);
        printf("%f, %f\n", recvbuf[i], sendbuf[i]); // for debugging
    }
 
    // TODO: add benchmarking here
 
    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
