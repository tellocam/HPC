#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <mpi.h>

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

int MY_Reduce_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node)
{
    int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbufL = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufR = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufL_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufR_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;

    // rank = node (= id)
    if (node == 0) // node == root
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b+j] = tmpbufL[j] > sendbuf[b+j] ? tmpbufL[j] : sendbuf[b+j];
                recvbuf[b+j] = tmpbufR[j] > recvbuf[b+j] ? tmpbufR[j] : recvbuf[b+j];
            }
        }
        if(count % blockSize != 0){
             
            MPI_Recv(tmpbufL_trunc, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(tmpbufR_trunc, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbufL_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbufL_trunc[j] : sendbuf[trunc_chunk_idx +j];
                recvbuf[trunc_chunk_idx +j] = tmpbufR_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbufR_trunc[j] : recvbuf[trunc_chunk_idx +j];
            }
        }
    }
    else if (node >= (size-1)/2) // node == leaf
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
        }
        if(count % blockSize != 0){
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
        }
    }
    else // node == interior
    {
        for (int j = 0; j < blockSize; j++) 
        {
            MPI_Recv(tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[j] = tmpbufL[j] > sendbuf[j] ? tmpbufL[j] : sendbuf[j];
                recvbuf[j] = tmpbufR[j] > recvbuf[j] ? tmpbufR[j] : recvbuf[j];
            }
            MPI_Send(recvbuf + j*blockSize, blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
        }

        // MPI_Recv(tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (int j = 0; j < blockSize; j++)
        // {
        //     recvbuf[j] = tmpbufL[j] > sendbuf[j] ? tmpbufL[j] : sendbuf[j];
        //     recvbuf[j] = tmpbufR[j] > recvbuf[j] ? tmpbufR[j] : recvbuf[j];
        // }
        // for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
        // {
        //     MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, get_parent(node), 0, tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, get_parent(node), 0, tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     for (int j = 0; j < blockSize; j++)
        //     {
        //         recvbuf[b+j] = tmpbufL[j] > sendbuf[b+j] ? tmpbufL[j] : sendbuf[b+j];
        //         recvbuf[b+j] = tmpbufR[j] > recvbuf[b+j] ? tmpbufR[j] : recvbuf[b+j];
        //     }
        // }
        // MPI_Send(recvbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
        if(count%blockSize != 0)
        {
            MPI_Recv(tmpbufL_trunc, count%blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(tmpbufR_trunc, count%blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbufL_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbufL_trunc[j] : sendbuf[trunc_chunk_idx +j];
                recvbuf[trunc_chunk_idx +j] = tmpbufR_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbufR_trunc[j] : recvbuf[trunc_chunk_idx +j];
            }
            MPI_Send(recvbuf + trunc_chunk_idx, count%blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
        }
    }  
    return MPI_SUCCESS;
}

int MY_Bcast_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node)
{
    int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbufL = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufR = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufL_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbufR_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;
    
    if (node == 0) // node == root
    {
        memcpy(&recvbuf[0], sendbuf, count * sizeof(tuwtype_t)); 
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0)
        {
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
        }
    }
    else if (node >= (size-1)/2) // node == leaf
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf+b, blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
    }
    else // node == interior
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf+b, blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
            MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
        }
    } 

    return MPI_SUCCESS;
}

int main(int argc, char *argv[])
{
    int rank, size;
    int count;
    tuwtype_t *sendbuf, *recvbuf, *testbuf, *testbuf1;

    double start, stop;
    double runtime[REPEAT];

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) printf("----- \n Number of processors must be of size 2^n-1. \n----- \n");

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    count = 7;
   
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);
    testbuf1 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf1 != NULL);

    for (int i = 0; i < count; i++)
        // sendbuf[i] = (tuwtype_t)i;
        sendbuf[i] = (tuwtype_t)i*(rank+1);
        // sendbuf[i] = (tuwtype_t)i*(size-rank);
        // sendbuf[i] = (tuwtype_t)((i*(rank+1))%3);
    for (int i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (int i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function
    MY_Reduce_T(sendbuf, testbuf1, count, size, rank);
    MY_Bcast_T(testbuf1, testbuf, count, size, rank);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("node: %d -> %f, %f\n", rank, recvbuf[i], testbuf[i]); // for debugging
    }
 
    // TODO: add benchmarking here

    // start timing
    // MY_Reduce_T(sendbuf, testbuf1, count, size, rank);
    // MY_Bcast_T(testbuf1, testbuf, count, size, rank);
    // end timing
 
    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);
    free(testbuf1);

    return 0;
}
