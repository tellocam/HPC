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
#ifndef ALLTOALL
#define ALLTOALL Alltoall_fully
// #define ALLTOALL Alltoall_telephone
// #define ALLTOALL Alltoall_factor
// #define ALLTOALL MPI_Alltoall
#endif

#define ALLTOALLTAG 7777

// Benchmarking parameters
#define WARMUP 8
#define REPEAT 40
#define MICRO 1000000.0

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

#define FLIP(_fac, _mul1, _mul2) (_fac = ((_fac == _mul1) ? (_mul2) : (_mul1)))

#define blockSize 4

// Apparently we do bcast of maximum.
int MY_Allreduce(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int rank)
{
    int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b+j] = tmpbuf[j] > recvbuf[b+j] ? tmpbuf[j] : recvbuf[b+j];
            }
            MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            
        }
        if(count % blockSize != 0){
             
            MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbuf_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbuf_trunc[j] : recvbuf[trunc_chunk_idx +j];
            }
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
    }
    else if (rank == size -1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        }
        if(count % blockSize != 0){
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        }
    }
    else 
    {
        MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int j = 0; j < blockSize; j++)
        {
            recvbuf[j] = tmpbuf[j] > recvbuf[j] ? tmpbuf[j] : recvbuf[j];
        }

        for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Sendrecv(sendbuf + b - blockSize, blockSize, MPI_DOUBLE, rank - 1, 0, tmpbuf, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b+j] = tmpbuf[j] > recvbuf[b+j] ? tmpbuf[j] : recvbuf[b+j];
            }

        }
        MPI_Send(sendbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        if(count%blockSize != 0){
            MPI_Recv(tmpbuf_trunc, count%blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbuf_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbuf_trunc[j] : recvbuf[trunc_chunk_idx +j];
            }
            MPI_Send(sendbuf + trunc_chunk_idx, count%blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pipelined Bcast [traeff_lecturenotes]
    /*
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0)
        {
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
    } else */
    if (rank == size -1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf+b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
    }
    else 
    {
        // this line can probably be moved aswell
        MPI_Recv(recvbuf, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Sendrecv(sendbuf + b - blockSize, blockSize, MPI_DOUBLE, rank + 1, 0, recvbuf + b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Send(sendbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    
    // warum funktioniert free() nicht???
    // free(tmpbuf);
    // free(tmpbuf_trunc);
    return MPI_SUCCESS;
}




int main(int argc, char *argv[])
{
    int rank, size;
    int count, c;
    tuwtype_t *sendbuf, *recvbuf, *testbuf;

    double start, stop;
    double runtime[REPEAT];

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    count = 10;
   
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
    c = count - 1;
    // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    // MPI_Reduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, 0, MPI_COMM_WORLD);
    // reduce_BCast(sendbuf, testbuf, count, size, rank);
    MY_Allreduce(sendbuf, testbuf, count, size, rank);
    // MY_Bcast_P(sendbuf, testbuf, count, size, rank);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("%f, %f\n", recvbuf[i], testbuf[i]); // for debugging
    }
 
    // int r, t;
    // int f = 5;
    // for (c = 1; c <= count; c *= f, FLIP(f, 2, 5))
    // {
    //     if (c * size > size * count)
    //         break;

    //     for (r = 0, t = 0; r < WARMUP + REPEAT; r++)
    //     {
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         MPI_Barrier(MPI_COMM_WORLD);

    //         start = MPI_Wtime();

    //         // ALLTOALL(sendbuf, c, TUW_TYPE, testbuf, c, TUW_TYPE, MPI_COMM_WORLD);
    //         // int reduce_BCast(void *sendbuf, void *recvbuf, int count, int size, int rank, int comm)
    //         // reduce_BCast(sendbuf, recvbuf, count, size, rank, MPI_COMM_WORLD);
    //         reduce_BCast(sendbuf, testbuf, count, size, rank);


    //         stop = MPI_Wtime();

    //         MPI_Barrier(MPI_COMM_WORLD);
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         if (r < WARMUP)
    //             continue;
    //         runtime[t++] = stop - start;
    //     }
    //     MPI_Allreduce(MPI_IN_PLACE, runtime, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    //     if (rank == 0)
    //     {
    //         double tuwavg, tuwmin;

    //         tuwavg = 0.0;
    //         tuwmin = runtime[0];
    //         for (t = 0; t < REPEAT; t++)
    //         {
    //             tuwavg += runtime[t];
    //             if (runtime[t] < tuwmin)
    //                 tuwmin = runtime[t];
    //         }
    //         tuwavg /= REPEAT;

    //         fprintf(stderr, "%d & %d & %ld & %.2f & %.2f \\\\\n",
    //                 c, c , c * sizeof(tuwtype_t), tuwavg * MICRO, tuwmin * MICRO);
    //     }
    // }
 
    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
