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

#define blockSize 1


int reduce_BCast(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int rank, int comm)
{
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;
    if (rank == 0)
    {
        for (int b = 0; b < count; b++)
        {
            //MPI_Recv(tmpbuf,blockSize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);
            MPI_Send(sendbuf, blockSize, MPI_DOUBLE, 1, 0 ,comm);
            tmp = tmpbuf[0];
            for(int i = 1; i < blockSize; i++) {
                tmp = tmpbuf[i] > tmp ? tmpbuf[i] : tmp;
            }
            recvbuf[b] = tmp > recvbuf[b] ? tmp : recvbuf[b];
        }
    }
    else if (rank == size -1){
        for (int b = 0; b < count; b++)
        {
            MPI_Send(sendbuf, blockSize, MPI_DOUBLE, 0, 0 ,comm);
        }
    }
    else 
    {
        MPI_Recv(tmpbuf,blockSize, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE );
        tmp = tmpbuf[0];
        for(int i = 1; i < blockSize; i++) {
            tmp = tmpbuf[i] > tmp ? tmpbuf[i] : tmp;
        }
        recvbuf[0] = tmp > recvbuf[0] ? tmp : recvbuf[0];
        for (int b = 1; b < count; b++)
        {
            MPI_Sendrecv(sendbuf + b-1, blockSize, MPI_DOUBLE, rank + 1, 0, tmpbuf,blockSize, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE);
            tmp = tmpbuf[0];
            for(int i = 1; i < blockSize; i++) {
              tmp = tmpbuf[i] > tmp ? tmpbuf[i] : tmp;
            }
            recvbuf[b] = tmp > recvbuf[b] ? tmp : recvbuf[b];
        }
        MPI_Send(sendbuf + count-1, blockSize, MPI_DOUBLE, rank + 1, 0 ,comm);
    }
    printf("at bcast%d\n",rank);
    for(int i = 0; i < count; i++) {
        printf("%f",recvbuf[i]);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pipelined Bcast [traeff_lecturenotes]
    if (rank == 0)
    {
        for (int b = 0; b < count; b++)
        {
            MPI_Send(sendbuf+b*blockSize, blockSize, MPI_DOUBLE, rank + 1, 0 ,comm);
        }
    }
    else if (rank == size -1){
        for (int b = 0; b < count; b++)
        {
            MPI_Recv(recvbuf+b,blockSize, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE );
        }
    }
    else 
    // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
    // int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //              int dest, int sendtag,
    //              void *recvbuf, int recvcount, MPI_Datatype recvtype,
    //              int source, int recvtag, MPI_Comm comm, MPI_Status * status)

    {
        MPI_Recv(recvbuf,blockSize, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE );
        // printf("%d\n",rank);
        for (int b = 1; b < count; b++)
        {
            // printf("before sendrc%d\n",rank);
            MPI_Sendrecv(sendbuf + b-1, blockSize, MPI_DOUBLE, rank + 1, 0, recvbuf + b,blockSize, MPI_DOUBLE, rank - 1, 0, comm, MPI_STATUS_IGNORE);
            // printf("after sendrc%d\n",rank);
        }
        MPI_Send(sendbuf + count-1, blockSize, MPI_DOUBLE, rank + 1, 0 ,comm);
    }
    return MPI_SUCCESS;
}

int main(int argc, char *argv[])
{
    int rank, size;
    int count, c;

    tuwtype_t *sendbuf, *recvbuf, *testbuf;

    int i;

    int r, t;

    double start, stop;
    double runtime[REPEAT];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    count = 4;

   
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);

    for (i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)i;
    for (i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function
    c = count - 1;
    // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    //MPI_Reduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, 0, MPI_COMM_WORLD);
    reduce_BCast(sendbuf, testbuf, count, size, rank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("%f, %f\n", recvbuf[i], testbuf[i]);
    }
/* 
    int f = 5;
    for (c = 1; c <= count; c *= f, FLIP(f, 2, 5))
    {
        if (c * size > size * count)
            break;

        for (r = 0, t = 0; r < WARMUP + REPEAT; r++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            start = MPI_Wtime();

            // ALLTOALL(sendbuf, c, TUW_TYPE, testbuf, c, TUW_TYPE, MPI_COMM_WORLD);
            // int reduce_BCast(void *sendbuf, void *recvbuf, int count, int size, int rank, int comm)
            // reduce_BCast(sendbuf, recvbuf, count, size, rank, MPI_COMM_WORLD);


            stop = MPI_Wtime();

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (r < WARMUP)
                continue;
            runtime[t++] = stop - start;
        }
        MPI_Allreduce(MPI_IN_PLACE, runtime, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (rank == 0)
        {
            double tuwavg, tuwmin;

            tuwavg = 0.0;
            tuwmin = runtime[0];
            for (t = 0; t < REPEAT; t++)
            {
                tuwavg += runtime[t];
                if (runtime[t] < tuwmin)
                    tuwmin = runtime[t];
            }
            tuwavg /= REPEAT;

            fprintf(stderr, "%d & %d & %ld & %.2f & %.2f \\\\\n",
                    c, c , c * sizeof(tuwtype_t), tuwavg * MICRO, tuwmin * MICRO);
        }
    }
 */
    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
