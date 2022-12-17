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


int reduce_BCast(void *sendbuf, void *recvbuf, int count, int size, int rank, int comm)
{

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
        for (int b = 1; b < count; b++)
        {
            MPI_Sendrecv(sendbuf + b-1, blockSize, MPI_DOUBLE, rank - 1, 0, recvbuf + b,blockSize, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE );
        }
        MPI_Send(sendbuf + count-1, blockSize, MPI_DOUBLE, rank + 1, 0 ,comm);
    }
}

int main(int argc, char *argv[])
{
  printf("ayo");

    return 0;
}
