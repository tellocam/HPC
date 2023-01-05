/* (C) Jesper Larsson Traff, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <mpi.h>


// /* Algorithm selection done at compile time */
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

//#define blockSize 4

#define indent

int MY_Allreduce(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int rank, int blockSize)
{
    int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count%blockSize * sizeof(tuwtype_t));
    int bctag = 1;
    //tuwtype_t tmp;
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b+j] = tmpbuf[j] > sendbuf[b+j] ? tmpbuf[j] : sendbuf[b+j];
            }
            // from bcast, should change a bit for small number of communication rounds
            MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, 1, bctag, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0){
             
            MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbuf_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbuf_trunc[j] : sendbuf[trunc_chunk_idx +j];
            }
            // from bcast, shouldn't change much
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, 1, bctag, MPI_COMM_WORLD);
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
            recvbuf[j] = tmpbuf[j] > sendbuf[j] ? tmpbuf[j] : sendbuf[j];
        }
        for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, rank - 1, 0, tmpbuf, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b+j] = tmpbuf[j] > sendbuf[b+j] ? tmpbuf[j] : sendbuf[b+j];
            }
        }
        MPI_Send(recvbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        if(count%blockSize != 0){
            MPI_Recv(tmpbuf_trunc, count%blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx +j] = tmpbuf_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbuf_trunc[j] : sendbuf[trunc_chunk_idx +j];
            }
            MPI_Send(recvbuf + trunc_chunk_idx, count%blockSize, MPI_DOUBLE, rank - 1, 0 , MPI_COMM_WORLD);
        }
        // from bcast, shouldn't rly matter
        MPI_Recv(recvbuf, blockSize, MPI_DOUBLE, rank - 1, bctag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }  

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pipelined Bcast [traeff_lecturenotes]
    /*
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, 1, bctag, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0)
        {
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, 1, bctag, MPI_COMM_WORLD);
        }
    }
    else */
    // this part can stay
    if (rank == size -1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf+b, blockSize, MPI_DOUBLE, rank - 1, bctag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, bctag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
    }
    else 
    {
        for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, rank + 1, bctag, recvbuf + b, blockSize, MPI_DOUBLE, rank - 1, bctag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Send(recvbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, rank + 1, bctag, MPI_COMM_WORLD);
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, bctag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank + 1, bctag, MPI_COMM_WORLD);
        }
    }
    return MPI_SUCCESS;
}

// qsort helper function to find median in compute_stats
int cmp_tuwtype(const void *a, const void *b) {
  tuwtype_t da = *(tuwtype_t *)a;
  tuwtype_t db = *(tuwtype_t *)b;
  return (da > db) - (da < db);
}

// Find median function
tuwtype_t find_median(tuwtype_t *runtime, size_t N){
  tuwtype_t *runtime_copy = malloc(N * sizeof(tuwtype_t));
  memcpy(runtime_copy, runtime, N * sizeof(tuwtype_t));
  // sort runtime_copy
  qsort(runtime_copy, N, sizeof(tuwtype_t), cmp_tuwtype);
  // compute median
    if (N % 2 == 1) {
       return runtime_copy[N / 2];
    } else {
       return (runtime_copy[N / 2 - 1] + runtime_copy[N / 2]) / 2;
      }
  free(runtime_copy);
}

tuwtype_t find_stddev(tuwtype_t *runtime, size_t N, tuwtype_t avg){
  tuwtype_t var = 0;
  for(size_t i = 0; i < N; i++){
    var += (runtime[i] - avg) * (runtime[i] - avg);
  }
  var /= N-1;
  return sqrt(var);
}

tuwtype_t find_CI_MOR(tuwtype_t stddev, size_t N){
  tuwtype_t t = 1.96; // t-value for 95% confidence interval with N-1 degrees of freedom -> wikipedia :-)
  return t * stddev / sqrt(N);
}

int main(int argc, char *argv[])
{
    int rank, size;
    int power, count, c, gentxt, blockSize;
    tuwtype_t *sendbuf, *recvbuf, *testbuf;
    tuwtype_t start, stop;
    int i, r, t;
    tuwtype_t *runtime;
    runtime = (tuwtype_t*)malloc(REPEAT * sizeof(tuwtype_t));


    runtime = (tuwtype_t*)malloc(REPEAT * sizeof(tuwtype_t));
    assert(runtime!=0);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    count = 10;
    power = 2;
    gentxt = 0;
    blockSize = 4;
    for (i=1; i<argc&&argv[i][0]=='-'; i++) {
    if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count); // commandline arg -c for adjusting max. count, if none given count = 10
    if (argv[i][1]=='p') i++, sscanf(argv[i], "%d", &power); // commandline arg -p for adjusting the power, if none given power = 2
    if (argv[i][1]=='b') i++, sscanf(argv[i], "%d", &blockSize); // commandline arg. -b for adjusting blocksize, if none given blockSize = 4
    if (argv[i][1]=='g') i++, sscanf(argv[i], "%d", &gentxt); // commandline arg. -g for generating a txt, if none given, no .txt
    }

    /* 
    mpirun -np 8 ./Ex2 -c 50 -p 2 -b 4 -g 1
    here a commandline example which executes the binary with 8 MPI Processes, max. count of 50, iterators raised to the power of 2
    pipelining blocksize of 2 and will generate a .txt file with commas as value separators.
    */

    FILE *fp; // file pointer
    if(rank==0){

    /*
    If Commandline argument -g is not 0, this will generate a txt file of following form: EX2_8_2_4.txt
    EX2 obviously stands for exercise 2, the first number after that 8 in this case, is the number
    of MPI Processes. The second number, 2 in this case, is the power to which every iteration number is raised.
    The last number, 4 in this case, is the blockSize which is the size of the pipelining packages!
    */

    if (gentxt!=0){
        char file_name[127];
        char file_suffix[64];
        char pow_char[8];
        char blockSize_char[8];
        char uline[8] = "_";
        sprintf(file_suffix, "%X", size);
        snprintf(pow_char, sizeof(pow_char), "%d", power);
        snprintf(blockSize_char, sizeof(blockSize_char), "%d", blockSize);
        strcat(file_suffix, uline);
        strcat(file_suffix, pow_char);
        strcat(file_suffix, uline);
        strcat(file_suffix, blockSize_char);
        sprintf(file_name, "EX2_%s.txt", file_suffix);
        fp = fopen(file_name, "w");
        if (fp == NULL) {
        printf("Error opening file!\n");
        }
    }
  }
    // allocate and initialize data for "correctness tests":
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
    // MY_Reduce_P(sendbuf, testbuf, count, size, rank, blockSize);
    // MY_Bcast_P(testbuf, testbuf, count, size, rank, blockSize);

    MY_Allreduce(sendbuf, testbuf, count, size, rank, blockSize);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("%f, %f\n", recvbuf[i], testbuf[i]); // for debugging
    }
 
  // prepare headline for benchmarking results
  if (rank==0) {
    fprintf(stderr,"Results for MY_Allreduce(), MPI Processes=%d \n",size); 
    fprintf(stderr,"count, m (Bytes), avg, min, median, stddev,  CIMOR \n"); 
  }

  // TIME MEASURE FOR POWERS OF power , can be command line arg, otherwise its 2!
  for (c = 1; c <= count; c *= power)
  {
    if (c > count) break;
    
    for (r = 0, t = 0; r < WARMUP+REPEAT; r++) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
      
    start = MPI_Wtime();
      
    // MY_Reduce_P(sendbuf, testbuf, count, size, rank, blockSize);
    // MY_Bcast_P(testbuf, testbuf, count, size, rank, blockSize);

    MY_Allreduce(sendbuf, testbuf, count, size, rank, blockSize);
      
    stop = MPI_Wtime();
      
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (r < WARMUP) continue;
    runtime[t++] = stop-start;
    }
    MPI_Allreduce(MPI_IN_PLACE, runtime, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    

    // compute mean, minimum (slowest process)
    // compute median, CI, standard deviation
    if (rank==0) {
      tuwtype_t tuwavg, tuwmin;
      tuwtype_t tuwmed, tuwstddev, CI_MOR;
      
      tuwavg = 0.0; 
      tuwmin = runtime[0]; 
      for (t = 0; t < REPEAT; t++) {
        tuwavg += runtime[t];
        if (runtime[t]<tuwmin) tuwmin = runtime[t];
      }
      tuwavg /= REPEAT;
      tuwmed = find_median(runtime, (size_t)REPEAT);
      tuwstddev = find_stddev(runtime, (size_t)REPEAT, tuwavg);
      CI_MOR = find_CI_MOR(tuwstddev, (size_t)REPEAT);

      fprintf(stderr, "%d, %ld, %.2f, %.2f, %.2f, %.2f, %.2f \n", 
	      c, c*sizeof(tuwtype_t), tuwavg*MICRO, tuwmin*MICRO, tuwmed*MICRO, tuwstddev*MICRO, CI_MOR*MICRO);
      if (gentxt!=0){
        fprintf(fp, "%d, %ld, %.2f, %.2f, %.2f, %.2f, %.2f \n",
          c, c*sizeof(tuwtype_t), tuwavg*MICRO, tuwmin*MICRO, tuwmed*MICRO, tuwstddev*MICRO, CI_MOR*MICRO);
      }
      
    }
  
  }

  if(gentxt!=0){
    if(rank==0){
      fclose(fp);
    }
  }
 
    MPI_Finalize();
    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
