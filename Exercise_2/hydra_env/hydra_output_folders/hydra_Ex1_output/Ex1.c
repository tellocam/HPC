/* (C) Manuela Maria Magdalena Raidl Avg. Git Enthusiast und Camilo Tello Fachin aka BananaMan, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#pragma GCC diagnostic ignored "-Wformat-truncation"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
// #include "tuw_alltoall.h"

/* Algorithm selection done at compile time */
// #ifndef ALLTOALL
// #define ALLTOALL Alltoall_fully
// #define ALLTOALL Alltoall_telephone
// #define ALLTOALL Alltoall_factor
// #define ALLTOALL MPI_Alltoall
// #endif

// #define ALLTOALLTAG 7777

// Benchmarking parameters

#define WARMUP 10
#define REPEAT 100
#define MICRO  1000000.0

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

// TODO: probably not needed?
// #define FLIP(_fac,_mul1,_mul2) (_fac = ((_fac==_mul1) ? (_mul2):(_mul1)) )

int MY_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int rank, size;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Reduce(sendbuf, recvbuf, count, datatype, op, 0, comm);
  MPI_Bcast(recvbuf, count, datatype, 0, comm);

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
  int power, count, c, hydra_nodes, gentxt;
  
  tuwtype_t *sendbuf, *recvbuf, *testbuf;

  int i;
  int r, t;
  
  tuwtype_t start, stop;
  tuwtype_t *runtime;
  runtime = (tuwtype_t*)malloc(REPEAT * sizeof(tuwtype_t));
  assert(runtime!=0);


 
  count = 1;
  power = 2;
  hydra_nodes = 1;
  gentxt = 0;
  for (i=1; i<argc&&argv[i][0]=='-'; i++) {
    if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count);
    if (argv[i][1]=='p') i++, sscanf(argv[i], "%d", &power);
    if (argv[i][1]=='h') i++, sscanf(argv[i], "%d", &hydra_nodes);
    if (argv[i][1]=='g') i++, sscanf(argv[i], "%d", &gentxt);
  }

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  FILE *fp; // file pointer
  if(rank==0){
      if (gentxt!=0){
        char file_name[127];
        char file_suffix[64];
        char sizepn_char[8];
        char pow_char[8];
        char uline[8] = "_";
        sprintf(file_suffix, "%d", hydra_nodes);
        snprintf(pow_char, sizeof(pow_char), "%d", power);
        snprintf(sizepn_char, sizeof(sizepn_char), "%d", size/hydra_nodes);
        strcat(file_suffix, uline);
        strcat(file_suffix, sizepn_char);
        strcat(file_suffix, uline);
        strcat(file_suffix, pow_char);
        sprintf(file_name, "EX1_%s.txt", file_suffix);  
        // filename is gonna be like "EX1_<nodes>_<processesPerNode>_<Power>.txt"
        // be careful to always hand over the -h commandline argument when running on multiple nodes, its just gonna be used for .txt file suffix!
        // mpirun -np 8 ./Ex1 -c 50 -p 2 -h 1 -g 1
        fp = fopen(file_name, "w");
        if (fp == NULL) {
        printf("Error opening file!\n");
        }
      }
  }



  // allocate and initialize data for "correctness tests":
  sendbuf = (tuwtype_t*)malloc(count*sizeof(tuwtype_t));
  assert(sendbuf!=NULL);
  recvbuf = (tuwtype_t*)malloc(count*sizeof(tuwtype_t));
  assert(recvbuf!=NULL);
  testbuf = (tuwtype_t*)malloc(count*sizeof(tuwtype_t));
  assert(testbuf!=NULL);
  
  for (i=0; i<count; i++) sendbuf[i] = (tuwtype_t)(rank+i); 
  for (i=0; i<count; i++) recvbuf[i] = (tuwtype_t)-1;
  for (i=0; i<count; i++) testbuf[i] = (tuwtype_t)-1;

  // "correctness test": compare against result from library function
  MPI_Allreduce(sendbuf, testbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
  MY_Allreduce (sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
  for (i=0; i<count; i++) 
  {
    assert(recvbuf[i]==testbuf[i]);
    //printf("Correctness test: %f, %f\n", recvbuf[i], testbuf[i]);
  }

  // prepare headline for benchmarking results
  if (rank==0) {
    fprintf(stderr,"Results for MY_Allreduce(), MPI Processes=%d \n",size); 
    fprintf(stderr,"count, m (Bytes), avg, min, median, stddev,  CIMOE \n"); 
  }

  // TIME MEASURE FOR POWERS OF power , can be command line arg, otherwise its 2!
  for (c = 1; c <= count; c *= power)
  {
    if (c > count) break;
    
    for (r = 0, t = 0; r < WARMUP+REPEAT; r++) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      
      start = MPI_Wtime();
      
      // MY_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
      MY_Allreduce(sendbuf, recvbuf, c, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
      
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
  free(runtime);

  
  return 0;
}
