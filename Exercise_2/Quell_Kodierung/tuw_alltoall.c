/* (C) Jesper Larsson Traff, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <mpi.h>

#include "tuw_alltoall.h"

/* Algorithm selection done at compile time */
#ifndef ALLTOALL
#define ALLTOALL Alltoall_fully
//#define ALLTOALL Alltoall_telephone
//#define ALLTOALL Alltoall_factor
//#define ALLTOALL MPI_Alltoall
#endif

#define ALLTOALLTAG 7777

// Benchmarking parameters
#define WARMUP 8
#define REPEAT 40
#define MICRO  1000000.0

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

#define FLIP(_fac,_mul1,_mul2) (_fac = ((_fac==_mul1) ? (_mul2):(_mul1)) )
 
int Alltoall_fully(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype,
		   MPI_Comm comm)
{
  int rank, size;
  int i;

  MPI_Aint lb, sendextent, recvextent;
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Type_get_extent(sendtype,&lb,&sendextent);
  MPI_Type_get_extent(recvtype,&lb,&recvextent);

  // self-copy done in round 0
  for (i=0; i<size; i++) {
    int t, f;
    t = (rank+i)%size;
    f = (rank-i+size)%size;
    
    MPI_Sendrecv((char*)sendbuf+t*sendcount*sendextent,sendcount,sendtype,
		 t,ALLTOALLTAG,
		 (char*)recvbuf+f*recvcount*recvextent,recvcount,recvtype,
		 f,ALLTOALLTAG,
		 comm,MPI_STATUS_IGNORE);
  }

  return MPI_SUCCESS;
}

int Alltoall_telephone(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		       void *recvbuf, int recvcount, MPI_Datatype recvtype,
		       MPI_Comm comm)
{
  int rank, size;
  int i;

  MPI_Aint lb, sendextent, recvextent;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Type_get_extent(sendtype,&lb,&sendextent);
  MPI_Type_get_extent(recvtype,&lb,&recvextent);

  // self-copy done when rank is paired with itself
  for (i=0; i<size; i++) {
    int tf;
    tf = (i-rank+size)%size;

    MPI_Sendrecv((char*)sendbuf+tf*sendcount*sendextent,sendcount,sendtype,
		 tf,ALLTOALLTAG,
		 (char*)recvbuf+tf*recvcount*recvextent,recvcount,recvtype,
		 tf,ALLTOALLTAG,
		 comm,MPI_STATUS_IGNORE);
  }
  
  return MPI_SUCCESS;
}

int Alltoall_factor(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int recvcount, MPI_Datatype recvtype,
		    MPI_Comm comm)
{
  int rank, size;
  int i;

  MPI_Aint lb, sendextent, recvextent;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Type_get_extent(sendtype,&lb,&sendextent);
  MPI_Type_get_extent(recvtype,&lb,&recvextent);

  if ((size&0x1)==0x1) {
    // odd size
    for (i=0; i<size; i++) {
      int tf;
      tf = (i-rank+size)%size;
      
      MPI_Sendrecv((char*)sendbuf+tf*sendcount*sendextent,sendcount,sendtype,
		   tf,ALLTOALLTAG,
		   (char*)recvbuf+tf*recvcount*recvextent,recvcount,recvtype,
		   tf,ALLTOALLTAG,
		   comm,MPI_STATUS_IGNORE);
    }
  } else {
    // even size

    // self-copy before (or after) communication
    MPI_Sendrecv((char*)sendbuf+rank*sendcount*sendextent,sendcount,sendtype,
		 rank,ALLTOALLTAG,
		 (char*)recvbuf+rank*recvcount*recvextent,recvcount,recvtype,
		 rank,ALLTOALLTAG,
		 comm,MPI_STATUS_IGNORE);

    for (i=0; i<size-1; i++) {
      int tf;
      if (rank==size-1) {
	tf = i/2+(i%2)*(size/2);
      } else if ((2*rank)%(size-1)==i) {
	tf = size-1;
      } else {
	tf = (i-rank+(size-1))%(size-1);
      }
      
      MPI_Sendrecv((char*)sendbuf+tf*sendcount*sendextent,sendcount,sendtype,
		   tf,ALLTOALLTAG,
		   (char*)recvbuf+tf*recvcount*recvextent,recvcount,recvtype,
		   tf,ALLTOALLTAG,
		   comm,MPI_STATUS_IGNORE);
    }
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

  MPI_Init(&argc,&argv);
  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  count = 1;
  
  for (i=1; i<argc&&argv[i][0]=='-'; i++) {
    if (argv[i][1]=='c') i++,sscanf(argv[i],"%d",&count);
  }
  
  sendbuf = (tuwtype_t*)malloc(size*count*sizeof(tuwtype_t));
  assert(sendbuf!=NULL);
  recvbuf = (tuwtype_t*)malloc(size*count*sizeof(tuwtype_t));
  assert(recvbuf!=NULL);
  testbuf = (tuwtype_t*)malloc(size*count*sizeof(tuwtype_t));
  assert(testbuf!=NULL);
  
  for (i=0; i<size*count; i++) sendbuf[i] = (tuwtype_t)(rank+i); 
  for (i=0; i<size*count; i++) recvbuf[i] = (tuwtype_t)-1;
  for (i=0; i<size*count; i++) testbuf[i] = (tuwtype_t)-1;

  // "correctness test": compare against result from library function
  c = count-1;
  MPI_Alltoall(sendbuf,c,TUW_TYPE,recvbuf,c,TUW_TYPE,MPI_COMM_WORLD);
  ALLTOALL(sendbuf,c,TUW_TYPE,testbuf,c,TUW_TYPE,MPI_COMM_WORLD);
  for (i=0; i<size*c; i++) assert(recvbuf[i]==testbuf[i]);

  if (rank==0) {
    fprintf(stderr,"Alltoall, size=%d\\\\\n",size); 
    fprintf(stderr,"count & m (count) & m (Bytes) & avg & min \\\\\n"); 
  }

  int f = 5;
  for (c=1; c<=count; c*=f,FLIP(f,2,5)) {
    if (c*size>size*count) break;
    
    for (r=0,t=0; r<WARMUP+REPEAT; r++) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      
      start = MPI_Wtime();
      
      ALLTOALL(sendbuf,c,TUW_TYPE,testbuf,c,TUW_TYPE,MPI_COMM_WORLD);
      
      stop = MPI_Wtime();
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      if (r<WARMUP) continue;
      runtime[t++] = stop-start;
    }
    MPI_Allreduce(MPI_IN_PLACE,runtime,t,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    if (rank==0) {
      double tuwavg, tuwmin;
      
      tuwavg = 0.0; tuwmin = runtime[0]; 
      for (t=0; t<REPEAT; t++) {
	tuwavg += runtime[t];
	if (runtime[t]<tuwmin) tuwmin = runtime[t];
      }
      tuwavg /= REPEAT;
    
      fprintf(stderr,"%d & %d & %d & %.2f & %.2f \\\\\n", 
	      c,c*size,c*size*sizeof(tuwtype_t),tuwavg*MICRO,tuwmin*MICRO);
    }
  }
  
  MPI_Finalize();
  
  free(sendbuf);
  free(recvbuf);
  free(testbuf);
  
  return 0;
}
