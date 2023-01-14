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

// Character arrays for filenaming
char file_name[128];
char file_suffix[64];
char node_char[16];
char sizepn_char[16];
char pow_char[16];
char bs_char[16];
char uline[16] = "_";
char NCHAR[16] = "N";
char TCHAR[16] = "T";
char PCHAR[16] = "P";
char BSCHAR[16] = "B";

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

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

int MY_Reduce_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node, int blockSize)
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
            if (get_childL(node) < size)
            {
                MPI_Recv(tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < blockSize; j++)
                {
                    recvbuf[b+j] = tmpbufL[j] > sendbuf[b+j] ? tmpbufL[j] : sendbuf[b+j];
                }
                if (get_childR(node) < size)
                {
                    MPI_Recv(tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int j = 0; j < blockSize; j++)
                    {
                        recvbuf[b+j] = tmpbufR[j] > recvbuf[b+j] ? tmpbufR[j] : recvbuf[b+j];
                    }
                } 
            }
            
        }
        if(count % blockSize != 0){
             
            if (get_childL(node) < size)
            {
                MPI_Recv(tmpbufL_trunc, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < count % blockSize; j++)
                {
                    recvbuf[trunc_chunk_idx +j] = tmpbufL_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbufL_trunc[j] : sendbuf[trunc_chunk_idx +j];
                }
                if (get_childR(node) < size)
                {
                    MPI_Recv(tmpbufR_trunc, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int j = 0; j < count % blockSize; j++)
                    {
                        recvbuf[trunc_chunk_idx +j] = tmpbufR_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbufR_trunc[j] : recvbuf[trunc_chunk_idx +j];
                    }
                }
            }   
        }
    }
    else if (get_childL(node) >= size) // node == leaf (if no children <=> if no left child)
    {
        // printf("leaf: %d \n", node);
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD);
        }
        if(count % blockSize != 0){
            MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD);
        }
    }
    else // node == interior
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            if (get_childL(node) < size) // immer true
            {
                // printf("interior with left child: %d , b = %d\n", node, b);
                MPI_Recv(tmpbufL, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < blockSize; j++)
                {
                    recvbuf[b+j] = tmpbufL[j] > sendbuf[b+j] ? tmpbufL[j] : sendbuf[b+j];
                }
                if (get_childR(node) < size)
                {
                    // printf("interior with right child: %d , b = %d\n", node, b);
                    MPI_Recv(tmpbufR, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int j = 0; j < blockSize; j++)
                    {
                        recvbuf[b+j] = tmpbufR[j] > recvbuf[b+j] ? tmpbufR[j] : recvbuf[b+j];
                    }
                } 
                MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
            }
            // MPI_Send(recvbuf + (int)(floor(count/blockSize)-1)*blockSize, blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);
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
            if (get_childL(node) < size)
            {
                // printf("interior with left child: %d , b = trunk\n", node);
                MPI_Recv(tmpbufL_trunc, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < count % blockSize; j++)
                {
                    recvbuf[trunc_chunk_idx +j] = tmpbufL_trunc[j] > sendbuf[trunc_chunk_idx +j] ? tmpbufL_trunc[j] : sendbuf[trunc_chunk_idx +j];
                }
                if (get_childR(node) < size)
                {
                    // printf("interior with left child: %d , b = trunk\n", node);
                    MPI_Recv(tmpbufR_trunc, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int j = 0; j < count % blockSize; j++)
                    {
                        recvbuf[trunc_chunk_idx +j] = tmpbufR_trunc[j] > recvbuf[trunc_chunk_idx +j] ? tmpbufR_trunc[j] : recvbuf[trunc_chunk_idx +j];
                    }
                }
            }
            MPI_Send(recvbuf + trunc_chunk_idx, count%blockSize, MPI_DOUBLE, get_parent(node), 0 , MPI_COMM_WORLD);   
        }
    }  
    return MPI_SUCCESS;
}

int MY_Bcast_T(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node, int blockSize)
{
    int trunc_chunk_idx = floor(count/blockSize)*blockSize; // index of first truncated chunk
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
            if (get_childL(node) < size)
            {
                MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
                if (get_childR(node) < size)
                {
                    MPI_Send(sendbuf+b, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
                }
            }       
        }
        if(count % blockSize != 0)
        {
            if (get_childL(node) < size) 
            {
                MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
                if (get_childR(node) < size)
                {
                    MPI_Send(sendbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
                }
            }   
        }
    }
    else if (get_childL(node) >= size) // node == leaf (if no children <=> if no left child)
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
            if (get_childR(node) < size)
            {
                MPI_Send(recvbuf+b, blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
            }
            
        }
        if(count % blockSize != 0)
        {
            MPI_Recv(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_parent(node), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childL(node), 0, MPI_COMM_WORLD);
            if (get_childR(node) < size)
            {
                MPI_Send(recvbuf+trunc_chunk_idx, count % blockSize, MPI_DOUBLE, get_childR(node), 0, MPI_COMM_WORLD);
            }
            
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
    int power, count, c, hydra_nodes, gentxt, blockSize;
    tuwtype_t *sendbuf, *recvbuf, *testbuf, *testbuf1;
    tuwtype_t start, stop;
    int i, r, t;
    tuwtype_t *runtime;
    runtime = (tuwtype_t*)malloc(REPEAT * sizeof(tuwtype_t));
    assert(runtime!=0);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    count = 10;
    power = 2;
    gentxt = 0;
    blockSize = 4;
    hydra_nodes = 1;

    for (i=1; i<argc&&argv[i][0]=='-'; i++) {
        if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count); // commandline arg -c for adjusting max. count, if none given count = 10
        if (argv[i][1]=='p') i++, sscanf(argv[i], "%d", &power); // commandline arg -p for adjusting the power, if none given power = 2
        if (argv[i][1]=='b') i++, sscanf(argv[i], "%d", &blockSize); // commandline arg. -b for adjusting blocksize, if none given blockSize = 4
        if (argv[i][1]=='h') i++, sscanf(argv[i], "%d", &hydra_nodes); // commandline arg. -b for adjusting blocksize, if none given blockSize = 4
        if (argv[i][1]=='g') i++, sscanf(argv[i], "%d", &gentxt); // commandline arg. -g for generating a txt, if none given, no .txt
    }

    if(rank == 0){
        //printf("----- \n Number of processors must be of size 2^n-1. \n----- \n");
        fprintf(stderr,"Results for Ex4 with %d Processes, on %d Nodes with Blocksize %d and powers of %d \n",size/hydra_nodes, hydra_nodes, blockSize, power); 
        fprintf(stderr,"count, m (Bytes), avg, min, median, stddev,  CIMOE \n");
    }
    // mpirun -np 8 ./Ex4 -c 50 -p 2 -b 4 -h 1 -g 1

    FILE *fp; // file pointer
    // Filenaming is: "EX1_N36_T32_P2.txt" Where N36 = 36 Nodes, T32 = 32 Task per Node, P2 = Powers of 2
    if(rank==0){
        if (gentxt!=0){
            sprintf(file_suffix, "%s", NCHAR);
            snprintf(node_char, sizeof(node_char), "%d", hydra_nodes);
            snprintf(pow_char, sizeof(pow_char), "%d", power);
            snprintf(bs_char, sizeof(bs_char), "%d", blockSize);
            snprintf(sizepn_char, sizeof(sizepn_char), "%d", size/hydra_nodes);
            strcat(file_suffix, node_char);
            strcat(file_suffix, uline);
            strcat(file_suffix, TCHAR);
            strcat(file_suffix, sizepn_char);
            strcat(file_suffix, uline);
            strcat(file_suffix, PCHAR);
            strcat(file_suffix, pow_char);
            strcat(file_suffix, uline);
            strcat(file_suffix, BSCHAR);
            strcat(file_suffix, bs_char);
            sprintf(file_name, "EX2_%s.txt", file_suffix);  
            // mpicc -o Ex4 Ex4.c -lm -O3
            // mpirun -np 8 ./Ex4 -c 50 -h 1 -p 2 -b 4 -g 1
            fp = fopen(file_name, "w");
            if (fp == NULL) {
                printf("Error opening file!\n");
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // allocate and initialize data for "correctness tests":
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
    MY_Reduce_T(sendbuf, testbuf1, count, size, rank, blockSize);
    MY_Bcast_T(testbuf1, testbuf, count, size, rank, blockSize);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        // printf("node: %d, i = %d -> %f, %f, %f\n", rank, i, recvbuf[i], testbuf1[i], testbuf[i]); // for debugging
    }
 
        // TIME MEASURE FOR POWERS OF power , can be command line arg, otherwise its 2!
    for (c = 1; c <= count; c *= power)
    {
    if (c > count) break;
    
    for (r = 0, t = 0; r < WARMUP+REPEAT; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        start = MPI_Wtime();
        // start timing
        MY_Reduce_T(sendbuf, testbuf1, count, size, rank, blockSize);
        MY_Bcast_T(testbuf1, testbuf, count, size, rank, blockSize);
        // end timing
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

    // start timing
    // MY_Reduce_T(sendbuf, testbuf1, count, size, rank, blockSize);
    // MY_Bcast_T(testbuf1, testbuf, count, size, rank, blockSize);
    // end timing
 
    MPI_Finalize();
    free(sendbuf);
    free(recvbuf);
    free(testbuf);
    free(testbuf1);
    return 0;
}
