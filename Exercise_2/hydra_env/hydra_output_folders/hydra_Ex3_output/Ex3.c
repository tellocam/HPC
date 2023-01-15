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

int get_blockSize(int count){
    if (count < 1e4){
        return 1e2;
    }
    else if (count < 1e6){
        return 1e3;
    }
    else return 1e6;
}

int get_parent(int id)
{
    return id-1;
}

int get_childL(int id)
{
    return id+1;
}

int MY_Allreduce_P(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int size, int node, int blockSize)
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
            sendSize = (j-d-1+1)*blockSize <= count ? blockSize : (j-d-1)*blockSize < count ? count%blockSize : 0;
            sendSize = (j-(d+1)<0 || j-(d+1)>=b) ? 0 : sendSize;
            
            MPI_Sendrecv(workbuf + (j-(d+1))*blockSize, sendSize, MPI_DOUBLE, get_childL(node), j, tmpbuf, blockSize, MPI_DOUBLE, get_childL(node), j, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            for (int i = 0; i < recvcount; i++)
            {
                workbuf[j*blockSize+i] = tmpbuf[i] > workbuf[j*blockSize+i] ? tmpbuf[i] : workbuf[j*blockSize+i];
            }
        }
        if (node != 0) // node != root
        {
            sendSize = (j+1)*blockSize <= count ? blockSize : j*blockSize < count ? count%blockSize : 0;
            sendSize = (j<0 || j>=b) ? 0 : sendSize;
            // if(node == size-1) printf("sendSize = %d, node = %d, j = %d \n", sendSize, node, j);
            MPI_Sendrecv(workbuf + j*blockSize, sendSize, MPI_DOUBLE, get_parent(node), j, tmpbuf, blockSize, MPI_DOUBLE, get_parent(node), j, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvcount);
            if (j-d >= 0 && j-d < b && recvcount > 0)
            {
                memcpy(&workbuf[(j-d)*blockSize], tmpbuf, sizeof(tuwtype_t)*recvcount);
            }
        }
    }

    // End:  recvbuf = workbuf
    memcpy(&recvbuf[0], workbuf, count * sizeof(tuwtype_t));  
    free(workbuf);
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
    tuwtype_t *sendbuf, *recvbuf, *testbuf;
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
    blockSize = 100;
    hydra_nodes = 1;

    for (i=1; i<argc&&argv[i][0]=='-'; i++) {
        if (argv[i][1]=='c') i++, sscanf(argv[i],"%d",&count); // commandline arg -c for adjusting max. count, if none given count = 10
        if (argv[i][1]=='p') i++, sscanf(argv[i], "%d", &power); // commandline arg -p for adjusting the power, if none given power = 2
        if (argv[i][1]=='b') i++, sscanf(argv[i], "%d", &blockSize); // commandline arg. -b for adjusting blocksize, if none given blockSize = 4
        if (argv[i][1]=='h') i++, sscanf(argv[i], "%d", &hydra_nodes); // commandline arg. -b for adjusting blocksize, if none given blockSize = 4
        if (argv[i][1]=='g') i++, sscanf(argv[i], "%d", &gentxt); // commandline arg. -g for generating a txt, if none given, no .txt
    }

    // the blockSize argument -b is not needed anymore from this point on
    // to prevent mistakes its gonna be kept in place also for Ex4 and Ex5..

    if(rank == 0){
        //printf("----- \n Number of processors must be of size 2^n-1. \n----- \n");
        fprintf(stderr,"Results for Ex3 with %d Processes, on %d Nodes with variable blockSize and powers of %d \n",size/hydra_nodes, hydra_nodes, power); 
        fprintf(stderr,"count, m (Bytes), avg, min, median, stddev,  CIMOE \n");
    }
    // mpirun -np 8 ./Ex3 -c 50 -p 2 -b 4 -h 1 -g 1

    FILE *fp; // file pointer
    // Filenaming is: "EX3_N36_T32_P2.txt" Where N36 = 36 Nodes, T32 = 32 Task per Node, P2 = Powers of 2
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
            //strcat(file_suffix, uline);
            //strcat(file_suffix, BSCHAR);    // commentend out bcs blocksize is variable now!
            //strcat(file_suffix, bs_char);
            sprintf(file_name, "EX3_%s.txt", file_suffix); 
            // mpicc -o Ex3 Ex3.c -lm O3 
            // mpirun -np 8 ./Ex3 -c 50 -p 2 -b 4 -h 1 -g 1
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

    for (int i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)i*(rank+1);
        // sendbuf[i] = (tuwtype_t)i*(size-rank);
        // sendbuf[i] = (tuwtype_t)((i*(rank+1))%3);

    for (int i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;

    for (int i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function
    MY_Allreduce_P(sendbuf, testbuf, count, size, rank, get_blockSize(count)); // get_blockSize redundant but for completeness :)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < count; i++) {
        assert(recvbuf[i] == testbuf[i]);
        //if(recvbuf[i] != testbuf[i]) printf("rank: %d -> %f, %f\n", rank, recvbuf[i], testbuf[i]); // for debugging
    }
 
    // TIME MEASURE FOR POWERS OF power , can be command line arg, otherwise its 2!
    for (c = 1; c <= count; c *= power)
    {
    if (c > count) break;
    
    for (r = 0, t = 0; r < WARMUP+REPEAT; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        start = MPI_Wtime();
        MY_Allreduce_P(sendbuf, testbuf, c, size, rank, get_blockSize(c));
        
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

 
    MPI_Finalize();
    free(sendbuf);
    free(recvbuf);
    free(testbuf);
    return 0;
}
