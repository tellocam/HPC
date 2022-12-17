#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NUM_PROBLEM_SIZES 10

int main(int argc, char** argv) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Create a communicator with all processes
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the rank and size of the communicator
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Define the range of problem sizes to test
  int problem_sizes[NUM_PROBLEM_SIZES] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

  // Allocate space to store the measured times
  double* times = malloc(NUM_PROBLEM_SIZES * size * sizeof(double));

  // Loop over the problem sizes
  for (int i = 0; i < NUM_PROBLEM_SIZES; i++) {
    // Allocate an array of the specified size and initialize it with random numbers
    int n = problem_sizes[i];
    float* array = malloc(n * sizeof(float));
    for (int j = 0; j < n; j++) {
      array[j] = rand() / (float) RAND_MAX;
    }

    // Synchronize all processes
    MPI_Barrier(comm);

    // Start the timer
    double t_start = MPI_Wtime();

    // Invoke the MPI_Allreduce operation using the MPI_MAX operator on the array
    MPI_Allreduce(MPI_IN_PLACE, array, n, MPI_FLOAT, MPI_MAX, comm);

    // Stop the timer and record the elapsed time
    double t_elapsed = MPI_Wtime() - t_start;
    times[i * size + rank] = t_elapsed;

    // Free the array
    free(array);
  }

  // Gather the measured times from all processes
  MPI_Gather(times + rank, NUM_PROBLEM_SIZES, MPI_DOUBLE, times, NUM_PROBLEM_SIZES, MPI_DOUBLE, 0, comm);

  // On the root process, compute the mean, median, and best seen time for each problem size, as well as confidence intervals and standard deviations
  if (rank == 0) {
    for (int i = 0; i < NUM_PROBLEM_SIZES; i++) {
      // Compute the mean time
      double mean = 0.0;
      for (int j = 0; j < size; j++) {
        mean += times[i * size + j];
      }
      mean /= size;

      // Compute the median time
      double* times_sorted = malloc(size * sizeof(double));
     
    }
  }
}