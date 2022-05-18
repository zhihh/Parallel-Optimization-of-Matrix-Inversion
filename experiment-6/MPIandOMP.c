#include "stdio.h"
#include "mpi.h"
#include "omp.h"
#define NUM_THREADS 8
int main(int argc,char*argv[])
{
    int my_rank,numprocs,thread_id,nthreads;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
     #pragma omp parallel num_threads(NUM_THREADS) private(thread_id, nthreads)
    {
        thread_id=omp_get_thread_num();
        nthreads=omp_get_num_threads();
        printf("Hello,The World!from thread number %d (on %d)for the MPI process number %d (on %d)\n",thread_id, nthreads, my_rank, numprocs);
    }
    MPI_Finalize();
    return 0;
}
// mpicc demo2.c -o demo2 -fopenmp
// mpirun ./demo2 -mp 6

