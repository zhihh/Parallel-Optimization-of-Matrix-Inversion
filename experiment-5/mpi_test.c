//测试mpi并行正优化的数量级

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <unistd.h>

double Trap(int n, int rank);

int main(int argc, char* argv[]){
    int n= strtol(argv[1], NULL, 10);
   

    int my_rank, comm_sz, local_n;
    double local_int, total_int;
    int source;
   
    double t1, t2;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);



    //并行
    t1 = MPI_Wtime();
    local_n = n/comm_sz;
 
    local_int = Trap(local_n, my_rank);

    if (my_rank != 0){
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&total_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }else{
        total_int = local_int;
        for (source = 1; source < comm_sz; source++)
        {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
        }
        t2 = MPI_Wtime();
        printf("并行:\nfinished in %fseconds\n", t2-t1);
    
        for(int i = 1; i < comm_sz; i++){
            MPI_Send(&total_int, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }

    printf("rank: %d, total_int: %lf\n", my_rank, total_int);

    
    //串行

    MPI_Barrier(MPI_COMM_WORLD);
    t1=MPI_Wtime();
    if (my_rank == 0){
        total_int = Trap(n, 0);
        t2=MPI_Wtime();
        printf("串行:\nfinished in %fseconds\n", t2-t1);
    }

   

    MPI_Finalize();
  
    return 0;
    
}

double Trap(int n, int rank){
    double temp = 0;
    for(int i = rank*n+1; i < (rank+1)*n+1; i++){
        temp += i;
    }
    return temp;

}



/*

With 4 ores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999309e+00
   finished in : 2.968043 seconds
With 1 cores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999517e+00
   finished in : 11.484140 seconds



With 8 ores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999162e+00
   finished in : 2.471632 seconds
With 1 cores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999517e+00
   finished in : 11.534453 seconds


With 2 ores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999922e+00
   finished in : 5.656407 seconds
With 1 cores:
   With n = 1000000000 trapezoids, 
   our estimate of the integral from 0.000000 to 3.000000 = 8.999999999999517e+00
   finished in : 11.368958 seconds
*/