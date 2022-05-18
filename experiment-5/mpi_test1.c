//测试mpi是否能传递二维数组

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <unistd.h>
#include <time.h>

#define n 10
typedef struct {
	double SquareMatrix[n][n];/*方阵定义*/
	int order; /*矩阵阶数*/
}matrix;

int main(int argc, char* argv[]){
   
   
    int my_rank, comm_sz;

    double t1, t2;


    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    matrix o;
    matrix* m = &o;
 
    
    if(my_rank == 0){
        time_t t;
	    /* 初始化随机数发生器 */
    	srand((unsigned) time(&t));

	    for(int i = 0; i < n; i++){
	    	for(int j = 0; j < n; j++){
	    		m->SquareMatrix[i][j] = rand()%100;
    		}
    	}

        printf("rank %d 矩阵为:\n", my_rank);
        
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                printf("%.2lf ", m->SquareMatrix[i][j]);
            }
            printf("\n");
        }
        
        MPI_Send(m->SquareMatrix, n*n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }else{
        t1 = MPI_Wtime();
        MPI_Recv(m->SquareMatrix, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        t2 = MPI_Wtime();
        printf("rank %d 接受用时%lf 矩阵为:\n ", my_rank, t2-t1);
        
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                printf("%.2lf ", m->SquareMatrix[i][j]);
            }
            printf("\n");
        }
        
        
    }



    MPI_Finalize();
  
    return 0;
    
}



