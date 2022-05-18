//omp优化测试代码

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char* argv[]){
    int max_num = strtol(argv[1], NULL, 10);
    int num_t = strtol(argv[2], NULL, 10);
    
    int sum = 0;
    
    //串行
    double t0 = omp_get_wtime();

    for(int i = 1; i <= max_num; i++)
        sum+=i*99;
    
    double t1 = omp_get_wtime();
    printf("串行:\nfinished in %f seconds\n\n", t1-t0);


    //并行
    t0 = omp_get_wtime();

    #pragma omp parallel for num_threads(num_t) reduction(+:sum)
    for(int i = 1; i <= max_num; i++)
        sum+=i*99;
    
    t1 = omp_get_wtime();
    printf("并行:\nfinished in %f seconds\n\n", t1-t0);

}