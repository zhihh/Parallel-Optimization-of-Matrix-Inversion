


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>




int num_t; /* num of threads */
typedef struct {
	double** SquareMatrix;/*方阵定义*/
	int order; /*矩阵阶数*/
}matrix;


//求逆算法
void Inverse_Matrix(matrix* m);	//求方阵O的逆矩阵
double FindElement_U(matrix* m, int row, int column);	//求U矩阵元素U(row, column)
double FindElement_L(matrix* m, int row, int column);	//求L矩阵的元素L(row, column)
double FindElement_Inverse_U(matrix* m, int row, int column);	//求U矩阵的逆矩阵的元素Inverse_U(row, column)
double FindElement_Inverse_L(matrix* m, int row, int column);	//求L矩阵的逆矩阵的元素Inverse_L(row, column)
double FindElement_Inverse_O(matrix* m, int row, int column);   //求原矩阵(O矩阵)的逆矩阵的元素Inverse_O(row, column)

//输入输出
void Print_Matrix(matrix* m);//打印矩阵
void Mode_0(matrix* m);//手动输入矩阵
void Mode_1(matrix* m);//自动生成矩阵

double Loop_Compute(matrix* m,double entity, int lower, int upper, int row, int column);//循环计算过程

void IO_Matrix(matrix* m);//矩阵的输入输出


//第一次参数：矩阵的阶 第二个参数：并行线程数
int main(int argc, char* argv[]){
    
	matrix o;
	matrix* m = &o;
	m->order = strtol(argv[1], NULL, 10);
    num_t = strtol(argv[2], NULL, 10);

	//申请内存空间
	m->SquareMatrix = (double **)malloc(sizeof(double*)*m->order);
	for(int i = 0; i < m->order; i++){
		m->SquareMatrix[i] = (double*)malloc(sizeof(double)*m->order);
	}

	printf("输入0选择手动输入矩阵, 输入1选择自动输入矩阵:\n");
	int option;
	scanf("%d",&option);

	switch (option)
	{
	case 0:
		Mode_0(m);
		break;

	case 1:
		Mode_1(m);
		break;

	default:
		break;
	}
	



}

//矩阵的输入输出
void IO_Matrix(matrix* m){


	printf("原始矩阵:\n");
//	Print_Matrix(m);

	//并行
	printf("并行:\n");

	double t0 =  omp_get_wtime();//开始时间
	
	Inverse_Matrix(m);
	double t1 = omp_get_wtime();//结束时间
//	Print_Matrix(m);

	printf("MatrixInverse_Parallel finished in : %f seconds\n", t1 - t0);
	




	
}

//手动输入矩阵
void Mode_0(matrix* m){
	
	
	for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			scanf("%lf", &m->SquareMatrix[i][j]);
		}
		
	}
	IO_Matrix(m);



}

//自动生成矩阵
void Mode_1(matrix* m){
	time_t t;

	

	/* 初始化随机数发生器 */
	srand((unsigned) time(&t));

	for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			m->SquareMatrix[i][j] = rand()%100;
		}
	}

	IO_Matrix(m);

	
}


//打印矩阵
void Print_Matrix(matrix* m){
	for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			printf("%.4lf ", m->SquareMatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void Inverse_Matrix(matrix* m){

	/*
		矩阵相乘是左矩阵的行与右矩阵的列相乘；
		而对于三角矩阵，有上三角与下三角之分；
		在进行上述乘法时可以优化掉行列中0的成分，因此
			1. 可以优化算法
			2. 可以使得待求元素只对已经求出的元素和尚且“存在”的元素有依赖，故可以将其存于一处，以优化存储


		1. 由于L矩阵和U矩阵的每一个元素都依赖左上元素，而顺序遍历，左上元素已经求出，
		故，这里将L矩阵和U矩阵都存放在m矩阵（原始矩阵中），以节省存储开销；

		2. L*矩阵对L矩阵的左上分块依赖，故将L*矩阵放在m矩阵中；
		U*矩阵对U矩阵的右下分块依赖，故将求解U*矩阵的遍历顺序完全倒置，以便将U*矩阵放在m矩阵中；

		最后U*矩阵与L*矩阵相乘的出目标逆矩阵，根据元素的依赖性分析得出，可将目标矩阵存放在m矩阵中；
	*/


	//求LU矩阵
    for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			if(i <= j){
				/* 求U矩阵 */
				m->SquareMatrix[i][j] = FindElement_U(m, i, j);
			}else if (i > j)
			{
				/* 求L矩阵 */
				m->SquareMatrix[i][j] = FindElement_L(m, i, j);
			}
			
		}
	}
	//Print_Matrix(m);

	//分别求LU矩阵的逆矩阵
	for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			if(i > j){
				/* 求L矩阵的逆矩阵 */
				m->SquareMatrix[i][j] = FindElement_Inverse_L(m, i, j);
				//这里或可与上面求L矩阵合并

				/* 求U矩阵的逆矩阵 */
				//U矩阵的逆矩阵需要完全逆序遍历
				m->SquareMatrix[m->order-1 - i][m->order-1 - j] = FindElement_Inverse_U(m, m->order-1 - i, m->order-1 - j);
			}else if(i == j){
				m->SquareMatrix[m->order-1 - i][m->order-1 - j] = FindElement_Inverse_U(m, m->order-1 - i, m->order-1 - j);
			}
		}
	}
	//Print_Matrix(m);

	//求原(O)矩阵的逆矩阵
	for(int i = 0; i < m->order; i++){
		for(int j = 0; j < m->order; j++){
			m->SquareMatrix[i][j] = FindElement_Inverse_O(m, i, j);
		}
	}
	

}


//求U矩阵元素U(row, column)
double FindElement_U(matrix* m, int row, int column){
	/*
		row:行
		column:列
		(行列从0开始计数)
	*/
	double entity_U = 0; 

	entity_U = Loop_Compute(m, entity_U, 0, row, row, column);

//	for(int i = 0; i < row; i++){
//		entity_U += m->SquareMatrix[row][i]*m->SquareMatrix[i][column];
			/*
				求U矩阵算法：
					U矩阵，L矩阵，O矩阵（原始矩阵）
				m->SquareMatrix[i][column] = U[i][column];
				m->SquareMatrix[row][i] = L[row][i];
				L[row][0]*U[0][column] + L[row][1]*U[1][column] + ... + L[row][row-1]*U[row-1][column] + L[row][row]*U[row][column]= O[row][column];
				
				其中：
				L[row][row] = 1;


				技巧：
					求矩阵的某个元素总是 左边的行 与 右边的列 相乘；
						总是 Left[row][i] * Right[i][column]
				且 左矩阵为下三角矩阵 右矩阵为上三角矩阵
				则 变量i 上限受左矩阵束（即做矩阵到row列就为零了）
				故 变量i应从0到row;


				抑或着说，这里求的是上三角区域的值，上三角区域row<column，所以受到较小的row的束缚

			*/
//	}


	entity_U = m->SquareMatrix[row][column] - entity_U;

	return entity_U;
}


//求L矩阵的元素L(row, column)
double FindElement_L(matrix* m, int row, int column){

	double entity_L = 0;

	entity_L = Loop_Compute(m, entity_L, 0, column, row, column);	

//	for(int i = 0; i < column; i++){
//		entity_L += m->SquareMatrix[row][i]*m->SquareMatrix[i][column];
		/*
			求L矩阵算法：
			m->SquareMatrix[row][j] = L[row][j];
			m->SquareMatrix[j][column] = U[j][column];

			L[row][0]*U[0][column] + L[row][1]*U[0]*[column] + ... + L[row][column-1]*U[column-1][column] + L[row][column]*U[column][column] = O[row][column];
			
			技巧：
				求矩阵的某个元素总是 左边的行 与 右边的列 相乘；
					总是 Left[row][i] * Right[i][column]
				且 左矩阵为下三角矩阵 右矩阵为上三角矩阵
				则 变量i 上限受右矩阵束（即右矩阵到column行就为零了）
				故 变量i应从0到column;


			抑或着说，这里求的是下三角区域的值，下三角区域row>column，所以受到较小的column的束缚
				
		*/
		 
//	}
	entity_L = (m->SquareMatrix[row][column] - entity_L)/m->SquareMatrix[column][column];

	return entity_L;
}


//求U矩阵的逆矩阵的元素Inverse_U(row, column)
double FindElement_Inverse_U(matrix* m, int row, int column){

	double  entity_Inverse_U;

	if(row == column){
		entity_Inverse_U = 1;

		return entity_Inverse_U/m->SquareMatrix[row][column];
		/*
			此时，计算的是对角线上的值；
			U[row][column] * Inverse_U[row][column] = 1;
		*/

	}else{
		entity_Inverse_U = 0;

	}

	entity_Inverse_U = Loop_Compute(m, entity_Inverse_U, row+1, column+1, row, column);

//	for(int i = row+1; i <= column; i++){
//		entity_Inverse_U += m->SquareMatrix[row][i] * m->SquareMatrix[i][column];
	
		/*
			m->SquareMatrix[row][i] = U[row][i];
			m->SquareMatrix[i][column] = Inverse_U[i][column];

			U[row][row] * Inverse_U[row][column] + U[row][row+1] * Inverse_U[row+1][column] + U[row][row+2] * Inverse_U[row+2][column] + ... + U[row][column] * Inverse_U[column][column] = 0;
		
			技巧：
			求矩阵的某个元素总是 左边的行 与 右边的列 相乘；
				总是 Left[row][i] * Right[i][column]
			这里求Inverse_U矩阵，即右边的矩阵 
			且 左矩阵为上三角矩阵 右矩阵为上三角矩阵
			则 变量i 下限受到左矩阵束 上限受右矩阵束
			故 变量i应从row到column;

			【注】因为所求是上三角区域所以column>row

		*/
//	}

	entity_Inverse_U = -entity_Inverse_U/m->SquareMatrix[row][row];

	return entity_Inverse_U;

}


//求L矩阵的逆矩阵的元素Inverse_L(row, column)
double FindElement_Inverse_L(matrix* m, int row, int column){

	double entity_Inverse_L = m->SquareMatrix[row][column]; // entity_Inverse_L = L[row][column]

	entity_Inverse_L = Loop_Compute(m, entity_Inverse_L, column+1, row, row, column);

//	for(int i = column+1; i < row; i++){
//		entity_Inverse_L += m->SquareMatrix[row][i] * m->SquareMatrix[i][column];
		/*
			求L的逆矩阵算法：
			m->SquareMatrix[row][i] = L[row][i]
			m->SquareMatrix[i][column] = Inverse_L[i][column]

			满足方程:

			L[row][column]*Inverse_L[column][column] + L[row][column+1]*Inverse_L[column+1][column] + ... + L[row][row]*Inverse_L[row][column] = 0
			其中：
				Inverse_L[column][column] == 1
				L[row][row] == 1


			技巧：
			求矩阵的某个元素总是 左边的行 与 右边的列 相乘；
				总是 Left[row][i] * Right[i][column]
			这里求L*矩阵，即右边的矩阵 
			且 左矩阵为下三角矩阵 右矩阵为下三角矩阵
			则 变量i 下限受到右矩阵束 上限受左矩阵束
			故 变量i应从column到row;
			
			因为所求是下三角区域所以column<row
		*/
//	}

	entity_Inverse_L = - entity_Inverse_L;

	return entity_Inverse_L;
}


//求原矩阵(O矩阵)的逆矩阵的元素Inverse_O(row, column)
double FindElement_Inverse_O(matrix* m, int row, int column){
	double entity_Inverse_O ;//所求Oij实体
	int i ;//循环变量

	//根据左右矩阵的特点对循环变量的初值优化
	if(row > column){
		i = row;
		entity_Inverse_O = 0;
	}else{
		i = column+1; //这里取column+1，因为L[column][column] = 1，故U[row][column]*L[column][column] = U[row][column]，直接赋值如下
		entity_Inverse_O = m->SquareMatrix[row][column];
	}

	entity_Inverse_O = Loop_Compute(m, entity_Inverse_O, i, m->order, row, column);

	//计算
//	for(; i < m->order; i++){
//		entity_Inverse_O += m->SquareMatrix[row][i] * m->SquareMatrix[i][column];

	/*
			求U的逆✖️L的逆算法：
			左边：m->SquareMatrix[row][i] = U的逆[row][i]；
			右边：m->SquareMatrix[i][column] = L的逆[i][column]；

 
			满足方程:(方便起见，方程里用U和L表示二者的逆矩阵)
			
			O[row][column] = U[row][0]*L[0][column] + U[row][1]*L[1][column] + ... + U[row][order-1]*L[order-1][column] 

			其中：
				U[row][i<row] = 0;
				L[i<column][column] = 0;
				L[column][column] = 1;

				故
				当 row > column 时， 取i=row
					0[row][column] = U[row][row]*L[row][column] + U[row][row+1]*L[row+1][column] + ... + U[row][order-1]*L[order-1][column]
				当 row < column 时， 取i=column
					0[row][column] = U[row][column]*L[column][column] + U[row][column+1]*L[column+1][column] + ... + U[row][order-1]*L[order-1][column]
					又，L[column][column] = 1;
					化简为
					0[row][column] = U[row][column] + U[row][column+1]*L[column+1][column] + ... + U[row][order-1]*L[order-1][column]
			
			
			技巧：
			求矩阵的某个元素总是 左边的行 与 右边的列 相乘；
				总是 Left[row][i] * Right[i][column]
			这里求L*矩阵，即右边的矩阵 
			且 左矩阵为上三角矩阵 右矩阵为下三角矩阵
			则 变量i下限受到左右矩阵束 
			故 变量i应从(row>column?row:column)开始
			
			
		*/
//	}

	return entity_Inverse_O;
}


//循环计算过程
double Loop_Compute(matrix* m,double entity, int lower, int upper, int row, int column){
	/* entity是过程变量
		lower 循环下限
		upper 循环上限
	*/

	

	#pragma omp parallel for num_threads(num_t) reduction(+: entity)
		for(int i = lower; i < upper; i++)
			entity += m->SquareMatrix[row][i] * m->SquareMatrix[i][column];
	
	return entity;
}