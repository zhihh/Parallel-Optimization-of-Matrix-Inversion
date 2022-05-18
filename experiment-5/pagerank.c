/*原则：1.所有矩阵向量，必须初始化，默认零矩阵 2.结果存放在第一个变量中*/

#define _CRT_SECURE_NO_WARNINGS 1
#include<stdio.h>
#include<math.h>




/*邻接矩阵存储结构代码*/
typedef char VertexType; /*顶点类型为char型*/
typedef double EdgeType; /*边上的权值类型定义为double型,因为后面会计算概率*/
#define MAXVEX 100 /*最大顶点数*/
#define INFINITY 0 /*给无穷大赋具体值，表示不区分无穷大和零*/
#define d 0.85 /*阻尼因子*/
#define e 0.00001 /*计算精度*/
typedef struct {
	VertexType ver[MAXVEX][MAXVEX]; /*顶点集，Vertex*/
	EdgeType adjM[MAXVEX][MAXVEX]; /*邻接矩阵-Adjacency Matrix*/
	int numNodes, numEdges; /*当前的顶点数和边数*/
}MGraph;

/*矩阵*/
//EdgeType(*M)[MAXVEX]; /*概率转移矩阵*/
//EdgeType(*A)[MAXVEX]; /*一般转移矩阵*/
//EdgeType(*E)[MAXVEX]; /*幺矩阵*/
//EdgeType(*T)[MAXVEX]; /*一般矩阵*/
//EdgeType P[MAXVEX]; /*PageRank向量*/
//int n; /*矩阵的阶*/
//int a, b, t;/*矩阵的数乘因子*/
//int t = 0; /*时刻*/



/*函数声明*/
/*图*/
void CreateMGraph(MGraph* G);//建立图
void PGraph(MGraph G);//打印图

/*向量*/
void InitialNullVector(EdgeType P[MAXVEX], int n);//零向量初始化（元素均为0）Initial NULL Vector
void InitialUnitVector(EdgeType P[MAXVEX], int n);//单位向量初始化（元素均为1）Initial Unit Vector
void PVector(EdgeType P[MAXVEX], int n);//打印向量 Print Vector
void InitialPRVector(EdgeType P[MAXVEX], int n);//PangeRank向量初始化 Vector Initial
void VSdd(EdgeType P[MAXVEX], int n);//PangeRank向量规范化 Vector Standard
/*向量运算*/
void VCopy(EdgeType P[MAXVEX], EdgeType T[MAXVEX], int n);//向量复制 Vector Copy  T->P
void VNMultiply(EdgeType T[MAXVEX], double num, int n);//向量数乘 Vector Number Multiply
void VVAdd(EdgeType T1[MAXVEX], EdgeType T2[MAXVEX], int n);//向量加法 Vector Addition
double VModule(EdgeType T[MAXVEX], int n);//向量求模 Vector Module

/*矩阵*/
void InitialUnitMatrix(EdgeType(*T)[MAXVEX], int n);//1矩阵初始化 Initial Unit Matrix (元素全为1的矩阵)
void InitialNullMatrix(EdgeType(*T)[MAXVEX], int n);//零矩阵初始化 Initial Null Matrix (元素全为0的矩阵)
void InitialYaoMatrix(EdgeType(*T)[MAXVEX], int n);//幺矩阵初始化 Initial Yao Matrix (元素对角线为1的矩阵，其它全为零)
void MatrixTrans(EdgeType(*T)[MAXVEX], int n);//矩阵转置 Matrix Transposition
void MNMultiply(EdgeType(*T)[MAXVEX], double t, int n);//矩阵数乘 Matrix Number Multiply
void MVMultiply(EdgeType(*T)[MAXVEX], EdgeType P[MAXVEX], int n);//矩阵×向量 Matrix Vector Multiply
void MMAdd(EdgeType(*T1)[MAXVEX], EdgeType(*T2)[MAXVEX], int n);//矩阵加法 Matrix Matrix Addition
void MMMultiply(EdgeType(*R)[MAXVEX], EdgeType(*l)[MAXVEX], EdgeType(*u)[MAXVEX], int n);//同阶方阵乘法	Matrix Matrix Multiply
void PMatrix(EdgeType(*T)[MAXVEX], int n);//打印矩阵 Print Matrix
void MCopy(EdgeType(*T1)[MAXVEX], EdgeType(*T2)[MAXVEX], int n);//矩阵复制 Matrix Copy

/*两个问题*/
void DeadEnds(EdgeType(*T)[MAXVEX], int n);//Dead Ends检测与修正
int SpiderTraps(EdgeType(*T)[MAXVEX], int n);//Spider Traps问题检测与修正
void RandomTeleport(EdgeType(*T)[MAXVEX], int n);//Spider Traps问题修正方法 

/*算法基础*/
int MSdd(EdgeType(*T)[MAXVEX], int n);//概率转移矩阵规范化,并检测Dead ends 和 Spider Traps问题 MatrixStandardization	
void MSddToP(EdgeType(*T)[MAXVEX], int n);//矩阵的初步规范化，使其表示概率 Matrix Standard To Probability
int JudgeAlgorithm(EdgeType Tb[MAXVEX], EdgeType Ta[MAXVEX], int n);//判断是否达到计算精度 Judge Algorithm

/*算法*/
int InterativeAlgorithm(EdgeType(*M)[MAXVEX], EdgeType P[MAXVEX], int n); //迭代算法 Iterative Algorithm 对应公式 P = (dM+(1-d)E/n)P
int PowerAlgorithm(EdgeType(*A)[MAXVEX], EdgeType P[MAXVEX], int n); //幂法 PowerAlgorithm

void AlgebraicAlgotithm(EdgeType(*M)[MAXVEX], EdgeType P[MAXVEX], int n);//代数算法 Algebraic Algotithm
double LUFind_Sum(EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], int max, int min, int left, int right, int n); //求LU矩阵过程中的求和公式 LU Find Sum
void LUFind(EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], EdgeType(*A)[MAXVEX], int n); //1.求LU矩阵 LU Find
void luFind(EdgeType(*l)[MAXVEX], EdgeType(*u)[MAXVEX], EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], int n);//2.求L和U矩阵的逆矩阵l和u lu Find
void InverseMFind(EdgeType(*TA)[MAXVEX], EdgeType(*TB)[MAXVEX], int n); //3.求逆矩阵 Find Inverse Matrix 


void Pend(VertexType(*v)[MAXVEX], EdgeType P[MAXVEX], int n);


/*----------------------------------------------------------------------------*/
/*图*/
/*邻接矩阵表示*/
void CreateMGraph(MGraph* G) {
	printf("输入顶点数和边数：\n");
	scanf("%d %d", &G->numNodes, &G->numEdges); /*输入顶点数和边数*/
	printf("请输入顶点(当中不能有空格)：\n");
	
	
	getchar();//读掉\n，避免读入下一个scanf
	/*读入顶点信息，建立顶点表*/
	for (int i = 0; i < G->numNodes; i++) {
		scanf("%s", G->ver[i]);
		
	} 
	
	/*邻接矩阵初始化*/
	for (int i = 0; i < G->numNodes; i++) {
		for (int j = 0; j < G->numNodes; j++) {
			G->adjM[i][j] = INFINITY;
		}
	} 
	
	/*输入数据*/
	int i, j;
	printf("请输入边(vi,vj)的下标i j：\n");
	for (int k = 0; k < G->numEdges; k++) {	
		scanf("%d %d", &i, &j);
		G->adjM[i][j] = 1;
	}
}

/*打印图-PrintGraph*/
void PGraph(MGraph G) {
	printf("图如下：\n");
	for (int i = 0; i < G.numNodes; i++) {
		for (int j = 0; j < G.numNodes; j++) {
			printf("   %.2f", G.adjM[i][j]);
		}
		printf("\n");
	}
}


/*-------------------------------------------------------------------------------------------------*/
/*向量*/

//打印向量 Print Vector
void PVector(EdgeType P[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		printf("%.4f\n", P[i]);
	}

}


//PangeRank向量初始化 Vector Initial
void InitialPRVector(EdgeType P[MAXVEX], int n) {
	//初始化
	InitialUnitVector(P, n);
	//规范化
	VSdd(P, n);
}

//PangeRank向量规范化 Vector Standard
void VSdd(EdgeType P[MAXVEX], int n) {
	double num = 0;

	for (int i = 0; i < n; i++) {
		num += P[i];/*求得所有元素的和*/
	}

	for (int i = 0; i < n; i++) {
		P[i] /= num;/*规范化，使其表示概率*/
	}
}

//单位向量初始化（元素均为1）Initial Unit Vector
void InitialUnitVector(EdgeType P[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		P[i] = 1;
	}
}

//零向量初始化（元素均为0）Initial NULL Vector
void InitialNullVector(EdgeType P[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		P[i] = 0;
	}
}

//向量复制 Vector Copy  T->P
void VCopy(EdgeType P[MAXVEX], EdgeType T[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		P[i] = T[i];
	}
}

//向量数乘 Vector Number Multiply
void VNMultiply(EdgeType T[MAXVEX], double num, int n) {
	for (int i = 0; i < n; i++) {
		T[i] *= num;
	}
}

//向量加法 Vector Addition
void VVAdd(EdgeType T1[MAXVEX], EdgeType T2[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		T1[i] += T2[i];/*结果放在T1里面*/
	}
}

//向量求模 Vector Module
double VModule(EdgeType T[MAXVEX], int n) {/*即求向量元素的平方根*/
	double mod = 0;
	for (int i = 0; i < n; i++) {
		mod += T[i] * T[i];
	}
	mod = sqrt(mod);
	return mod;
}



/*------------------------------------------------------------------------------------------------*/
/*矩阵*/

//打印矩阵 Print Matrix
void PMatrix(EdgeType (*T)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("   %.5f", T[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

//矩阵复制 Matrix Copy
void MCopy(EdgeType(*T1)[MAXVEX], EdgeType(*T2)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T1[i][j] = T2[i][j];
		}
	}
}

//1矩阵初始化 Initial Unit Matrix (元素全为1的矩阵)
void InitialUnitMatrix(EdgeType(*T)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[i][j] = 1;
		}
	}
}

//零矩阵初始化 Initial Null Matrix (元素全为0的矩阵)
void InitialNullMatrix(EdgeType(*T)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[i][j] = 0;
		}
	}
}

//幺矩阵初始化 Initial Yao Matrix (元素对角线为1的矩阵，其它全为零)
void InitialYaoMatrix(EdgeType(*T)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				T[i][j] = 1;
			}
			else {
				T[i][j] = 0;
			}
		}
	}
}

//矩阵转置 Matrix Transposition
void MatrixTrans(EdgeType(*T)[MAXVEX], int n) {
	EdgeType temp[MAXVEX][MAXVEX];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp[i][j] = T[j][i];
		}
	}
	MCopy(T, temp, n);
	
}

//矩阵数乘 Matrix Number Multiply
void MNMultiply(EdgeType(*T)[MAXVEX], double t, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[i][j] *= t;
		}
	}
}

//矩阵×向量 Matrix Vector Multiply
void MVMultiply(EdgeType(*T)[MAXVEX], EdgeType P[MAXVEX], int n) {
	EdgeType temp[MAXVEX];
	InitialNullVector(temp, n); /*临时零向量初始化*/
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp[i] += T[i][j] * P[j];
		}
	}
	VCopy(P, temp, n);/*将结果复制到原P向量*/

}

//矩阵加法 Matrix Matrix Addition
void MMAdd(EdgeType(*T1)[MAXVEX], EdgeType(*T2)[MAXVEX], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T1[i][j] += T2[i][j];
		}
	}
}

//同阶方阵乘法	Matrix Matrix Multiply
void MMMultiply(EdgeType(*R)[MAXVEX], EdgeType(*l)[MAXVEX], EdgeType(*u)[MAXVEX], int n) {/*R是Result结果矩阵，l是左矩阵，u是右矩阵*/
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//行列相乘
			for (int k = 0; k < n; k++) {
				R[i][j] += l[i][k] * u[k][j];
			}
			
		}
	}
}





/*-------------------------------------------------------------------------------------------------------------------*/
/*两大问题*/


//Dead Ends检测与修正
void DeadEnds(EdgeType(*T)[MAXVEX], int n) {
	//检测是否某一列全为零
	double flag;
	for (int i = 0; i < n; i++) {
		flag = 0;
		for (int j = 0; j < n; j++) {
			flag += T[j][i]; //若第i列全为零，则flag一直为零
		}

		//Teleport修正
		if (flag == 0) {/*若第i列全为零，则进行修正*/
			for (int k = 0; k < n; k++) {
				T[k][i] = (double)1.0 / n; 
			}
		}
	}
}

//Spider Traps问题检测与修正
int SpiderTraps(EdgeType(*T)[MAXVEX], int n) {
	/*检测是否某一列只有此列下标对应的行为1，即出现自环*/
	/*在SpiderTraps之前一定要进行列规范化*/
	/*因为每一列的和都为1，所以只要此列的对应行出现1，则该列其他行全为零*/
	for (int i = 0; i < n; i++) {
			if (T[i][i] == 1) {

				//Random Teleport修正
				RandomTeleport(T, n);
				return 1;//若进行了修正则返回true
				/*应为这里牵扯到对于逆矩阵方法，不论如何都进行了RandomTeleport修正，而其他两种方法只有在出现
				SpiderTraps时才会进行RandomTeleport修正，所以为了保证结果一致性，和避免重复修正，我们需要知道
				是否真的出现了SpiderTraps*/
			}
	}
	return 0;
}

//Spider Traps问题修正方法 
void RandomTeleport(EdgeType(*T)[MAXVEX], int n) {
	//Random Teleport修正
	EdgeType E[MAXVEX][MAXVEX];/*定义一个1矩阵*/
	InitialUnitMatrix(E, n);/*1矩阵初始化*/
	MSddToP(E, n);/*1矩阵列规范化*/

	/*1. 概率转移矩阵M修正*/
	MNMultiply(T, d, n);/*d是阻尼因子*/
	/*2. 1矩阵修正*/
	MNMultiply(E, 1 - d, n);
	/*3. 最终修正*/
	MMAdd(T, E, n);
}




/*-----------------------------------------------------------------------------------------*/
/*算法基础*/

//概率转移矩阵规范化,并检测Dead ends 和 Spider Traps问题 MatrixStandardization	
int MSdd(EdgeType(*T)[MAXVEX], int n) {
	/*Dead ends检测*/
	DeadEnds(T, n);
	//这里DeadEnds要在MSddToP之前运行，因为如果矩阵出现DeadEnds问题，MSddToP中将会出现0做分母的情况
	/*列规范化*/
	MSddToP(T, n);
	/*Spider Traps检测*/
	int flag = SpiderTraps(T, n);
	return flag;
}

//矩阵的初步规范化，使其表示概率 Matrix Standard To Probability
void MSddToP(EdgeType(*T)[MAXVEX], int n) {
	/*这是对列的规范化，以使其表示概率*/
	double sum[MAXVEX];/*用于记录每列不为零的行个数，即出度的个数*/
	for (int i = 0; i < n; i++) {
		sum[i] = 0;
	}/*sum数组初始化*/
	/*列规范化*/
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum[i] += T[j][i];/*因为从邻接矩阵转置到概率转移矩阵时，矩阵元素都是1或0的整数*/
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[j][i] /= sum[i];/*将每列的值进行规范化，使得矩阵每列之和为1*/
		}
	}
}

//判断是否达到计算精度 Judge Algorithm
int JudgeAlgorithm(EdgeType Tb[MAXVEX], EdgeType Ta[MAXVEX], int n) {
	/*Tb是此次计算之前的PR值，Ta是此次计算之后的PR值*/
	EdgeType temp[MAXVEX]; /*暂存Ta-Tb的结果*/
	/*计算Ta-Tb，并暂存结果与temp*/
	VCopy(temp, Ta, n);/*Ta备份到temp*/
	VNMultiply(Tb, -1, n);/*Tb取负值*/
	VVAdd(temp, Tb, n);/*temp与Tb的负值相加，相当于Ta-Tb*/
	VNMultiply(Tb, -1, n);/*恢复P*/
	/*求模*/
	double mod = VModule(temp, n);
	/*判断PageRank差值是否小于计算精度e*/
	if (mod < e) {
		return 1;
	}
	else {
		return 0;
	}
}





/*-----------------------------------------------------------------------------------------------*/
//算法

//迭代算法 Iterative Algorithm
//对应公式 P = (dM+(1-d)E/n)P
int InterativeAlgorithm(EdgeType(*M)[MAXVEX], EdgeType P[MAXVEX], int n) {
	int t = 0;//记录迭代次数
	/*这里认为传入的M是完全修正过的M*/
	EdgeType N[MAXVEX]; /*记录P的下一个状态*/
	int flag = 0; /*用于判断是否达到计算精度*/
	VCopy(N, P, n);/*初始时刻，将N与P置为相同项*/
	for (int i = 0; ; i++) {
		t++;
		MVMultiply(M, N, n);/*计算一次，结果存储在N中*/
		flag = JudgeAlgorithm(P, N, n);
		if (flag == 1) {
			/*进行规范化*/
			VSdd(P, n);
			
			return t;//返回迭代次数
		}
		else {
			VCopy(P, N, n);/*若没有达到计算精度，则将此次运算结果N赋值给P*/
		}
	}

	/*至此P中的值即为各个节点的PangRank值*/
}




//幂法 PowerAlgorithm
//幂法不同之处在于，每次迭代，会将运算结果进行单位化
int PowerAlgorithm(EdgeType(*A)[MAXVEX], EdgeType P[MAXVEX], int n) {
	int t = 0;//存储迭代次数
	/*这里认为传入的M是完全修正过的M*/
	EdgeType N[MAXVEX]; /*记录P的下一个状态*/
	int flag = 0; /*用于判断是否达到计算精度*/
	double mod = 0; /*用于每次的结果的单位化*/
	VCopy(N, P, n);/*初始时刻，将N与P置为相同项*/
	for (int i = 0; ; i++) {
		t++;
		MVMultiply(A, N, n);/*计算一次，结果存储在N中*/

		mod = VModule(N, n);/*结果向量求模*/
		VNMultiply(N, 1 / mod, n);/*结果向量单位化*/

		flag = JudgeAlgorithm(P, N, n);
		if (flag == 1) {
			/*进行规范化*/
			VSdd(P, n);
			return t;//返回迭代次数
		}
		else {
			VCopy(P, N, n);/*若没有达到计算精度，则将此次运算结果N赋值给P*/
		}
	}

	/*至此P中的值即为各个节点的PangRank值*/
}



//代数算法

//求逆矩阵算法--LU算法

//求LU矩阵过程中的求和公式 LU Find Sum
double LUFind_Sum(EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], int max, int min, int left, int right, int n) {
	/*max是循环上界，min是循环下界，left是矩阵左参即行，right是矩阵右参即列*/
	double sum = 0; /*用于存储求和结果*/
	for (int k = min; k <= max; k++) {
		sum += L[left][k] * U[k][right];
	}
	return sum;
}
//1.求LU矩阵 LU Find
void LUFind(EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], EdgeType(*A)[MAXVEX], int n) {

	/*1.行列大循环*/
	for (int r = 0; r < n; r++) {
		/*U先于L一轮小循环*/
		/*求Urj*/
		for (int j = r; j < n; j++) {
			U[r][j] = A[r][j] - LUFind_Sum(L, U, r - 1, 0, r, j, n);/*运用公式求出Urj*/
		}
		/*求Lir*/
		for (int i = r; i < n; i++) {
			L[i][r] = (A[i][r] - LUFind_Sum(L, U, r - 1, 0, i, r, n)) / U[r][r];/*运用公式求出Lir*/
		}
	}
}

//2.求L和U矩阵的逆矩阵l和u lu Find
void luFind(EdgeType(*l)[MAXVEX], EdgeType(*u)[MAXVEX], EdgeType(*L)[MAXVEX], EdgeType(*U)[MAXVEX], int n) {
	//l矩阵
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				l[i][j] = 1 / L[i][j];
			}
			else if (i < j) {
				//行递增
				l[j][i] = -LUFind_Sum(L, l, j - 1, i, j, i, n) * l[i][i];
			}
		}
	}
	//u矩阵
	for (int i = 0; i < n; i++) {
		for (int j = i; j >= 0; j--) {
			if (i == j) {
				u[j][i] = 1 / U[j][i];
			}
			else if (i > j) {
				//行递减
				u[j][i] = -LUFind_Sum(U, u, i, j+1, j, i, n) / U[j][j];
			}
		}
	}
}

//3.求逆矩阵 Find Inverse Matrix
void InverseMFind(EdgeType(*TA)[MAXVEX], EdgeType(*TB)[MAXVEX], int n) {
	/*TA is the one that after inversing*/
	/*TB is the one that before inversing*/

	EdgeType L[MAXVEX][MAXVEX];
	EdgeType U[MAXVEX][MAXVEX];
	EdgeType l[MAXVEX][MAXVEX];
	EdgeType u[MAXVEX][MAXVEX];

	//初始化
	InitialNullMatrix(L, n);
	InitialNullMatrix(U, n);
	InitialNullMatrix(l, n);
	InitialNullMatrix(u, n);

	/*求LU矩阵*/
	LUFind(L, U, TB, n);
	/*求lu矩阵*/
	luFind(l, u, L, U, n);
	/*求TA矩阵*/
	MMMultiply(TA, u, l, n);
}


//代数算法 Algebraic Algotithm
void AlgebraicAlgotithm( EdgeType(*M)[MAXVEX], EdgeType P[MAXVEX], int n) {
	EdgeType MI[MAXVEX][MAXVEX];/*逆矩阵 Matrix Inverse*/
	EdgeType I[MAXVEX][MAXVEX];/*幺矩阵*/
	EdgeType E[MAXVEX];/*单位向量*/
	/*逆矩阵初始化*/
	InitialNullMatrix(MI, n);
	/*幺矩阵初始化*/
	InitialYaoMatrix(I, n);
	/*单位向量初始化*/
	InitialUnitVector(E, n);


	/*求系数逆矩阵*/
	MNMultiply(M, -d, n);
	MMAdd(M, I, n);
	InverseMFind(MI, M, n);
	/*求修正向量*/
	VNMultiply(E, (double)(1 - d) / n, n);
	/*求PangRank值*/
	MVMultiply(MI, E, n);
	VCopy(P, E, n);


}
/*-------------------------------------------------------------------------------------------------------------*/
//主函数
int main() {
	/*建立图*/
	MGraph G;

	/*建立概率转移矩阵*/
	int n;//概率转移矩阵阶数
	EdgeType M[MAXVEX][MAXVEX];//概率转移矩阵
	EdgeType M1[MAXVEX][MAXVEX];//For InterativeAlgorithm
	EdgeType M2[MAXVEX][MAXVEX];//For PowerAlgorithm
	EdgeType M3[MAXVEX][MAXVEX];//For AlgebraicAlgorithm

	/*建立初始PageRank向量*/
	EdgeType P[MAXVEX];//PageRank向量
	EdgeType P1[MAXVEX];//For InterativeAlgorithm
	EdgeType P2[MAXVEX];//For PowerAlgorithm
	EdgeType P3[MAXVEX];//For AlgebraicAlgorithm


	int t;//记录迭代次数

	/*算法*/
	printf("                  PageRank - ( ゜- ゜)つロ 淦~                          \n");
	printf("--------------------------------------------------------------------\n");
	printf("        -------------1  输入图              ---------------\n");
	printf("        -------------2  展示图              ---------------\n");
	printf("        -------------3  InterativeAlgorithm---------------\n");
	printf("        -------------4  PowerAlgorithm     ---------------\n");
	printf("        -------------5  AlgebraicAlgorithm ---------------\n");
	printf("        -------------0  退出系统            ---------------\n");
	printf("--------------------------------------------------------------------\n");
	while (1) { 
		printf("请选择：\n");
		int cho;
		scanf("%d", &cho);
		switch (cho) {
		case 0:
			return 0;
		case 1:

			CreateMGraph(&G);//建立图
			n = G.numNodes;

			/*矩阵*/
			MCopy(M, G.adjM, n);//从邻接矩阵复制到概率转移矩阵
			MatrixTrans(M, n);//矩阵转置
			MSdd(M, n);//矩阵规范化
			
			
			/*向量*/
			InitialPRVector(P, n);//PageRank向量初始化
			
			break;
		case 2:
			PGraph(G);
			break;
		case 3:
			MCopy(M1, M, n);
			VCopy(P1, P, n);
			t = InterativeAlgorithm(M1, P1, n);
			printf("This is InterativeAlgorithm:\n");
			printf("  迭代次数：%d\n", t);
			printf("  PageRank值：\n");
			PVector(P1, n);
			printf("排名：\n");
			Pend(G.ver, P1, n);
			printf("\n");
			break;
		case 4:
			MCopy(M2, M, n);
			VCopy(P2, P, n);

			t = PowerAlgorithm(M2, P2, n);
			printf("This is PowerAlgorithm:\n");
			printf("  迭代次数：%d\n", t);
			printf("  PageRank值：\n");
			PVector(P2, n);
			printf("排名：\n");
			Pend(G.ver, P2, n);
			printf("\n");
			break;
		case 5:
			MCopy(M3, M, n);
			VCopy(P3, P, n);

			AlgebraicAlgotithm(M3, P3, n);
			printf("This is AlgebraicAlgotithm:\n");
			printf("  PageRank值：\n");
			PVector(P3, n);
			printf("排名：\n");
			Pend(G.ver, P3, n);
			printf("\n");
			break;
		}
	
	}




}void Pend(VertexType (*v)[MAXVEX], EdgeType P[MAXVEX], int n) {
	EdgeType t[6];
	for (int i = 0; i < n; i++) {
		t[i] = P[i];
	}
	EdgeType temp;
	int max;
	for (int i = 0; i < n; i++) {
		max = 0;
		for (int j = 0; j < n - i; j++) {
			if (t[j] > t[max]) {
				max = j;
			}
		}
		temp = t[max];
		t[max] = t[n - i - 1];
		t[n - i - 1] = temp;
	}
	for (int i = n-1; i >=0; i--) {
		for (int j = 0; j < n; j++) {
			if (P[j] == t[i]) {
				P[j] = 9999;
				printf("%s\n", v[j]);
				break;
			}
		}
	}
}