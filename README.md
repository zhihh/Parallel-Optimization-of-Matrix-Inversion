### 介绍
矩阵求逆的算法是LU算法，优化使用OpenMP和MPI

### 矩阵求逆在线计算网站
http://www.yunsuan.info/matrixcomputations/solvematrixinverse.htmlCancel changes

### 不同的几个阶段
* experiment-5  最基本的算法，尝试进行并行
* experiment-6  并行优化 方法包括: omp 和 mpi，主要优化点：通用计算过程、LU矩阵分别求逆
* experiment-7  两个比较完善的版本（但结果是错的，暂时不想改了，等考完研再说吧）
* experiment-8  应老师要求添加文件读写和对逗号分隔矩阵的兼容

### 注意点
Mac系统下，
对于omp，使用gcc-11编译
对于MPI，使用export OMPI_CC=gcc-11
参考：https://blog.csdn.net/baimafujinji/article/details/52769930

### 总结
* 完全独立地做一个小项目，收获的是满满的自信（答辩的时候很硬气啊）
* 做出贡献、付出精力，和队友友好相处，我真棒！


