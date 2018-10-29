#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "Point.h"
#include "PDE.h"
#include "TmpEle.h"
#include "Mesh.h"

extern double Dt;           // 时间剖分的步长

class Parabolic
{
private:
	// 数据结构
	PDE pde;	
	Mesh mesh;     
	TmpEle tmpEle;  
	std::vector<double> u_h;
    // 私有函数
	// gauss-seidel迭代
	void GaussSeidel(std::vector<std::vector<double>>&,std::vector<double>&, std::vector<double>&);   
    // 利用Eigen库求解稀疏代数方程组
    void SolveLinearASparse(Eigen::SparseMatrix<double>&,std::vector<double>&, std::vector<double>&, int);
	// 基函数的构造
	void basisValue(Point&, std::vector<Point>&, std::vector<double>&);
	void basisValue(std::vector<Point>&, std::vector<Point>&, std::vector<std::vector<double>>&);
	// 梯度算子的计算
	void basisGrad(Point&, std::vector<Point>&, std::vector<std::vector<double>>&);
	void basisGrad(std::vector<Point>&, std::vector<Point>&, 
                    std::vector<std::vector<std::vector<double>>>&);
	// 单元刚度矩阵
	void buildEleMatrix_RHS(int, std::vector<std::vector<double>>&, 
                            std::vector<std::vector<double>>&, std::vector<double>&, double);    
    // 总刚度矩阵(稀疏存储)
    void buildMatrix_RHS_Sparse(Eigen::SparseMatrix<double>&, 
                                Eigen::SparseMatrix<double>&, std::vector<double>&, double); 
    // 稀疏矩阵的边界条件的处理
    void DealBoundary_Sparse(Eigen::SparseMatrix<double>&,
                             std::vector<double>&, double);
    // 求解抛物方程的全离散显格式
    void fullDiscreteFormat(Eigen::SparseMatrix<double> &,
                            const Eigen::SparseMatrix<double>&,
                            const std::vector<double> &, std::vector<double> &);
    // 求解抛物方程的全离散隐格式
    void fullDiscreteImplicitFormat(Eigen::SparseMatrix<double> &,
                                    Eigen::SparseMatrix<double> &,
                                    const std::vector<double> &, std::vector<double> &);
	// 计算L2误差
	double ComputerL2Error(const std::vector<double> &, double); 
    // 计算H1误差
    double ComputerH1Error(const std::vector<double> &, double);
    // 输出边界处的值
    void checkBoundary(const std::vector<double> &, double);
    // 输出解析解和数值解
    void checkValue(const std::vector<double> &, double);
public:
    Parabolic() = default;
	Parabolic(const std::string&, int);
    // 求解方程组(全离散显格式)
    void runFullDiscreteFormat(const std::string &);
    // 求解方程组(全离散隐格式)
    void runFullDiscreteImplicitFormat(const std::string &);
};

#endif
