#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>

class Point
{
	friend std::ostream& operator<< (std::ostream&, const Point&);
	friend std::istream& operator>> (std::istream&, Point&);
	friend Point operator+ (const Point&, const Point&);
	friend Point operator- (const Point&, const Point&);
	friend Point midpoint(const Point&, const Point&);
	friend double distance(const Point&, const Point&);
	friend Point barycenter(const std::vector<Point>&, const double *); // 求出三角形的质心(可以求加权的质心)
	friend double AREA(const Point&, const Point&, const Point&);

	private:
		double x[2];   
	
	public:
		Point();  
		Point(const double*);  
 		Point(const Point&);
 		Point(const double&, const double&);  

 		Point& operator= (const Point&);
 		double& operator[] (int);
 		const double& operator[] (int) const;
 		double length() const;
 		Point& operator+= (const Point&);		
 		Point& operator-= (const Point&);
 		Point& operator*= (const double&);
 		Point& operator/= (const double&);
};

class TmpEle
{
	friend std::ostream& operator<< (std::ostream&, const TmpEle&);

	private:
		std::vector<Point> Pnt;
		std::vector<Point> GaussPnt;
		std::vector<double> GaussWeight;

	public:
		// 构造函数
		TmpEle();
		TmpEle(int);
		void buildTE(int);

		Point getPnt(int);
		double getArea();
		std::vector<Point> getGaussPnt();
		std::vector<double> getGaussWeight();

		// 标准单元到普通单元
		Point Local_to_Global(const Point&, const std::vector<Point> &) const;
		// 普通单元到标准单元
		Point Global_to_Local(const Point&, const std::vector<Point> &) const;
		// 标准单元到普通单元
		std::vector<Point> Local_to_Global(const std::vector<Point> &, const std::vector<Point> & ) const;
		// 普通单元到标准单元
		std::vector<Point> Global_to_Local(const std::vector<Point>&, const std::vector<Point>&) const;
		// 从标准单元到普通单元的雅克比矩阵
		double Local_to_Global_jacobian(const Point&, const std::vector<Point>&) const;
		// 从普通单元到标准单元的雅克比矩阵
		double Global_to_Local_jacobian(const Point&, const std::vector<Point>&) const;
		// 从标准单元到普通单元的雅克比矩阵
		std::vector<double> Local_to_Global_jacobian(const std::vector<Point>&, const std::vector<Point>&) const;
		// 从普通单元到标准单元的雅克比矩阵
		std::vector<double> Global_to_Local_jacobian(const std::vector<Point>&, const std::vector<Point>&) const; 
};

class Mesh
{
	friend std::ostream& operator<< (std::ostream&, const Mesh&);

	private:
		std::vector<Point> Pnt;   // 三角单元的点
		std::vector<int> BndPnt;  // 三角单元的边界点
		std::vector<std::vector<int>> Ele;  // 三角单元的边

	public:
		Mesh() = default;  // 默认构造函数

		void readData(const std::string&);
		int getEleVtx(int , int);
		std::vector<int> getEleVtx(int);
		Point getPnt(int);
		std::vector<Point> getPnt(std::vector<int>&);
		std::vector<int> getBndPnt();
		int n_element(); // 得到三角单元的数量
		int n_point();  // 得到点的个数
		int n_boundaryPoint(); // 得到边界点的个数
};

class PDE
{	
	public:
		// 偏微分方程的边界条件
		double u_boundary(const Point&);

		// 偏微分方程的右端项
		double f(const Point&);

		// 偏微分方程解析解
		double u_exact(const Point&);
};


class Matrix
{
	public:
		// 数据结构
		PDE pde;	
		Mesh mesh;     
		TmpEle tmpEle;  
		std::vector<double> u_h;

		Matrix() = default;
		Matrix(const std::string&, int);  
		// gauss-seidel迭代
		void GaussSeidel(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);

		// 基函数的构造
		void basisValue(Point&, std::vector<Point>&, std::vector<double>&);
		void basisValue(std::vector<Point>&, std::vector<Point>&, std::vector<std::vector<double>>&);

		// 梯度算子的计算
		void basisGrad(Point&, std::vector<Point>&, std::vector<std::vector<double>>&);
		void basisGrad(std::vector<Point>&, std::vector<Point>&, std::vector<std::vector<std::vector<double>>>&);

		// 单元刚度矩阵
		void buildEleMatrix_RHS(int, std::vector<std::vector<double>>&, std::vector<double>&);

		// 总刚度矩阵
		void buildMatrix_RHS(std::vector<std::vector<double>>&, std::vector<double>&);

		
		// 边界条件的处理
		void DealBoundary(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);

		// 计算L2误差
		double ComputerL2Error(std::vector<double> &);
	
};
