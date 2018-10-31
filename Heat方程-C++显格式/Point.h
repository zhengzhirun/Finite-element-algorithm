/*************************************************************************
	> File Name: Point.h
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 08时51分59秒
 ************************************************************************/

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>  // 线性代数库
#include <Eigen/Sparse> // 稀疏存储矩阵库,以及稀疏存储的线性方程组求解库
#include <Eigen/IterativeLinearSolvers> // 迭代法求解线性方程组的求解库
#include <unsupported/Eigen/IterativeSolvers>

class Point
{
	friend std::ostream& operator<< (std::ostream&, const Point&);
	friend std::istream& operator>> (std::istream&, Point&);
	friend Point operator+ (const Point&, const Point&);
	friend Point operator- (const Point&, const Point&);
	friend Point midpoint(const Point&, const Point&);
	friend double distance(const Point&, const Point&);
    // 求三角形的质心,可以求加权质心
	friend Point barycenter(const std::vector<Point>&, const double *);
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

#endif
