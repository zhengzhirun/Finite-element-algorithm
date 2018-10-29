/*************************************************************************
	> File Name: PDE.h
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时42分05秒
 ************************************************************************/

#ifndef PDE_H
#define PDE_H

#include "Point.h"

const double Pi = 4.0 * atan(1.0);

class PDE
{	
public:
	// 偏微分方程的边界条件(第一类边界条件)
	double u_boundary(const Point &, double);
	// 偏微分方程的右端项
	double f(const Point &, double);
	// 偏微分方程解析解
	double u_exact(const Point&, double);

    // 偏微分方程解析解的梯度
    std::vector<double> u_exact_grad(const Point &, double);
};

#endif
