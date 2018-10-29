/*************************************************************************
	> File Name: PDE.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时45分37秒
 ************************************************************************/

#include "PDE.h"

// 边界条件(第一类边界条件)
double PDE::u_boundary(const Point &p, double t)
{
	double value;
    value = exp(0.5 * (p[0] + p[1]) - t);
    return value;
}

double PDE::f(const Point &p, double t)
{
	double value;
	value = -3/2.0 * exp(0.5 * (p[0] + p[1]) - t);
    return value; 
}

double PDE::u_exact(const Point& p, double t)
{
	double value;
    value = exp(0.5 * (p[0] + p[1]) - t);
    return value;
}

std::vector<double> PDE::u_exact_grad(const Point &p, double t)
{
    std::vector<double> val(2);
    val[0] = 0.5 * exp(0.5 * (p[0] + p[1]) - t);
    val[1] = 0.5 * exp(0.5 * (p[0] + p[1]) - t);
    return val;
}
