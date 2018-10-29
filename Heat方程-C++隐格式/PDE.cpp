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
    value = 0.0;
    return value;
}

double PDE::f(const Point &p, double t)
{
	double value;
    value = exp(t) * ( (p[0] * p[0] - p[0]) * (p[1] * p[1] - p[1]) - 
                    2 * (p[0] * p[0] - p[0] + p[1] * p[1] - p[1]) ); 
	return value; 
}

double PDE::u_exact(const Point& p, double t)
{
	double value;
	value = (p[0] * p[0] - p[0]) * (p[1] * p[1] - p[1]) * exp(t);
    return value;
}

std::vector<double> PDE::u_exact_grad(const Point &p, double t)
{
    std::vector<double> val(2);
    val[0] = -exp(t) * (-p[1] * p[1] + p[1]) * (2 * p[0] - 1);
    val[1] = -exp(t) * (-p[0] * p[0] + p[0]) * (2 * p[1] - 1);
    return val;
}
