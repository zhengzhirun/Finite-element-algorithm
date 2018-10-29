/*************************************************************************
	> File Name: TmpEle.h
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时06分26秒
 ************************************************************************/

#ifndef TMPELE_H
#define TMPELE_H

#include "Point.h"
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
	std::vector<Point> Local_to_Global(const std::vector<Point> &,const std::vector<Point> & ) const;
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

#endif
