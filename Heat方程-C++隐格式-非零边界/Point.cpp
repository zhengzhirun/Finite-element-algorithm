/*************************************************************************
	> File Name: Point.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时01分06秒
 ************************************************************************/
#include "Point.h"

//非成员接口函数
std::ostream& operator<< (std::ostream& os, const Point& p)
{
	for (auto i : p.x)
		os << i << '\t';
	return os;
}

std::istream& operator>> (std::istream& is, Point& p)
{
	for (auto &i : p.x)
		is >> i;
	return is;
}

Point operator+ (const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = p1.x[i] + p2.x[i];
	return p;
}

Point operator- (const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = p1.x[i] - p2.x[i];
	return p;
}

Point midpoint(const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = (p1.x[i] + p2.x[i]) / 2.0;
	return p;
}

double distance(const Point& p1, const Point& p2)
{
	double p = 0.;
	for (int i = 0; i != 2; ++i)
		p += (p1.x[i] - p2.x[i]) * (p1.x[i] - p2.x[i]);
	return sqrt(p);
}

Point barycenter(const std::vector<Point>& p, const double * w)
{
	double bc[2] = {0,0};
	int k = p.size();
	if (w == NULL){
		for (int i = 0; i < k; i++){
			bc[0] += p[i][0];
			bc[1] += p[i][1];
		}
		bc[0] /= k;
		bc[1] /= k;
	}
	else{
		double sw = 0;
		for (int i = 0; i < k; i++) sw += w[i];
		for (int i = 0; i < k; i++){
			bc[0] += w[i] * p[i][0];
			bc[1] += w[i] * p[i][1];
		}
		bc[0] /= sw;
		bc[1] /= sw;
	}
	return bc;
}

double AREA(const Point& p1, const Point& p2, const Point& p3)
{
	double area = ((p2.x[0] - p1.x[0]) * (p3.x[1] - p1.x[1]) - 
		(p3.x[0] - p1.x[0]) * (p2.x[1] - p1.x[1])) * 0.5;
	if (area < 0)
		area = -area;
	return area;
}
// 构造函数
Point::Point()
{
	for (int i = 0; i != 2; ++i)
		x[i] = 0;

}

Point::Point(const double* data)
{
	for (int i = 0; i != 2; ++i)
		x[i] = data[i];
}

Point::Point(const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] = p.x[i];
}

Point::Point(const double& v1, const double& v2)
{
	x[0] = v1;
	x[1] = v2;
}
// 成员函数
Point& Point::operator= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] = p.x[i];
	return *this;
}

double& Point::operator[] (int i)
{
	return x[i];
}

const double& Point::operator[] (int i) const
{
	return x[i];
}

double Point::length() const
{
	double v = 0.0;
	for (auto i : x)
		v += i * i;
	return sqrt(v);
}

Point& Point::operator+= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] += p.x[i];
	return *this;
}

Point& Point::operator-= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] -= p.x[i];
	return *this;
}

Point& Point::operator*= (const double& s)
{
	for (auto &i : x)
		i *= s;
	return *this;
}

Point& Point::operator/= (const double& s)
{
	for (auto &i : x)
		i /= s;
	return *this;
}
