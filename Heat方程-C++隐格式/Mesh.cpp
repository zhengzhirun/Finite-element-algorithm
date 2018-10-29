/*************************************************************************
	> File Name: Mesh.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时39分19秒
 ************************************************************************/

#include "Mesh.h"

// 非成员函数
std::ostream& operator<< (std::ostream& os, const Mesh& M)
{
	for (auto i : M.Pnt)
		os << i << '\n';

	for (auto &i : M.Ele){
		for (auto j : i)
			os << j << '\t';
		os << '\n';
	}

	for (auto i : M.BndPnt)
		os << i << '\n';
		os << '\n';
	return os;
}

// 成员函数
void Mesh::readData(const std::string& f)
{
	int i;
	std::ifstream is(f,std::ifstream::in); // 以只读的模式打开
	is >> i;
	Pnt.resize(i);
	for (int j = 0; j < i; j++)
		is >> Pnt[j][0] >> Pnt[j][1];

	int n;
	is >> n;
	Ele.resize(n);
	for (int j = 0; j < n; j++){
		Ele[j].resize(3);
		is >> Ele[j][0] >> Ele[j][1] >> Ele[j][2];
	}

	int bn;
	is >> bn;
	BndPnt.resize(bn);
	for (int j = 0; j < bn; j++)
		is >> BndPnt[j];

	is.close();
}

int Mesh::getEleVtx(int i, int j)
{
	return Ele[i][j];
}

std::vector<int> Mesh::getEleVtx(int i)
{
	return Ele[i];
}

Point Mesh::getPnt(int i)
{
	return Pnt[i];
}

std::vector<Point> Mesh::getPnt(std::vector<int>& vt)
{
	std::vector<Point> vec;
	for (int x : vt){
		Point point = Pnt[x];
		vec.push_back(point);
	}
	return vec;
}

std::vector<int> Mesh::getBndPnt()
{
	return BndPnt;
}

int Mesh::n_element()
{
	return Ele.size();
}

int Mesh::n_point()
{
	return Pnt.size();
}

int Mesh::n_boundaryPoint()
{
	return BndPnt.size();
}
