/*************************************************************************
	> File Name: Mesh.h
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时35分13秒
 ************************************************************************/

#ifndef MESH_H
#define MESH_H

#include "Point.h"
#include <fstream>
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

#endif
