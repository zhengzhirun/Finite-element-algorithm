/*************************************************************************
	> File Name: Possion.h
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年08月19日 星期日 22时04分12秒
 ************************************************************************/

#ifndef _POSSION_H
#define _POSSION_H

#include "Matrix.h"

double Possion(const std::string, const std::string, std::vector<double>& ,int);    // 求解Possion方程
void nodesf2dat(const std::string, std::vector<Point>&, std::vector<int>&);    // 读取nodes.dat文档
void trigsf2dat(const std::string, std::vector<std::vector<int>>&);     // 读取trigs.dat文档
void writeTriangle(std::vector<Point>&, std::vector<std::vector<int>>&, std::vector<int>&, 
        const std::string);

#endif
