/*************************************************************************
	> File Name: density.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年08月20日 星期一 22时33分51秒
 ************************************************************************/

#include "density.h"

void density(const std::vector<Point>& node, const std::vector<std::vector<int>>& ele, const std::vector<double>& eta)
{
    auto N = node.size();   // 网格节点数目
    auto NT = ele.size();   // 单元数目
    // 生成文件"dens_nodes.dat"
    std::ofstream out("dens_nodes.dat",std::ofstream::out);     // 以只读的模式打开

    auto radi = 6371000;
    std::vector<std::vector<double>> sp_crd;
    // 初始化
    sp_crd.resize(N);
    for (decltype(sp_crd.size()) i = 0; i != N; ++i)
        sp_crd[i].resize(3);
    
    for (decltype(sp_crd.size()) i = 0; i != N; ++i){
        auto x = 0.5 * node[i][0] / radi;
        auto y = 0.5 * node[i][1] / radi;
        auto xy2 = x * x + y * y;
        sp_crd[i][0] = radi * 2 * x / (1+xy2);
        sp_crd[i][1] = radi * 2 * y / (1+xy2);
        sp_crd[i][2] = radi * (1-xy2) / (1+xy2);
    }
    // 输出到文档中
    for (decltype(sp_crd.size()) i = 0; i != sp_crd.size(); ++i){
        out << sp_crd[i][0] << "\t" << sp_crd[i][1] << "\t" << sp_crd[i][2] << "\n";
    }
    out.close();
    
    // 生成文件"dens_valus.dat"
    // 检查节点i是否在单元j内,是为1,不是为0
    std::vector<std::vector<int>> p2t;
    p2t.resize(N);
    for (decltype(p2t.size()) i = 0; i != N; ++i)
        p2t[i].resize(NT);

    for (decltype(p2t.size()) i = 0; i != N; ++i){
        for (decltype(p2t.size()) j = 0; j != NT; ++j){
            bool direction_0 = false;
            bool direction_1 = false;
            bool direction_2 = false;
            if (i == ele[j][0]) direction_0 = true;
            if (i == ele[j][1]) direction_1 = true;
            if (i == ele[j][2]) direction_2 = true;
            if (direction_0 || direction_1 || direction_2)
                p2t[i][j] = 1;
            else
                p2t[i][j] = 0;
        }
    }
    
    // 生成单元尺寸
    std::vector<double> rho;    // 点的尺寸N
    std::vector<double> area;
    area.resize(NT);    // 初始化
    rho.resize(N);    // 初始化

    std::vector<double> temporary_rho;
    temporary_rho.resize(NT);   // 边的尺寸

    for (size_t i = 0; i != NT; ++i){
        Point point_0 = node[ele[i][0]];
        Point point_1 = node[ele[i][1]];
        Point point_2 = node[ele[i][2]];
        area[i] = AREA(point_0,point_1,point_2);
        auto distance_0 = distance(point_0,point_1);
        auto distance_1 = distance(point_0,point_2);
        auto distance_2 = distance(point_1,point_2);
        auto htri = (distance_0 + distance_1 + distance_2) / 3;
        // 计算暂时的rho
        temporary_rho[i] = pow(eta[i],2.0) / pow(htri,4.0);
    }
    // 转换为点的(这里有一个尺寸的变化)
    for (decltype(p2t.size()) i = 0; i != N; ++i){
        double temporary_0 = 0.0;
        double temporary_1 = 0.0;
        for (decltype(p2t[i].size()) j = 0; j != NT; ++j){
            temporary_0 += p2t[i][j] * area[j] * temporary_rho[j];
            temporary_1 += p2t[i][j] * area[j];
        }
        rho[i] = temporary_0 / temporary_1;
    }

    for (size_t i = 0; i != 5; ++i){
        for (size_t j = 0; j != NT; ++j){
            double temporary = 0.0;
            for (size_t k = 0; k != N; ++k){
                temporary += p2t[k][j] * rho[k] / 3.0;
            }
            temporary_rho[j] = temporary;
        }
        for (size_t j = 0; j != N; ++j){
            double temporary_0 = 0.0;
            double temporary_1 = 0.0;
            for (size_t k = 0; k != NT; ++k){
                temporary_0 += p2t[j][k] * temporary_rho[k] * area[k];
                temporary_1 += p2t[j][k] * area[k];
            }
            rho[j] = temporary_0 / temporary_1; 
        }
    }

    // 写入文件
    out.open("dens_valus.dat",std::ofstream::out);
    for (decltype(rho.size()) i = 0; i != rho.size(); ++i){
        out << rho[i] << "\n";
    }
    out.close();
}


