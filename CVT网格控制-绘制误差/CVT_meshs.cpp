/*************************************************************************
	> File Name: CVT_Adaptive_finite.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年08月19日 星期日 22时24分34秒
 ************************************************************************/

#include "Possion.h"
#include "density.h"
int main()
{
    // 参数设置
    double maxN = 2e3, theta = 0.5;
    int maxIt = 50;
    // 初始化网格
    system("./mesh_init -f square.poly -q 1 -a 0.01");
    // 把初始网格信息写入文档中
    std::vector<Point> node;
    std::vector<std::vector<int>> ele;
    std::vector<int> boundary_node;
    writeTriangle(node,ele,boundary_node,"Trimesh.txt");
    // 定义记录L2误差
    std::vector<double> L2Errors;
    // 定义记录点的个数
    std::vector<int> numbers_node;
    numbers_node.push_back(node.size());
    // 开始自适应的过程
    for (int index = 0; index != maxIt; ++index){
        // 求解Possion方程并且得到残量型误差指示子
        std::vector<double> eta;    // 存放误差指示子
        double L2Error = Possion("Trimesh.txt","Results.txt",eta,7);
        L2Errors.push_back(L2Error);
        if (node.size() > maxN) break;  // 当网格点的数量大于maxN时候退出循环
        // 优化和细化网格
        density(node,ele,eta);
        system("./mesh_refi -b 1 -p 0.2");
        system("./mesh_opti -b 1");
        std::vector<int> boundary_node1;
        writeTriangle(node,ele,boundary_node1,"Trimesh.txt");
        numbers_node.push_back(node.size());
    }
    // 将L2范误差输出到文档"L2Error.txt"中
    std::ofstream out("L2Error.txt",std::ofstream::out);   
    for (decltype(L2Errors.size()) i = 0; i != L2Errors.size(); ++i)
        out << L2Errors[i] << "\t" << numbers_node[i] << "\n";
    out.close();
    return 0;
}

