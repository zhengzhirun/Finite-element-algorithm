/*************************************************************************
	> File Name: CVT_Error.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年08月23日 星期四 17时13分50秒
 ************************************************************************/
#include "Possion.h"
int main()
{
    // 网格初始化
    system("./mesh_init -f square.poly -q 1 -a 0.01");
    // 把初始网格信息写入文档中
    std::vector<Point> node;
    std::vector<std::vector<int>> ele;
    std::vector<int> boundary_node;
    writeTriangle(node,ele,boundary_node,"Trimesh_0.txt");  // 初始网格信息
    std::vector<double> Error;  // Error[0]为L2范误差,Error[1]为H1范误差
    std::vector<std::vector<double>> Errors; // 存储每一次的误差
    std::vector<double> eta;
    Error = Possion("Trimesh_0.txt","Results.txt",eta,7);
    // 记录点的个数
    std::vector<int> nodes; 
    system("./dens_info");
    Errors.push_back(Error);
    nodes.push_back(node.size());
    system("./mesh_opti -b 1");

    // 文件处理
    for (int i = 0; i != 200; ++i){
        std::string filename = "Trimesh_";
        filename = filename + std::to_string(i+1) + ".txt";
        writeTriangle(node,ele,boundary_node,filename,i);
        nodes.push_back(node.size());
        std::cout << "i = " << i << std::endl;
        Error = Possion(filename,"Results.txt",eta,7);
        Errors.push_back(Error);
    }
    // 输出文档
    std::ofstream out("Error.txt",std::ofstream::out);    
    for (decltype(Errors.size()) i = 0; i != Errors.size(); ++i)
        out << i+1 << "\t" << Errors[i][0] << "\t\t" << Errors[i][1] << "\t\t" << nodes[i] << "\n";
    out.close();
    
    return 0;
}
