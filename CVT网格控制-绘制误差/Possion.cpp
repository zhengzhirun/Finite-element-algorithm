/*************************************************************************
	> File Name: Possion.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年08月19日 星期日 21时23分39秒
 ************************************************************************/
#include "Possion.h"

std::vector<double> Possion(const std::string filename_input,const std::string filename_output, 
        std::vector<double>& eta,int i)
{
    Matrix possion(filename_input,i);
    PDE pde;

    int n_pnt = possion.mesh.n_point();     // 得到网格点的个数
    possion.u_h.resize(n_pnt);  
    
    std::vector<std::vector<double>> A_matrix;
    std::vector<double> Rhs;
    
    // 定义一个用于存储A_matrix的稀疏矩阵A_matrix_sparse
    Eigen::SparseMatrix<double> A_matrix_sparse(A_matrix.size(),A_matrix.size());
    possion.buildMatrix_RHS_Sparse(A_matrix_sparse,Rhs);

    possion.DealBoundary_Sparse(A_matrix_sparse,possion.u_h,Rhs);
    possion.SolveLinearASparse(A_matrix_sparse,possion.u_h,Rhs,3);  // k 可以取1-5;这里取的3:gamres方法
    double L2Error = possion.ComputerL2Error(possion.u_h);
    double H1Error = possion.ComputerH1Error(possion.u_h);

    std::vector<double> Error(2);
    Error[0] = L2Error;
    Error[1] = H1Error;
    // 计算eta误差指示子
    possion.Estimateresidual(possion.u_h,eta);
    // 把数据u_h写入结果文档
    std::ofstream out(filename_output);
    for (auto i : possion.u_h)
        out << i << "\n";
    out.close();
    return Error;
}

std::vector<double> Possion(const std::string filename_input, const std::string filename_output, int i)
{
    Matrix possion(filename_input,i);
    PDE pde;

    int n_pnt = possion.mesh.n_point();     // 得到网格点的个数
    possion.u_h.resize(n_pnt);  
    
    std::vector<std::vector<double>> A_matrix;
    std::vector<double> Rhs;
    
    // 定义一个用于存储A_matrix的稀疏矩阵A_matrix_sparse
    Eigen::SparseMatrix<double> A_matrix_sparse(A_matrix.size(),A_matrix.size());
    possion.buildMatrix_RHS_Sparse(A_matrix_sparse,Rhs);

    possion.DealBoundary_Sparse(A_matrix_sparse,possion.u_h,Rhs);
    possion.SolveLinearASparse(A_matrix_sparse,possion.u_h,Rhs,3);  // k 可以取1-5;这里取的3:gamres方法
    double L2Error = possion.ComputerL2Error(possion.u_h);
    double H1Error = possion.ComputerH1Error(possion.u_h);

    std::vector<double> Error(2);
    Error[0] = L2Error;
    Error[1] = H1Error;
    // 把数据u_h写入结果文档
    std::ofstream out(filename_output);
    for (auto i : possion.u_h)
        out << i << "\n";
    out.close();
    return Error;
}

void nodesf2dat(const std::string datafile, std::vector<Point>& node, std::vector<int>& boundary_node)
{
    std::ifstream is(datafile,std::ifstream::in);   // 以只读模式打开
    int numbers; 
    is >> numbers;
    is.close();
    
    std::ifstream read(datafile,std::ifstream::in);     // 以只读模式打开
    std::vector<std::vector<double>> data;
    // 初始化
    data.resize(numbers + 1);
    for (decltype(data.size()) i = 0; i != numbers+1; ++i){
        data[i].resize(4);
    }
    for (decltype(data.size()) i = 0; i != numbers+1; ++i){
        read >> data[i][0] >> data[i][1] >> data[i][2] >> data[i][3];
    }
    read.close();

    // 初始化
    node.resize(numbers);
    for (decltype(node.size()) i = 0; i != numbers; ++i){
        node[i][0] = data[i+1][1];
        node[i][1] = data[i+1][2];
        if (data[i+1][3] != 0)
            boundary_node.push_back(data[i+1][0]-1);
    }
}
 
void trigsf2dat(const std::string datafile, std::vector<std::vector<int>>& ele)
{
    std::ifstream is(datafile,std::ifstream::in);   // 以只读模式打开
    int numbers;
    is >> numbers;
    is.close();

    std::ifstream read(datafile,std::ifstream::in);     // 以只读模式打开
    std::vector<std::vector<double>> infor;
    // 初始化
    infor.resize(numbers + 1);
    for (decltype(infor.size()) i = 0; i != numbers+1; ++i){
        infor[i].resize(4);
    }
    read >> infor[0][0] >> infor[0][1] >> infor[0][2];
    for (decltype(infor.size()) i = 0; i != numbers; ++i){
        read >> infor[i+1][0] >> infor[i+1][1] >> infor[i+1][2] >> infor[i+1][3];
    }
    read.close();

    // 初始化
    ele.resize(numbers);
    for (decltype(ele.size()) i = 0; i != numbers; ++i){
        ele[i].resize(3);
    }
    for (decltype(ele.size()) i = 0; i != numbers; ++i){
        ele[i][0] = infor[i+1][1] - 1;  // triangle 输出边的编号从1开始
        ele[i][1] = infor[i+1][2] - 1;  
        ele[i][2] = infor[i+1][3] - 1;
    }
}

void writeTriangle(std::vector<Point>& node, std::vector<std::vector<int>>& ele, 
                   std::vector<int>& boundary_node, const std::string output_file)
{
    nodesf2dat("nodes.dat",node,boundary_node);
    trigsf2dat("trigs.dat",ele);
    // 输出内容
    std::ofstream out(output_file,std::ofstream::out);
    out << node.size() << "\n";
    for(decltype(node.size()) i = 0; i != node.size(); ++i)
        out << node[i][0] << "\t" << node[i][1] << "\n";
    out << "\n";
   
    out << ele.size() << "\n";
    for (decltype(ele.size()) i = 0; i != ele.size(); ++i)
        out << ele[i][0] << "\t" << ele[i][1] << "\t" << ele[i][2] << "\n";
    out << "\n";

    out << boundary_node.size() << "\n";
    for (decltype(boundary_node.size()) i = 0; i != boundary_node.size(); ++i)
        out << boundary_node[i] << "\n";
    out << "\n";
    out.close();
}

void writeTriangle(std::vector<Point>& node, std::vector<std::vector<int>>& ele, 
        std::vector<int>& boundary_node, const std::string output_file, int iterator)
{
    std::string filename0 = "nodes_";
    filename0 = filename0 + std::to_string(iterator) + ".dat";
    
    std::string filename1 = "trigs_";
    filename1 = filename1 + std::to_string(iterator) + ".dat";

    nodesf2dat(filename0,node,boundary_node);
    trigsf2dat(filename1,ele);
    // 输出内容
    std::ofstream out(output_file,std::ofstream::out);
    out.precision(12);      // 将输出精度设置为12
    out << node.size() << "\n";
    for(decltype(node.size()) i = 0; i != node.size(); ++i)
        out << node[i][0] << "\t" << node[i][1] << "\n";
    out << "\n";
   
    out << ele.size() << "\n";
    for (decltype(ele.size()) i = 0; i != ele.size(); ++i)
        out << ele[i][0] << "\t" << ele[i][1] << "\t" << ele[i][2] << "\n";
    out << "\n";

    out << boundary_node.size() << "\n";
    for (decltype(boundary_node.size()) i = 0; i != boundary_node.size(); ++i)
        out << boundary_node[i] << "\n";
    out << "\n";
    out.close();
}

