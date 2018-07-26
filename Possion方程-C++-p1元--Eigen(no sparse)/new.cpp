/*************************************************************************
	> File Name: new.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年07月25日 星期三 19时27分41秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<Eigen/Dense>

int main()
{
    std::vector<std::vector<double>> A_matrix = {{0.68,0.597},{-0.211,0.823},{0.566,-0.605}};
    std::vector<double> b_matrix = {-0.33,0.563,-0.444};
    
    auto A_size = A_matrix.size();
    std::cout << A_size << std::endl;
    
    Eigen::MatrixXd A;
    A.resize(A_size,A_size);
    std::cout << A << std::endl;
    
    
    for (size_t i = 0; i != A_size; ++i){
        for (size_t j = 0; j != A_size; ++j){
            double value = A_matrix[i][j];
            A(i,j) = value;
        }
    }

    std::cout << A << std::endl;

    return 0;
}



