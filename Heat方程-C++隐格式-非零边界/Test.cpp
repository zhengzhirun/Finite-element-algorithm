/*************************************************************************
	> File Name: Test.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 15时07分53秒
 ************************************************************************/
#include "Parabolic.h"

int main()
{
    int index = 0;
    for (int i = 0; i != 5; ++i){
        std::cout << "*********************"<<i<<"******************************"<<std::endl;
        std::string inputFile = "Trimesh_";
        std::string outputFile = "Error_";
        inputFile = inputFile + std::to_string(index + i) + ".txt";
        outputFile = outputFile + std::to_string(index + i) + ".txt";
        Parabolic solve(inputFile, 5);
        solve.runFullDiscreteImplicitFormat(outputFile);        // 全离散显格式
        Dt = Dt / 4.0;
    }
    return 0;
}
