#include "Matrix.h"

using std::cout;
using std::endl;
using std::cin;

int main()
{
	Matrix Possion("Trimesh.txt",7);
	PDE pde;

	int n_pnt = Possion.mesh.n_point();  // 得到网格点的个数
	Possion.u_h.resize(n_pnt);         

	std::vector<std::vector<double>> A_matrix;
	std::vector<double> Rhs;

	Possion.buildMatrix_RHS(A_matrix,Rhs);   // 得到总的左边刚度矩阵和右端刚度矩阵

	Possion.DealBoundary(A_matrix, Possion.u_h, Rhs);  // 边界处理

	Possion.GaussSeidel(A_matrix,Possion.u_h, Rhs);  // 高斯赛德尔迭代

	double L2Error=Possion.ComputerL2Error(Possion.u_h);  // 求解L2误差
   	std::cout<<"L2Error="<<L2Error<<std::endl;

   	// 把数值解u_h和误差L2写入到文件中
   	std::ofstream out("Results.txt");
   	for (auto i : Possion.u_h)
   		out << i << "\n";
   	out.close();

   	return 0;
}