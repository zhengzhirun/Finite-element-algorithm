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

	//Possion.buildMatrix_RHS(A_matrix,Rhs);   // 得到总的左边刚度矩阵和右端刚度矩阵
   
    // 定义一个用于存储A_matrix的稀疏矩阵A_matrix_sparse
    Eigen::SparseMatrix<double> A_matrix_sparse(A_matrix.size(),A_matrix.size());
    Possion.buildMatrix_RHS_Sparse(A_matrix_sparse,Rhs);
    
    Possion.DealBoundary_Sparse(A_matrix_sparse, Possion.u_h, Rhs);
        
    // k 可以取1-5 ; gmres:case 3
    Possion.SolveLinearASparse(A_matrix_sparse,Possion.u_h,Rhs,3);

    double L2Error=Possion.ComputerL2Error(Possion.u_h);  // 求解L2误差
   	std::cout<<"L2Error="<<L2Error<<std::endl;

    // 把数值解u_h和误差L2写入到文件中
   	std::ofstream out("Results.txt");
   	for (auto i : Possion.u_h)
   		out << i << "\n";
   	out.close();

    return 0;
}
