#include "Parabolic.h"

double Dt = 0.2;

void Parabolic::GaussSeidel(std::vector<std::vector<double> > &A_matrix, std::vector<double> &u, std::vector<double> &rhs)
{
   std::vector<double> last_u(A_matrix.size());  
   do{
      last_u=u; 
      double error=0.0;     
      for(int i=0;i<A_matrix.size();i++){
         double temp=rhs[i];
         for(int j=0; j<A_matrix[i].size(); j++){
             if(j != i){
                temp-=A_matrix[i][j]*u[j];
             }   
         }
         u[i]=temp/A_matrix[i][i];
      }
      u_int Vsize=u.size();
      for(int k=0;k<Vsize;k++){
         error+=(last_u[k]-u[k])*(last_u[k]-u[k]);
      }
      error=sqrt(error);
      
      if(error<1.0e-10)
          break;
   }while(1); 
}

void Parabolic::SolveLinearASparse(Eigen::SparseMatrix<double>& A_matrix_sparse,std::vector<double>& u,
        std::vector<double>& rhs, int k)
{
    auto b_size = rhs.size();
    Eigen::VectorXd b;
    b.resize(b_size);
    for (size_t i = 0; i != b_size; ++i)
        b(i) = rhs[i];
    Eigen::VectorXd u_h_now;
    // 求解器声明
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> bicg;
    Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> gmres;
    Eigen::DGMRES<Eigen::SparseMatrix<double>,Eigen::IdentityPreconditioner> dgmres;
    Eigen::MINRES<Eigen::SparseMatrix<double>,Eigen::Lower|Eigen::Upper,Eigen::IdentityPreconditioner> minres;
    switch (k){
        // 迭代法
        case 1:
            cg.compute(A_matrix_sparse);
            u_h_now = cg.solve(b);
            std::cout << "CG:\t #iterations: " << cg.iterations() << ", estimated error: " << cg.error() 
                    << std::endl;
            break;
        case 2:
            bicg.compute(A_matrix_sparse);
            u_h_now = bicg.solve(b);
            std::cout << "BiCGSTAB:\t #iterations: " << bicg.iterations() << ", estimated error: " 
                    << bicg.error() << std::endl;
            break;
        case 3:
            gmres.compute(A_matrix_sparse);
            u_h_now = gmres.solve(b);
            std::cout << "GMRES:\t #iterator: " << gmres.iterations() << ", estimated error: "
                    << gmres.error() << std::endl;
            break;
        case 4:
            dgmres.compute(A_matrix_sparse);
            u_h_now = dgmres.solve(b);
            std::cout << "DGMRES:\t #iterator: " << dgmres.iterations() << ", estimate error: "
                    << dgmres.error() << std::endl;
            break;
        case 5:
            minres.compute(A_matrix_sparse);
            u_h_now = minres.solve(b);
            std::cout << "MINRES:\t #iterator: " << minres.iterations() << ", estimate error: "
                    << minres.error() << std::endl;
            break;
    }

    // 将求解结果存储到u中
    for (size_t i = 0; i != u.size(); ++i)
        u[i] = u_h_now(i);
}

Parabolic::Parabolic(const std::string& file, int i)
{
	mesh.readData(file);
	tmpEle.buildTE(i);
}

// 成员函数
// 基函数的构造(面积坐标)
void Parabolic::basisValue(Point& p, std::vector<Point>& v, std::vector<double>& val)
{
	val.resize(3); // 分片线性单元
	double area = AREA(v[0], v[1], v[2]);
	val[0] = AREA(p,v[1],v[2]);
	val[1] = AREA(v[0],p,v[2]);
	val[2] = AREA(v[0],v[1],p);
	val[0] /= area;
	val[1] /= area;
	val[2] /= area;
}

void Parabolic::basisValue(std::vector<Point>& p, std::vector<Point>& v, 
	std::vector<std::vector<double>>& val)
{
	int n = p.size();
	val.resize(3);
	for(int i = 0; i < 3; i++)
		val[i].resize(n);
	double area = AREA(v[0],v[1],v[2]);

	for (int i=0; i < n; i++){
		val[0][i] = AREA(p[i],v[1],v[2]);
		val[1][i] = AREA(v[0],p[i],v[2]);
		val[2][i] = AREA(v[0],v[1],p[i]);
		val[0][i] /= area;
		val[1][i] /= area;
		val[2][i] /= area;
	}
}

void Parabolic::basisGrad(Point&, std::vector<Point>& v, 
	std::vector<std::vector<double>>& val)
{
	val.resize(3);
	val[0].resize(2);
	val[1].resize(2);
	val[2].resize(2);

	double area = AREA(v[0],v[1],v[2]) * 2;  // 面积坐标到直角坐标的雅克比为2s

	val[0][0] = (v[1][1] - v[2][1]) / area;
	val[0][1] = (v[2][0] - v[1][0]) / area;

	val[1][0] = (v[2][1] - v[0][1]) / area;
	val[1][1] = (v[0][0] - v[2][0]) / area;

	val[2][0] = (v[0][1] - v[1][1]) / area;
	val[2][1] = (v[1][0] - v[0][0]) / area;
}

void Parabolic::basisGrad(std::vector<Point>& p, std::vector<Point>& v, 
	std::vector<std::vector<std::vector<double>>>& val)
{
	int n = p.size();
	val.resize(3);
	for (int i = 0; i < 3; i++){
		val[i].resize(n);
	}
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < n; j++){
			val[i][j].resize(2);
		}
	}
	double area = AREA(v[0],v[1],v[2]) * 2; 
	for (int i = 0; i < n; i++){

		val[0][i][0] = (v[1][1] - v[2][1]) / area;
		val[0][i][1] = (v[2][0] - v[1][0]) / area;

		val[1][i][0] = (v[2][1] - v[0][1]) / area;
		val[1][i][1] = (v[0][0] - v[2][0]) / area;

		val[2][i][0] = (v[0][1] - v[1][1]) / area;
		val[2][i][1] = (v[1][0] - v[0][0]) / area;
	}
}

void Parabolic::buildEleMatrix_RHS(int n, std::vector<std::vector<double>> &Ele_Amatrix, 
                                   std::vector<std::vector<double>> &Ele_Mmatrix,
                                   std::vector<double> &Ele_Rhs, double t)
{
	Ele_Amatrix.resize(3);  // 单元刚度矩阵(左端a(u,v)),维度是3*3的
    Ele_Mmatrix.resize(3);  // 单元刚度矩阵(左端(u,v)),维度是3*3的
	Ele_Rhs.resize(3, 0.0);   // 单元刚度矩阵(右端),维度是3*1的
	for(int i=0; i<3; i++){
		Ele_Amatrix[i].resize(3, 0.0);
        Ele_Mmatrix[i].resize(3, 0.0);
	}

	// 得到第n个单元的边
	std::vector<int> EV = mesh.getEleVtx(n);
	// 得到单元的三个点
	std::vector<Point> EP(3);
	EP[0] = mesh.getPnt(EV[0]);
	EP[1] = mesh.getPnt(EV[1]);
	EP[2] = mesh.getPnt(EV[2]);
	// 标准三角单元的面积
	double volume = tmpEle.getArea();
	// 得到标准三角单元的Gauss点
	std::vector<Point> GP = tmpEle.getGaussPnt();
	// 得到标准三角单元的Gauss权重
	std::vector<double> GW = tmpEle.getGaussWeight();
	// 把标准单元的Gauss点投影到普通单元上
	std::vector<Point> q_point = tmpEle.Local_to_Global(GP, EP); 
	// 得到标准单元到普通单元的雅克比矩阵
	std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
	// 得到三角单元上基函数的值
	std::vector<std::vector<double>> basis_value;
	basisValue(q_point, EP, basis_value);
	// 得到三角单元上梯度算子的值
	std::vector<std::vector<std::vector<double>>> basis_grad;
	basisGrad(q_point, EP, basis_grad);
	// 计算右端的单元刚度矩阵
	for(int i = 0; i < q_point.size(); i++){
		double Jxw = GW[i] * jacobian[i] * volume;
		double f_value = pde.f(q_point[i],t);
		for(int j = 0; j < basis_value.size(); j++){
			Ele_Rhs[j] += Jxw * f_value * basis_value[j][i];
		}
	}
	for(int i = 0; i < q_point.size(); i++){
		double Jxw = GW[i] * jacobian[i] * volume;
		for(int k = 0; k<basis_grad.size(); k++){
			for(int j = 0; j<basis_grad.size(); j++){
				// 左端的单元刚度矩阵(a(u,v))
                Ele_Amatrix[j][k] += Jxw * (basis_grad[j][i][0] * 
                                            basis_grad[k][i][0] + basis_grad[j][i][1] * basis_grad[k][i][1]);
            }
		}
        for (int k = 0; k < basis_value.size(); ++k){
            for (int j = 0; j < basis_value.size(); ++j){
                // 左端的单元刚度矩阵((u,v))
                Ele_Mmatrix[j][k] += Jxw * basis_value[j][i] * basis_value[k][i];
            }
        }
	}
}

void Parabolic::buildMatrix_RHS_Sparse(Eigen::SparseMatrix<double> &A_matrix_sparse,
                                       Eigen::SparseMatrix<double> &M_matrix_sparse, 
                                       std::vector<double>& Rhs, double time)
{
    auto n_ele = mesh.n_element();
    auto n_pnt = mesh.n_point();
    A_matrix_sparse.resize(n_pnt,n_pnt);
    M_matrix_sparse.resize(n_pnt,n_pnt);
    Rhs.resize(n_pnt);
    // 定义一个用于存储稀疏矩阵的三元组
    std::vector<Eigen::Triplet<double>> Atriple, Mtriple;
    for (size_t i = 0; i != n_ele; ++i){
    	std::vector<int> NV = mesh.getEleVtx(i);
        std::vector<std::vector<double>> Ele_Amatrix(3);
        std::vector<std::vector<double>> Ele_Mmatrix(3);
        for (size_t j = 0; j != 3; ++j){
            Ele_Amatrix[j].resize(3,0.0);
            Ele_Mmatrix[j].resize(3,0.0);
        }
        std::vector<double> Ele_Rhs(3,0.0);
        buildEleMatrix_RHS(i,Ele_Amatrix,Ele_Mmatrix,Ele_Rhs,time);
        // 稀疏矩阵存储技术
        // 利用经典的三元组插入方式来存储A_matrix
        Atriple.push_back(Eigen::Triplet<double>(NV[0],NV[0],Ele_Amatrix[0][0]));
        Atriple.push_back(Eigen::Triplet<double>(NV[0],NV[1],Ele_Amatrix[0][1]));
        Atriple.push_back(Eigen::Triplet<double>(NV[0],NV[2],Ele_Amatrix[0][2]));
        Atriple.push_back(Eigen::Triplet<double>(NV[1],NV[0],Ele_Amatrix[1][0]));
        Atriple.push_back(Eigen::Triplet<double>(NV[1],NV[1],Ele_Amatrix[1][1]));
        Atriple.push_back(Eigen::Triplet<double>(NV[1],NV[2],Ele_Amatrix[1][2]));
        Atriple.push_back(Eigen::Triplet<double>(NV[2],NV[0],Ele_Amatrix[2][0]));
        Atriple.push_back(Eigen::Triplet<double>(NV[2],NV[1],Ele_Amatrix[2][1]));
        Atriple.push_back(Eigen::Triplet<double>(NV[2],NV[2],Ele_Amatrix[2][2]));
        // 利用经典的三元组插入方式来存储M_matrix
        Mtriple.push_back(Eigen::Triplet<double>(NV[0],NV[0],Ele_Mmatrix[0][0]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[0],NV[1],Ele_Mmatrix[0][1]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[0],NV[2],Ele_Mmatrix[0][2]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[1],NV[0],Ele_Mmatrix[1][0]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[1],NV[1],Ele_Mmatrix[1][1]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[1],NV[2],Ele_Mmatrix[1][2]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[2],NV[0],Ele_Mmatrix[2][0]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[2],NV[1],Ele_Mmatrix[2][1]));
        Mtriple.push_back(Eigen::Triplet<double>(NV[2],NV[2],Ele_Mmatrix[2][2]));

        // 不利用稀疏矩阵存储右端项Rhs
        Rhs[NV[0]] += Ele_Rhs[0];
        Rhs[NV[1]] += Ele_Rhs[1];
        Rhs[NV[2]] += Ele_Rhs[2];
    }
    // 把三元组转换成稀疏矩阵
    A_matrix_sparse.setFromTriplets(Atriple.begin(),Atriple.end());
    M_matrix_sparse.setFromTriplets(Mtriple.begin(),Mtriple.end());
}

// 对于稀疏矩阵A_matrix_sparse的边界条件处理(dirichlet边界条件)
void Parabolic::DealBoundary_Sparse(Eigen::SparseMatrix<double> &M_matrix_sparse, 
                                    std::vector<double>& Rhs, double t)
{
    auto n_pnt = mesh.n_point();
    auto n_bp = mesh.n_boundaryPoint();
    std::vector<int> BV = mesh.getBndPnt(); // 得到边界点
    for (size_t i = 0; i != n_bp; ++i){
        Point Pnt_i = mesh.getPnt(BV[i]); // 得到第i个单元的点
        double val = pde.u_boundary(Pnt_i,t);
        Rhs[BV[i]] = M_matrix_sparse.coeffRef(BV[i],BV[i]) * val;
        for (size_t j = 0; j != n_pnt; ++j){
            if (j != BV[i]){
                M_matrix_sparse.coeffRef(BV[i],j) = 0.0;
            }
        }
        for (size_t j = 0; j != n_pnt; ++j){
            if (j != BV[i]){
                Rhs[j] -= M_matrix_sparse.coeffRef(j, BV[i]) * val;
                M_matrix_sparse.coeffRef(j,BV[i]) = 0.0;
            }
        }
    }
    // 去除稀疏矩阵中的零元素
    M_matrix_sparse = M_matrix_sparse.pruned();
    // 压缩剩余空间
    M_matrix_sparse.makeCompressed();
}

// 全离散显格式
void Parabolic::fullDiscreteFormat(Eigen::SparseMatrix<double> &A_matrix_sparse,
                        const Eigen::SparseMatrix<double> &M_matrix_sparse,
                        const std::vector<double> &u_h, std::vector<double> &Rhs)
{
    // 去除矩阵A_matrix_sparse中的零元素
    A_matrix_sparse = A_matrix_sparse.pruned();
    // 压缩剩余空间
    A_matrix_sparse.makeCompressed();
    // 利用Eigen中的Matrix存储u_h,Rhs
    Eigen::MatrixXd u_h_matrix(u_h.size(),1);
    Eigen::MatrixXd Rhs_matrix(Rhs.size(),1);
    for (decltype(u_h.size()) i = 0; i != u_h.size(); ++i)
        u_h_matrix(i,0) = u_h[i];
    for (decltype(Rhs.size()) i = 0; i != Rhs.size(); ++i)
        Rhs_matrix(i,0) = Rhs[i];
    // Eigen中二元操作符支持稀疏矩阵和密集矩阵的混合操作
    // 实现全离散显格式
    Rhs_matrix = (M_matrix_sparse - Dt * A_matrix_sparse) * u_h_matrix + Rhs_matrix * Dt;
    // 把得到的Eigen密集型矩阵转换为vector类型
    for (decltype(Rhs.size()) i = 0; i != Rhs.size(); ++i)
        Rhs[i] = Rhs_matrix(i,0);
}

// 全离散隐格式
void Parabolic::fullDiscreteImplicitFormat(Eigen::SparseMatrix<double> &A_matrix_sparse,
                        Eigen::SparseMatrix<double> &M_matrix_sparse,
                        const std::vector<double> &u_h, std::vector<double> &Rhs)
{
    // 去除矩阵A_matrix_sparse中的零元素
    A_matrix_sparse = A_matrix_sparse.pruned();
    // 压缩剩余空间
    A_matrix_sparse.makeCompressed();
    // 利用Eigen中的Matrix存储u_h,Rhs
    Eigen::MatrixXd u_h_matrix(u_h.size(),1);
    Eigen::MatrixXd Rhs_matrix(Rhs.size(),1);
    for (decltype(u_h.size()) i = 0; i != u_h.size(); ++i)
        u_h_matrix(i,0) = u_h[i];
    for (decltype(Rhs.size()) i = 0; i != Rhs.size(); ++i)
        Rhs_matrix(i,0) = Rhs[i];
    // 实现全离散隐格式
    Rhs_matrix = Dt * Rhs_matrix + M_matrix_sparse * u_h_matrix;
    M_matrix_sparse = M_matrix_sparse + Dt * A_matrix_sparse;
    // 把得到的Eigen密集型矩阵转换为vector类型
    for (decltype(Rhs.size()) i = 0; i != Rhs.size(); ++i)
        Rhs[i] = Rhs_matrix(i,0);
}

double Parabolic::ComputerL2Error(const std::vector<double> &f, double time)
{
	double err=0.0;
	int n_ele=mesh.n_element();
	for(int j=0; j<n_ele; j++){
		std::vector<int> NV = mesh.getEleVtx(j);
		std::vector<Point> EP(3);
		EP[0]=mesh.getPnt(NV[0]);
		EP[1]=mesh.getPnt(NV[1]);
		EP[2]=mesh.getPnt(NV[2]);
		double volume = tmpEle.getArea();
		std::vector<Point> GP = tmpEle.getGaussPnt();

		std::vector<double> GW = tmpEle.getGaussWeight();
		std::vector<Point> q_point = tmpEle.Local_to_Global(GP, EP); 
		std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
		std::vector<std::vector<double> > basis_value;
		basisValue(q_point, EP, basis_value);

		for(int i = 0; i < q_point.size(); i++){
			double Jxw = GW[i]*jacobian[i]*volume;
			double f_value = pde.u_exact(q_point[i],time);
			double u_h_val = f[NV[0]]*basis_value[0][i]+f[NV[1]]*basis_value[1][i]+f[NV[2]]*basis_value[2][i];
			double df_val = f_value-u_h_val;
			err += Jxw*df_val*df_val;
		}
	}
	err=sqrt(err);
	return err;
}

double Parabolic::ComputerH1Error(const std::vector<double> &f, double time)
{
   double err=0.0;
   int n_ele=mesh.n_element();
   for(int j=0; j<n_ele; j++){
        std::vector<int> NV=mesh.getEleVtx(j);
        std::vector<Point> EP(3);
        EP[0]=mesh.getPnt(NV[0]);
        EP[1]=mesh.getPnt(NV[1]);
        EP[2]=mesh.getPnt(NV[2]);
        std::vector<Point> GP=tmpEle.getGaussPnt();

        std::vector<double> GW=tmpEle.getGaussWeight();
        std::vector<Point> q_point = tmpEle.Local_to_Global(GP, EP); 
        std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP, EP);
        std::vector<std::vector<std::vector<double>>> basis_grad;
        basisGrad(q_point, EP, basis_grad);
	  
        for(int i = 0; i < q_point.size(); i++){
            double Jxw = GW[i]*jacobian[i];
            std::vector<double> f_value = pde.u_exact_grad(q_point[i], time);
            double u_h_grad1=0.0;
            double u_h_grad2=0.0;

            u_h_grad1 = f[NV[0]] * basis_grad[0][i][0] + f[NV[1]] * basis_grad[1][i][0] +
                f[NV[2]] * basis_grad[2][i][0];
            u_h_grad2 = f[NV[0]] * basis_grad[0][i][1] + f[NV[1]] * basis_grad[1][i][1] +
                f[NV[2]] * basis_grad[2][i][1];
            
            double df_val1 = f_value[0] - u_h_grad1;
            double df_val2 = f_value[1] - u_h_grad2;
            err += Jxw * df_val1 * df_val1;
            err += Jxw * df_val2 * df_val2;
        }
   }
   err = sqrt(err);
   return err;
}

// 输出u_h在边界处的值
void Parabolic::checkBoundary(const std::vector<double> &u_h, double t)
{
    std::vector<int> boundaryPointIndex = mesh.getBndPnt();
    for (size_t i = 0; i != mesh.n_boundaryPoint(); ++i){
        std::cout << u_h[boundaryPointIndex[i]] << "\t" << pde.u_boundary(mesh.getPnt(boundaryPointIndex[i]),t) << "\t"
                    << pde.u_boundary(mesh.getPnt(boundaryPointIndex[i]),t) - u_h[boundaryPointIndex[i]] << "\n";
    }
    std::cout << "\n";
}

void Parabolic::checkValue(const std::vector<double> &u_h, double t)
{
    for (size_t i = 0; i != mesh.n_point(); ++i){
        std::cout << u_h[i] << "\t" << pde.u_exact(mesh.getPnt(i), t)
                    <<"\t" << pde.u_exact(mesh.getPnt(i),t) - u_h[i] << "\n";
    }
    std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// 类的接口函数
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Parabolic::runFullDiscreteFormat(const std::string &outputFile)
{
    int n_pnt = mesh.n_point();             // 得到网格点的个数
    u_h.resize(n_pnt);
    // 给u_h赋初始值
    for (size_t i = 0; i != n_pnt; ++i)
        u_h[i] = pde.u_exact(mesh.getPnt(i),0);  
    double L2Error, H1Error;
    std::ofstream out(outputFile);
    out << "time\tL2Error\tH1Error" << "\n";  
    for (double time = Dt; time <= 1.000009; time = time + Dt){
        // 定义一个用于存储A_matrix的稀疏矩阵A_matrix_sparse
        Eigen::SparseMatrix<double> A_matrix_sparse, M_matrix_sparse;
        std::vector<double> Rhs;
        // 求解刚度矩阵A(a(u,v)),M((u,v)),Rhs((f,v))
        buildMatrix_RHS_Sparse(A_matrix_sparse, M_matrix_sparse, Rhs, time);
        // 求解抛物方程的全离散显格式
        fullDiscreteFormat(A_matrix_sparse, M_matrix_sparse, u_h, Rhs);
        // 边界条件的处理
        DealBoundary_Sparse(M_matrix_sparse, Rhs, time);
        // k可以取1-5; gmres:case 3(求解PDE离散后形成的代数系统)
        std::cout << time << std::endl;
        SolveLinearASparse(M_matrix_sparse, u_h, Rhs, 1);
        L2Error = ComputerL2Error(u_h, time);      // 求解L2误差
        H1Error = ComputerH1Error(u_h, time);      // 求解H1误差 
        out << time << "\t" << L2Error << "\t" << H1Error << "\n";
    }
    out.close();
}

void Parabolic::runFullDiscreteImplicitFormat(const std::string &outputFile)
{
    int n_pnt = mesh.n_point();             // 得到网格点的个数
    u_h.resize(n_pnt);
    double L2Error, H1Error;
    // 给u_h赋初始值
    for (size_t i = 0; i != n_pnt; ++i)
        u_h[i] = pde.u_exact(mesh.getPnt(i),0); 
    L2Error = ComputerL2Error(u_h,0);
    H1Error = ComputerH1Error(u_h,0);
    std::ofstream out(outputFile);
    out << "time\tL2Error\tH1Error" << "\n";  
    for (double time = Dt; time <= 1.00009; time += Dt){
        // 定义一个用于存储A_matrix的稀疏矩阵A_matrix_sparse
        Eigen::SparseMatrix<double> A_matrix_sparse, M_matrix_sparse;
        std::vector<double> Rhs;
        // 求解刚度矩阵A(a(u,v)),M((u,v)),Rhs((f,v))
        buildMatrix_RHS_Sparse(A_matrix_sparse, M_matrix_sparse, Rhs, time);
        // 求解抛物方程的全离散隐格式
        fullDiscreteImplicitFormat(A_matrix_sparse, M_matrix_sparse, u_h, Rhs);
        // 边界条件的处理
        DealBoundary_Sparse(M_matrix_sparse, Rhs, time);
        // k可以取1-5; gmres:case 3(求解PDE离散后形成的代数系统)
        std::cout << time << std::endl;
        SolveLinearASparse(M_matrix_sparse, u_h, Rhs, 3);
        L2Error = ComputerL2Error(u_h, time);      // 求解L2误差
        H1Error = ComputerH1Error(u_h, time);      // 求解H1误差 
        out << time << "\t" << L2Error << "\t" << H1Error << "\n";
        //checkBoundary(u_h,time);
        //checkValue(u_h,time);
    }
    out.close();
}
