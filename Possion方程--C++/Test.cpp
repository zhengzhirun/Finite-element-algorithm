#include "Matrix.h"

using std::cout;
using std::endl;
using std::cin;

void PointTest()
{
	////////////////////////////////////////////////////////////////////////////
	// 测试Point类
	////////////////////////////////////////////////////////////////////////////
	// 第一种初始化的方式
	cout << "*********************************************************" << endl;
	cout << "\n\t\t类Point的测试\n\n\n";
	cout << "*********************************************************" << endl;
	Point p;
	cout << "默认初始化:\n" << p << endl;

	//第二种初始化的方式
	double a[] = {1,1};
	Point p1(a);
	cout << "初始化方式p(数组):\n" << p1 << endl;

	//第三种初始化的方式
	Point p2(p1);
	cout << "初始化方式p(Point类型):\n" << p2 << endl;

	//第四种初值化方式
	Point p3 = p2;
	cout << "初始化方式p = Point类型:\n" << p3 << endl;

	Point pp1(0.5,0.5);
	cout << "初始化方式p(double, double):\n" << pp1 << endl;

	cout << "重载运算符[]:\n" << p3[0] << '\t' << p3[1] << endl;

	cout << "输出Point类型的长度:\n " << p3.length() << endl;

	p3 += p2;
	cout << "重载运算符+=:\n" << p3 << endl;

	p3 -= p2;
	cout << "重载运算符-=:\n" << p3 << endl;

	p3 *= 2;
	cout << "重载运算符*=:\n" << p3 << endl;

	p3 /= 2;
	cout << "重载运算符/=:\n" << p3 << endl;

	p3 = p2 + p1;
	cout << "重载运算符+:\n " << p3 << endl;

	p3 = p2 - p1;
	cout << "重载运算符-:\n" << p3 << endl;
	cout << "输入两个数:" << endl;
	cin >> p3;
	cout << "重载输入流 >>:\n" <<p3 << endl;

	cout << "输出两个点的中点:\n" << midpoint(p1,p3) << endl;

	cout << "输出两个点的距离:\n" << distance(p1,p3) << endl;

	Point pp2(0,0), pp3(0.5,0), pp4(0,0.5);
	cout << "输出三角形的面积:\n" << AREA(pp2,pp3,pp4) << endl;

	std::vector<Point> pp{pp2,pp3,pp4};
	Point barycenter0 = barycenter(pp, NULL);
	cout << "输出三角形的质心:\n" << barycenter0 << endl;
}

void TmpEleTest()
{
	////////////////////////////////////////////////////////////////////////////
	// 测试TmpEle类
	////////////////////////////////////////////////////////////////////////////
	cout << "*********************************************************" << endl;
	cout << "\n\t\t类TmpEle的测试\n\n\n";
	cout << "*********************************************************" << endl;
	TmpEle te;
	cout << "默认初始化以及输出流重载的测试:\n" << te <<endl;

	TmpEle te1(1);
	cout << "给私有成员赋予一阶的Gauss点和权重:\n" << te1 << endl;

	TmpEle te2(2);
	cout << "给私有成员赋予二阶的Gauss点和权重:\n" << te2 << endl;

	TmpEle te3(3);
	cout << "给私有成员赋予三阶的Gauss点和权重:\n" << te3 << endl;

	TmpEle te4(4);
	cout << "给私有成员赋予四阶的Gauss点和权重:\n" << te4 << endl;

	TmpEle te5(5);
	cout << "给私有成员赋予五阶的Gauss点和权重:\n" << te5 << endl;

	TmpEle te6(6);
	cout << "给私有成员赋予六阶的Gauss点和权重:\n" << te6 << endl;

	TmpEle te7(7);
	cout << "给私有成员赋予七阶的Gauss点和权重:\n" << te7 << endl;

	TmpEle te10(10);
	cout << "给私有成员赋予十阶的Gauss点和权重:\n" << te10 << endl;

	TmpEle te12(12);
	cout << "给私有成员赋予十二阶的Gauss点和权重：\n" << te12 << endl;


	Point p4 = te12.getPnt(1);
	cout << "成员函数getPnt(int):\n" << p4 << endl;

	double Area = te12.getArea();
	cout << "成员函数getArea():\n" << Area << endl;

	cout << "成员函数getGaussPnt():" << endl;
	std::vector<Point> GP = te12.getGaussPnt();
	for (auto i : GP)
		cout << i << endl;

	cout << "成员函数getGaussWeight():" << endl; 
	std::vector<double> GW = te12.getGaussWeight();
	for (auto i : GW)
		cout << i << endl;

	Point p5(0,0), p6(0.5,0), p7(0,0.5);
	std::vector<Point> v1{p5,p6,p7};

	Point ppp1(0,0);
	Point p8 = te.Local_to_Global(ppp1,v1);
	cout << "测试函数Local_to_Global:\n" << "(0,0)->>\t"<< p8 << endl;
	
	Point ppp2(0,1);
	p8 = te.Local_to_Global(ppp2,v1);
	cout << "(0,1)->>\t" << p8 << endl;

	Point ppp3(1,0);
	p8 = te.Local_to_Global(ppp3,v1);
	cout << "(1,0)->>\t" << p8 << endl;

	Point ppp4(0,0);
	Point p9 = te.Global_to_Local(ppp4,v1);
	cout << "测试函数Global_to_Local:\n" << "(0,0)->> \t" << p9 << endl;

	Point ppp5(0,0.5);
	p9 = te.Global_to_Local(ppp5,v1);
	cout << "(0,0.5)->>\t" << p9 << endl;

	Point ppp6(0.5,0);
	p9 = te.Global_to_Local(ppp6,v1);
	cout << "(0.5,0)->>\t" << p9 << endl;

	std::vector<Point> v2{ppp1,ppp2,ppp3};
	std::vector<Point> vv = te.Local_to_Global(v2,v1);
	cout << "测试函数Local_to_Global:\n";
	for (int i = 0; i < v2.size(); i++)
		cout << v2[i] << "-->>" <<vv[i] << endl;

	std::vector<Point> vv1 = te.Global_to_Local(v1,v1);
	cout << "测试函数Global_to_Local:\n";
	for (int i = 0; i < v1.size(); i++)
		cout << v1[i] << "-->>" << vv1[i] << endl;

	double jacobian = te.Local_to_Global_jacobian(ppp1,v1);
	cout << "从标准单元到普通单元的雅克比矩阵:\n" << jacobian << endl;

	jacobian = te.Global_to_Local_jacobian(ppp1,v1);
	cout << "从普通单元到标准单元的雅克比矩阵:\n" << jacobian << endl;

	std::vector<double> jacobian1 = te.Local_to_Global_jacobian(v1,v1);
	cout << "从标准单元到普通单元的雅克比矩阵:\n";
	for (auto i : jacobian1)
		cout << i << "\t";
	cout << endl;

	jacobian1 = te.Global_to_Local_jacobian(v1,v1);
	cout << "从普通单元到标准单元的雅克比矩阵:\n";
	for (auto i : jacobian1)
		cout << i << "\t";
	cout << endl;
}

void MeshTest()
{
	 ////////////////////////////////////////////////////////////////////////////
    //测试类Mesh
    ////////////////////////////////////////////////////////////////////////////
	cout << "*********************************************************" << endl;
	cout << "\n\t\t类Mesh的测试\n\n\n";
	cout << "*********************************************************" << endl;

	Mesh mesh;
	cout << "测试默认构造函数:\n" << mesh << endl;

	mesh.readData("Trimesh.txt");
	cout << "测试函数readData\n" << mesh << endl;

	int ele = mesh.getEleVtx(2,2);  // 选取第三行边的第三个边
	cout << "测试函数getEleData\n" << ele << endl;

	std::vector<int> ele1 = mesh.getEleVtx(2);
	cout << "测试函数getEleData:\n";
	for (auto i : ele1)
		cout << i << "\t";
	cout << endl;

	Point p = mesh.getPnt(2);
	cout << "测试函数getPnt:\n" << p << endl;

	int n_ele = mesh.n_element();
	cout << "得到三角单元的个数:\n" << n_ele << endl;

	int n_point = mesh.n_point();
	cout << "得到剖分单元点的个数:\n" << n_point << endl;

	int n_boundaryPoint = mesh.n_boundaryPoint();
	cout << "得到剖分单元边界点的个数:\n" << n_boundaryPoint << endl;

	std::vector<int> num{1,2,3,4,5,6,7,8,9};
	std::vector<Point> point = mesh.getPnt(num);
	cout << "测试函数getPnt,输入一系列整数:\n";
	for (auto i : point)
		cout << i << "\n";
	cout << endl;
}

void PDETest()
{
	cout << "*********************************************************" << endl;
	cout << "\n\t\t类PDE的测试\n\n\n";
	cout << "*********************************************************" << endl;

	PDE pde;
	double p[2] = {1,1};
	cout << "PDE中成员函数u_boundary的测试:\n" << pde.u_boundary(p) << endl;
	Point p1(1,1);
	cout << pde.u_boundary(p) << endl;

	cout << "PDE中成员函数u_exact的测试:\n" << pde.u_exact(p) << endl;
	cout << pde.u_exact(p1) << endl;

	cout << "PDE中成员函数f的测试:\n" << pde.f(p) << endl;
	cout << pde.f(p1) << endl;
}


void MatrixTest()
{
	////////////////////////////////////////////////////////////////////////////
	//测试类Possion
	////////////////////////////////////////////////////////////////////////////
	cout << "*********************************************************" << endl;
	cout << "\n\t\t类Matrix的测试\n\n\n";
	cout << "*********************************************************" << endl;

	Matrix possion;
	cout << "****基函数构造测试****" << endl;
	Point p0(0,0);
	Point p1(1,0);
	Point p2(0,1);
	std::vector<Point> Te{p0,p1,p2};
	Point center = barycenter(Te,NULL);
	cout << "输出三角单元的质心:\n" << center << endl;
	std::vector<double> basis;
	cout << "标准单元上基函数的值:\n" << endl;
	possion.basisValue(p0,Te,basis);
	for (auto i : basis)
		cout << i << "\t";
	cout << endl;

	possion.basisValue(p1,Te,basis);
	for (auto i : basis)
		cout << i << "\t";
	cout << endl;

	possion.basisValue(p2,Te,basis);
	for (auto i : basis)
		cout << i << "\t";
	cout << endl;

	cout << "测试标准单元上的基函数的值,(重载basisValue函数):\n";
	std::vector<std::vector<double>> basis1;
	std::vector<Point> Te1{p0,p1,p2,p0,p1,p2};
	possion.basisValue(Te1,Te,basis1);
	for (auto& i : basis1){
		for (auto j : i)
			cout << j << "\t";
		cout << "\n";
	}
	cout << endl;

	cout << "****梯度算子计算测试****" << endl;
	std::vector<std::vector<double>> basisGrad1;
	possion.basisGrad(p0,Te,basisGrad1);
	for (auto &i : basisGrad1){
		for (auto j : i)
			cout << j << "\t";
		cout << "\n";
	}
	cout << endl;

	cout << "重载梯度算子的函数测试:\n";
	std::vector<std::vector<std::vector<double>>> basisGrad2;
	possion.basisGrad(Te,Te,basisGrad2);
	for (auto& i : basisGrad2){
		for (auto& j : i){
			for (auto k : j){
				cout << k << "\t";
			}
			cout << "\n";
		}
	} 
	cout <<endl;

	////////////////////////////////////////////////////////////////////////////
	//左端刚度矩阵和右端矩阵的测试(未经过边界处理)
	////////////////////////////////////////////////////////////////////////////

	Matrix Possion_matrix("Trimesh.txt",7);
	PDE pde;

	int n_pnt = Possion_matrix.mesh.n_point();
	Possion_matrix.u_h.resize(n_pnt);
	
	std::vector<std::vector<double>> A_matrix;
	std::vector<double> Rhs;

	Possion_matrix.buildMatrix_RHS(A_matrix, Rhs);

	cout << "输出(未经过边界处理)左端刚度矩阵和右端项:\n" << endl;
	// 输出左端刚度矩阵和右端项
	for (auto &i : A_matrix){
		for (auto j : i)
			cout << j << "\t";
		cout << "\n";
	}
	cout << endl;

	cout << "\n\n\n";
	for (auto i : Rhs)
		cout << i << "\n";
	cout << endl;

	Possion_matrix.DealBoundary(A_matrix, Possion_matrix.u_h, Rhs);

	cout << "输出(经过边界处理)左端刚度矩阵和右端项:\n" << endl;
	// 输出左端刚度矩阵和右端项
	for (auto &i : A_matrix){
		for (auto j : i)
			cout << j << "\t";
		cout << "\n";
	}
	cout << endl;

	cout << "\n\n\n";
	for (auto i : Rhs)
		cout << i << "\n";
	cout << endl;

 	Possion_matrix.GaussSeidel(A_matrix,Possion_matrix.u_h, Rhs);
	cout << "数值解:\n" << endl;

	for (auto i : Possion_matrix.u_h)
		cout << i << "\n";
	cout << endl;


 	double L2Error=Possion_matrix.ComputerL2Error(Possion_matrix.u_h);
   	std::cout<<"L2Error="<<L2Error<<std::endl;


}

int main()
{
	PointTest();

	TmpEleTest();
	
	MeshTest();	

	PDETest();

	MatrixTest();

	return 0;
}