/*************************************************************************
	> File Name: TmpEle.cpp
	> Author: ZhirunZheng 
	> Mail: jiangxizhengzhirun@163.com 
	> Created Time: 2018年10月22日 星期一 09时10分27秒
 ************************************************************************/

#include "TmpEle.h"
// 非成员函数
std::ostream& operator<< (std::ostream& os, const TmpEle& te)
{	
	for (auto i : te.Pnt)
		os << i << '\n';
	for (auto i : te.GaussPnt)
		os << i << '\n';
	for (auto i : te.GaussWeight)
		os << i << '\n';
	return os;
}

// 构造函数
TmpEle::TmpEle()
{
	Pnt.resize(3);
	Pnt[0][0]=0.0;
	Pnt[0][1]=0.0; 
	Pnt[1][0]=1.0;
	Pnt[1][1]=0.0; 
	Pnt[2][0]=0.0;
	Pnt[2][1]=1.0; 
	
	GaussPnt.resize(7);
	GaussWeight.resize(7);
	GaussPnt[0][0]=0.333333333333;
	GaussPnt[1][0]=0.470142064105;
	GaussPnt[2][0]=0.470142064105;
	GaussPnt[3][0]=0.059715871790;
	GaussPnt[4][0]=0.101286507323;
	GaussPnt[5][0]=0.101286507323;
	GaussPnt[6][0]=0.797426985353;
	GaussPnt[0][1]=0.333333333333;
	GaussPnt[1][1]=0.470142064105;
	GaussPnt[2][1]=0.059715871790;
	GaussPnt[3][1]=0.470142064105;
	GaussPnt[4][1]=0.101286507323;
	GaussPnt[5][1]=0.797426985353;
	GaussPnt[6][1]=0.101286507323;
	GaussWeight[0]=0.225000000000;
	GaussWeight[1]=0.132394152788;
	GaussWeight[2]=0.132394152789;
	GaussWeight[3]=0.132394152788;
	GaussWeight[4]=0.125939180545;
	GaussWeight[5]=0.125939180545; 
	GaussWeight[6]=0.125939180545; 
}

TmpEle::TmpEle(int i)
{
  Pnt.resize(3);
  Pnt[0][0]=0.0;
  Pnt[0][1]=0.0; 
  Pnt[1][0]=1.0;
  Pnt[1][1]=0.0; 
  Pnt[2][0]=0.0;
  Pnt[2][1]=1.0; 
  switch(i){
     case 1:
        GaussPnt.resize(1);
		GaussWeight.resize(1);
		GaussPnt[0][0]=0.3333333333333333;
		GaussPnt[0][1]=0.3333333333333333;
		GaussWeight[0]=1.0;
        break;
     case 2:
        GaussPnt.resize(3);
		GaussWeight.resize(3);
		GaussPnt[0][0]=0.166666666667;
		GaussPnt[1][0]=0.166666666667;
		GaussPnt[2][0]=0.666666666667;
		GaussPnt[0][1]=0.166666666667;
		GaussPnt[1][1]=0.666666666667;
		GaussPnt[2][1]=0.166666666667;
		GaussWeight[0]=0.333333333333;
		GaussWeight[1]=0.333333333334;
		GaussWeight[2]=0.333333333333;
        break;
     case 3:
        GaussPnt.resize(4);
		GaussWeight.resize(4);
		GaussPnt[0][0]=0.333333333333;
		GaussPnt[1][0]=0.200000000000;
		GaussPnt[2][0]=0.200000000000;
		GaussPnt[3][0]=0.600000000000;
		GaussPnt[0][1]=0.333333333333;
		GaussPnt[1][1]=0.200000000000;
		GaussPnt[2][1]=0.600000000000;
		GaussPnt[3][1]=0.200000000000;
		GaussWeight[0]=-0.562500000000;
		GaussWeight[1]=0.520833333333;
		GaussWeight[2]=0.520833333334;
		GaussWeight[3]=0.520833333333;
        break;
     case 4:
        GaussPnt.resize(6);
		GaussWeight.resize(6);
		GaussPnt[0][0]=0.445948490916;
		GaussPnt[1][0]=0.445948490916;
		GaussPnt[2][0]=0.108103018168;
		GaussPnt[3][0]=0.091576213510;
		GaussPnt[4][0]=0.091576213510;
		GaussPnt[5][0]=0.816847572980;
		GaussPnt[0][1]=0.445948490916;
		GaussPnt[1][1]=0.108103018168;
		GaussPnt[2][1]=0.445948490916;
		GaussPnt[3][1]=0.091576213510;
		GaussPnt[4][1]=0.816847572980;
		GaussPnt[5][1]=0.091576213510;
		GaussWeight[0]=0.223381589678;
		GaussWeight[1]=0.223381589678;
		GaussWeight[2]=0.223381589678;
		GaussWeight[3]=0.109951743655;
		GaussWeight[4]=0.109951743656;
		GaussWeight[5]=0.109951743655;
        break;
     case 5:
        GaussPnt.resize(7);
		GaussWeight.resize(7);
		GaussPnt[0][0]=0.333333333333;
		GaussPnt[1][0]=0.470142064105;
		GaussPnt[2][0]=0.470142064105;
		GaussPnt[3][0]=0.059715871790;
		GaussPnt[4][0]=0.101286507323;
		GaussPnt[5][0]=0.101286507323;
		GaussPnt[6][0]=0.797426985353;
		GaussPnt[0][1]=0.333333333333;
		GaussPnt[1][1]=0.470142064105;
		GaussPnt[2][1]=0.059715871790;
		GaussPnt[3][1]=0.470142064105;
		GaussPnt[4][1]=0.101286507323;
		GaussPnt[5][1]=0.797426985353;
		GaussPnt[6][1]=0.101286507323;
		GaussWeight[0]=0.225000000000;
		GaussWeight[1]=0.132394152788;
		GaussWeight[2]=0.132394152789;
		GaussWeight[3]=0.132394152788;
		GaussWeight[4]=0.125939180545;
		GaussWeight[5]=0.125939180545; 
		GaussWeight[6]=0.125939180545; 
        break;
     case 6:
        GaussPnt.resize(12);
		GaussWeight.resize(12);
		GaussPnt[0][0]=0.249286745171;
		GaussPnt[1][0]=0.249286745171;
		GaussPnt[2][0]=0.501426509658;
		GaussPnt[3][0]=0.063089014492;
		GaussPnt[4][0]=0.063089014492;
		GaussPnt[5][0]=0.873821971017;
		GaussPnt[6][0]=0.310352451034;
		GaussPnt[7][0]=0.636502499121;
		GaussPnt[8][0]=0.053145049845;
		GaussPnt[9][0]=0.636502499121;
		GaussPnt[10][0]=0.310352451034;
		GaussPnt[11][0]=0.053145049845;
		GaussPnt[0][1]=0.249286745171;
		GaussPnt[1][1]=0.501426509658;
		GaussPnt[2][1]=0.249286745171;
		GaussPnt[3][1]=0.063089014492;
		GaussPnt[4][1]=0.873821971017;
		GaussPnt[5][1]=0.063089014492;
		GaussPnt[6][1]=0.636502499121;
		GaussPnt[7][1]=0.053145049845;
		GaussPnt[8][1]=0.310352451034;
		GaussPnt[9][1]=0.310352451034;
		GaussPnt[10][1]=0.053145049845;
		GaussPnt[11][1]=0.636502499121;
		GaussWeight[0]=0.116786275726;
		GaussWeight[1]=0.116786275726;
		GaussWeight[2]=0.116786275726;
		GaussWeight[3]=0.050844906370;
		GaussWeight[4]=0.050844906370;
		GaussWeight[5]=0.050844906370; 
		GaussWeight[6]=0.082851075618;
		GaussWeight[7]=0.082851075618; 
		GaussWeight[8]=0.082851075618; 
		GaussWeight[9]=0.082851075618; 
		GaussWeight[10]=0.082851075618; 
		GaussWeight[11]=0.082851075618; 
        break;
     case 7:
        GaussPnt.resize(13);
		GaussWeight.resize(13);
		GaussPnt[0][0]=0.333333333333;
		GaussPnt[1][0]=0.260345966079;
		GaussPnt[2][0]=0.260345966079;
		GaussPnt[3][0]=0.479308067842;
		GaussPnt[4][0]=0.065130102902;
		GaussPnt[5][0]=0.065130102902;
		GaussPnt[6][0]=0.869739794196;
		GaussPnt[7][0]=0.312865496005;
		GaussPnt[8][0]=0.638444188570;
		GaussPnt[9][0]=0.048690315425;
		GaussPnt[10][0]=0.638444188570;
		GaussPnt[11][0]=0.312865496005;
		GaussPnt[12][0]=0.048690315425;
		GaussPnt[0][1]=0.333333333333;
		GaussPnt[1][1]=0.260345966079;
		GaussPnt[2][1]=0.479308067842;
		GaussPnt[3][1]=0.260345966079;
		GaussPnt[4][1]=0.065130102902;
		GaussPnt[5][1]=0.869739794196;
		GaussPnt[6][1]=0.065130102902;
		GaussPnt[7][1]=0.638444188570;
		GaussPnt[8][1]=0.048690315425;
		GaussPnt[9][1]=0.312865496005;
		GaussPnt[10][1]=0.312865496005;
		GaussPnt[11][1]=0.048690315425;
		GaussPnt[12][1]=0.638444188570;
		GaussWeight[0]=-0.149570044468;
		GaussWeight[1]=0.175615257433;
		GaussWeight[2]=0.175615257433;
		GaussWeight[3]=0.175615257433;
		GaussWeight[4]=0.053347235609;
		GaussWeight[5]=0.053347235609; 
		GaussWeight[6]=0.053347235609;
		GaussWeight[7]=0.077113760890; 
		GaussWeight[8]=0.077113760890; 
		GaussWeight[9]=0.077113760890; 
		GaussWeight[10]=0.077113760890; 
		GaussWeight[11]=0.077113760890; 
		GaussWeight[12]=0.077113760890;  
        break;
     case 10:
        GaussPnt.resize(25);
		GaussWeight.resize(25);
		GaussPnt[0][0]=0.333333333333;
		GaussPnt[1][0]=0.485577633384;
		GaussPnt[2][0]=0.485577633384;
		GaussPnt[3][0]=0.028844733233;
		GaussPnt[4][0]=0.109481575485;
		GaussPnt[5][0]=0.109481575485;
		GaussPnt[6][0]=0.781036849030;
		GaussPnt[7][0]=0.307939838764;
		GaussPnt[8][0]=0.550352941821;
		GaussPnt[9][0]=0.141707219415;
		GaussPnt[10][0]=0.550352941821;
		GaussPnt[11][0]=0.307939838764;
		GaussPnt[12][0]=0.141707219415;
		GaussPnt[13][0]=0.246672560640;
		GaussPnt[14][0]=0.728323904597;
		GaussPnt[15][0]=0.025003534763;
		GaussPnt[16][0]=0.728323904597;
		GaussPnt[17][0]=0.246672560640;
		GaussPnt[18][0]=0.025003534763;
		GaussPnt[19][0]=0.066803251012;
		GaussPnt[20][0]=0.923655933587;
		GaussPnt[21][0]=0.009540815400;
		GaussPnt[22][0]=0.923655933587;
		GaussPnt[23][0]=0.066803251012;
		GaussPnt[24][0]=0.009540815400;
		GaussPnt[0][1]=0.333333333333;
		GaussPnt[1][1]=0.485577633384;
		GaussPnt[2][1]=0.028844733233;
		GaussPnt[3][1]=0.485577633384;
		GaussPnt[4][1]=0.109481575485;
		GaussPnt[5][1]=0.781036849030;
		GaussPnt[6][1]=0.109481575485;
		GaussPnt[7][1]=0.550352941821;
		GaussPnt[8][1]=0.141707219415;
		GaussPnt[9][1]=0.307939838764;
		GaussPnt[10][1]=0.307939838764;
		GaussPnt[11][1]=0.141707219415;
		GaussPnt[12][1]=0.550352941821;
		GaussPnt[13][1]=0.728323904597;
		GaussPnt[14][1]=0.025003534763;
		GaussPnt[15][1]=0.246672560640;
		GaussPnt[16][1]=0.246672560640;
		GaussPnt[17][1]=0.025003534763;
		GaussPnt[18][1]=0.728323904597;
		GaussPnt[19][1]=0.923655933587;
		GaussPnt[20][1]=0.009540815400;
		GaussPnt[21][1]=0.066803251012;
		GaussPnt[22][1]=0.066803251012;
		GaussPnt[23][1]=0.009540815400;
		GaussPnt[24][1]=0.923655933587;
		GaussWeight[0]=0.090817990383;
		GaussWeight[1]=0.036725957756;
		GaussWeight[2]=0.036725957756;
		GaussWeight[3]=0.036725957756;
		GaussWeight[4]=0.045321059436;
		GaussWeight[5]=0.045321059436; 
		GaussWeight[6]=0.045321059436;
		GaussWeight[7]=0.072757916845; 
		GaussWeight[8]=0.072757916845; 
		GaussWeight[9]=0.072757916845; 
		GaussWeight[10]=0.072757916845; 
		GaussWeight[11]=0.072757916845; 
		GaussWeight[12]=0.072757916845; 
		GaussWeight[13]=0.028327242531; 
		GaussWeight[14]=0.028327242531; 
		GaussWeight[15]=0.028327242531; 
		GaussWeight[16]=0.028327242531; 
		GaussWeight[17]=0.028327242531; 
		GaussWeight[18]=0.028327242531; 
		GaussWeight[19]=0.009421666964; 
		GaussWeight[20]=0.009421666964; 
		GaussWeight[21]=0.009421666964; 
		GaussWeight[22]=0.009421666964; 
		GaussWeight[23]=0.009421666964; 
		GaussWeight[24]=0.009421666964; 
        break;
     case 12:
        GaussPnt.resize(33);
		GaussWeight.resize(33);
		GaussPnt[0][0]=0.488217389774;
		GaussPnt[1][0]=0.488217389774;
		GaussPnt[2][0]=0.023565220452;
		GaussPnt[3][0]=0.439724392294;
		GaussPnt[4][0]=0.439724392294;
		GaussPnt[5][0]=0.120551215411;
		GaussPnt[6][0]=0.271210385012;
		GaussPnt[7][0]=0.271210385012;
		GaussPnt[8][0]=0.457579229976;
		GaussPnt[9][0]=0.127576145542;
		GaussPnt[10][0]=0.127576145542;
		GaussPnt[11][0]=0.744847708917;
		GaussPnt[12][0]=0.021317350453;
		GaussPnt[13][0]=0.021317350453;
		GaussPnt[14][0]=0.957365299094;
		GaussPnt[15][0]=0.275713269686;
		GaussPnt[16][0]=0.608943235780;
		GaussPnt[17][0]=0.115343494535;
		GaussPnt[18][0]=0.608943235780;
		GaussPnt[19][0]=0.275713269686;
		GaussPnt[20][0]=0.115343494535;
		GaussPnt[21][0]=0.281325580990;
		GaussPnt[22][0]=0.695836086788;
		GaussPnt[23][0]=0.022838332222;
		GaussPnt[24][0]=0.695836086788;
		GaussPnt[25][0]=0.281325580990;
		GaussPnt[26][0]=0.022838332222;
		GaussPnt[27][0]=0.116251915908;
		GaussPnt[28][0]=0.858014033544;
		GaussPnt[29][0]=0.025734050548;
		GaussPnt[30][0]=0.858014033544;
		GaussPnt[31][0]=0.116251915908;
		GaussPnt[32][0]=0.025734050548;
		GaussPnt[0][1]=0.488217389774;
		GaussPnt[1][1]=0.023565220452;
		GaussPnt[2][1]=0.488217389774;
		GaussPnt[3][1]=0.439724392294;
		GaussPnt[4][1]=0.120551215411;
		GaussPnt[5][1]=0.439724392294;
		GaussPnt[6][1]=0.271210385012;
		GaussPnt[7][1]=0.457579229976;
		GaussPnt[8][1]=0.271210385012;
		GaussPnt[9][1]=0.127576145542;
		GaussPnt[10][1]=0.744847708917;
		GaussPnt[11][1]=0.127576145542;
		GaussPnt[12][1]=0.021317350453;
		GaussPnt[13][1]=0.957365299094;
		GaussPnt[14][1]=0.021317350453;
		GaussPnt[15][1]=0.608943235780;
		GaussPnt[16][1]=0.115343494535;
		GaussPnt[17][1]=0.275713269686;
		GaussPnt[18][1]=0.275713269686;
		GaussPnt[19][1]=0.115343494535;
		GaussPnt[20][1]=0.608943235780;
		GaussPnt[21][1]=0.695836086788;
		GaussPnt[22][1]=0.022838332222;
		GaussPnt[23][1]=0.281325580990;
		GaussPnt[24][1]=0.281325580990;
		GaussPnt[25][1]=0.022838332222;
		GaussPnt[26][1]=0.695836086788;
		GaussPnt[27][1]=0.858014033544;
		GaussPnt[28][1]=0.025734050548;
		GaussPnt[29][1]=0.116251915908;
		GaussPnt[30][1]=0.116251915908;
		GaussPnt[31][1]=0.025734050548;
		GaussPnt[32][1]=0.858014033544;
		GaussWeight[0]=0.025731066440;
		GaussWeight[1]=0.025731066440;
		GaussWeight[2]=0.025731066440;
		GaussWeight[3]=0.043692544538;
		GaussWeight[4]=0.043692544538;
		GaussWeight[5]=0.043692544538; 
		GaussWeight[6]=0.062858224218;
		GaussWeight[7]=0.062858224218; 
		GaussWeight[8]=0.062858224218; 
		GaussWeight[9]=0.034796112931; 
		GaussWeight[10]=0.034796112931; 
		GaussWeight[11]=0.034796112931; 
		GaussWeight[12]=0.006166261052; 
		GaussWeight[13]=0.006166261052; 
		GaussWeight[14]=0.006166261052; 
		GaussWeight[15]=0.040371557766; 
		GaussWeight[16]=0.040371557766; 
		GaussWeight[17]=0.040371557766; 
		GaussWeight[18]=0.040371557766; 
		GaussWeight[19]=0.040371557766; 
		GaussWeight[20]=0.040371557766; 
		GaussWeight[21]=0.022356773202; 
		GaussWeight[22]=0.022356773202; 
		GaussWeight[23]=0.022356773202; 
		GaussWeight[24]=0.022356773202; 
		GaussWeight[25]=0.022356773202; 
		GaussWeight[26]=0.022356773202; 
		GaussWeight[27]=0.017316231109; 
		GaussWeight[28]=0.017316231109; 
		GaussWeight[29]=0.017316231109; 
		GaussWeight[30]=0.017316231109; 
		GaussWeight[31]=0.017316231109; 
		GaussWeight[32]=0.017316231109; 
        break;
     default:
        std::cout<<"There are no Gauss quadrature points with accuracy of "<<i<<"-th order"<<std::endl;
  }
}

// 成员函数

void TmpEle::buildTE(int i)
{
	TmpEle tmpEle(i);
	for (auto x : tmpEle.Pnt)
		Pnt.push_back(x);
	for (auto x : tmpEle.GaussPnt)
		GaussPnt.push_back(x);
	for (auto x : tmpEle.GaussWeight)
		GaussWeight.push_back(x);
}

Point TmpEle::getPnt(int i)
{
	return Pnt[i];
}

std::vector<Point> TmpEle::getGaussPnt()
{
	return GaussPnt;
}

std::vector<double> TmpEle::getGaussWeight()
{
	return GaussWeight;
}

double TmpEle::getArea()
{
	return 0.5*((Pnt[1][0] - Pnt[0][0]) * (Pnt[2][1] - Pnt[0][1])
				- (Pnt[1][1] - Pnt[0][1]) * (Pnt[2][0] - Pnt[0][0]));
}

// 标准单元到普通单元
Point TmpEle::Local_to_Global(const Point & lp, const std::vector<Point> & gv) const
{
	Point gp;
	double lambda[3];
	double area = AREA(Pnt[0], Pnt[1], Pnt[2]);
	lambda[0] = AREA(lp, Pnt[1], Pnt[2]) / area;
	lambda[1] = AREA(lp, Pnt[2], Pnt[0]) / area;
	lambda[2] = AREA(lp, Pnt[0], Pnt[1]) / area;
	gp[0] = lambda[0] * gv[0][0] + lambda[1] * gv[1][0] + lambda[2] * gv[2][0];
	gp[1] = lambda[0] * gv[0][1] + lambda[1] * gv[1][1] + lambda[2] * gv[2][1];
	return gp;
}

// 普通单元到标准单元
Point TmpEle::Global_to_Local(const Point& gp, const std::vector<Point>& gv) const
{
	Point lp;
	double area = AREA(gv[0],gv[1],gv[2]);

	lp[0] = ((gv[2][1] - gv[0][1]) * (gp[0] - gv[0][0]) - 
		(gv[2][0] - gv[0][0]) * (gp[1] - gv[0][1])) / (2 * area);
	lp[1] = (-(gv[1][1] - gv[0][1]) * (gp[0] - gv[0][0]) +
		(gv[1][0] - gv[0][0]) * (gp[1] - gv[0][1])) / (2 * area);

	return lp;
}

// 标准单元到普通单元
std::vector<Point> TmpEle::Local_to_Global(const std::vector<Point>& lp, const std::vector<Point> & gv) const
{
	std::vector<Point> gp(lp.size());
	double area = AREA(Pnt[0], Pnt[1], Pnt[2]);
	for(int i=0; i<gp.size(); i++){
		double lambda[3];
		lambda[0] = AREA(lp[i], Pnt[1], Pnt[2]) / area;
		lambda[1] = AREA(lp[i], Pnt[2], Pnt[0]) / area;
		lambda[2] = AREA(lp[i], Pnt[0], Pnt[1]) / area;
		gp[i][0] = lambda[0] * gv[0][0] + lambda[1] * gv[1][0] + lambda[2] * gv[2][0];
		gp[i][1] = lambda[0] * gv[0][1] + lambda[1] * gv[1][1] + lambda[2] * gv[2][1];
	}
	return gp;
}

// 普通单元到标准单元
std::vector<Point> TmpEle::Global_to_Local(const std::vector<Point>& gp, 
	const std::vector<Point>& gv) const
{
	std::vector<Point> lp;
	double area = AREA(gv[0],gv[1],gv[2]);

	for (auto point : gp){
		Point p;
		p[0] = ((gv[2][1] - gv[0][1]) * (point[0] - gv[0][0]) - 
				(gv[2][0] - gv[0][0]) * (point[1] - gv[0][1])) / (2 * area);
		p[1] = (-(gv[1][1] - gv[0][1]) * (point[0] - gv[0][0]) + 
				(gv[1][0] - gv[0][0]) * (point[1] - gv[0][1])) / (2 * area);
		lp.push_back(p);
	}
	return lp;
}

// 从标准单元到普通单元的雅克比矩阵
double TmpEle::Local_to_Global_jacobian(const Point& lp, const std::vector<Point>& gv) const
{
	// lp只是提供维度
	double jacobian =  AREA(gv[0],gv[1],gv[2]) / AREA(Pnt[0], Pnt[1], Pnt[2]);
	return jacobian;
}

// 从普通单元到标准单元的雅克比矩阵
double TmpEle::Global_to_Local_jacobian(const Point& lp, const std::vector<Point>& gv) const
{
	// lp只是提供维度
	double jacobian = AREA(Pnt[0],Pnt[1],Pnt[2]) / AREA(gv[0],gv[1],gv[2]);
	return jacobian;
}

// 从标准单元到普通单元的雅克比矩阵
std::vector<double> TmpEle::Local_to_Global_jacobian(const std::vector<Point>& lp, const std::vector<Point>& gv) const
{
	std::vector<double> gj(lp.size());
	double larea = AREA(Pnt[0], Pnt[1], Pnt[2]);
	double garea = AREA(gv[0], gv[1], gv[2]);
	for (int i = 0; i < gj.size(); i++){
		gj[i] = garea / larea;
	}
	return gj;
}

// 从普通单元到标准单元的雅克比矩阵
std::vector<double> TmpEle::Global_to_Local_jacobian(const std::vector<Point>& lp, const std::vector<Point>& gv) const
{
	std::vector<double> gj(lp.size());
	double larea = AREA(Pnt[0], Pnt[1], Pnt[2]);
	double garea = AREA(gv[0], gv[1], gv[2]);
	for (int i = 0; i < gj.size(); i++){
		gj[i] = larea / garea;
	}
	return gj;
}
