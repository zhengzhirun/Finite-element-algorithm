#include "Matrix.h"

#define PI (4.0*atan(1.0))

////////////////////////////////////////////////////////////////////////////////
// 类Point
////////////////////////////////////////////////////////////////////////////////

//非成员接口函数
std::ostream& operator<< (std::ostream& os, const Point& p)
{
	for (auto i : p.x)
		os << i << '\t';
	return os;
}

std::istream& operator>> (std::istream& is, Point& p)
{
	for (auto &i : p.x)
		is >> i;
	return is;
}

Point operator+ (const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = p1.x[i] + p2.x[i];
	return p;
}

Point operator- (const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = p1.x[i] - p2.x[i];
	return p;
}

Point midpoint(const Point& p1, const Point& p2)
{
	Point p;
	for (int i = 0; i != 2; ++i)
		p.x[i] = (p1.x[i] + p2.x[i]) / 2.0;
	return p;
}

double distance(const Point& p1, const Point& p2)
{
	double p = 0.;
	for (int i = 0; i != 2; ++i)
		p += (p1.x[i] - p2.x[i]) * (p1.x[i] - p2.x[i]);
	return sqrt(p);
}

Point barycenter(const std::vector<Point>& p, const double * w)
{
	double bc[2] = {0,0};
	int k = p.size();
	if (w == NULL){
		for (int i = 0; i < k; i++){
			bc[0] += p[i][0];
			bc[1] += p[i][1];
		}
		bc[0] /= k;
		bc[1] /= k;
	}
	else{
		double sw = 0;
		for (int i = 0; i < k; i++) sw += w[i];
		for (int i = 0; i < k; i++){
			bc[0] += w[i] * p[i][0];
			bc[1] += w[i] * p[i][1];
		}
		bc[0] /= sw;
		bc[1] /= sw;
	}
	return bc;
}

double AREA(const Point& p1, const Point& p2, const Point& p3)
{
	double area = ((p2.x[0] - p1.x[0]) * (p3.x[1] - p1.x[1]) - 
		(p3.x[0] - p1.x[0]) * (p2.x[1] - p1.x[1])) * 0.5;
	if (area < 0)
		area = -area;
	return area;
}
// 构造函数
Point::Point()
{
	for (int i = 0; i != 2; ++i)
		x[i] = 0;

}

Point::Point(const double* data)
{
	for (int i = 0; i != 2; ++i)
		x[i] = data[i];
}

Point::Point(const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] = p.x[i];
}

Point::Point(const double& v1, const double& v2)
{
	x[0] = v1;
	x[1] = v2;
}
// 成员函数
Point& Point::operator= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] = p.x[i];
	return *this;
}

double& Point::operator[] (int i)
{
	return x[i];
}

const double& Point::operator[] (int i) const
{
	return x[i];
}

double Point::length() const
{
	double v = 0.0;
	for (auto i : x)
		v += i * i;
	return sqrt(v);
}

Point& Point::operator+= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] += p.x[i];
	return *this;
}

Point& Point::operator-= (const Point& p)
{
	for (int i = 0; i != 2; ++i)
		x[i] -= p.x[i];
	return *this;
}

Point& Point::operator*= (const double& s)
{
	for (auto &i : x)
		i *= s;
	return *this;
}

Point& Point::operator/= (const double& s)
{
	for (auto &i : x)
		i /= s;
	return *this;
}

////////////////////////////////////////////////////////////////////////////////
//类TmpEle
////////////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////////////
//类Mesh
////////////////////////////////////////////////////////////////////////////////

// 非成员函数
std::ostream& operator<< (std::ostream& os, const Mesh& M)
{
	for (auto i : M.Pnt)
		os << i << '\n';

	for (auto &i : M.Ele){
		for (auto j : i)
			os << j << '\t';
		os << '\n';
	}

	for (auto i : M.BndPnt)
		os << i << '\n';
		os << '\n';
	return os;
}

// 成员函数
void Mesh::readData(const std::string& f)
{
	int i;
	std::ifstream is(f,std::ifstream::in); // 以只读的模式打开
	is >> i;
	Pnt.resize(i);
	for (int j = 0; j < i; j++)
		is >> Pnt[j][0] >> Pnt[j][1];

	int n;
	is >> n;
	Ele.resize(n);
	for (int j = 0; j < n; j++){
		Ele[j].resize(3);
		is >> Ele[j][0] >> Ele[j][1] >> Ele[j][2];
	}

	int bn;
	is >> bn;
	BndPnt.resize(bn);
	for (int j = 0; j < bn; j++)
		is >> BndPnt[j];

	is.close();
}

int Mesh::getEleVtx(int i, int j)
{
	return Ele[i][j];
}

std::vector<int> Mesh::getEleVtx(int i)
{
	return Ele[i];
}

Point Mesh::getPnt(int i)
{
	return Pnt[i];
}

std::vector<Point> Mesh::getPnt(std::vector<int>& vt)
{
	std::vector<Point> vec;
	for (int x : vt){
		Point point = Pnt[x];
		vec.push_back(point);
	}
	return vec;
}

std::vector<int> Mesh::getBndPnt()
{
	return BndPnt;
}

int Mesh::n_element()
{
	return Ele.size();
}

int Mesh::n_point()
{
	return Pnt.size();
}

int Mesh::n_boundaryPoint()
{
	return BndPnt.size();
}

void Mesh::neighbor(std::vector<std::vector<int>>& neighbor)
{
    auto n_ele = n_element();   // 得到剖分小三角形的数量
    // 初始化
    neighbor.resize(n_ele);
    for (decltype(neighbor.size()) i = 0; i != n_ele; ++i){
        neighbor[i].resize(3);
    }
    for (decltype(neighbor.size()) i = 0; i != n_ele; ++i){
        // 检测顶点所对应的边是不是边界边
        bool boundary_mark_00 = false, boundary_mark_01 = false;
        bool boundary_mark_10 = false, boundary_mark_11 = false;
        bool boundary_mark_20 = false, boundary_mark_21 = false;
        for (auto boundary : BndPnt){
            // 第一个顶点所对应的边是不是边界边
            if (boundary == Ele[i][1]) boundary_mark_00 = true;
            if (boundary == Ele[i][2]) boundary_mark_01 = true;
            // 第二个顶点所对应的边是不是边界边
            if (boundary == Ele[i][0]) boundary_mark_10 = true;
            if (boundary == Ele[i][2]) boundary_mark_11 = true;
            // 第三个顶点所对应的边是不是边界边
            if (boundary == Ele[i][0]) boundary_mark_20 = true;
            if (boundary == Ele[i][1]) boundary_mark_21 = true;
        }
        // 如果是边界边的话邻居单元就是自己
        if (boundary_mark_00 && boundary_mark_01) neighbor[i][0] = i;
        if (boundary_mark_10 && boundary_mark_11) neighbor[i][1] = i;
        if (boundary_mark_20 && boundary_mark_21) neighbor[i][2] = i;
        
        for (decltype(neighbor.size()) j = 0; j != n_ele; ++j){
            if (i == j) continue;
            bool mark_00 = false, mark_01 = false;
            bool mark_10 = false, mark_11 = false;
            bool mark_20 = false, mark_21 = false;
            for (decltype(neighbor.size()) k = 0; k != 3; ++k){
                // 第一个顶点对应边的邻居单元(如果是边界边的话邻居单元就是自己)
                if (Ele[j][k] == Ele[i][1]) mark_00 = true;
                if (Ele[j][k] == Ele[i][2]) mark_01 = true;
                // 第二个顶点对应边的邻居单元(如果是边界边的话邻居单元就是自己)
                if (Ele[j][k] == Ele[i][0]) mark_10 = true;
                if (Ele[j][k] == Ele[i][2]) mark_11 = true;
                // 第三个顶点对应边的邻居单元(如果是边界边的话邻居单元就是自己)
                if (Ele[j][k] == Ele[i][0]) mark_20 = true;
                if (Ele[j][k] == Ele[i][1]) mark_21 = true;
            }
            if (mark_00 && mark_01) neighbor[i][0] = j;
            if (mark_10 && mark_11) neighbor[i][1] = j;
            if (mark_20 && mark_21) neighbor[i][2] = j; 
        }
        
    }
    
}


////////////////////////////////////////////////////////////////////////////////
//偏微分方程,边界条件和初值条件结构体,PDE
//
////////////////////////////////////////////////////////////////////////////////

// 边界条件(第一类边界条件)

double PDE::u_boundary(const Point& p)
{
	return 0;
}

double PDE::f(const Point& p)
{
	double value;
	value = - 2 * PI * PI * exp(PI * (p[0] + p[1])) * (sin(PI * p[0]) * cos(PI * p[1]) + cos(PI * p[0]) * sin(PI * p[1]));
	return value; 
}

double PDE::u_exact(const Point& p)
{
	double value;
	value = exp(PI * (p[0] + p[1])) * sin(PI * p[0]) * sin(PI * p[1]);
	return value;
}

////////////////////////////////////////////////////////////////////////////////
//类Matrix
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


void Matrix::GaussSeidel(const std::vector<std::vector<double>>& A_matrix, std::vector<double>& u, 
       const std::vector<double>& rhs)
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

void Matrix::SolveLinearASparse(const Eigen::SparseMatrix<double>& A_matrix_sparse,std::vector<double>& u,
       const std::vector<double>& rhs, int k)
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

Matrix::Matrix(const std::string& file, int i)
{
	mesh.readData(file);
	tmpEle.buildTE(i);
}

// 成员函数
// 基函数的构造(面积坐标)
void Matrix::basisValue(const Point& p, const std::vector<Point>& v, std::vector<double>& val)
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

void Matrix::basisValue(const std::vector<Point>& p, const std::vector<Point>& v, 
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

void Matrix::basisGrad(const Point&, const std::vector<Point>& v, 
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

void Matrix::basisGrad(const std::vector<Point>& p, const std::vector<Point>& v, 
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

void Matrix::buildEleMatrix_RHS(int n, std::vector<std::vector<double>> &Ele_matrix, std::vector<double> &Ele_Rhs)
{
	Ele_matrix.resize(3);  // 单元刚度矩阵(左端),维度是3*3的
	Ele_Rhs.resize(3);   // 单元刚度矩阵(右端),维度是3*1的
	for(int i=0; i<3; i++){
		Ele_matrix[i].resize(3);
	}

	// 得到第n个单元的边
	std::vector<int> EV=mesh.getEleVtx(n);
	// 得到单元的三个点
	std::vector<Point> EP(3);
	EP[0]=mesh.getPnt(EV[0]);
	EP[1]=mesh.getPnt(EV[1]);
	EP[2]=mesh.getPnt(EV[2]);
	// 标准三角单元的面积
	double volume = tmpEle.getArea();
	// 得到标准三角单元的Gauss点
	std::vector<Point> GP=tmpEle.getGaussPnt();
	// 得到标准三角单元的Gauss权重
	std::vector<double> GW=tmpEle.getGaussWeight();
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
		double f_value = pde.f(q_point[i]);
		for(int j = 0; j < basis_value.size(); j++){
			Ele_Rhs[j] += Jxw*f_value*basis_value[j][i];
		}
	}

	for(int i = 0; i < q_point.size(); i++){
		double Jxw = GW[i] * jacobian[i] * volume;
		for(int k=0; k<basis_grad.size(); k++){
			for(int j=0; j<basis_grad.size(); j++){
				Ele_matrix[j][k] += Jxw * (basis_grad[j][i][0] * basis_grad[k][i][0] + basis_grad[j][i][1] *
                               basis_grad[k][i][1]);
			}
		}
	}
}

void Matrix::buildMatrix_RHS_Sparse(Eigen::SparseMatrix<double>& A_matrix_sparse,std::vector<double>& Rhs)
{
    auto n_ele = mesh.n_element();
    auto n_pnt = mesh.n_point();
    A_matrix_sparse.resize(n_pnt,n_pnt);
    Rhs.resize(n_pnt);
    // 定义一个用于存储稀疏矩阵的三元组
    std::vector<Eigen::Triplet<double>> triple;
    for (size_t i = 0; i != n_ele; ++i){
    	std::vector<int> NV = mesh.getEleVtx(i);

        std::vector<std::vector<double>> Ele_matrix(3);
        for (size_t j = 0; j != 3; ++j){
            Ele_matrix[j].resize(3);
        }
        std::vector<double> Ele_Rhs(3);
        buildEleMatrix_RHS(i,Ele_matrix,Ele_Rhs);
        // 稀疏矩阵存储技术
        // 利用经典的三元组插入方式来存储A_matrix
        triple.push_back(Eigen::Triplet<double>(NV[0],NV[0],Ele_matrix[0][0]));
        triple.push_back(Eigen::Triplet<double>(NV[0],NV[1],Ele_matrix[0][1]));
        triple.push_back(Eigen::Triplet<double>(NV[0],NV[2],Ele_matrix[0][2]));
        triple.push_back(Eigen::Triplet<double>(NV[1],NV[0],Ele_matrix[1][0]));
        triple.push_back(Eigen::Triplet<double>(NV[1],NV[1],Ele_matrix[1][1]));
        triple.push_back(Eigen::Triplet<double>(NV[1],NV[2],Ele_matrix[1][2]));
        triple.push_back(Eigen::Triplet<double>(NV[2],NV[0],Ele_matrix[2][0]));
        triple.push_back(Eigen::Triplet<double>(NV[2],NV[1],Ele_matrix[2][1]));
        triple.push_back(Eigen::Triplet<double>(NV[2],NV[2],Ele_matrix[2][2]));

        // 不利用稀疏矩阵存储右端项Rhs
        Rhs[NV[0]] += Ele_Rhs[0];
        Rhs[NV[1]] += Ele_Rhs[1];
        Rhs[NV[2]] += Ele_Rhs[2];
    }
    // 把三元组转换成稀疏矩阵
    A_matrix_sparse.setFromTriplets(triple.begin(),triple.end());
}

// 对于稀疏矩阵A_matrix_sparse的边界条件处理(dirichlet边界条件)
void Matrix::DealBoundary_Sparse(Eigen::SparseMatrix<double>& A_matrix_sparse, const std::vector<double>& uh,
        std::vector<double>& Rhs)
{
    auto n_pnt = mesh.n_point();
    auto n_bp = mesh.n_boundaryPoint();
    std::vector<int> BV = mesh.getBndPnt(); // 得到边界点
    for (size_t i = 0; i != n_bp; ++i){
        Point Pnt_i = mesh.getPnt(BV[i]); // 得到第i个单元的点
        double val = pde.u_boundary(Pnt_i);
        Rhs[BV[i]] = A_matrix_sparse.coeffRef(BV[i],BV[i]) * val;
        for (size_t j = 0; j != n_pnt; ++j){
            if (j != BV[i]){
                A_matrix_sparse.coeffRef(BV[i],j) = 0.0;
            }
        }
        for (size_t j = 0; j != n_pnt; ++j){
            if (j != BV[i]){
                Rhs[j] -= A_matrix_sparse.coeffRef(j,BV[i]) * val;
                A_matrix_sparse.coeffRef(j,BV[i]) = 0.0;
            }
        }
    }
}

double Matrix::ComputerL2Error(const std::vector<double> &f)
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
			double f_value = pde.u_exact(q_point[i]);
			double u_h_val = f[NV[0]]*basis_value[0][i]+f[NV[1]]*basis_value[1][i]+f[NV[2]]*basis_value[2][i];
			double df_val = f_value-u_h_val;
			err += Jxw*df_val*df_val;
		}
	}
	err=sqrt(err);
	return err;
}

void Matrix::GradBasis(std::vector<std::vector<std::vector<double>>>& gradbasis)
{
    auto n_ele = mesh.n_element();  // 三角单元的数量
    auto n_point =  mesh.n_point();    // 三角单元点的个数
    gradbasis.resize(3);
    for (decltype(gradbasis.size()) i = 0; i != 3; ++i)
        gradbasis[i].resize(n_ele);
    for (decltype(gradbasis.size()) i = 0; i != 3; ++i)
        for (decltype(gradbasis.size()) j = 0; j != n_ele; ++j)
            gradbasis[i][j].resize(2);
    
    std::vector<int> EV;    // 存储三角单元的边
    std::vector<Point> EP(3);   // 存储三角单元的三个顶点
    for (decltype(gradbasis.size()) i = 0; i != n_ele; ++i){
        // 得到第i个单元的边
        EV = mesh.getEleVtx(i);
        // 得到单元的三个顶点
        EP[0] = mesh.getPnt(EV[0]);
        EP[1] = mesh.getPnt(EV[1]);
        EP[2] = mesh.getPnt(EV[2]);
       
        gradbasis[0][i][0] = (EP[1][1] - EP[2][1]) / (2 * AREA(EP[0],EP[1],EP[2]));
        gradbasis[0][i][1] = (EP[2][0] - EP[1][0]) / (2 * AREA(EP[0],EP[1],EP[2]));

        gradbasis[1][i][0] = (EP[2][1] - EP[0][1]) / (2 * AREA(EP[0],EP[1],EP[2]));
        gradbasis[1][i][1] = (EP[0][0] - EP[2][0]) / (2 * AREA(EP[0],EP[1],EP[2]));

        gradbasis[2][i][0] = (EP[0][1] - EP[1][1]) / (2 * AREA(EP[0],EP[1],EP[2]));
        gradbasis[2][i][1] = (EP[1][0] - EP[0][0]) / (2 * AREA(EP[0],EP[1],EP[2]));
    }
    
}

void Matrix::Gradu(const std::vector<std::vector<std::vector<double>>>& gradbasis, 
        const std::vector<double>& u_h, std::vector<std::vector<double>>& Du)
{
    auto n_ele = mesh.n_element();  // 三角形单元的数量
    Du.resize(n_ele);
    for (decltype(Du.size()) i = 0; i != n_ele; ++i)
        Du[i].resize(2);

    std::vector<int> EV;    // 存储三角单元的边
    std::vector<Point> EP(3);   // 存储三角单元的三个顶点
    
    for (decltype(Du.size()) i =0; i != n_ele; ++i){
        // 得到第i个单元的边
        EV = mesh.getEleVtx(i);
        // 得到单元的三个顶点
        EP[0] = mesh.getPnt(EV[0]);
        EP[1] = mesh.getPnt(EV[1]);
        EP[2] = mesh.getPnt(EV[2]);

        Du[i][0] = u_h[EV[0]] * gradbasis[0][i][0] + u_h[EV[1]] * gradbasis[1][i][0] + u_h[EV[2]] * gradbasis[2][i][0];
        Du[i][1] = u_h[EV[0]] * gradbasis[0][i][1] + u_h[EV[1]] * gradbasis[1][i][1] + u_h[EV[2]] * gradbasis[2][i][1];
    }
    
}

void Matrix::eleVector(std::vector<std::vector<std::vector<double>>>& ve)
{
    // 初始化
    auto n_ele = mesh.n_element();
    ve.resize(n_ele);
    for (decltype(ve.size()) i = 0; i != n_ele; ++i){
        ve[i].resize(2);
        for (decltype(ve.size()) j = 0; j != 2; ++j){
            ve[i][j].resize(3);
        }
    }
    // 边的向量
    for (decltype(ve.size()) i = 0; i != n_ele; ++i){
        for (decltype(ve.size()) j = 0; j != 2; ++j){
            // 得到第i个单元的边
            std::vector<int> EV = mesh.getEleVtx(i);
            std::vector<Point> EP(3);
            EP[0] = mesh.getPnt(EV[0]);
            EP[1] = mesh.getPnt(EV[1]);
            EP[2] = mesh.getPnt(EV[2]);
            ve[i][j][0] = EP[2][j] - EP[1][j];
            ve[i][j][1] = EP[0][j] - EP[2][j];
            ve[i][j][2] = EP[1][j] - EP[0][j];
        }
    }
}

void Matrix::normalVector(const std::vector<std::vector<std::vector<double>>>& ve, 
        std::vector<std::vector<std::vector<double>>>& ne)
{
    // 初始化
    auto n_ele = mesh.n_element();
    ne.resize(n_ele);
    for (decltype(ne.size()) i = 0; i != n_ele; ++i){
        ne[i].resize(2);
        for (decltype(ne.size()) j = 0; j != 2; ++j){
            ne[i][j].resize(3);
        }
    }
    // 把ve旋转90度得到ne
    for (decltype(ne.size()) i = 0; i != n_ele; ++i){
        ne[i][0][0] = ve[i][1][0]; ne[i][1][0] = -ve[i][0][0];
        ne[i][0][1] = ve[i][1][1]; ne[i][1][1] = -ve[i][0][1];
        ne[i][0][2] = ve[i][1][2]; ne[i][1][2] = -ve[i][0][2];
    }
}

double Matrix::diameter(const std::vector<Point>& EP)
{
    // 求出三角单元三条边的边长
    auto distance_0 = distance(EP[0],EP[1]);
    auto distance_1 = distance(EP[1],EP[2]);
    auto distance_2 = distance(EP[0],EP[2]);

    double max = distance_0;
    if (max < distance_1)
        max = distance_1;
    if (max < distance_2);
        max = distance_2;
    return max;
}


void Matrix::Estimateresidual(std::vector<double>& u_h,std::vector<double>& eta)
{
    auto n_ele = mesh.n_element();  // 三角形单元的数量
    auto n_point = mesh.n_point();  // 三角形点的个数
    // 得到剖分上每一个单元的梯度
    std::vector<std::vector<std::vector<double>>> gradbasis;
    GradBasis(gradbasis);
    // 得到u_h的梯度
    std::vector<std::vector<double>> Du;
    Gradu(gradbasis,u_h,Du);
    // 得到邻居单元的索引
    std::vector<std::vector<int>> neighbor_index;
    mesh.neighbor(neighbor_index);
    // 得到剖分单元上每一条边的向量
    std::vector<std::vector<std::vector<double>>> ve;     
    eleVector(ve);
    // 得到向量的法向量
    std::vector<std::vector<std::vector<double>>> ne;
    normalVector(ve,ne);
    // 对于边界边的跳量为0
    std::vector<double> edgeJump;
    edgeJump.resize(n_ele);
    for (decltype(edgeJump.size()) i = 0; i != n_ele; ++i){
        auto index_0 = neighbor_index[i][0];
        auto index_1 = neighbor_index[i][1];
        auto index_2 = neighbor_index[i][2];

        std::vector<int> EV = mesh.getEleVtx(i);
        std::vector<Point> EP(3);
        EP[0] = mesh.getPnt(EV[0]);
        EP[1] = mesh.getPnt(EV[1]);
        EP[2] = mesh.getPnt(EV[2]);

        edgeJump[i] = pow(((Du[i][0] - Du[index_0][0]) * ne[i][0][0] + (Du[i][1] - Du[index_0][1]) * ne[i][1][0]),2.0) 
            +  pow(((Du[i][0] - Du[index_1][0]) * ne[i][0][1] + (Du[i][1] - Du[index_1][1]) * ne[i][1][1]),2.0)  + 
            pow(((Du[i][0] - Du[index_2][0]) * ne[i][0][2] + (Du[i][1] - Du[index_2][1]) * ne[i][1][2]),2.0);
    }
    // 初始化
    eta.resize(n_ele);
    // 计算残差||r||_{L^2}^2
    for (decltype(eta.size()) i = 0; i != n_ele; ++i ){
        std::vector<int> NV = mesh.getEleVtx(i);
        std::vector<Point> EP(3);
        EP[0] = mesh.getPnt(NV[0]);
        EP[1] = mesh.getPnt(NV[1]);
        EP[2] = mesh.getPnt(NV[2]);
        double volume = tmpEle.getArea();
        // 得到高斯点
        std::vector<Point> GP = tmpEle.getGaussPnt();
        // 得到高斯权重
        std::vector<double> GW = tmpEle.getGaussWeight();
        // 把高斯点投影到普通单元
        std::vector<Point> q_point = tmpEle.Local_to_Global(GP,EP);
        std::vector<double> jacobian = tmpEle.Local_to_Global_jacobian(GP,EP);
        // 得到改三角单元的直径
        double dm = diameter(EP);
        for (decltype(eta.size()) j = 0; j != q_point.size(); ++j){
            double Jxw = GW[j] * jacobian[j] * volume;
            double f_value = pde.f(q_point[j]);
            eta[i] += Jxw * f_value * f_value * dm * dm;
        }
        // 加入跳量
        eta[i] += 0.5 * dm * edgeJump[i];
    }
}









