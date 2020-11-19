#pragma once
#include"MAT.h"
#include<fstream>
#include<iostream>
#include<string>

using namespace std;

/*
平差模型：V=BX-L
数据格式：
观测值数 未知参数个数
系数阵B
常数阵L
权阵P
*/

class JJPCH_zhikui {
public:
	bool readdata(string filename);	//读入数据
	bool caculate();	//平差计算 模型为V=BX-L
	bool writedata(string filename);	//以文件形式输出数据
	bool deal(string ifilename, string ofilename);//将输入、计算、输出封装成一个函数，方便使用
private:
	int n;		//观测值个数
	int u;		//未知参数个数
	int t;		//必要观测个数=R(B)
	double m0;	//单位权中误差
	MAT B, L, P, N, Nm_, W, X, V, Qxx, Qvv, Dxx;	//平差要用到的矩阵
};

//读入数据
bool JJPCH_zhikui::readdata(string filename) 
{
	ifstream in(filename, ios::_Nocreate);
	if (!in) return false;
	in >> n >> u;
	B = MAT(n, u);
	L = MAT(n, 1);
	P = MAT(n, n);
	in >> B >> L >> P;
	in.close();
	cout << "观测值个数：" << n << "    未知参数个数：" << u << endl<<endl;
	cout << "系数阵B:" << endl<<B << endl;
	cout << "常数阵L:" << endl << L << endl;
	cout << "权阵P:" <<endl<< P << endl;
	return true;
}

//平差计算
bool JJPCH_zhikui::caculate()
{
	N = B.T()*P*B;
	W = B.T()*P*L;
	if (N.hl() == 0) {
		Nm_ = N*(N*N).G_inv1();
		X = Nm_*W;
		Qxx = N.P_inv();
	}
	else {
		X = N.inverse()*W;
		Qxx = N.inverse();
	}
	V = B*X - L;
	Qvv = P.inverse() - B*Qxx*B.T();
	t = B.R();
	m0 = sqrt((V.T()*P*V).hl() / (n - t));
	Dxx = m0*m0*Qxx;
	return true;
}

//以文件形式输出数据
bool JJPCH_zhikui::writedata(string filename)
{
	ofstream out(filename);
	if (!out)	return false;
	out << "观测值个数：" << n << "    未知参数个数：" << u << endl << endl;
	out << "系数阵B:" << endl << B << endl;
	out << "常数阵L:" << endl << L << endl;
	out << "权阵P:" << endl << P << endl;

	out << "N:" << endl << N << endl;
	out << "W:" << endl << W << endl;
	out << "X:" << endl << X << endl;
	out << "V:" << endl << V << endl;
	out << "单位权中误差m0:" << endl<< "+-" << m0 << endl;
	out << "Dxx:" << endl << Dxx << endl;
	out.close();

	cout << endl;
	cout << "N:" << endl << N << endl;
	cout << "W:" << endl << W << endl;
	cout << "X:" << endl << X << endl;
	cout << "V:" << endl << V << endl;
	cout << "单位权中误差m0:" << endl<< "+-" << m0 << endl;
	cout << "Dxx:" << endl << Dxx << endl;
	return true;
}

//将输入、计算、输出封装成一个函数，方便使用
bool JJPCH_zhikui::deal(string ifilename, string ofilename)
{
	if (!readdata(ifilename)) {
		cout << "错误：读入文件失败！" << endl;
		return false;
	}
	if (!caculate()) {
		cout << "错误：平差计算失败！" << endl;
		return false;
	}
	if (!writedata(ofilename)) {
		cout << "错误：输出数据文件失败！" << endl;
		return false;
	}
	return true;
}