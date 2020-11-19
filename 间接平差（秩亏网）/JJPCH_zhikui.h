#pragma once
#include"MAT.h"
#include<fstream>
#include<iostream>
#include<string>

using namespace std;

/*
ƽ��ģ�ͣ�V=BX-L
���ݸ�ʽ��
�۲�ֵ�� δ֪��������
ϵ����B
������L
Ȩ��P
*/

class JJPCH_zhikui {
public:
	bool readdata(string filename);	//��������
	bool caculate();	//ƽ����� ģ��ΪV=BX-L
	bool writedata(string filename);	//���ļ���ʽ�������
	bool deal(string ifilename, string ofilename);//�����롢���㡢�����װ��һ������������ʹ��
private:
	int n;		//�۲�ֵ����
	int u;		//δ֪��������
	int t;		//��Ҫ�۲����=R(B)
	double m0;	//��λȨ�����
	MAT B, L, P, N, Nm_, W, X, V, Qxx, Qvv, Dxx;	//ƽ��Ҫ�õ��ľ���
};

//��������
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
	cout << "�۲�ֵ������" << n << "    δ֪����������" << u << endl<<endl;
	cout << "ϵ����B:" << endl<<B << endl;
	cout << "������L:" << endl << L << endl;
	cout << "Ȩ��P:" <<endl<< P << endl;
	return true;
}

//ƽ�����
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

//���ļ���ʽ�������
bool JJPCH_zhikui::writedata(string filename)
{
	ofstream out(filename);
	if (!out)	return false;
	out << "�۲�ֵ������" << n << "    δ֪����������" << u << endl << endl;
	out << "ϵ����B:" << endl << B << endl;
	out << "������L:" << endl << L << endl;
	out << "Ȩ��P:" << endl << P << endl;

	out << "N:" << endl << N << endl;
	out << "W:" << endl << W << endl;
	out << "X:" << endl << X << endl;
	out << "V:" << endl << V << endl;
	out << "��λȨ�����m0:" << endl<< "+-" << m0 << endl;
	out << "Dxx:" << endl << Dxx << endl;
	out.close();

	cout << endl;
	cout << "N:" << endl << N << endl;
	cout << "W:" << endl << W << endl;
	cout << "X:" << endl << X << endl;
	cout << "V:" << endl << V << endl;
	cout << "��λȨ�����m0:" << endl<< "+-" << m0 << endl;
	cout << "Dxx:" << endl << Dxx << endl;
	return true;
}

//�����롢���㡢�����װ��һ������������ʹ��
bool JJPCH_zhikui::deal(string ifilename, string ofilename)
{
	if (!readdata(ifilename)) {
		cout << "���󣺶����ļ�ʧ�ܣ�" << endl;
		return false;
	}
	if (!caculate()) {
		cout << "����ƽ�����ʧ�ܣ�" << endl;
		return false;
	}
	if (!writedata(ofilename)) {
		cout << "������������ļ�ʧ�ܣ�" << endl;
		return false;
	}
	return true;
}