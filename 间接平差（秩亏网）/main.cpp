#include"MAT.h"
#include"JJPCH_zhikui.h"
#include<iostream>

using namespace std;

int main()
{
	/*
	�����������������ȷ�ԣ�
	A��	1 2 1
		2 4 5
		3 6 4
	��A�Ĺ��������
	*/

	MAT A(3, 3);
	A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 1;
	A(1, 0) = 2; A(1, 1) = 4; A(1, 2) = 5;
	A(2, 0) = 3; A(2, 1) = 6; A(2, 2) = 4;
	cout << A << endl;
	MAT B = A.G_inv1();//���ȷ�
	cout << B << endl;
	cout << A*B*A << endl;
	MAT C = A.G_inv2();
	cout << C << endl;
	cout << A*C*A << endl;//���ȷֽⷨ
	cout << endl << endl;

	//ʹ��JJPCH_zhikui������ȿ����������ƽ�����
	JJPCH_zhikui aa;
	aa.deal("data.txt", "result.txt");
	return 0;
}