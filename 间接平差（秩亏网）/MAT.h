#if !defined(CHDADJ_MAT_H__INCLUDED_)                    // MAT.h �ļ�
#define CHDADJ_MAT_H__INCLUDED_     //VS: 2017:11 Version Code:2017.1
//********************************************************************
//***************************�����ඨ��*****************   (lines:743)
//********************************************************************
#include<iostream>                                      // ����/�����
#include<fstream>                                       // �ļ���
using namespace std;
#if !defined PI                                    // ����PI����
#define PI (3.14159265358979312)                   // ����PIΪԲ���ʦ�
#endif                                             // ��������PI��if
#if !defined rou                                   // ����rou����
#define rou (180.0/PI*3600)                        // ����1���ȵ���ֵ
#endif                                             // ��������rou��if
#if !defined UNK                                   // ����UNK����
#define UNK (-PI*PI*PI/rou)                        // ����UNK
#endif                                             // ��������UNK��if
//**********************************************************************
double setf(double a,int t)                      // ����a����tλС����ֵ
{                             // ��Ҫ�������ݵ����λ�����ƣ����ı�a��ֵ
	double b=fabs(a);                              // ����a�ľ���ֵ 
    if(b<0.5*pow(10.0,-t)) return 0;               // ������С�����λ��ֵ
	b*=pow(10.0,t);                                // ��b����10��t�η�
    if(b-floor(b)>=0.5) b=floor(b)+1;              // ��С�����ֵ�ֵ>=0.5
	else  b=floor(b);                              // ��b��С��ȡ��
	b/=pow(10.0,t);                                // ��b����10��t�η�
    if(a<0) b=-b;                                  // ��aΪ�����Ĵ���
    return b;                                      // ����b��ֵ
}
//******************************************************************** 
class MAT                                                    // ������        
{                                // �����������ݽṹ����������㹦��
  public:                                       //���г�Ա����
	MAT();                                      //��ʼ��δ֪���еľ���	
	MAT(int hang,int lie);                      //��ʼ����֪���еľ���
    MAT( MAT &S);                               //�������캯��
	MAT(int hl);                                //���������ξ���  
	~MAT();                                     //��������	
	void SetRow(int h);                         //���þ�������   
	void SetRank(int l);                        //���þ������� 
	void SetRR(int h,int l);                    //���þ����С�����
	int  GetRow();                              //���ؾ�������   
	int  GetRank();                             //���ؾ�������   
    void GetRR(int &row,int &rank);             //���ؾ����������
	void SetElem(int h,int l,double m);         //���þ���h��lλ��Ԫ�ص�ֵ
	double& GetElem(int h, int l);              //���ؾ���h��lλ�õ�Ԫ��
	void SetALL();                              //�����������Ԫ��
	void SetI();                                //���þ���Ϊ��λ��
    void Set1();                                //���þ�������Ԫ��Ϊ1
	void Set0();                                //���þ�������Ԫ��Ϊ0
	bool CheckI();                              //�жϾ����Ƿ�Ϊ��λ��
	MAT & operator = ( MAT &other);             //����Ծ���ֵ	
	bool operator ==(MAT &A);                   //�ж����������Ƿ����
    MAT & operator +=(MAT &other);              //�������ĺ͵��ڲ���
	MAT operator + (MAT &B);                    //���������������
	MAT operator - (MAT &B);                    //���������������
    MAT operator * (MAT &B);                    //���������������
    friend MAT operator * (double A,MAT &B);    //���˾�������
    MAT operator * (double B);                  //�����������
    MAT operator / (double B);                  //�������������
    friend ostream & operator <<(ostream &ou,MAT &A);//�������
    friend istream & operator >>(istream &in,MAT &A);//�������
    double & operator () (int h,int l);              //����h��l��λ�õ�Ԫ��	
    bool Is_Valid();              //�жϾ����Ƿ���Ч(�Ƿ�������)
	MAT inverse();                //���㷽���������ʧ�ܣ��򷵻�0��0�ľ���
 	MAT T();                      //����ת�ü��㣬����ת�þ���
    int R();                      //����������          
    double hl();                  //�����������ʽ��ֵ   
	double Trace();               //�������ļ�
	MAT ChildMAT(int row1,int rank1,int row2,int rank2);  
	                              //����row1-row2~rank1-rank2���Ӿ���
	MAT ChildMAT(int row,int rank);
	                              //����(row,rank)Ԫ�ص��Ӿ���
    MAT _row(int i);              //ȥ����i�еõ��¾���
    MAT _rank(int i);             //ȥ����i�еõ��¾���
	MAT _row_rank(int i,int j);   //ȥ����i�к͵�j�еõ��¾���
	static double Get_Precision();              //���ؼ��㾫��
  	static void Set_Precision(double Pre=1E-12);//���þ�����㾫��
	//****** ���� *************
	MAT G_inv1();	//���ȷ��������
	MAT G_inv2();	//���ȷֽⷨ�������
	MAT P_inv();	//��α��
	//*************************

  private:                        //˽�г�Ա���� 
  	int row;                      //���������������ʱ�±꣺0��row-1��
	int rank;                     //���������������ʱ�±꣺0��rank-1��
    double *elem;                 //����Ԫ�أ���ֱһά��̬���飩
	static double Precision;      //���㾫��
	MAT inverse1();               //��������ֱ�ӱ任,����>5ʱӦ�ã�
	void SetRR();     //��row��rank��ֵ���þ������������Ԫ�ظ�ֵΪ0
	bool exrow(int,int);          //�������������Ԫ��
	bool exrank(int,int);	      //�������������Ԫ��
};
double MAT::Precision=1E-12;
//********************************************************************
 inline bool MAT::Is_Valid()                //�жϾ����Ƿ���Ч(������)
 {
     if(row<1 || rank<1) return false;         // ��Ч��������������
     else return true;                         // ��Ч����                      
 }
//********************************************************************
MAT & MAT:: operator +=(MAT &other)               //�������͵��ڲ���
{
      *this=*this+other;               // ��ǰ����=��ǰ����+other����
	  return *this;                    // ���ص�ǰ����
}
//********************************************************************
MAT MAT::_row(int ii)                           //ȥ����ii�й����¾���
{
	if(ii<0 || ii>row-1) 
	{
	  cout<<"MAT::_row(int ii)���󣺸���������ii����ԭ���������"<<endl;
	  return MAT(0,0);
	}
    MAT C(row-1,rank);                                 //�¾���������1
	for(int i=0;i<row;i++)
	  for(int j=0;j<rank;j++)
	    if(i<ii) C(i,j)=(*this)(i,j);
		else if(i==ii) continue;                       //ɾ����ii��Ԫ�� 
		else C(i-1,j)=(*this)(i,j);
	return C;
}
//********************************************************************
MAT MAT::_rank(int ii)                          //ȥ����ii�й����¾���
{
	if(ii<0 || ii>rank-1) 
	{
	  cout<<"MAT::_rank(int ii)���󣺸�������ii����ԭ���������"<<endl;
	  return MAT(0,0);
	}
    MAT C(row,rank-1);                                 //���巵�ؾ���
	for(int i=0;i<rank;i++)
	 for(int j=0;j<row;j++)
	    if(i<ii) C(j,i)=(*this)(j,i);                  //���ؾ���ֵ
		else if(i==ii) continue;                       //ɾ����ii��Ԫ��                      
		else C(j,i-1)=(*this)(j,i);                    //���ؾ���ֵ
    return C;                                          //���ؽ������
}
//********************************************************************
 MAT MAT::_row_rank(int i,int j)          //ȥ����i�С���j�й����¾���
{
     MAT C=_row(i);                    // ɾ����i��Ԫ�أ��γɷ��ؾ���C                     
     return C._rank(j);                // ɾ��C��j��Ԫ�أ�����C
}
//********************************************************************
ostream & operator<<(ostream &ou,MAT& A)// �����������壨�����������
{   
	for(int i=0;i<A.row;i++)                        // ��ѭ��
	{
		 ou<<"     ";                               // ������ո�
		 for(int j=0;j<A.rank;j++)                  // ��ѭ��
		    ou<<" "<<A(i,j);                       // �����Ԫ��
		 ou<<endl;                                  // ת����һ��  
	}
	return ou;                                      // �������������
 }
//********************************************************************
 istream & operator >>(istream &in,MAT &A)               // ������ȡ��
 {
      for(int i=0;i<A.row;i++)
       for(int j=0;j<A.rank;j++)  
	      in>>A(i,j);                                  // �������Ԫ��
	  return in;                                         // ����������
 }
//********************************************************************
inline MAT::MAT()                                 //����δ֪�����ľ���	
{
       row=rank=0;                                   //������������Ϊ0
       elem=new double;                              //����Ϊ��̬����
	   *elem=0;                                      //�ñ�����ֵΪ0
}
//********************************************************************
inline MAT::~MAT()                                          //��������	
{  
     if(row+rank>1)  delete [] elem;           // �ͷž���Ԫ�ض�̬����   
     else delete elem;
     elem=NULL;
}
//********************************************************************
bool MAT::operator==(MAT &A)                   // �������Ƿ���ȵ��ж� 
{
    if(row!=A.row || rank!=A.rank)	
	   return false;                       // �������������������
    for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	  if(fabs(GetElem(i,j)-A(i,j))>Precision)
	     return false;                   // ����Ԫ�ز������������
    return true;                                          //���������
}
//********************************************************************
inline int MAT::GetRow()                                //���ؾ�������   
{
	 return row;                                   // ���ص�ǰ��������
}
//********************************************************************
inline int MAT::GetRank()                               //���ؾ�������   
{
     return rank;                                  // ���ص�ǰ��������
}
//********************************************************************
inline void MAT::GetRR(int &row,int &rank)         // ���ؾ����С�����
{
     row=this->row;                                    // ��ǰ��������
	 rank=this->rank;                                  // ��ǰ��������
} 
//********************************************************************
double MAT::Trace()                                    // �������ļ�
{
    if( row!=rank ) return UNK;                  // �Ƿ��󲻼����伣
	double S=0;                                  // ���巵������
	for(int i=0;i<row;i++)
	   S+=GetElem(i,i);                          // ����Խ���Ԫ��֮��
	return S;                                    // ����S
}
//********************************************************************
void MAT::SetRR()                                      // ���þ������
{
    delete []elem; elem=NULL;                      // �ͷ�Ԫ�ض�̬����
    elem = new double[row*rank];                     // �¿����ڴ�ռ�
    Set0();                                    // �¾���Ԫ��ȫ����Ϊ��
}  
//********************************************************************
inline void MAT::SetRow(int h)                          //���þ�������   
{
	 row=h;                                                 //������ֵ
	 SetRR();                                       //���·����ڴ�ռ�
}
//********************************************************************
void MAT::SetRank(int l)                                //���þ������� 
{
	 rank=l;                                                //������ֵ
     SetRR();                                       //���·����ڴ�ռ�
}
//********************************************************************
void MAT::SetRR(int row,int rank)                   //���þ����С�����
{
	this->row=row; 
	this->rank=rank;                                   // �����С�����
	SetRR();                                           // �����ڴ�ռ�
}
//********************************************************************	
class node{              //����ࣨ�ھ���������㺯��inverse1()��ʹ�ã�
     public:  
        int dat1;                                      //����������1 
		int dat2;                                      //����������2
   		node(int d1=0,int d2=0);
};
// ********node���Ա��������*****************************************
inline node::node(int d1,int d2)
{
	 dat1=d1;
	 dat2=d2;
}
//********************************************************************	
class stack{               //ջ�ࣨ�ھ���������㺯��inverse1()��ʹ�ã�
   public:
	 node *SZ;                                      //����ջ�����ָ��                               
	 stack(){i=0;SZ=new node[1];}
	 ~stack(){delete []SZ;SZ=NULL;}
	 void push(int,int);                                //����ѹ�뺯��
     node pop();                                        //���嵯������
	 int i;
};
// ********stack���Ա��������****************************************
void stack::push(int t1,int t2)                     //ѹ��һ������ջ 
{   
	 node New_node(t1,t2);        // �����½��
     node *LS=new node[i];
     for(int j=0;j<i;j++)
		LS[j]=SZ[j];
	 delete []SZ;                 // ɾ��ԭ�������
	 SZ=NULL;
     SZ=new node[i+1];            // �����½������
	 for(int j=0;j<i;j++)
	   SZ[j]=LS[j];	              // ����ԭ�������
	 delete []LS;
	 LS=NULL;
	 SZ[i]=New_node;              // �����½��
     i++;
}
//********************************************************************
node stack::pop()                           // ����ջ������inverse1()��
{
    if(i) 
	{
		i--;
		return SZ[i];
	}
    return node(0,0);   
}
//********************************************************************	
inline double & MAT::GetElem(int h,int l)       //���ؾ���h�С�l��Ԫ��
{ 
    if(h<0 || h>row-1 || l<0 || l>rank-1)
	{
	    cout<<" MAT::GetElem(int h, int l)����: �������������ޣ�"<<endl;
	    exit(1);	
	}
	return elem[l+h*rank];         //���ص�l�к�h��Ԫ��(��ʼ��������Ϊ0) 
}
//**********************************************************************	
double & MAT:: operator () (int h,int l)          //���ؾ���h�С�l��Ԫ�� 
{
    return GetElem(h,l);           //���ص�l�к�h�е�Ԫ��(��ʼ������Ϊ0)
}
//**********************************************************************
void MAT::SetElem(int h,int l,double m)                   //���þ���Ԫ��
{   
     GetElem(h,l)=m;                             // h��l��Ӧλ��Ԫ�ظ�ֵ
}
//**********************************************************************	
void MAT::SetALL()                                //���������������Ԫ��
{
     cout<<"������("<<row<<"��"<<rank<<")�����Ԫ��:"<<endl;
     for(int i=0;i<row*rank;i++)
       cin>>elem[i];                            //����Ļ��һ�������Ԫ��
}
//**********************************************************************	
bool MAT::CheckI()                                  //�жϾ����Ƿ�λ�� 
{ 
	if(row!=rank||!Is_Valid()) return false;//������Ч��Ƿ�����ǵ�λ��
	for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	   if(i==j && fabs(GetElem(i,j)-1)>Precision || 
		  i!=j && fabs(GetElem(i,j))>Precision)              // Ԫ���ж�
		   return false;
    return true;
}
//**********************************************************************	
MAT & MAT::operator = ( MAT &other)                    // ����Ծ���ֵ
{	
//1 ����Ƿ����Ը�ֵ�������磺A=A��
   if(this == &other)
	  return *this;
//2 �ͷ�ԭ�е��ڴ���Դ
   if(row+rank>1) delete [] elem;                // �ͷž���Ԫ�ض�̬����
   else delete elem;                                   // �ͷ�ԭ�ڴ�ռ�
   elem=NULL;
   row=other.row;
   rank=other.rank;
//3 �����µ��ڴ���Դ������������
   elem = new double[row*rank];
   for(int i=0;i<row*rank;i++)
	 elem[i]=other.elem[i];                              // ��ӦԪ�ظ�ֵ		
//4 ���ص�ǰ����
   return *this;
}	
//**********************************************************************	
void MAT::Set1()                                   //���þ�������Ԫ��Ϊ1
{
	for(int i=0;i<row*rank;i++)
	  elem[i]=1;                                   //����Ԫ��ֵ������Ϊ1              
}
//**********************************************************************	
void MAT::Set0()                                   //���þ�������Ԫ��Ϊ0
{
	for(int i=0;i<row*rank;i++)
	  elem[i]=0;
}
//**********************************************************************	
void MAT::SetI()                                      //���þ���Ϊ��λ��
{
	if(row!=rank)                                     // �Ƿ��󲻲���
	{
		cout<<"MAT::SetI()���󣺵�ǰ����Ƿ���"<<endl;
	    return;
	}
	for(int i=0;i<row;i++)
	  for(int j=0;j<rank;j++)
		if(i==j) SetElem(i,j,1);                       //�Խ���Ԫ����Ϊ1
		else SetElem(i,j,0);                         //�ǶԽ���Ԫ����Ϊ0
}
//**********************************************************************	
MAT::MAT( MAT &S)                                   //������ʼ�����캯��
{
     row=S.row;
     rank=S.rank;                                      // ������ȸ�ֵ
     elem=new double[row*rank];                        // �������ڴ�ռ�
     for(int i=0;i<row*rank;i++)
	    elem[i]=S.elem[i];                             // ����Ԫ�ظ�ֵ
}		
//**********************************************************************	
MAT::MAT(int hang,int lie)                        //����ȷ���������ľ���
{
   	 row=hang;
   	 rank=lie;                                               // ȷ������
	 elem=new double[hang*lie];                      // ����Ԫ�ض�̬����
     Set0();                                         // ����Ԫ�ظ�ֵΪ0
}
//**********************************************************************	
MAT::MAT(int hl)                                             // ���巽��
{ 
     elem=new double[hl*hl];                        //�����ڴ�ռ� 
	 SetRR(hl,hl);                                  //����������
}
//**********************************************************************	
bool MAT::exrow(int row1,int row2)                  //�������������Ԫ��
{   
	if(row1<0 || row1>row-1 || row2<0 || row2>row-1)
    { 
	    cout<<"MAT::exrow(int row1,int row2) ����::�����������ޣ�"<<endl;
	    return false;
	}
  	if(row1==row2) return true;
	for(int i=0;i<rank;i++) 
	{ 
	    double ex=GetElem(row1,i); 
        SetElem(row1,i,GetElem(row2,i));      
        SetElem(row2,i,ex);
   }
   return true;
}
//**********************************************************************	
bool MAT::exrank(int rank1,int rank2)               //�������������Ԫ��
{
	if(rank1<0 || rank1>rank-1 || rank2<0 || rank2>rank-1)
	{
	    cout<<"MAT::exrank() ����::�����������ޣ�"<<endl;
	    exit(1);
	}
 	if(rank1==rank2) return true;
    for(int i=0;i<row;i++) 
    {
	   double ex=GetElem(i,rank1);
       SetElem(i,rank1,GetElem(i,rank2));      
       SetElem(i,rank2,ex);
    }
    return true;
}
//**********************************************************************	
MAT MAT::T()                             //����ת�ã�����ֵΪת�ú����
{
	MAT C(rank,row);                                 // ����ת�ý������
	for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	   C(j,i)=GetElem(i,j);                          // ����ת�þ����ֵ
	return C;                                        // ����ת�þ���
}
//**********************************************************************	
MAT MAT::operator + (MAT &B)                      //������ͣ����غ;���
{   
    if(row!=B.row || rank!=B.rank)
	{ 
		cout<<"MAT::������Ӽ�����������������ͬ!"<<endl;
		return MAT(0,0);
	}
	MAT C(row,rank);                             //�������������
	for(int i=0;i<row*rank;i++)
		C.elem[i]=this->elem[i]+B.elem[i];       //�����ӦԪ�����
	return C;                                    //���ؽ������
}
//**********************************************************************	
MAT MAT::operator - (MAT &B)              // ����������Ĳ���ز����
{
	if(row!=B.row || rank!=B.rank)
	{ 
		cout<<"MAT::�������������������������ͬ!"<<endl;
		return MAT(0,0);
	}
	MAT C(row,rank);                              //�������������
	for(int i=0;i<row*rank;i++)
		C.elem[i]=this->elem[i]-B.elem[i];        //��ӦԪ�����
	return C;                                     //���ؼ�����
}
//**********************************************************************	
MAT MAT::operator * (MAT &B)               // ������˼��㣬���س˻�����
{ 
   if(this->rank!=B.row) 
   {
	   cout<<"MAT::����˻��������:��������������������! "<<endl;
       return MAT(0,0);
   }
   int Crow=this->row;
   int Crank=B.rank;
   int Brow=B.row;
   MAT C(Crow,Crank);                                 //�������������
   for(int i=0;i<Crow;i++)
	 for(int j=0;j<Crank;j++)
	   for(int k=0;k<Brow;k++)
	     C(i,j)+=this->GetElem(i,k)*B(k,j);           //����˻�����Ԫ��
   return C;                                          //���ؽ������
}
//**********************************************************************	
MAT operator * (double A,MAT &B)           // ���˾�����㣬���ؽ������
{ 
	int Crow=B.row;
	int Crank=B.rank;
	MAT C(Crow,Crank); 
	for(int i=0;i<(Crow*Crank);i++) 
		C.elem[i]=A*B.elem[i];
	return C;
}
//**********************************************************************	
MAT MAT::operator * (double B)              //����������㣬���ؽ������
{   
	 MAT A(row,rank);
  	 A=B*(*this); 
	 return A;
}
//**********************************************************************	
MAT MAT::operator / (double B)           // ������������㣬���ؽ������
{   
	if(B==0) 
	{
		cout<<"MAT::������������� :����Ϊ0��"<<endl;
	    return MAT(0,0);
	}
 	MAT C(row,rank); 
	for(int i=0;i<(row*rank);i++) 
		C.elem[i]=this->elem[i]/B;
	return C;
}
//**********************************************************************	
MAT MAT::inverse()                      //���б任����������㣬��������
{
    if(row!=rank || row==0)                           //�жϾ����Ƿ���
    { 
	   cout<<" MAT::inverse() ����: ��������Ƿ��� "<<endl;
	   return MAT(0,0);           // �������棬����0��
    }
    if(row>=5)                         
      return inverse1();          // �������ϴ󣬲��þ���ֱ�ӱ任������
    MAT C(row,row);             
    C.SetI();	                  // ���嵥λ��
	MAT A=*this;                  // ����任����
	double MAXX,b;
	int h,m, i,j,k;
	m=row;
	for(i=0;i<m;i++)
	{  
		MAXX=0;  h=0;
		if(fabs(A(i,i))<1E-5 && i<m-1)     // �Ծ�����ԪΪ���΢С�Ĵ���
		{
		   for(int l=i+1;l<m;l++)
			 if(fabs(A(l,i))>MAXX) 
			 { 
				 MAXX=fabs(A(l,i));       
                 h=l;                      // �ҵ���Ԫ��Ӧ�����Ԫ������
			 }
			 if(h!=i)
			  for(int k1=0;k1<m;k1++)
			  {
                A(i,k1)+=A(h,k1);
			    C(i,k1)+=C(h,k1);
			  }
		}
	    if(fabs(A(i,i))<Precision)                // ��ԪΪ��
		{
		     cout<<" MAT::inverse()����: ����������� ! "<<endl;
		     return MAT(0,0);                        // �������棬����0��
		}
		b=A(i,i);	                       
		for(j=0;j<m;j++)                             // �����б任������
		{
			A(i,j)/=b;
			C(i,j)/=b;
        }
		for(j=0;j<m;j++)	   
		{  
			b=A(j,i);
			for(k=0;k<m;k++)
			if(i!=j) 
			{
				A(j,k)-=b*A(i,k);
				C(j,k)-=b*C(i,k);
			}
		}
	} 
    return C; 
}
//**********************************************************************	
MAT MAT::inverse1()                  // ����ֱ�ӱ任�����棬�����ٶȽϿ�
{
    MAT C=*this;                                     // ��������������
	stack exrow;
	for(int i=0;i<row;i++)
	{ 
	    int h=i;
	    double b=fabs(C(i,i));
		if(i<row-1)
		{
		  for(int j=i;j<row;j++)
		   if(fabs(C(j,i))>b)
		   {
		       h=j;
			   b=fabs(C(j,i));
		   }
		  if(h!=i)
		  {
			  C.exrow(h,i);
			  exrow.push(h,i);
		  }
		}  
    	if(fabs(C(i,i))<Precision) 
		{
			cout<<"MAT::inverse1()����: ����������� !  \n";
			return MAT(0,0);          // �������棬����0��
		}
		double a=1/C(i,i);
		C(i,i)=a;
	   for(int j=0;j<row;j++)
        if(j!=i)
	    	C(j,i)=-a*C(j,i);
	   for(int j=0;j<row;j++)
		for(int k=0;k<row;k++)
		 if(i!=k && j!=i)
		   C(j,k)=C(j,k)+C(j,i)*C(i,k);			
	   for(int j=0;j<row;j++)
         if(j!=i)
			C(i,j)=a*C(i,j);
	} 
    while(exrow.i)
    {
	    node p=exrow.pop();
		C.exrank(p.dat1,p.dat2);
	};
 	return C;
}
//**********************************************************************
int MAT::R()                                  //���������ȣ������ȵ�ֵ          
{    // ������������ͨ���б任��Ϊ�Խ�����ʽ��Ȼ�����Խ���Ԫ�ط������
   MAT C;
   if(row>rank) C=this->T();
   else C=*this;
   int row=C.row;   
   int rank=C.rank;
   double b,js(0);
   int h,Z(row),i,j,k;
	for(i=0;i<row;i++)
     {
		h=i;
		b=fabs(C(i,i));
   	    if(b<1E-5)                                 //����ԪΪ����С�Ĵ���
		    for(int l=i+1;l<row;l++)
			 if(fabs(C(l,i))>b)
			 {
		        h=l;
				b=fabs(C(l,i));
             }
			if(h!=i) C.exrow(h,i);
	       	if(fabs(C(i,i))<Precision) 
			{ 
				Z--;
				continue;
			} 
	        for(j=i+1;j<row;j++)
			{
			   b=C(j,i)/C(i,i);
	           for(k=0;k<rank;k++) 
			     C(j,k)-=b*C(i,k);                         // �����б任
			}
	  }
	return Z;	 
}	  
//**********************************************************************	
double MAT::hl()                        //�����������ʽ����������ʽ��ֵ  
{      // ������������ͨ���б任��Ϊ�Խ�����ʽ��Ȼ�����Խ���Ԫ�صĳ˻�
  if(row!=rank) 
  { 
	  cout<<"����ʽ���㺯��MAT::hl(MAT &A)���󣺼������Ƿ���"<<endl;
      return UNK;
  }
  MAT C=(*this);
  double b;
  int i,j,k;
	for(i=0;i<row;i++)
      for(j=i+1;j<row;j++)
      { 
		  if(fabs(C(i,i))<1E-5)               //����ԪΪ��Ĵ���
		    for(int l=i+1;l<row;l++)          //��ԪΪ��ʱ����������Ԫ��
			 if(C(l,i)!=0)                    //����Ԫ�ض�Ӧ�м�����Ԫ��
			 {
				 for(int k1=0;k1<rank;k1++)
                   C(i,k1)+=C(l,k1);
			     break;
			 }
		   if(fabs(C(i,i))<Precision)         //�������Ԫ��Ϊ��
			   return 0;				      //������ʽֵΪ��			  
		    b=C.GetElem(j,i)/C(i,i);          //����j�ж���Ԫ�ı���
	        for(k=0;k<rank;k++) 
			 if(i!=j) 
		        C(j,k)-=b*C(i,k);             //�����ݷ���Ԫλ�ñ�Ϊ��
	  }
	  b=1.0;
	  for(i=0;i<row;i++) b*=C(i,i);          //����Խ���Ԫ�س˻�
      return b;
}
//**********************************************************************
MAT MAT::ChildMAT(int row1,int rank1,int row2,int rank2)    
                            // ��(row1,rank1)~(row2,rank2)Ԫ�ع����¾���
{
    if(row1<0 || row1>row-1 || row2<0 || row2>row-1 ||
		rank1<0 || rank1>rank-1||rank2<0 || rank2>rank-1)
	{
		cout<<"MAT::ChildMAT() ���� : ѡȡ�����������ޣ� "<<endl;
	    return MAT(0,0);
	}
    if(row1>row2 || rank1>rank2) 
	{
		cout<<"MAT::ChildMAT() ���� : ��ʼ������>��ֹ��������"<<endl;
		return MAT(0,0);
	}
   	int row;
	row=row2-row1+1;
	int rank;
	rank=rank2-rank1+1;
	MAT C(row,rank);
	for(int i=row1;i<row+row1;i++)
	  for(int j=rank1;j<rank+rank1;j++)
		C.SetElem(i-row1,j-rank1,GetElem(i,j));
	return C;
}
//**********************************************************************	
MAT MAT::ChildMAT(int row,int rank)        // ����(row,rank)Ԫ�ص��Ӿ���
{
	MAT B=*this;
    if(row<0 || row>B.row-1 || rank<0 || rank>B.rank-1) 
	{
	   cout<<"MAT::ChildMAT(int row,int rank)�������������ޣ�"<<endl;
	   return MAT(0,0);
	}
  	MAT C(B.row-1,B.rank-1); 
	for(int i=0;i<row;i++)
	 for(int j=0;j<rank;j++)
	   C.SetElem(i,j,B(i,j));
	for(int i=row+1;i<B.row;i++)
	 for(int j=rank+1;j<B.rank;j++)
	  C.SetElem(i-1,j-1,B(i,j));
    for(int i=0;i<row;i++)
	 for(int j=rank+1;j<B.rank;j++)
	  C.SetElem(i,j-1,B(i,j));
	for(int i=row+1;i<B.row;i++)
	 for(int j=0;j<rank;j++)
	  C.SetElem(i-1,j,B(i,j));
	return C;		
}
//**********************************************************************
inline double MAT::Get_Precision()                        //���ؼ��㾫��
{
     return Precision;
}
//**********************************************************************
inline void MAT::Set_Precision(double precision)          //���ü��㾫�� 
{
     Precision=precision;
}
//**********************************************************************


MAT MAT::G_inv1()	//���ȷ��������
{
	MAT T = *this;
	MAT P(row, rank);
	MAT Q(row, rank);
	P.SetI(); Q.SetI();
	double* temple = new double[row];
	for (int i = 0; i < row; i++) {
		if (T(i, i) == 0) {
			for (int j = i + 1; j < row; j++) {
				if (T(j, i) != 0) {
					for (int x = 0; x < rank; x++)
						temple[x] = T(i, x);
					for (int x = 0; x < rank; x++)
						T(i, x) =T(j,x) ;
					for (int x = 0; x < rank; x++)
						T(j, x) = temple[x];

					for (int x = 0; x < rank; x++)
						temple[x] = P(i, x);
					for (int x = 0; x < rank; x++)
						P(i, x) = P(j, x);
					for (int x = 0; x < rank; x++)
						P(j, x) = temple[x];
					break;
				}
			}
		}
		if (T(i, i) == 0) continue;
		for (int j = i + 1; j < row; j++) {
			double yzi = -T(j, i) / T(i, i);
			for (int x = 0; x < rank; x++)
				T(j, x) += T(i, x)*yzi;
		}
	}

	for (int i = 0; i < rank; i++) {
		if (T(i, i) == 0) {
			for (int j = i + 1; j < rank; j++) {
				if (T(i, j) != 0) {
					for (int x = 0; x < row; x++)
						temple[x] = T(x, j);
					for (int x = 0; x < row; x++)
						T(x, i) = T(x, j);
					for (int x = 0; x < row; x++)
						T(x, j) = temple[x];

					for (int x = 0; x < row; x++)
						temple[x] = Q(x, j);
					for (int x = 0; x < row; x++)
						Q(x, i) = Q(x, j);
					for (int x = 0; x < row; x++)
						Q(x, j) = temple[x];
					break;
				}
			}
		}
	}
	delete[] temple;
	T= P*(*this)*Q;
	MAT A(row, rank);
	int r = R();
	MAT B(r, r);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < r; j++)
			B(i, j) = T.GetElem(i, j);
	B = B.inverse();
	for (int i = 0; i < r; i++)
		for (int j = 0; j < r; j++)
			A(i, j) = B(i, j);
	return Q*A*P;
}

MAT MAT::G_inv2()	//���ȷֽⷨ���������
{
	MAT T = *this;
	MAT P(row, rank);
	MAT Q(row, rank);
	P.SetI(); Q.SetI();
	double* temple = new double[row];
	for (int i = 0; i < row; i++) {
		if (T(i, i) == 0) {
			for (int j = i + 1; j < row; j++) {
				if (T(j, i) != 0) {
					for (int x = 0; x < rank; x++)
						temple[x] = T(i, x);
					for (int x = 0; x < rank; x++)
						T(i, x) = T(j, x);
					for (int x = 0; x < rank; x++)
						T(j, x) = temple[x];

					for (int x = 0; x < rank; x++)
						temple[x] = P(i, x);
					for (int x = 0; x < rank; x++)
						P(i, x) = P(j, x);
					for (int x = 0; x < rank; x++)
						P(j, x) = temple[x];
					break;
				}
			}
		}
		if (T(i, i) == 0) continue;
		for (int j = i + 1; j < row; j++) {
			double yzi = -T(j, i) / T(i, i);
			for (int x = 0; x < rank; x++)
				T(j, x) += T(i, x)*yzi;
		}
	}

	for (int i = 0; i < rank; i++) {
		if (T(i, i) == 0) {
			for (int j = i + 1; j < rank; j++) {
				if (T(i, j) != 0) {
					for (int x = 0; x < row; x++)
						temple[x] = T(x, j);
					for (int x = 0; x < row; x++)
						T(x, i) = T(x, j);
					for (int x = 0; x < row; x++)
						T(x, j) = temple[x];

					for (int x = 0; x < row; x++)
						temple[x] = Q(x, j);
					for (int x = 0; x < row; x++)
						Q(x, i) = Q(x, j);
					for (int x = 0; x < row; x++)
						Q(x, j) = temple[x];
					break;
				}
			}
		}
	}
	delete[] temple;
	T = P*(*this)*Q;
	MAT A(row, rank);
	int r = R();
	MAT A1(row, r), A2(row, rank - r);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < r; j++)
			A1(i, j) = T.GetElem(i, j);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < rank - r; j++)
			A2(i, j) = T.GetElem(i, j+r);
	MAT B = A1;
	MAT BL_1 = (B.T()*B).inverse()*B.T();
	MAT C = BL_1*T;
	MAT CR_1 = C.T()*(C*C.T()).inverse();
	A = CR_1*BL_1;
	return Q*A*P;
}

MAT MAT::P_inv()	//��α��
{
	MAT A=(*this);
	MAT B = (A*A.T()).G_inv1();
	MAT C = (A.T()*A).G_inv1();
	MAT D = A.T()*B*A*C*A.T();
	return D;
}


#endif                             //��������CHDADJ_MAT_H__INCLUDED_��if
