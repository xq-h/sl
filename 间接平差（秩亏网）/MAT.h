#if !defined(CHDADJ_MAT_H__INCLUDED_)                    // MAT.h 文件
#define CHDADJ_MAT_H__INCLUDED_     //VS: 2017:11 Version Code:2017.1
//********************************************************************
//***************************矩阵类定义*****************   (lines:743)
//********************************************************************
#include<iostream>                                      // 输入/输出流
#include<fstream>                                       // 文件流
using namespace std;
#if !defined PI                                    // 若无PI定义
#define PI (3.14159265358979312)                   // 定义PI为圆周率π
#endif                                             // 结束关于PI的if
#if !defined rou                                   // 若无rou定义
#define rou (180.0/PI*3600)                        // 定义1弧度的秒值
#endif                                             // 结束关于rou的if
#if !defined UNK                                   // 若无UNK定义
#define UNK (-PI*PI*PI/rou)                        // 定义UNK
#endif                                             // 结束关于UNK的if
//**********************************************************************
double setf(double a,int t)                      // 返回a保留t位小数的值
{                             // 主要用于数据的输出位数控制，不改变a的值
	double b=fabs(a);                              // 计算a的绝对值 
    if(b<0.5*pow(10.0,-t)) return 0;               // 若数据小于输出位数值
	b*=pow(10.0,t);                                // 对b乘以10的t次方
    if(b-floor(b)>=0.5) b=floor(b)+1;              // 若小数部分的值>=0.5
	else  b=floor(b);                              // 对b的小数取整
	b/=pow(10.0,t);                                // 对b除以10的t次方
    if(a<0) b=-b;                                  // 对a为负数的处理
    return b;                                      // 返回b的值
}
//******************************************************************** 
class MAT                                                    // 矩阵类        
{                                // 定义矩阵的数据结构及其相关运算功能
  public:                                       //公有成员声明
	MAT();                                      //初始化未知行列的矩阵	
	MAT(int hang,int lie);                      //初始化已知行列的矩阵
    MAT( MAT &S);                               //拷贝构造函数
	MAT(int hl);                                //定义正方形矩阵  
	~MAT();                                     //析构函数	
	void SetRow(int h);                         //设置矩阵行数   
	void SetRank(int l);                        //设置矩阵列数 
	void SetRR(int h,int l);                    //设置矩阵行、列数
	int  GetRow();                              //返回矩阵行数   
	int  GetRank();                             //返回矩阵列数   
    void GetRR(int &row,int &rank);             //返回矩阵的行列数
	void SetElem(int h,int l,double m);         //设置矩阵h、l位置元素的值
	double& GetElem(int h, int l);              //返回矩阵h、l位置的元素
	void SetALL();                              //键盘输入矩阵元素
	void SetI();                                //设置矩阵为单位阵
    void Set1();                                //设置矩阵所有元素为1
	void Set0();                                //设置矩阵所有元素为0
	bool CheckI();                              //判断矩阵是否为单位阵
	MAT & operator = ( MAT &other);             //矩阵对矩阵赋值	
	bool operator ==(MAT &A);                   //判断两个矩阵是否相等
    MAT & operator +=(MAT &other);              //定义矩阵的和等于操作
	MAT operator + (MAT &B);                    //两个矩阵相加运算
	MAT operator - (MAT &B);                    //两个矩阵相减运算
    MAT operator * (MAT &B);                    //两个矩阵相乘运算
    friend MAT operator * (double A,MAT &B);    //数乘矩阵运算
    MAT operator * (double B);                  //矩阵乘数运算
    MAT operator / (double B);                  //矩阵除以数运算
    friend ostream & operator <<(ostream &ou,MAT &A);//输出矩阵
    friend istream & operator >>(istream &in,MAT &A);//输入矩阵
    double & operator () (int h,int l);              //返回h行l列位置的元素	
    bool Is_Valid();              //判断矩阵是否有效(是否有数据)
	MAT inverse();                //计算方阵的逆阵（若失败，则返回0×0的矩阵）
 	MAT T();                      //矩阵转置计算，返回转置矩阵
    int R();                      //计算矩阵的秩          
    double hl();                  //计算矩阵行列式的值   
	double Trace();               //计算矩阵的迹
	MAT ChildMAT(int row1,int rank1,int row2,int rank2);  
	                              //构造row1-row2~rank1-rank2的子矩阵
	MAT ChildMAT(int row,int rank);
	                              //构成(row,rank)元素的子矩阵
    MAT _row(int i);              //去掉第i行得到新矩阵
    MAT _rank(int i);             //去掉第i列得到新矩阵
	MAT _row_rank(int i,int j);   //去掉第i行和第j列得到新矩阵
	static double Get_Precision();              //返回计算精度
  	static void Set_Precision(double Pre=1E-12);//设置矩阵计算精度
	//****** 新增 *************
	MAT G_inv1();	//降秩法求广义逆
	MAT G_inv2();	//满秩分解法求广义逆
	MAT P_inv();	//求伪逆
	//*************************

  private:                        //私有成员声明 
  	int row;                      //矩阵的行数（访问时下标：0到row-1）
	int rank;                     //矩阵的列数（访问时下标：0到rank-1）
    double *elem;                 //矩阵元素（拉直一维动态数组）
	static double Precision;      //计算精度
	MAT inverse1();               //计算逆阵（直接变换,阶数>5时应用）
	void SetRR();     //用row和rank的值设置矩阵阶数，并将元素赋值为0
	bool exrow(int,int);          //交换矩阵的两行元素
	bool exrank(int,int);	      //交换矩阵的两列元素
};
double MAT::Precision=1E-12;
//********************************************************************
 inline bool MAT::Is_Valid()                //判断矩阵是否有效(有数据)
 {
     if(row<1 || rank<1) return false;         // 无效矩阵（无行列数）
     else return true;                         // 有效矩阵                      
 }
//********************************************************************
MAT & MAT:: operator +=(MAT &other)               //定义矩阵和等于操作
{
      *this=*this+other;               // 当前矩阵=当前矩阵+other矩阵
	  return *this;                    // 返回当前矩阵
}
//********************************************************************
MAT MAT::_row(int ii)                           //去掉第ii行构造新矩阵
{
	if(ii<0 || ii>row-1) 
	{
	  cout<<"MAT::_row(int ii)错误：给定的行数ii超出原矩阵阶数！"<<endl;
	  return MAT(0,0);
	}
    MAT C(row-1,rank);                                 //新矩阵行数减1
	for(int i=0;i<row;i++)
	  for(int j=0;j<rank;j++)
	    if(i<ii) C(i,j)=(*this)(i,j);
		else if(i==ii) continue;                       //删除第ii行元素 
		else C(i-1,j)=(*this)(i,j);
	return C;
}
//********************************************************************
MAT MAT::_rank(int ii)                          //去掉第ii列构造新矩阵
{
	if(ii<0 || ii>rank-1) 
	{
	  cout<<"MAT::_rank(int ii)错误：给定列数ii超出原矩阵阶数！"<<endl;
	  return MAT(0,0);
	}
    MAT C(row,rank-1);                                 //定义返回矩阵
	for(int i=0;i<rank;i++)
	 for(int j=0;j<row;j++)
	    if(i<ii) C(j,i)=(*this)(j,i);                  //返回矩阵赋值
		else if(i==ii) continue;                       //删除第ii列元素                      
		else C(j,i-1)=(*this)(j,i);                    //返回矩阵赋值
    return C;                                          //返回结果矩阵
}
//********************************************************************
 MAT MAT::_row_rank(int i,int j)          //去掉第i行、第j列构造新矩阵
{
     MAT C=_row(i);                    // 删除第i行元素，形成返回矩阵C                     
     return C._rank(j);                // 删除C第j列元素，返回C
}
//********************************************************************
ostream & operator<<(ostream &ou,MAT& A)// 矩阵插入符定义（矩阵流输出）
{   
	for(int i=0;i<A.row;i++)                        // 行循环
	{
		 ou<<"     ";                               // 先输出空格
		 for(int j=0;j<A.rank;j++)                  // 列循环
		    ou<<" "<<A(i,j);                       // 输出该元素
		 ou<<endl;                                  // 转到下一行  
	}
	return ou;                                      // 返回输出流对象
 }
//********************************************************************
 istream & operator >>(istream &in,MAT &A)               // 矩阵提取符
 {
      for(int i=0;i<A.row;i++)
       for(int j=0;j<A.rank;j++)  
	      in>>A(i,j);                                  // 输入矩阵元素
	  return in;                                         // 返回输入流
 }
//********************************************************************
inline MAT::MAT()                                 //构造未知阶数的矩阵	
{
       row=rank=0;                                   //行列数均设置为0
       elem=new double;                              //定义为动态变量
	   *elem=0;                                      //该变量赋值为0
}
//********************************************************************
inline MAT::~MAT()                                          //析构函数	
{  
     if(row+rank>1)  delete [] elem;           // 释放矩阵元素动态数组   
     else delete elem;
     elem=NULL;
}
//********************************************************************
bool MAT::operator==(MAT &A)                   // 两矩阵是否相等的判断 
{
    if(row!=A.row || rank!=A.rank)	
	   return false;                       // 行列数不相等则矩阵不相等
    for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	  if(fabs(GetElem(i,j)-A(i,j))>Precision)
	     return false;                   // 若有元素不相等则矩阵不相等
    return true;                                          //两矩阵相等
}
//********************************************************************
inline int MAT::GetRow()                                //返回矩阵行数   
{
	 return row;                                   // 返回当前矩阵行数
}
//********************************************************************
inline int MAT::GetRank()                               //返回矩阵列数   
{
     return rank;                                  // 返回当前矩阵列数
}
//********************************************************************
inline void MAT::GetRR(int &row,int &rank)         // 返回矩阵行、列数
{
     row=this->row;                                    // 当前矩阵行数
	 rank=this->rank;                                  // 当前矩阵列数
} 
//********************************************************************
double MAT::Trace()                                    // 计算矩阵的迹
{
    if( row!=rank ) return UNK;                  // 非方阵不计算其迹
	double S=0;                                  // 定义返回数据
	for(int i=0;i<row;i++)
	   S+=GetElem(i,i);                          // 计算对角线元素之和
	return S;                                    // 返回S
}
//********************************************************************
void MAT::SetRR()                                      // 设置矩阵阶数
{
    delete []elem; elem=NULL;                      // 释放元素动态数组
    elem = new double[row*rank];                     // 新开辟内存空间
    Set0();                                    // 新矩阵元素全部设为零
}  
//********************************************************************
inline void MAT::SetRow(int h)                          //设置矩阵行数   
{
	 row=h;                                                 //行数赋值
	 SetRR();                                       //重新分配内存空间
}
//********************************************************************
void MAT::SetRank(int l)                                //设置矩阵列数 
{
	 rank=l;                                                //列数赋值
     SetRR();                                       //重新分配内存空间
}
//********************************************************************
void MAT::SetRR(int row,int rank)                   //设置矩阵行、列数
{
	this->row=row; 
	this->rank=rank;                                   // 设置行、列数
	SetRR();                                           // 分配内存空间
}
//********************************************************************	
class node{              //结点类（在矩阵逆阵计算函数inverse1()中使用）
     public:  
        int dat1;                                      //定义结点数据1 
		int dat2;                                      //定义结点数据2
   		node(int d1=0,int d2=0);
};
// ********node类成员函数定义*****************************************
inline node::node(int d1,int d2)
{
	 dat1=d1;
	 dat2=d2;
}
//********************************************************************	
class stack{               //栈类（在矩阵逆阵计算函数inverse1()中使用）
   public:
	 node *SZ;                                      //定义栈顶结点指针                               
	 stack(){i=0;SZ=new node[1];}
	 ~stack(){delete []SZ;SZ=NULL;}
	 void push(int,int);                                //定义压入函数
     node pop();                                        //定义弹出函数
	 int i;
};
// ********stack类成员函数定义****************************************
void stack::push(int t1,int t2)                     //压入一个结点进栈 
{   
	 node New_node(t1,t2);        // 定义新结点
     node *LS=new node[i];
     for(int j=0;j<i;j++)
		LS[j]=SZ[j];
	 delete []SZ;                 // 删除原结点数组
	 SZ=NULL;
     SZ=new node[i+1];            // 定义新结点数组
	 for(int j=0;j<i;j++)
	   SZ[j]=LS[j];	              // 保存原结点数组
	 delete []LS;
	 LS=NULL;
	 SZ[i]=New_node;              // 插入新结点
     i++;
}
//********************************************************************
node stack::pop()                           // 弹出栈（用于inverse1()）
{
    if(i) 
	{
		i--;
		return SZ[i];
	}
    return node(0,0);   
}
//********************************************************************	
inline double & MAT::GetElem(int h,int l)       //返回矩阵h行、l列元素
{ 
    if(h<0 || h>row-1 || l<0 || l>rank-1)
	{
	    cout<<" MAT::GetElem(int h, int l)错误: 矩阵行列数超限！"<<endl;
	    exit(1);	
	}
	return elem[l+h*rank];         //返回第l行和h列元素(起始行列数均为0) 
}
//**********************************************************************	
double & MAT:: operator () (int h,int l)          //返回矩阵h行、l列元素 
{
    return GetElem(h,l);           //返回第l行和h列的元素(起始行列数为0)
}
//**********************************************************************
void MAT::SetElem(int h,int l,double m)                   //设置矩阵元素
{   
     GetElem(h,l)=m;                             // h、l对应位置元素赋值
}
//**********************************************************************	
void MAT::SetALL()                                //键盘输入矩阵所有元素
{
     cout<<"请输入("<<row<<"×"<<rank<<")矩阵的元素:"<<endl;
     for(int i=0;i<row*rank;i++)
       cin>>elem[i];                            //从屏幕逐一输入矩阵元素
}
//**********************************************************************	
bool MAT::CheckI()                                  //判断矩阵是否单位阵 
{ 
	if(row!=rank||!Is_Valid()) return false;//矩阵无效或非方阵则非单位阵
	for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	   if(i==j && fabs(GetElem(i,j)-1)>Precision || 
		  i!=j && fabs(GetElem(i,j))>Precision)              // 元素判断
		   return false;
    return true;
}
//**********************************************************************	
MAT & MAT::operator = ( MAT &other)                    // 矩阵对矩阵赋值
{	
//1 检查是否是自赋值操作（如：A=A）
   if(this == &other)
	  return *this;
//2 释放原有的内存资源
   if(row+rank>1) delete [] elem;                // 释放矩阵元素动态数组
   else delete elem;                                   // 释放原内存空间
   elem=NULL;
   row=other.row;
   rank=other.rank;
//3 分配新的内存资源，并复制数据
   elem = new double[row*rank];
   for(int i=0;i<row*rank;i++)
	 elem[i]=other.elem[i];                              // 对应元素赋值		
//4 返回当前矩阵
   return *this;
}	
//**********************************************************************	
void MAT::Set1()                                   //设置矩阵所有元素为1
{
	for(int i=0;i<row*rank;i++)
	  elem[i]=1;                                   //所有元素值均设置为1              
}
//**********************************************************************	
void MAT::Set0()                                   //设置矩阵所有元素为0
{
	for(int i=0;i<row*rank;i++)
	  elem[i]=0;
}
//**********************************************************************	
void MAT::SetI()                                      //设置矩阵为单位阵
{
	if(row!=rank)                                     // 非方阵不操作
	{
		cout<<"MAT::SetI()错误：当前矩阵非方阵！"<<endl;
	    return;
	}
	for(int i=0;i<row;i++)
	  for(int j=0;j<rank;j++)
		if(i==j) SetElem(i,j,1);                       //对角线元素设为1
		else SetElem(i,j,0);                         //非对角线元素设为0
}
//**********************************************************************	
MAT::MAT( MAT &S)                                   //拷贝初始化构造函数
{
     row=S.row;
     rank=S.rank;                                      // 阶数相等赋值
     elem=new double[row*rank];                        // 开辟新内存空间
     for(int i=0;i<row*rank;i++)
	    elem[i]=S.elem[i];                             // 矩阵元素赋值
}		
//**********************************************************************	
MAT::MAT(int hang,int lie)                        //构造确定行列数的矩阵
{
   	 row=hang;
   	 rank=lie;                                               // 确定阶数
	 elem=new double[hang*lie];                      // 定义元素动态数组
     Set0();                                         // 矩阵元素赋值为0
}
//**********************************************************************	
MAT::MAT(int hl)                                             // 定义方阵
{ 
     elem=new double[hl*hl];                        //分配内存空间 
	 SetRR(hl,hl);                                  //设置行列数
}
//**********************************************************************	
bool MAT::exrow(int row1,int row2)                  //交换矩阵的两行元素
{   
	if(row1<0 || row1>row-1 || row2<0 || row2>row-1)
    { 
	    cout<<"MAT::exrow(int row1,int row2) 错误::交换行数超限！"<<endl;
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
bool MAT::exrank(int rank1,int rank2)               //交换矩阵的两列元素
{
	if(rank1<0 || rank1>rank-1 || rank2<0 || rank2>rank-1)
	{
	    cout<<"MAT::exrank() 错误::交换列数超限！"<<endl;
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
MAT MAT::T()                             //矩阵转置（返回值为转置后矩阵）
{
	MAT C(rank,row);                                 // 构造转置结果矩阵
	for(int i=0;i<row;i++)
     for(int j=0;j<rank;j++)
	   C(j,i)=GetElem(i,j);                          // 设置转置矩阵的值
	return C;                                        // 返回转置矩阵
}
//**********************************************************************	
MAT MAT::operator + (MAT &B)                      //矩阵求和，返回和矩阵
{   
    if(row!=B.row || rank!=B.rank)
	{ 
		cout<<"MAT::矩阵相加计算错误：两矩阵阶数不同!"<<endl;
		return MAT(0,0);
	}
	MAT C(row,rank);                             //定义计算结果矩阵
	for(int i=0;i<row*rank;i++)
		C.elem[i]=this->elem[i]+B.elem[i];       //矩阵对应元素求和
	return C;                                    //返回结果矩阵
}
//**********************************************************************	
MAT MAT::operator - (MAT &B)              // 计算两矩阵的差（返回差矩阵）
{
	if(row!=B.row || rank!=B.rank)
	{ 
		cout<<"MAT::矩阵相减计算错误：两矩阵阶数不同!"<<endl;
		return MAT(0,0);
	}
	MAT C(row,rank);                              //定义计算结果矩阵
	for(int i=0;i<row*rank;i++)
		C.elem[i]=this->elem[i]-B.elem[i];        //对应元素求差
	return C;                                     //返回计算结果
}
//**********************************************************************	
MAT MAT::operator * (MAT &B)               // 矩阵相乘计算，返回乘积矩阵
{ 
   if(this->rank!=B.row) 
   {
	   cout<<"MAT::矩阵乘积计算错误:矩阵阶数不满足相乘条件! "<<endl;
       return MAT(0,0);
   }
   int Crow=this->row;
   int Crank=B.rank;
   int Brow=B.row;
   MAT C(Crow,Crank);                                 //定义计算结果矩阵
   for(int i=0;i<Crow;i++)
	 for(int j=0;j<Crank;j++)
	   for(int k=0;k<Brow;k++)
	     C(i,j)+=this->GetElem(i,k)*B(k,j);           //计算乘积矩阵元素
   return C;                                          //返回结果矩阵
}
//**********************************************************************	
MAT operator * (double A,MAT &B)           // 数乘矩阵计算，返回结果矩阵
{ 
	int Crow=B.row;
	int Crank=B.rank;
	MAT C(Crow,Crank); 
	for(int i=0;i<(Crow*Crank);i++) 
		C.elem[i]=A*B.elem[i];
	return C;
}
//**********************************************************************	
MAT MAT::operator * (double B)              //矩阵乘数计算，返回结果矩阵
{   
	 MAT A(row,rank);
  	 A=B*(*this); 
	 return A;
}
//**********************************************************************	
MAT MAT::operator / (double B)           // 矩阵除以数计算，返回结果矩阵
{   
	if(B==0) 
	{
		cout<<"MAT::矩阵除以数错误 :除数为0！"<<endl;
	    return MAT(0,0);
	}
 	MAT C(row,rank); 
	for(int i=0;i<(row*rank);i++) 
		C.elem[i]=this->elem[i]/B;
	return C;
}
//**********************************************************************	
MAT MAT::inverse()                      //行列变换矩阵求逆计算，返回逆阵
{
    if(row!=rank || row==0)                           //判断矩阵是否方阵
    { 
	   cout<<" MAT::inverse() 错误: 求逆矩阵不是方阵！ "<<endl;
	   return MAT(0,0);           // 不可求逆，返回0阵
    }
    if(row>=5)                         
      return inverse1();          // 如果矩阵较大，采用矩阵直接变换法计算
    MAT C(row,row);             
    C.SetI();	                  // 定义单位阵
	MAT A=*this;                  // 定义变换矩阵
	double MAXX,b;
	int h,m, i,j,k;
	m=row;
	for(i=0;i<m;i++)
	{  
		MAXX=0;  h=0;
		if(fabs(A(i,i))<1E-5 && i<m-1)     // 对矩阵主元为零或微小的处理
		{
		   for(int l=i+1;l<m;l++)
			 if(fabs(A(l,i))>MAXX) 
			 { 
				 MAXX=fabs(A(l,i));       
                 h=l;                      // 找到主元对应列最大元素行数
			 }
			 if(h!=i)
			  for(int k1=0;k1<m;k1++)
			  {
                A(i,k1)+=A(h,k1);
			    C(i,k1)+=C(h,k1);
			  }
		}
	    if(fabs(A(i,i))<Precision)                // 主元为零
		{
		     cout<<" MAT::inverse()错误: 求逆矩阵不满秩 ! "<<endl;
		     return MAT(0,0);                        // 不可求逆，返回0阵
		}
		b=A(i,i);	                       
		for(j=0;j<m;j++)                             // 矩阵行变换法求逆
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
MAT MAT::inverse1()                  // 矩阵直接变换法求逆，计算速度较快
{
    MAT C=*this;                                     // 定义求逆结果矩阵
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
			cout<<"MAT::inverse1()错误: 求逆矩阵不满秩 !  \n";
			return MAT(0,0);          // 不可求逆，返回0阵
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
int MAT::R()                                  //计算矩阵的秩，返回秩的值          
{    // 方法：将矩阵通过行变换变为对角线形式，然后计算对角线元素非零个数
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
   	    if(b<1E-5)                                 //对主元为零或很小的处理
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
			     C(j,k)-=b*C(i,k);                         // 矩阵行变换
			}
	  }
	return Z;	 
}	  
//**********************************************************************	
double MAT::hl()                        //计算矩阵行列式，返回行列式的值  
{      // 方法：将矩阵通过行变换变为对角线形式，然后计算对角线元素的乘积
  if(row!=rank) 
  { 
	  cout<<"行列式计算函数MAT::hl(MAT &A)错误：计算矩阵非方阵！"<<endl;
      return UNK;
  }
  MAT C=(*this);
  double b;
  int i,j,k;
	for(i=0;i<row;i++)
      for(j=i+1;j<row;j++)
      { 
		  if(fabs(C(i,i))<1E-5)               //对主元为零的处理：
		    for(int l=i+1;l<row;l++)          //主元为零时，将本列主元下
			 if(C(l,i)!=0)                    //非零元素对应行加于主元行
			 {
				 for(int k1=0;k1<rank;k1++)
                   C(i,k1)+=C(l,k1);
			     break;
			 }
		   if(fabs(C(i,i))<Precision)         //处理后主元仍为零
			   return 0;				      //则行列式值为零			  
		    b=C.GetElem(j,i)/C(i,i);          //计算j行对主元的倍数
	        for(k=0;k<rank;k++) 
			 if(i!=j) 
		        C(j,k)-=b*C(i,k);             //列数据非主元位置变为零
	  }
	  b=1.0;
	  for(i=0;i<row;i++) b*=C(i,i);          //计算对角线元素乘积
      return b;
}
//**********************************************************************
MAT MAT::ChildMAT(int row1,int rank1,int row2,int rank2)    
                            // 用(row1,rank1)~(row2,rank2)元素构成新矩阵
{
    if(row1<0 || row1>row-1 || row2<0 || row2>row-1 ||
		rank1<0 || rank1>rank-1||rank2<0 || rank2>rank-1)
	{
		cout<<"MAT::ChildMAT() 错误 : 选取起行列数超限！ "<<endl;
	    return MAT(0,0);
	}
    if(row1>row2 || rank1>rank2) 
	{
		cout<<"MAT::ChildMAT() 错误 : 起始行列数>终止行列数！"<<endl;
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
MAT MAT::ChildMAT(int row,int rank)        // 构成(row,rank)元素的子矩阵
{
	MAT B=*this;
    if(row<0 || row>B.row-1 || rank<0 || rank>B.rank-1) 
	{
	   cout<<"MAT::ChildMAT(int row,int rank)错误：行列数超限！"<<endl;
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
inline double MAT::Get_Precision()                        //返回计算精度
{
     return Precision;
}
//**********************************************************************
inline void MAT::Set_Precision(double precision)          //设置计算精度 
{
     Precision=precision;
}
//**********************************************************************


MAT MAT::G_inv1()	//降秩法求广义逆
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

MAT MAT::G_inv2()	//满秩分解法求最广义逆
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

MAT MAT::P_inv()	//求伪逆
{
	MAT A=(*this);
	MAT B = (A*A.T()).G_inv1();
	MAT C = (A.T()*A).G_inv1();
	MAT D = A.T()*B*A*C*A.T();
	return D;
}


#endif                             //结束关于CHDADJ_MAT_H__INCLUDED_的if
