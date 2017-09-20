// corr_0.cpp : 定义控制台应用程序的入口点。
//
//
//

#include "stdafx.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<time.h>
#include "corr.h"
#include "wavelet.h"

using namespace std;

float lorenz[LORENZ_LEN];	//当作全局数组，只需定义一次
int maxi_ref; //全局变量

/**
* 求输出数组的最大值
* 返回最大值下标maxi
**/
int max(float a[], int size)
{
	int maxi = 0;
	for (int i = 0; i < size; ++i)
	{
		if (abs(a[maxi]) < abs(a[i]))
			maxi = i;
	}
	return maxi; 
}

/**
* 求点积函数
**/
float dot_product(const float *x, const float *y, int length)
{
	float sum = 0;
	for (int i = 0; i < length; ++i)
	{
		sum += x[i] * y[i];
	}

	return sum;
}
/**
* corr函数，输入为每一列的数组，及其长度
* 内部调用DSPF_sp_convol_cn相关函数与标准洛伦兹相关
* 输出为峰值横坐标
**/
float corr(float *corr_input, float *lorenz, short input_len)
{
	float corr_output[CORR_OUTPUT_LEN];

	int maxi1, maxi2;	  //输出数组的最大值的下标值


	//调用dsp相关函数，得到output结果
	DSPF_sp_convol_cn(corr_input, lorenz, corr_output, input_len, LORENZ_LEN, CORR_OUTPUT_LEN, maxi_ref);

	//先把结果放到txt里，用其他软件，如matlab画出来看看
	// if (result.is_open())
	// {
	// for (int i = 0; i < CORR_OUTPUT_LEN; ++i)
	// {
	// result << fixed << setprecision(8) << corr_output[i] << endl;
	// }
	// }

	//求出输出序列最大的两个值maxi1 和 maxi2
	maxi1 = max(corr_output, CORR_OUTPUT_LEN);

	 if(corr_output[maxi1 + 1] > corr_output[maxi1 - 1])
	 {
		maxi2 = maxi1 +1;
	 }
	 else 
	 {
		maxi2 = maxi1;
		maxi1 = maxi1 - 1;
	 }

	// cout << "The output maxi1 and maxi12 : " << maxi1 <<"  "  << maxi2 << endl;

	//后续操作，插值处理
	int T, T1;	
	float maxi,frac;
	float c0_T, c0_T1;

	//先对齐采样信号和输出信号，得到采样信号对应的峰值横坐标maxi
	maxi = maxi1;    

	//找到采样信号峰值与参考信号峰值对齐时第一个点位置T
	T = maxi - (LORENZ_LEN - maxi_ref);
	T1 = T + 1;

	//越界的话，直接取frac = 0.5
	if (T < 0 || T + LORENZ_LEN > input_len)
		frac = 0.5;

	 //没越界则计算插值
	else
	{
		c0_T = corr_output[maxi1];
		c0_T1 = corr_output[maxi2];

		//求点积，信号长度要比参考信号长，避免越界,若洛伦兹取100，则已经越界，第一列取出T为23，23+100 > 121
		float cT_T, cT_T1, cT1_T1, denom;
		float *interpolation_input_T, *interpolation_input_T1;
		interpolation_input_T = &corr_input[T];
		interpolation_input_T1 = &corr_input[T1];

		//求数组,不需要这一步
		/*for (int i = 0; i < LORENZ_LEN; ++i)
		{
			interpolation_input_T[i] = *(&corr_input[T + i]);
		}

		for (int i = 0; i < LORENZ_LEN; ++i)
		{
			interpolation_input_T1[i] = *(&corr_input[T1 + i]);
		}*/

		cT_T = dot_product(interpolation_input_T, interpolation_input_T, LORENZ_LEN);
		cT_T1 = dot_product(interpolation_input_T, interpolation_input_T1, LORENZ_LEN);
		cT1_T1 = dot_product(interpolation_input_T1, interpolation_input_T1, LORENZ_LEN);

		denom = c0_T1*(cT_T - cT_T1) + c0_T*(cT1_T1 - cT_T1);   //分母
		//if (abs(denom) > 0.01)  这个分母>0.01不符合这里，改成先计算插值，插值如果大于1就定为0.5
		frac = (c0_T1*cT_T - c0_T*cT_T1) / denom;
		//frac大于1不取，小于0插值的时候有问题，插值算法使对的，这列数据有问题，可能后面有比较大的数，如第一列最后几个数莫名其妙很大
		if (frac > 1 || frac < 0)
			frac = 0.5;
	}
	maxi = maxi + frac;
	return maxi;
}



/**
* fitting函数: 参数collection_rows 表示采集了多少行， datacount_cols表示有多少列
* input指向输入的数据为一维数组，output保存输出结果,外部函数可以使用output
* 对输入数据的每一列先进行小波滤波，再进行相关运算，求出峰值横坐标
**/

 void fitting(int collection_rows, int datacount_cols, float *input, float *output)
{
	int space = 0;
	float **p = new float*[collection_rows];					//p指向二维数组
	float *corr_input = new float[collection_rows];
	//output = new float[datacount_cols];    //注意这个内存在哪释放？
	float maxi;

	//标准洛伦兹型的三个参数，a为系数， b为半高线宽，c为最高点的横坐标即中心频率
	float a = 0.1, b = 0.052, c = 0.06;
	//输入洛伦兹曲线数据，先得到横坐标F[]值，范围为-0.06~0.18，步长0.002 
	float F[LORENZ_LEN];
	for (int i = 0; i < LORENZ_LEN; i++)
	{
		F[i] = 0 + 0.002*i;
	}
	//计算标准洛伦兹曲线值
	for (int i = 0; i < LORENZ_LEN; i++)
	{
		lorenz[i] = a*(pow((b*0.5), 2)) / ((pow((b*0.5), 2)) + (pow((F[i] - c), 2)));
	}
	maxi_ref = max(lorenz, LORENZ_LEN);


	//给每一行初始化datacount_cols个数
	for (int i = 0; i < collection_rows; ++i)
		p[i] = new float[datacount_cols];

	//将输入的矩阵赋值,把input赋给p[][]
	for (int i = 0; i < collection_rows; ++i)
	{
		for (int j = 0; j < datacount_cols; ++j)
		{
			p[i][j] = input[j + space] ;
		}
		space = space + datacount_cols;
	}

	//对p[][]的每一列调用相关函数，得到每一列的峰值横坐标,添加到output
	for (int j = 0; j < datacount_cols; ++j)
	{
		for (int i = 0; i < collection_rows; ++i)
		{
			corr_input[i] = p[i][j]; //把一列数赋给corr_input[]作为相关的输入
		}

	    dwt_filter(corr_input, collection_rows); //调用小波滤波函数，先进行滤波
		maxi = corr(corr_input, lorenz, collection_rows); //调用相关函数，得到峰值横坐标
		output[j] = maxi;        //把结果添加到output	
	}

	for (int i = 0; i < collection_rows; ++i)
		delete[]p[i];
	
	delete[]p;
	delete[]corr_input;
}
