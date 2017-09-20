
#include "stdafx.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include <iomanip>
#include "wavelet.h"
#include "corr.h"
using namespace std;

int main()
{
	clock_t cstart, cend;
	ifstream input_file;
	input_file.open("mydata.txt");

	ofstream final_result;
	final_result.open("final_result.txt");


	int collection_rows = INPUT_LEN;
	int datacount_cols = DATA_COLS;
	float *gets = new float[collection_rows*datacount_cols];    //gets是一个一维数组
	float *output = new float[datacount_cols];


	//读取输入文件数据
	if (!input_file.is_open())
	{
		cout << "Error opening file" << endl;
		exit(1);
	}

	for (int i = 0; ; ++i)
	{
		input_file >> gets[i];
		if (input_file.eof())
			break;
	}

	cstart = clock();

	//调用fitting函数,每一列滤波后进行相关
	fitting(collection_rows, datacount_cols, gets, output);
	//先把结果放到txt里
	if (final_result.is_open())
	{
		for (int i = 0; i < datacount_cols; ++i)
		{
			//每一列下标从0开始，所以output[i]的下标即直接*0.002就是增量
			final_result << fixed << setprecision(8) << 10.8 + (0.002*output[i]) << endl;
			//final_result << fixed << setprecision(8) << output[i] << endl;
		}
	}


	delete[]gets;
	delete[]output;
	input_file.close();
	final_result.close();

	cend = clock();
	cout << "time is " << cend - cstart << endl;

	return 0;
}


