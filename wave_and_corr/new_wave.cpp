// new_wave.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include "wavelet.h"
#include<fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#define  FILTERLEN 16  //sym8的滤波器系数长度
#define  COFE_LENGTH 158 //分解出来的低频、高频系数长度都为158
#define DWT_STAGE	 4
using namespace std;


static const float sym8lowd_q15[] = {
	0.001890, -0.000303, -0.014952, 0.003809, 0.049137, -0.027219, -0.051946, 0.364442, 0.777186, 0.481360,-0.061273, -0.143294, 0.007607, 0.031695, -0.000542, -0.003382, };
static const float sym8highd_q15[] = {
	-0.003382, 0.000542, 0.031695, -0.007607, -0.143294, 0.061273, 0.481360, -0.777186, 0.364442, 0.051946,-0.027219, -0.049137, 0.003809, 0.014952, -0.000303, -0.001890, };
static const float sym8lowr_q15[] = {
	-0.003382, -0.000542, 0.031695, 0.007607, -0.143294, -0.061273, 0.481360, 0.777186, 0.364442, -0.051946,-0.027219, 0.049137, 0.003809, -0.014952, -0.000303, 0.001890, };
static const float sym8highr_q15[] = {
	-0.001890, -0.000303, 0.014952, 0.003809, -0.049137, -0.027219, 0.051946, 0.364442, -0.777186, 0.481360,0.061273, -0.143294, -0.007607, 0.031695, 0.000542, -0.003382, };



float dwt_filter(float *pdatain, int ndatain)
{
	//plowdataout[],phighdataout []的长度158是4层分解之后，输出的个数
	//对称延拓后每层系数长度(srcLen+filterlen-1)/2，即68+41+28+21
	float *phighdataout = new float[COFE_LENGTH];
	float *plowdataout = new float[COFE_LENGTH];

	decompose(pdatain, ndatain, plowdataout, phighdataout);

	//得到第一层高频系数，求得阈值thr
	int Det1len = 68;
	float *pDet1 = new float[Det1len];
	for (int i = 0; i < Det1len; ++i)
	{
		pDet1[i] = phighdataout[i];
	}
	//float thr = wavelet_getThr(pDet1, Det1len);
	////对高频系数进行阈值去噪
	//wthresh(phighdataout, thr, COFE_LENGTH);

	//对去噪后结果进行小波逆变换
	recompose(pdatain, ndatain, plowdataout, phighdataout);

	delete[]phighdataout;
	phighdataout = NULL;

	delete[]plowdataout;
	plowdataout = NULL;

	delete[]pDet1;
	pDet1 = NULL;
	return 1;
}


//输出H1 H2 H3 H4, L1 L2 L3 L4
float decompose(float *pdatain, int ndatain, float *plowdataout, float *phighdataout)
{
	int nlevel, decLen;
	decLen = ndatain;

	for (nlevel = 0; nlevel<DWT_STAGE; nlevel++)
	{
		decLen = DWT(pdatain, decLen, plowdataout, phighdataout); //调用DWT函数，返回分解系数的长度
		pdatain = plowdataout;      //低频部分再作为输入，继续分解       
		plowdataout += decLen;      //下一层的输出挨着上一层，最后得到一个数组，也可以分开来，先这样
		phighdataout += decLen;
	}
	return 1;
}


/**
 * @brief				多层重构函数
 * @param pdataout		输出数组
 * @param ndataout		输出长度
 * @param pLowDatain	低频分解系数
 * @param pHighDatain	高频分解系数
 * 简单过程：
 * 输入H1 H2 H3 H4, L1 L2 L3 L4，使用了H1 H2 H3 H4,  L4
 * L4+H4->L3, L3+H3->L2, L2+H2-> L1 , 最后L1+H1 -> out
 **/
float recompose(float *pdataout, int ndataout, float *pLowDatain, float *pHighDatain)
{
	int nlevel, cALength;
	float *pLow = new float[68]; //pLow最多到68,第一层个数
	float*pHigh = NULL;
	static int StageBack[] = { 68,41,28,21 };  //对称延拓后每层系数长度(srcLen+filterlen-1)/2

	//pLow = pLowDatain + StageBack[0] ; //先将pLow、pHigh定位到最后一层
	int init_plow = COFE_LENGTH - StageBack[3];
	for (int i = 0; i < StageBack[3]; ++i)
	{
		pLow[i] = pLowDatain[init_plow + i];
	}
	pHigh = pHighDatain + StageBack[0] + StageBack[1] + StageBack[2];
	for (nlevel = DWT_STAGE; nlevel >= 1; nlevel--)
	{
		cALength = StageBack[nlevel - 1];
		Idwt(pLow, pHigh, cALength, pdataout);
		if (nlevel > 1)
		{
			for (int i = 0; i < StageBack[nlevel - 2]; ++i) //pLow更新到上一层
			{
				pLow[i] = pdataout[i];
			}
			//pLow = pdataout; //此处不能直接赋指针，否则，pdataout值改变，则pLow值也立即改变，会影响后面的相关计算
			pHigh -= StageBack[nlevel - 2]; //高频系数直接往上走一层
		}		
	}
	delete[]pLow;
	pLow = NULL;
	return 1;
}


/**
* @brief 			小波变换之分解
* @param pSrcData 	分解的源信号
* @param srcLen 	源信号的长度
* @param cA 		分解后的近似部分序列-低频部分
* @param cD 		分解后的细节部分序列-高频部分
*/

int  DWT(float *pSrcData, int srcLen, float *cA, float *cD)
{
	//禁止出现这种情况，否则数据出错（对称拓延长度为filterLen-1，如果大于了signalLen将越界）

	//if (srcLen < m_dbFilter.filterLen - 1)    //filterLen为低通高通滤波系数长度
	if (srcLen < FILTERLEN - 1)
	{
		cout << "错误信息：滤波器长度大于信号!" << endl;
		exit(1);
	}
	int exLen = (srcLen + FILTERLEN - 1) / 2;//对称拓延后系数的长度
	int k = 0;
	float tmp = 0.0;
	for (int i = 0; i < exLen; i++)
	{

		cA[i] = 0.0;
		cD[i] = 0.0;
		for (int j = 0; j < FILTERLEN; j++)
		{
			k = 2 * i - j + 1;
			//信号边沿对称延拓
			if ((k<0) && (k >= -FILTERLEN + 1))//左边沿拓延
				tmp = pSrcData[-k - 1];
			else if ((k >= 0) && (k <= srcLen - 1))//保持不变
				tmp = pSrcData[k];
			else if ((k>srcLen - 1) && (k <= (srcLen + FILTERLEN - 2)))//右边沿拓延
				tmp = pSrcData[2 * srcLen - k - 1];
			else
				tmp = 0.0;
			cA[i] += sym8lowd_q15[j] * tmp;
			cD[i] += sym8highd_q15[j] * tmp;
		}
	}
	return  exLen;  //返回分解后系数的长度
}


/**
* @brief 			小波变换之重构
* @param cA 		分解后的近似部分序列-低频部分
* @param cD 		分解后的细节部分序列-高频部分
* @param cALength 	上述分解系数的数据长度
* @param recData 	重构后输出的数据
*/
void  Idwt(float *cA, float *cD, int cALength, float *recData)
{
	if ((NULL == cA) || (NULL == cD) || (NULL == recData))
		return;

	int n, k, p, recLen;

	//输入121，四层分解68,41,28,21，偶数的cALength重构时会少一个数
	if(cALength % 2  != 0)
		 recLen = 2 * cALength - FILTERLEN + 2;  
	else
		recLen = 2 * cALength - FILTERLEN + 1;
	//cout << "recLen=" << recLen << endl;

	for (n = 0; n<recLen; n++)
	{
		recData[n] = 0;
		for (k = 0; k<cALength; k++)
		{
			p = n - 2 * k + FILTERLEN - 2;

			//信号重构
			if ((p >= 0) && (p<FILTERLEN))
			{
				recData[n] += sym8lowr_q15[p] * cA[k] + sym8highr_q15[p] * cD[k];
				//cout<<"recData["<<n<<"]="<<recData[n]<<endl;
			}

		}
	}
}

void swap(float *data, int i, int j)
{
	float tmp = data[i];
	data[i] = data[j];
	data[j] = tmp;
}

/**
* @brief 			求取中位数
* @param data 		待求取数组
* @param start 		数组起始位置
* @param last 		数组最后一个元素的下一个位置
* @param nth 		寻找的第nth个数，当nth = (n+1)/2时，即为中位数
* return			由于我们的data是偶数68，所以返回中间两个数平均值
*/
float select_middle(float *data, int start, int last, int nth)
{
	int i = start + 1, j = last - 1;
	float pivot = data[start];
	int k;
	if (last - start < 2)
		return data[start];
	while (i <= j) {
		if (data[i] < pivot) {
			++i;
			continue;
		}
		if (data[j] >= pivot) {
			--j;
			continue;
		}
		swap(data, i, j);
	}
	swap(data, i - 1, start);
	k = i - start;
	if (k == nth)  return (data[i - 1] + data[i]) / 2;
	else if (k>nth)  return select_middle(data, start, i, nth);
	else return select_middle(data, i, last, nth - k);
}

/**
 * @brief			根据细节系数，及信号长度计算阈值
 * @param  pDetCoef 第一层细节系数
 * @param  detlen   其长度
 **/
float wavelet_getThr(float *pDetCoef, int detlen)
{
	float thr = 0.0;
	float sigma = 0.0;
	for (int i = 0; i < detlen; ++i)
	{
		pDetCoef[i] = abs(pDetCoef[i]);
	}

	sigma = select_middle(pDetCoef, 0, detlen, (detlen + 1) / 2) / 0.6745;
	float N = INPUT_LEN; 		//输入信号长度
	thr = sigma*sqrt(2.0*log(N)); //求得阈值

	return thr;

}

 
/**
* @brief				高频系数阈值去噪函数
* @param  pDstCoef		所有高频系数
* @param  high_length   高频系数总长度
**/
void wthresh(float *pDstCoef, float thr, int high_length)
{
	for (int i = 0; i < high_length; ++i)
	{
		if (abs(pDstCoef[i] < thr)) //小于阈值的置零
		{
			pDstCoef[i] = 0.0;
		}
		else 		//大于阈值收缩，软阈值处理
		{
			if (pDstCoef[i] < 0.0)
				pDstCoef[i] = thr - abs(pDstCoef[i]);
			else
				pDstCoef[i] = abs(pDstCoef[i]) - thr;
		}
	}
}

