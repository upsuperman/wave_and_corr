#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INPUT_LEN 121	

float decompose(float *pdatain, int ndatain, float *plowdataout, float *phighdataout);
float recompose(float *pdataout, int ndataout, float *pLowDatain, float *pHighDatain);
int  DWT(float *pSrcData, int srcLen, float *cA, float *cD);
void  Idwt(float *cA, float *cD, int cALength, float *recData);
float dwt_filter(float *pdatain, int ndatain);
float wavelet_getThr(float *pDetCoef, int detlen);
void wthresh(float *pDstCoef, float thr, int high_length);

