

#ifndef _DSPF_SP_CONVOL_CN_H_
#define _DSPF_SP_CONVOL_CN_H_

void DSPF_sp_convol_cn(const float *x, const float *h,
	float *y, const short nx, const short nh, const short ny, const short nmax_ref);

#endif /* _DSPF_SP_CONVOL_CN_H_ */

#define INPUT_LEN 121
#define LORENZ_LEN 60		
#define CORR_OUTPUT_LEN  121
#define DATA_COLS 10001

int max(float a[], int size);
float dot_product(const float *x, const float *y, int length);
float corr(float *corr_input, short input_len);
void fitting(int collection_rows, int datacount_cols, float *input, float *output);
