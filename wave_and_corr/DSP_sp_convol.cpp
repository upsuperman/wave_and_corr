//#pragma CODE_SECTION(DSPF_sp_convol_cn, ".text:ansi");
#include "stdafx.h"

void DSPF_sp_convol_cn(const float *x, const float *h,
    float *y, const short nx, const short nh, const short ny, const short nmax_ref)
{
	short i, j, ref_beg, ref_end, k;

	for (i = 0; i < ny; i++)
	{
		y[i] = 0;
		//前面nmax_ref次，x值都是从0开始计数
		if (i < nmax_ref)
		{
			ref_beg = nmax_ref - i;
			for (j = ref_beg, k = 0; j <= nh - 1; j++, k++)
				y[i] += x[k] * h[j];
		}

		//后面开始,x从0开始，随i递增往前走，到i = nx-nh+nmax_ref时，h的末尾开始收缩,进行后面一半的相关
		else
		{
			ref_end = i - (nx - nmax_ref);
			if (ref_end < 0) ref_end = 0;
			ref_end = nh - 1 - ref_end;

			for (j = 0, k = 0; j <= ref_end; j++, k++)
				y[i] += x[i - nmax_ref + k] * h[j];
		}
	}
} 

/* ======================================================================= */
/*  End of file:  DSPF_sp_convol_cn.c                                      */
/* ----------------------------------------------------------------------- */
/*            Copyright (c) 2011 Texas Instruments, Incorporated.          */
/*                           All Rights Reserved.                          */
/* ======================================================================= */

