// Power_Func_Api.cpp : DLL 응용 프로그램을 위해 내보낸 함수를 정의합니다.
//

#include "stdafx.h"

#include "Power_Func_Dll.h"
#include <stdio.h>
#include <stdarg.h>

#include "fftw\fftw3.h" 

//ver0.33
#include "DspFilters\Butterworth.h"
#include "DspFilters\ChebyshevI.h"
#include "DspFilters\ChebyshevII.h"
#include "DspFilters\Elliptic.h"
#include "DspFilters\Legendre.h"
#include "DspFilters\Filter.h"
#include "DspFilters\Common.h"
#include "DspFilters\Design.h"
#include "DspFilters\SmoothedFilter.h"


#include "Flicker\Filter_mat.h"
#include "Flicker\PQ_function.h"
//#pragma comment(lib, "\PQ_function_x64.lib")


#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _DEBUG
VOID _DebugString(LPCTSTR lpOutFmt, ...)
{
	{
		va_list ap;
		char szTemp[1024];

		va_start(ap, lpOutFmt);
		vsprintf((char *)szTemp, lpOutFmt, ap);
		OutputDebugString(szTemp);
		va_end(ap);
	}
}
#endif


BOOLEAN PWR_FUNC_EXPORT PwScaleRawData(
	void *pReadDataAI,
	void *pSelAiInfo,

	float **pRawBufAI,		//ver0.0.1.0

	int nChCountAi,
	int nStart,
	int nRecvBlockCount,
	int nScanSize,
	int nTotalDataCnt
	)
{
	PBOARD_READ_DATA_ADV_X64 pAI = (PBOARD_READ_DATA_ADV_X64)pReadDataAI;
	PSYS_SEL_AI_INFO pSelInfo = (PSYS_SEL_AI_INFO)pSelAiInfo;

	int i, k;
	int nChNo;
	int nBlockNo = 0;
	int *p32;

	int nCurPos = nStart;


#ifdef _DEBUG
	//for (i = 0; i < nChCountAi; i++)
	//	_DebugString("Neo -  (%ld), (%f)", pAI->BufferPosArr[i], pSelInfo->SelAi[i].SwFactor);
#endif

	/**/
	for (nBlockNo = 0; nBlockNo < nRecvBlockCount; nBlockNo++)
	{
		for (i = 0; i < nScanSize; i++)
		{
			for (k = 0; k < nChCountAi; k++)
			{
				//nChNo = pSelInfo->SelAi[k].ChNo;

				p32 = (int *)pAI->BufferPosArr[k];	// more stable ???
				p32 = (int *)(p32 + pAI->BufferIncArr[k] * (nBlockNo * nScanSize + i));
				
				if ((ULONGLONG)(p32) > pAI->BufferAbsEndPosArr[k])
				{
					p32 = (int *)(pAI->BufferAbsStartPosArr[k] + 
						(((ULONGLONG)p32 - pAI->BufferAbsEndPosArr[k]) - (ULONGLONG)(pAI->BufferIncArr[k] * sizeof(int))));
				}

				//ver0.30.1
				//pRawBufAI[k][nCurPos] = (double)(
				//	((*p32) * pSelInfo->SelAi[k].AdGain) * pSelInfo->SelAi[k].SwFactor + pSelInfo->SelAi[k].SwOffset);
				pRawBufAI[k][nCurPos] = ((*p32) * pSelInfo->SelAi[k].FinalGain) + pSelInfo->SelAi[k].FinalOffset;		//ver0.0.1.0
			} // channel


			//#######################################
			nCurPos++;
			if (nCurPos >= nTotalDataCnt)
				nCurPos = 0;
			//#######################################
		}// scan size
	}
	/**/

	return true;
}


void PWR_FUNC_EXPORT PwZeroCrossing(
	float **pRawBufAI,			//ver0.0.1.0
	//void	*pSelChAi,			// need not	-> //ver0.0.2.1.1 - del
	//int	 nSelChCount,		//ver0.0.2.1.1 - del
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	int		 nRefChIndex,
	int		 nCycle_Ref,
	float	 trigger_level,		//ver0.0.1.0

	int		*pZeroIndex,
	int		*pFirstIndex,
	int		*pLastIndex,
	int		*pLastZ,
	int		*pCalcSize
	)
{
	//PSYS_POWER_AI_CH_SEL_INFO pSelAiInfo = (PSYS_POWER_AI_CH_SEL_INFO)pSelChAi;	//ver0.0.2.1.1 - del
	memset(pZeroIndex, 0, nPwrScanSize);


#ifdef _DEBUG
	//_DebugString("Neo - Power_Func_Dll : ZeroCrossing() - Two mem?(%s), 1st pos(%d), 1st len(%d), Pwr Scan(%d)", 
	//	bIsSeparateMem ? "Yes" : "No", nRPos, nFirstLength, nPwrScanSize);
#endif

	int z = 0, z0 = 0;



	int i, j;
	bool U_trigger = false, L_trigger = false;
	float nZero, nTrigger_H, nTrigger_L;		//ver0.0.1.0

	int nPos1, nPos2;
	int nIndex1, nIndex2;
	for (i = 0; i < (nPwrScanSize - 1); i++)
	{
		nPos1 = nRPos + i;
		nPos2 = nRPos + i + 1;

		if (bIsSeparateMem && (i >= (nFirstLength - 1)))
		{
			//nPos1 = (nRPos + i) % (nDmaBlockCount*nDmaScanSize);
			//nPos2 = (nRPos + i + 1) % (nDmaBlockCount*nDmaScanSize);

			if ((nRPos + i) >= (nDmaBlockCount*nDmaScanSize))
				nPos1 = (nRPos + i) - (nDmaBlockCount*nDmaScanSize);

			if ((nRPos + i + 1) >= (nDmaBlockCount*nDmaScanSize))
				nPos2 = (nRPos + i + 1) - (nDmaBlockCount*nDmaScanSize);
		}



		//nZero = RawData[0][i] * RawData[0][i + 1];	// org
		//nZero = pRawBufAI[0][nPos1] * pRawBufAI[0][nPos2];
		nZero = pRawBufAI[nRefChIndex][nPos1] * pRawBufAI[nRefChIndex][nPos2];


		if (nZero <= 0) 
		{
			////////////////// Level Trigger /////////////////
			for (j = z0; j < i; j++) 
			{
				nIndex1 = nRPos + j;
				nIndex2 = nRPos + j + 1;

				if (bIsSeparateMem && (j >= (nFirstLength - 1)))
				{
					//nIndex1 = (nRPos + j) % (nDmaBlockCount*nDmaScanSize);
					//nIndex2 = (nRPos + j + 1) % (nDmaBlockCount*nDmaScanSize);

					if ((nRPos + j) >= (nDmaBlockCount*nDmaScanSize))
						nIndex1 = (nRPos + j) - (nDmaBlockCount*nDmaScanSize);

					if ((nRPos + j + 1) >= (nDmaBlockCount*nDmaScanSize))
						nIndex2 = (nRPos + j + 1) - (nDmaBlockCount*nDmaScanSize);
				}



				//nTrigger_H = (RawData[0][j] - trigger_level) * (RawData[0][j + 1] - trigger_level);
				//nTrigger_L = (RawData[0][j] + trigger_level) * (RawData[0][j + 1] + trigger_level);
				nTrigger_H = (pRawBufAI[nRefChIndex][nIndex1] - trigger_level) * (pRawBufAI[nRefChIndex][nIndex2] - trigger_level);
				nTrigger_L = (pRawBufAI[nRefChIndex][nIndex1] + trigger_level) * (pRawBufAI[nRefChIndex][nIndex2] + trigger_level);

				if (nTrigger_H <= 0) 
				{
					//if (RawData[0][j] > RawData[0][j + 1]) 
					if (pRawBufAI[nRefChIndex][nIndex1] > pRawBufAI[nRefChIndex][nIndex2])
						U_trigger = true;
					else 
						U_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}

				if (nTrigger_L <= 0) 
				{
					//if (RawData[0][j] < RawData[0][j + 1]) 
					if (pRawBufAI[nRefChIndex][nIndex1] < pRawBufAI[nRefChIndex][nIndex2])
						L_trigger = true;
					else 
						L_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}
			}
			////////////////////////////////////////////////////////				


			z0 = i;

			if (U_trigger != L_trigger) {
				//pZeroIndex[z] = i;
				pZeroIndex[z] = nRPos + i;
				z++;
				i += (nPwrScanSize / nCycle_Ref / 2 / 100);

				U_trigger = false;
				L_trigger = false;
			}

		} // end of if (nZero <= 0)

		if (z >= nPwrScanSize) {
			//i = nScanSize;
			break;
		}

	}

	
	*pFirstIndex = pZeroIndex[0];

	*pLastZ = 0;
	if (z > 0)
		*pLastZ = z - 1;

	*pLastIndex = pZeroIndex[*pLastZ];
	*pCalcSize = *pLastIndex - *pFirstIndex;

	if (*pCalcSize <= 0)
	{
		*pCalcSize = nPwrScanSize;
		//*pFirstIndex = 0;
		//*pLastIndex = nPwrScanSize - 1;
		*pFirstIndex = nRPos;
		*pLastIndex = nRPos + nPwrScanSize - 1;
	}
}


//#####################################################################
float *Raw_que[2];		//ver0.0.1.0
float *Raw_LPF;			//ver0.0.1.0
void PWR_FUNC_EXPORT PwAllocFilterMemory(int nCount)
{
	if (Raw_que[0] != 0)
		delete Raw_que[0];

	if (Raw_que[1] != 0)
		delete Raw_que[1];

	if (Raw_LPF != 0)
		delete Raw_LPF;

	Raw_que[0] = new float[nCount * 2];		//ver0.0.1.0
	Raw_que[1] = new float[nCount * 2];		//ver0.0.1.0

	Raw_LPF = new float[nCount];			//ver0.0.1.0
}

void PWR_FUNC_EXPORT PwFreeFilterMemory()
{
	delete Raw_que[0];
	delete Raw_que[1];

	delete Raw_LPF;
}


void PWR_FUNC_EXPORT PwZeroCrossingWithFilter(
	float **pRawBufAI,			//ver0.0.1.0
	//void	*pSelChAi,			// need not	-> //ver0.0.2.1.1 - del
	//int		 nSelChCount,	//ver0.0.2.1.1 - del
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,
	float	 fSampleRate,		//ver0.0.1.0

	int		 nLowOrder,
	float	 fLowCutOff,		//ver0.0.1.0
	int		 nHighOrder,
	float	 fHighCutOff,		//ver0.0.1.0


	int		 nRefChIndex,
	int		 nCycle_Ref,
	float	 trigger_level,		//ver0.0.1.0

	int		*pZeroIndex,
	int		*pFirstIndex,
	int		*pLastIndex,
	int		*pLastZ,
	int		*pCalcSize
	)
{
	//PSYS_POWER_AI_CH_SEL_INFO pSelAiInfo = (PSYS_POWER_AI_CH_SEL_INFO)pSelChAi;	//ver0.0.2.1.1 - del
	memset(pZeroIndex, 0, nPwrScanSize);


#ifdef _DEBUG
	//_DebugString("Neo - Power_Func_Dll : ZeroCrossing() - Two mem?(%s), 1st pos(%d), 1st len(%d), Pwr Scan(%d)", 
	//	bIsSeparateMem ? "Yes" : "No", nRPos, nFirstLength, nPwrScanSize);
#endif
	int i, j;
	int nPos1, nPos2, nPosPrevStart, nPosPrev;



	nPosPrevStart = nRPos - nPwrScanSize;
	//_DebugString("Neo - Power_Func_Dll : ZeroCrossingWithFilter() - (nPrevPos) [%d], (nRPos) [%d]", nPosPrevStart, nRPos);

	if (nPosPrevStart < 0)
	{
		nPosPrevStart = nPosPrevStart + (nDmaBlockCount*nDmaScanSize);
		//_DebugString("Neo - Power_Func_Dll : ZeroCrossingWithFilter() - nPrevPos < 0 : (nPrevPos) [%d], (nRPos) [%d]", nPosPrevStart, nRPos);
	}

	memcpy(&Raw_que[0][0], &pRawBufAI[nRefChIndex][nPosPrevStart], sizeof(float)* nPwrScanSize);		//ver0.0.1.0
	memcpy(&Raw_que[0][nPwrScanSize], &pRawBufAI[nRefChIndex][nRPos], sizeof(float)* nPwrScanSize);		//ver0.0.1.0



	//Dsp::SimpleFilter <Dsp::Butterworth::LowPass <5>, 1> lp;
	//lp.setup(5, dSampleRate, 65);
	//lp.process((nPwrScanSize * 2), Raw_que);

	//Dsp::SimpleFilter <Dsp::Butterworth::HighPass <5>, 1> hp;
	//hp.setup(5, dSampleRate, 45);
	//hp.process((nPwrScanSize * 2), Raw_que);

	Dsp::SimpleFilter <Dsp::Butterworth::LowPass <5>, 1> lp;
	lp.setup(nLowOrder, fSampleRate, fLowCutOff);				//ver0.0.1.0
	lp.process((nPwrScanSize * 2), Raw_que);

	Dsp::SimpleFilter <Dsp::Butterworth::HighPass <5>, 1> hp;
	hp.setup(nHighOrder, fSampleRate, fHighCutOff);				//ver0.0.1.0
	hp.process((nPwrScanSize * 2), Raw_que);



	memcpy(Raw_LPF, &Raw_que[0][nPwrScanSize], sizeof(float)* nPwrScanSize);		//ver0.0.1.0



	int z = 0, z0 = 0;



	bool U_trigger = false, L_trigger = false;
	float nZero, nTrigger_H, nTrigger_L;			//ver0.0.1.0

	for (i = 0; i < (nPwrScanSize - 1); i++)
	{
		nZero = Raw_LPF[i] * Raw_LPF[i + 1];

		if (nZero <= 0)
		{
			////////////////// Level Trigger /////////////////
			for (j = z0; j < i; j++)
			{
				nTrigger_H = (Raw_LPF[j] - trigger_level) * (Raw_LPF[j + 1] - trigger_level);
				nTrigger_L = (Raw_LPF[j] + trigger_level) * (Raw_LPF[j + 1] + trigger_level);


				if (nTrigger_H <= 0)
				{
					if (Raw_LPF[j] > Raw_LPF[j + 1])
						U_trigger = true;
					else
						U_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}

				if (nTrigger_L <= 0)
				{
					if (Raw_LPF[j] < Raw_LPF[j + 1])
						L_trigger = true;
					else
						L_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}
			}
			////////////////////////////////////////////////////////				


			z0 = i;

			if (U_trigger != L_trigger) {
				//pZeroIndex[z] = i;
				pZeroIndex[z] = nRPos + i;
				z++;
				i += (nPwrScanSize / nCycle_Ref / 2 / 100);

				U_trigger = false;
				L_trigger = false;
			}

		} // end of if (nZero <= 0)

		if (z >= nPwrScanSize) {
			//i = nScanSize;
			break;
		}

	}


	*pFirstIndex = pZeroIndex[0];

	*pLastZ = 0;
	if (z > 0)
		*pLastZ = z - 1;

	*pLastIndex = pZeroIndex[*pLastZ];
	*pCalcSize = *pLastIndex - *pFirstIndex;

	if (*pCalcSize <= 0)
	{
		*pCalcSize = nPwrScanSize;
		//*pFirstIndex = 0;
		//*pLastIndex = nPwrScanSize - 1;
		*pFirstIndex = nRPos;
		*pLastIndex = nRPos + nPwrScanSize - 1;
	}
#ifdef _DEBUG
	//if (nRefChIndex == 3)
	//	for (i = 0; i <= *pLastZ; i++)
	//		_DebugString("Neo - Power_Func_Dll : ZeroCrossing() ------------------> nZeroIndex[ %d ] : ( %d )", i, pZeroIndex[i]);
#endif
}

float *BPF_que[2];		//ver0.0.1.0
float *Prev_BPF;		//ver0.0.1.0

//ver0.0.1.7 - this is never used !!!!
void PWR_FUNC_EXPORT PwInitBPFilterMem(int nCount)
{
	if (BPF_que[0] != 0)
		delete BPF_que[0];

	if (BPF_que[1] != 0)
		delete BPF_que[1];

	BPF_que[0] = new float[nCount * 2];		//ver0.0.1.0
	BPF_que[1] = new float[nCount * 2];		//ver0.0.1.0



	if (Prev_BPF != 0)
		delete Prev_BPF;

	Prev_BPF = new float[nCount];					//ver0.0.1.0
	memset(Prev_BPF, 0, sizeof(float)* nCount);	//ver0.0.1.0
}

//ver0.0.1.7 - this is never used !!!!
void PWR_FUNC_EXPORT PwBPFilter(
	float	*fInput,			//ver0.0.1.0
	float	 fSampleRate,		//ver0.0.1.0
	int		 nScanSize,

	int		 nLowOrder,
	int		 nLowCutOff,

	int		 nHighOrder,
	int		 nHighCutOff,

	float	*fOutData			//ver0.0.1.0
	)
{
	memcpy(&BPF_que[0][0], Prev_BPF, sizeof(double) * nScanSize);			//ver0.0.1.0
	memcpy(&BPF_que[0][nScanSize], fInput, sizeof(double)* nScanSize);			//ver0.0.1.0


	//Dsp::SimpleFilter <Dsp::Butterworth::LowPass <5>, 1> lp;
	//lp.setup(5, dSampleRate, 65);
	//lp.process((nScanSize * 2), BPF_que);

	//Dsp::SimpleFilter <Dsp::Butterworth::HighPass <5>, 1> hp;
	//hp.setup(5, dSampleRate, 45);
	//hp.process((nScanSize * 2), BPF_que);

	Dsp::SimpleFilter <Dsp::Butterworth::LowPass <10>, 1> lp;
	lp.setup(nLowOrder, fSampleRate, nLowCutOff);			//ver0.0.1.0
	lp.process((nScanSize * 2), BPF_que);

	Dsp::SimpleFilter <Dsp::Butterworth::HighPass <10>, 1> hp;
	hp.setup(nHighOrder, fSampleRate, nHighCutOff);			//ver0.0.1.0
	hp.process((nScanSize * 2), BPF_que);


	memcpy(fOutData, &BPF_que[0][nScanSize], sizeof(float)* nScanSize);			//ver0.0.1.0
	memcpy(Prev_BPF, &BPF_que[0][nScanSize], sizeof(float)* nScanSize);			//ver0.0.1.0
}


//#####################################################################


void PWR_FUNC_EXPORT PwButterworth(
	float	*fInput, 			//ver0.0.1.0
	float	 fSamplerateHz,		//ver0.0.1.0
	float	 fCutOffFreqHz,		//ver0.0.1.0
	int		 iOrder, 
	float	*fOutData			//ver0.0.1.0
	)
{
	float m_fFilter_Coef[6][10];			//ver0.0.1.0
	float m_fMemory_Coef[3][10];			//ver0.0.1.0

	//int iScanSize = dSamplerateHz;//sizeof(fInput) / sizeof(double);
	int iScanSize = fSamplerateHz;			//ver0.0.1.0


	if (fCutOffFreqHz <= 0.0)					//ver0.0.1.0
	{
		memcpy(fInput, fOutData, iScanSize);	//ver0.0.1.0
		//return dOutData;
		return;
	}

	// 1. Calculate the filter coefficients
	//
	int iNs2, iModn;
	float Arg, Rep, Omega, OmegaSq, temp, W0, W1, m;	//ver0.0.1.0
	float Zero, ONE, TWO, HALF, Pi;						//ver0.0.1.0

	Zero = 0;
	ONE = 1;
	TWO = 2;
	HALF = 0.5;
	Pi = 3.1415926535;

	//Arg = Pi * Ts * Fc;
	float fTs = 1.0 / fSamplerateHz;	//ver0.0.1.0
	Arg = Pi * fTs * fCutOffFreqHz;		//ver0.0.1.0
	if (abs(Arg) > 2.0 * Pi)
	{
		m = (int)(Arg / 2.0 / Pi);
		Arg = Arg - (m * 2.0 * Pi);
	}

	Omega = tan(Arg);
	OmegaSq = Omega * Omega;
	iModn = (iOrder % 2);
	if (iModn == 0)
		temp = HALF;
	else
		temp = Zero;

	iNs2 = iOrder / 2;
	int iNSections = iNs2 + iModn;
	float fTg = Zero;				//ver0.0.1.0

	if (iOrder > 1)
	for (int i = 1; i < iNs2 + 1; i++)
	{
		Rep = Omega * cos(Pi * (i - temp) / iOrder);
		fTg = fTg + fTs * Rep / OmegaSq;					//ver0.0.1.0
		W0 = TWO * Rep;
		W1 = ONE + W0 + OmegaSq;
		m_fFilter_Coef[1][i] = -TWO * (OmegaSq - ONE) / W1;	//ver0.0.1.0
		m_fFilter_Coef[2][i] = -(ONE - W0 + OmegaSq) / W1;	//ver0.0.1.0
		m_fFilter_Coef[3][i] = TWO;							//ver0.0.1.0
		m_fFilter_Coef[4][i] = ONE;							//ver0.0.1.0
		m_fFilter_Coef[5][i] = OmegaSq / W1;				//ver0.0.1.0
	}
	if (temp == Zero)
	{
		m_fFilter_Coef[1][iNSections] = (ONE - Omega) / (ONE + Omega);	//ver0.0.1.0
		m_fFilter_Coef[2][iNSections] = Zero;					//ver0.0.1.0
		m_fFilter_Coef[3][iNSections] = ONE;					//ver0.0.1.0
		m_fFilter_Coef[4][iNSections] = Zero;					//ver0.0.1.0
		m_fFilter_Coef[5][iNSections] = Omega / (ONE + Omega);	//ver0.0.1.0
		fTg = fTg + fTs / (TWO * Omega);						//ver0.0.1.0
	}

	// 2. Initialize filter memory
	//
	float fSum;			//ver0.0.1.0
	float fDC = 0;		//ver0.0.1.0
	for (int j = 1; j < iNSections + 1; j++)
	{
		m_fMemory_Coef[2][j] = fDC / (1 - m_fFilter_Coef[1][j] - m_fFilter_Coef[2][j]);	//ver0.0.1.0
		m_fMemory_Coef[1][j] = m_fMemory_Coef[2][j];									//ver0.0.1.0
		fSum = 0;		//ver0.0.1.0

		for (int i = 1; i < 5; i++)
			fSum = fSum + m_fFilter_Coef[i][j];		//ver0.0.1.0

		fDC = m_fFilter_Coef[5][j] * (fDC + m_fMemory_Coef[2][j] * fSum);				//ver0.0.1.0
	}

	// 3. Recursively call Butterworth filter
	//
	float fData = 0, err = 0;	//ver0.0.1.0
	float sInput;				//ver0.0.1.0

	for (int j = 0; j < iScanSize; j++)
	{
		sInput = fInput[j];		//ver0.0.1.0

		for (int i = 1; i < iNSections + 1; i++)
		{
			err = sInput + m_fFilter_Coef[1][i] * m_fMemory_Coef[1][i] + m_fFilter_Coef[2][i] * m_fMemory_Coef[2][i];	//ver0.0.1.0
			fData = m_fFilter_Coef[5][i] * (err + m_fFilter_Coef[3][i] * m_fMemory_Coef[1][i] + m_fFilter_Coef[4][i] * m_fMemory_Coef[2][i]);	//ver0.0.1.0
			m_fMemory_Coef[2][i] = m_fMemory_Coef[1][i];	//ver0.0.1.0
			m_fMemory_Coef[1][i] = err;						//ver0.0.1.0
			sInput = fData;									//ver0.0.1.0
		}

		fOutData[j] = sInput;								//ver0.0.1.0
	}

	//return dOutData;
}

//Obsolete
void PWR_FUNC_EXPORT PwZeroCrossingForLPF(
	float   *pRaw_LPF,			//ver0.0.1.0

	int		 nRPos,
	int		 nPwrScanSize,

	int		 nCycle_Ref,
	float	 trigger_level,		//ver0.0.1.0

	int		*pZeroIndex,
	int		*pFirstIndex,
	int		*pLastIndex,
	int		*pLastZ,
	int		*pCalcSize
	)
{
	memset(pZeroIndex, 0, nPwrScanSize);


#ifdef _DEBUG
	//_DebugString("Neo - Power_Func_Dll : ZeroCrossing() - Two mem?(%s), 1st pos(%d), 1st len(%d), Pwr Scan(%d)", 
	//	bIsSeparateMem ? "Yes" : "No", nRPos, nFirstLength, nPwrScanSize);
#endif

	int z = 0, z0 = 0;



	int i, j;
	bool U_trigger = false, L_trigger = false;
	float nZero, nTrigger_H, nTrigger_L;		//ver0.0.1.0

	for (i = 0; i < (nPwrScanSize - 1); i++)
	{
		nZero = pRaw_LPF[i] * pRaw_LPF[i + 1];

		if (nZero <= 0)
		{
			////////////////// Level Trigger /////////////////
			for (j = z0; j < i; j++)
			{
				nTrigger_H = (pRaw_LPF[j] - trigger_level) * (pRaw_LPF[j + 1] - trigger_level);
				nTrigger_L = (pRaw_LPF[j] + trigger_level) * (pRaw_LPF[j + 1] + trigger_level);


				if (nTrigger_H <= 0)
				{
					if (pRaw_LPF[j] > pRaw_LPF[j + 1])
						U_trigger = true;
					else
						U_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}

				if (nTrigger_L <= 0)
				{
					if (pRaw_LPF[j] < pRaw_LPF[j + 1])
						L_trigger = true;
					else
						L_trigger = false;

					j += (nPwrScanSize / nCycle_Ref / 2 / 100);
				}
			}
			////////////////////////////////////////////////////////				


			z0 = i;

			if (U_trigger != L_trigger) {
				//pZeroIndex[z] = i;
				pZeroIndex[z] = nRPos + i;
				z++;
				i += (nPwrScanSize / nCycle_Ref / 2 / 100);

				U_trigger = false;
				L_trigger = false;
			}

		} // end of if (nZero <= 0)

		if (z >= nPwrScanSize) {
			//i = nScanSize;
			break;
		}

	}


	*pFirstIndex = pZeroIndex[0];

	*pLastZ = 0;
	if (z > 0)
		*pLastZ = z - 1;

	*pLastIndex = pZeroIndex[*pLastZ];
	*pCalcSize = *pLastIndex - *pFirstIndex;

	if (*pCalcSize <= 0)
	{
		*pCalcSize = nPwrScanSize;
		//*pFirstIndex = 0;
		//*pLastIndex = nPwrScanSize - 1;
		*pFirstIndex = nRPos;
		*pLastIndex = nRPos + nPwrScanSize - 1;
	}
}

void PWR_FUNC_EXPORT PwCalcWindowSize(
	int		nPwrScanSize,

	int		*pZeroIndex,

	int		 nCycle_Ref,
	int		 nFirstIndex,
	int		 nLastIndex,
	int		 nLastZ,
	int		 nCalcSize,

	int		*pIntIndex,
	int		*pIntCycleSamples,
	int		*pNumIntCycle
	)
{
	if ((nLastZ % 2) == 0) {
		*pIntIndex = nLastIndex;
		*pIntCycleSamples = nCalcSize;
		*pNumIntCycle = nLastZ / 2;
	}
	else if ((nLastZ % 2) != 0) {
		*pIntIndex = pZeroIndex[nLastZ - 1];
		*pIntCycleSamples = *pIntIndex - nFirstIndex;
		*pNumIntCycle = (nLastZ - 1) / 2;
	}

	if (*pNumIntCycle == 0) {
		*pNumIntCycle = nCycle_Ref;
		*pIntCycleSamples = nPwrScanSize;
	}

	if (*pIntCycleSamples == 0) {
		*pNumIntCycle = nCycle_Ref;
		*pIntCycleSamples = nPwrScanSize;
	}
}


void PWR_FUNC_EXPORT PwSetInFFTMem(
	int		 nChIndex,

	float  **pRawBufAI,			//ver0.0.1.0
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	int		 nFirstIndex,
	int		 nIntIndex,

	float   *pInFFTMem
	)
{
	int k;
	int nPos;

	//LARGE_INTEGER		g_nFreq1, g_nStartTime1, g_nEndTime1;
	//double				g_dPassTimeU;

	//QueryPerformanceFrequency(&g_nFreq1);
	//QueryPerformanceCounter(&g_nStartTime1);

#ifdef _DEBUG
	//_DebugString("Neo -  Dll :DMA scan size(%d), Blocks(%d)", nDmaScanSize, nDmaBlockCount);
	//_DebugString("Neo -  Dll :Separate?(%s), Pos(%d), FirstLen(%d), SCAN(%d), FirstIndex(%d), IntIndex(%d)",
	//	bIsSeparateMem ? "YES" : "NO", nRPos, nFirstLength, nPwrScanSize, nFirstIndex, nIntIndex);

	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll :(%d) (%d) Buf [%d][%d] = %f", nRPos, nFirstIndex, nChIndex, k, pRawBufAI[nChIndex][k]);
	//_DebugString("Neo -  Dll :--------------");
#endif


	for (k = nFirstIndex; k < (nRPos + nPwrScanSize); k++)
	{
		nPos = k;

		if (bIsSeparateMem && (k >= (nDmaBlockCount*nDmaScanSize)))
			nPos = (nRPos + k) - (nDmaBlockCount*nDmaScanSize);

		if (k < nIntIndex) {
		//if (nPos < nIntIndex) {
			pInFFTMem[(k - nFirstIndex) * 2] = pRawBufAI[nChIndex][nPos];		//ver0.0.1.0
			pInFFTMem[(k - nFirstIndex) * 2 + 1] = 0;
		}
		else {
			pInFFTMem[(k - nFirstIndex) * 2] = 0;
			pInFFTMem[(k - nFirstIndex) * 2 + 1] = 0;
		}
	}
#ifdef _DEBUG
	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll : After In  [%d][%d] = %f", nChIndex, k, pInFFTMem[k]);
	//_DebugString("Neo -  Dll :--------------");
#endif


	//QueryPerformanceCounter(&g_nEndTime1);
	//g_dPassTimeU = (double)(g_nEndTime1.QuadPart - g_nStartTime1.QuadPart) / ((double)(g_nFreq1.QuadPart));
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : ChNo[%d], Set In FFT(%10.7f ms)", nChIndex, g_dPassTimeU * 1000);
#endif

}


void PWR_FUNC_EXPORT PwSetInFFTMem4Down(
	int		 nChIndex,

	float  **pRawBufAI,			//ver0.0.1.0
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	int		 nFirstIndex,
	int		 nIntIndex,


	int		 nDownRatio,
	int		 nIntCycleSamples,
	float   *pDownRaw,
	float   *pInFFTMem
	)
{
	int k;
	int nPos;

	//LARGE_INTEGER		g_nFreq1, g_nStartTime1, g_nEndTime1;
	//double				g_dPassTimeU;

	//QueryPerformanceFrequency(&g_nFreq1);
	//QueryPerformanceCounter(&g_nStartTime1);

#ifdef _DEBUG
	//_DebugString("Neo -  Dll :DMA scan size(%d), Blocks(%d)", nDmaScanSize, nDmaBlockCount);
	//_DebugString("Neo -  Dll :Separate?(%s), Pos(%d), FirstLen(%d), SCAN(%d), FirstIndex(%d), IntIndex(%d)",
	//	bIsSeparateMem ? "YES" : "NO", nRPos, nFirstLength, nPwrScanSize, nFirstIndex, nIntIndex);

	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll :(%d) (%d) Buf [%d][%d] = %f", nRPos, nFirstIndex, nChIndex, k, pRawBufAI[nChIndex][k]);
	//_DebugString("Neo -  Dll :--------------");
#endif

	// down samples
	int DownSize = nPwrScanSize / nDownRatio;

	//ver0.0.3.2 - mod
	/**
	for (k = 0; k < DownSize; k++) 
	{
		nPos = nRPos + (k * nDownRatio);

		//ver0.0.2.1.1 - del - no need to consider separate memory
		//if (bIsSeparateMem && (nPos >= (nDmaBlockCount*nDmaScanSize)))
		//	nPos = nPos - (nDmaBlockCount*nDmaScanSize);
		//ver0.0.2.1.1 - del - no need to consider separate memory - end

		pDownRaw[k] = pRawBufAI[nChIndex][nPos];
	}


	for (k = 0; k < nPwrScanSize; k++)
	{
		nPos = nRPos + k;

		//ver0.0.2.1.1 - del - no need to consider separate memory
		//if (bIsSeparateMem && (nPos >= (nDmaBlockCount*nDmaScanSize)))
		//	nPos = nPos - (nDmaBlockCount*nDmaScanSize);
		//ver0.0.2.1.1 - del - no need to consider separate memory - end

		if (k < nIntCycleSamples) {
			pInFFTMem[k * 2] = pDownRaw[k];
			pInFFTMem[k * 2 + 1] = 0;
		}
		else {
			pInFFTMem[k * 2] = 0;
			pInFFTMem[k * 2 + 1] = 0;
		}
	}
	/**/

	for (k = 0; k < DownSize; k++)
	{
		//nPos = nRPos + (k * nDownRatio);
		nPos = nFirstIndex + (k * nDownRatio);			//ver0.0.3.5 - ????????????? think more ....

		if (k < nIntCycleSamples)
		{
			pInFFTMem[k * 2] = pRawBufAI[nChIndex][nPos];
			pInFFTMem[k * 2 + 1] = 0;
		}
		else
		{
			pInFFTMem[k * 2] = 0;
			pInFFTMem[k * 2 + 1] = 0;
		}
	}
	//ver0.0.3.2 - mod - end

#ifdef _DEBUG
	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll : After In  [%d][%d] = %f", nChIndex, k, pInFFTMem[k]);
	//_DebugString("Neo -  Dll :--------------");
#endif


	//QueryPerformanceCounter(&g_nEndTime1);
	//g_dPassTimeU = (double)(g_nEndTime1.QuadPart - g_nStartTime1.QuadPart) / ((double)(g_nFreq1.QuadPart));
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : ChNo[%d], Set In FFT(%10.7f ms)", nChIndex, g_dPassTimeU * 1000);
#endif

}


/**/
void PWR_FUNC_EXPORT PwCalcHarmonics(
	int		 nChIndex,

	float  **pRawBufAI,			//ver0.0.1.0
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	float   *pfft_out,
	int		 nIntCycleSamples,
	int		 nNumIntCycle,
	int		 nFirstIndex,
	int		 nLastIndex,
	int		 nIntIndex,
	int		 nCalcSize,
	float	 fft_p_ref,			//ver0.0.1.0

	float	*pfft_a,			//ver0.0.1.0
	float	*pfft_p,			//ver0.0.1.0
	float	*pfft_IH,			//ver0.0.1.0
	float	*pTHD,				//ver0.0.1.0	// ref

	float	*pFund_rms,			//ver0.0.1.0	// ref
	float	*pFund_angle,		//ver0.0.1.0	// ref
	float	*pDC_rms,			//ver0.0.1.0	// ref
	float	*prms				//ver0.0.1.0	// ref
	)
{
	int k;
	int nPos;

	//LARGE_INTEGER		g_nFreq1, g_nStartTime1, g_nEndTime1;
	//double				g_dPassTimeU;

	//QueryPerformanceFrequency(&g_nFreq1);
	//QueryPerformanceCounter(&g_nStartTime1);

#ifdef _DEBUG
	//_DebugString("Neo -  Dll :DMA scan size(%d), Blocks(%d)", nDmaScanSize, nDmaBlockCount);
	//_DebugString("Neo -  Dll :Separate?(%s), Pos(%d), FirstLen(%d), SCAN(%d), FirstIndex(%d), IntIndex(%d)",
	//	bIsSeparateMem ? "YES" : "NO", nRPos, nFirstLength, nPwrScanSize, nFirstIndex, nIntIndex);

	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll :(%d) (%d) Buf [%d][%d] = %f", nRPos, nFirstIndex, nChIndex, k, pRawBufAI[nChIndex][k]);
	//_DebugString("Neo -  Dll :--------------");
#endif
	//int		_nHarmRmsCnt = 0;
	int		_nHarmPhiCnt = 0;
	int		_nInterHPhiCnt = 0;

	int		ik = 0, r = 0;
	bool	bIsPI = false;
	float	nSum = 0;			//ver0.0.1.0



	//########################
	*pTHD = 0;
	//########################

	for (k = 0; k < (nIntCycleSamples / 2); k++)
	{
		/////////////// Magnitude //////////////
		pfft_a[k] = sqrt((pfft_out[k * 2] * pfft_out[k * 2]) + (pfft_out[k * 2 + 1] * pfft_out[k * 2 + 1])) / nIntCycleSamples / sqrt((float)2) * 2;	//ver0.0.1.0
		////////////////////////////////////////

		if ((k % nNumIntCycle) == 0)   // 정수배 일때만 계산
		{
			/////////////// Phase Angle //////////////////
			pfft_p[k] = (atan2(pfft_out[k * 2 + 1], pfft_out[k * 2]) * 180 / M_PI) - fft_p_ref;

			bIsPI = true;

			while (bIsPI)
			{
				if (pfft_p[k] > 180) {
					pfft_p[k] -= 360;
				}
				else if (pfft_p[k] < -180) {
					pfft_p[k] += 360;
				}
				else {
					bIsPI = false;
				}

			}

			// here !!!
			//	harmonic rms
			//	harmonic phi
			//	harmonics must be up to 50th. so need to check if count is < 50 or not
			//	increment harmonic count

			//_nHarmPhiCnt++;
			//////////////////////////////////////////////


			//////////////////// THD /////////////////////
			if (k > nNumIntCycle && k < nNumIntCycle * 51) {
				*pTHD = *pTHD + (pfft_a[k] * pfft_a[k]);
			}
			//////////////////////////////////////////////
		}

		//_nHarmRmsCnt++;

		// init for inter harmonic calculation
		//########################
		pfft_IH[k] = 0;
		//########################
	
	}// end of for (k = 0; k < (_nIntCycleSamples / 2); k++) 


	/////////////  Inter Harmonics ///////////////
	_nInterHPhiCnt = 0;
	for (ik = 0; ik < (nIntCycleSamples / 2); ik += nNumIntCycle)
	{
		for (r = 1; r < (nNumIntCycle); r++)
		{
			pfft_IH[ik] = pfft_IH[ik] + (pfft_a[ik + r] * pfft_a[ik + r]);
		}
		pfft_IH[ik] = sqrt(pfft_IH[ik]);

		// here !!!
		//	harmonic rms
		//	harmonic phi
		//	harmonics must be up to 50th. so need to check if count is < 50 or not
		//	increment harmonic count

		//_nInterHPhiCnt++;
	}
	//////////////////////////////////////////////



	/////////////////// THD //////////////////////
	*pTHD = sqrt(*pTHD) / pfft_a[nNumIntCycle] * 100;


	////// Fundamental RMS, PHI and DC RMS ///////
	// Fundamental RMS, PHI and DC RMS
	*pFund_rms = pfft_a[nNumIntCycle];
	*pFund_angle = pfft_p[nNumIntCycle];
	*pDC_rms = pfft_a[0] / sqrt((float)2) / 2;			//ver0.0.1.0


	/////////////////// RMS //////////////////////
	nSum = 0;
	for (k = nFirstIndex; k < nLastIndex; k++)
		nSum = nSum + (pRawBufAI[nChIndex][k] * pRawBufAI[nChIndex][k]);

	*prms = sqrt(nSum / nCalcSize);


#ifdef _DEBUG
	//if (nChIndex == 0)
	//	_DebugString("Neo -  Dll : Rms Count (%d), Phi Count (%d), IH Count (%d)", _nHarmRmsCnt, _nHarmPhiCnt, _nInterHPhiCnt);
	//_DebugString("Neo -  Dll :--------------");
#endif


	//QueryPerformanceCounter(&g_nEndTime1);
	//g_dPassTimeU = (double)(g_nEndTime1.QuadPart - g_nStartTime1.QuadPart) / ((double)(g_nFreq1.QuadPart));
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : ChNo[%d], Set In FFT(%10.7f ms)", nChIndex, g_dPassTimeU * 1000);
#endif
}
/**/

/**/
void PWR_FUNC_EXPORT PwCalcHarmonics2(
	int		 nChIndex,	//*

	float  **pRawBufAI,			//ver0.0.1.0
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	float   *pfft_out,	//*
	int		 nIntCycleSamples,	//*
	int		 nNumIntCycle,	//*
	int		 nFirstIndex,	//*
	int		 nLastIndex,	//*
	int		 nIntIndex,
	int		 nCalcSize,	//*
	float	 fft_p_ref,	//*			//ver0.0.1.0

	float	*pfft_a,			//ver0.0.1.0
	float	*pfft_p,			//ver0.0.1.0
	float	*pfft_IH,			//ver0.0.1.0
	float	*pTHD,				//ver0.0.1.0		// ref

	float	*pFund_rms,			//ver0.0.1.0		// ref
	float	*pFund_angle,		//ver0.0.1.0		// ref
	float	*pDC_rms,			//ver0.0.1.0		// ref
	float	*prms,				//ver0.0.1.0		// ref
	float	*pTotal_rms,		//ver0.0.1.0		// ref	//ver0.56 - for line voltage

	float	*pHarmonic_Rms,		//ver0.0.1.0
	float	*pHarmonic_Phi,		//ver0.0.1.0
	float	*pInterharmonic		//ver0.0.1.0
	)
{
	int k;
	int nPos;

	//LARGE_INTEGER		g_nFreq1, g_nStartTime1, g_nEndTime1;
	//double				g_dPassTimeU;

	//QueryPerformanceFrequency(&g_nFreq1);
	//QueryPerformanceCounter(&g_nStartTime1);

#ifdef _DEBUG
	//_DebugString("Neo -  Dll :DMA scan size(%d), Blocks(%d)", nDmaScanSize, nDmaBlockCount);
	//_DebugString("Neo -  Dll :Separate?(%s), Pos(%d), FirstLen(%d), SCAN(%d), FirstIndex(%d), IntIndex(%d)",
	//	bIsSeparateMem ? "YES" : "NO", nRPos, nFirstLength, nPwrScanSize, nFirstIndex, nIntIndex);

	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll :(%d) (%d) Buf [%d][%d] = %f", nRPos, nFirstIndex, nChIndex, k, pRawBufAI[nChIndex][k]);
	//_DebugString("Neo -  Dll :--------------");
#endif
	//int		_nHarmRmsCnt = 0;
	int		_nHarmPhiCnt = 0;
	int		_nInterHPhiCnt = 0;

	int		ik = 0, r = 0;
	bool	bIsPI = false;
	float	nSum = 0;			//ver0.0.1.0



	//########################
	*pTHD = 0;
	//########################

	*pTotal_rms = 0;	//ver0.56


	for (k = 0; k < (nIntCycleSamples / 2); k++)
	{
		/////////////// Magnitude //////////////
		pfft_a[k] = sqrt((pfft_out[k * 2] * pfft_out[k * 2]) + (pfft_out[k * 2 + 1] * pfft_out[k * 2 + 1])) / nIntCycleSamples / sqrt((float)2) * 2;	//ver0.0.1.0

		//ver0.56
		*pTotal_rms = *pTotal_rms + (pfft_a[k] * pfft_a[k]);
		////////////////////////////////////////




		if ((k % nNumIntCycle) == 0)   // 정수배 일때만 계산
		{
			/////////////// Phase Angle //////////////////
			pfft_p[k] = (atan2(pfft_out[k * 2 + 1], pfft_out[k * 2]) * 180 / M_PI) - fft_p_ref;

			bIsPI = true;

			while (bIsPI)
			{
				if (pfft_p[k] > 180) {
					pfft_p[k] -= 360;
				}
				else if (pfft_p[k] < -180) {
					pfft_p[k] += 360;
				}
				else {
					bIsPI = false;
				}

			}

			// here !!!
			//	harmonic rms
			//	harmonic phi
			//	harmonics must be up to 50th. so need to check if count is < 50 or not
			//	increment harmonic count
			if (_nHarmPhiCnt < 51)
			{
				pHarmonic_Rms[_nHarmPhiCnt] = pfft_a[k];
				pHarmonic_Phi[_nHarmPhiCnt] = pfft_p[k];
			}
			_nHarmPhiCnt++;
			//////////////////////////////////////////////


			//ver0.39 - del
			/**
			//////////////////// THD /////////////////////
			if (k > nNumIntCycle && k < nNumIntCycle * 51) {
				*pTHD = *pTHD + (pfft_a[k] * pfft_a[k]);
			}
			//////////////////////////////////////////////
			/**/
		}

		//_nHarmRmsCnt++;

		// init for inter harmonic calculation
		//########################
		pfft_IH[k] = 0;
		//########################
	
	}// end of for (k = 0; k < (_nIntCycleSamples / 2); k++) 


	/////////////  Inter Harmonics ///////////////
	_nInterHPhiCnt = 0;
	for (ik = 0; ik < (nIntCycleSamples / 2); ik += nNumIntCycle)
	{
		for (r = 1; r < (nNumIntCycle); r++)
		{
			pfft_IH[ik] = pfft_IH[ik] + (pfft_a[ik + r] * pfft_a[ik + r]);
		}
		pfft_IH[ik] = sqrt(pfft_IH[ik]);

		// here !!!
		//	harmonic rms
		//	harmonic phi
		//	harmonics must be up to 50th. so need to check if count is < 50 or not
		//	increment harmonic count
		if (_nInterHPhiCnt < 50)
			pInterharmonic[_nInterHPhiCnt] = pfft_IH[ik];

		_nInterHPhiCnt++;
	}
	//////////////////////////////////////////////


	//ver0.39 - del
	/**
	/////////////////// THD //////////////////////
	*pTHD = sqrt(*pTHD) / pfft_a[nNumIntCycle] * 100;
	/**/


	//ver0.56
	///////////// Total RMS for VLL ///////////////
	*pTotal_rms = sqrt(*pTotal_rms);


	/**
	////// Fundamental RMS, PHI and DC RMS ///////
	// Fundamental RMS, PHI and DC RMS
	*pFund_rms = pfft_a[nNumIntCycle];
	*pFund_angle = pfft_p[nNumIntCycle];
	*pDC_rms = pfft_a[0] / sqrt((double)2) / 2;
	/**/


	/////////////////// RMS //////////////////////
	nSum = 0;
	for (k = nFirstIndex; k < nLastIndex; k++)
		nSum = nSum + (pRawBufAI[nChIndex][k] * pRawBufAI[nChIndex][k]);

	*prms = sqrt(nSum / nCalcSize);


#ifdef _DEBUG
	//if (nChIndex == 0)
	//	_DebugString("Neo -  Dll : Phi Count (%d), IH Count (%d)", _nHarmPhiCnt, _nInterHPhiCnt);
	//_DebugString("Neo -  Dll :--------------");
#endif


	//QueryPerformanceCounter(&g_nEndTime1);
	//g_dPassTimeU = (double)(g_nEndTime1.QuadPart - g_nStartTime1.QuadPart) / ((double)(g_nFreq1.QuadPart));
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : ChNo[%d], Set In FFT(%10.7f ms)", nChIndex, g_dPassTimeU * 1000);
#endif
}
/**/

void PWR_FUNC_EXPORT PwCalcPower(
	float  **pRawBufAI,			//ver0.0.1.0
	void	*pMappedChannelAI,
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	int		 nFirstIndex,
	int		 nLastIndex,
	int		 nCalcSize,
	float	*pFund_angle,			//ver0.0.1.0	// for P/Q/S & unbalance
	float	*prms,					//ver0.0.1.0


	float	*pP_Active,				//ver0.0.1.0
	float	*pQ_Reactive,			//ver0.0.1.0
	float	*pS_Apparent,			//ver0.0.1.0
	float	*pPowerFactor,			//ver0.0.1.0

	float	*pP_accumulate,			//ver0.0.1.0
	float	*pQ_accumulate,			//ver0.0.1.0
	float	*pS_accumulate,			//ver0.0.1.0

	int		 nGridType,				//v0.0.1.1.a
	float	*pFirstHarmonicCalVal,	//v0.0.1.1.a
	float	*v_mag,					//v0.0.1.1.a
	float	*v_ang					//v0.0.1.1.a

	//int		*pEnergy_count,
	//int		*pEnergy_sec,
	//int		*pEnergy_min,
	//int		*pEnergy_hour
	)
{
	PSYS_POWER_MAPPED_AI pMappedAi = (PSYS_POWER_MAPPED_AI)pMappedChannelAI;


	//int nEnergy_count = *pEnergy_count;
	//int nEnergy_sec = *pEnergy_sec;
	//int	nEnergy_min = *pEnergy_min;
	//int	nEnergy_hour = *pEnergy_hour;



	int k;
	int PWR_THREE_PHASES_COUNT = 3;

	//----------------------------
	// Init for Total Power 
	//----------------------------
	pP_Active[3] = 0;
	pQ_Reactive[3] = 0;
	pS_Apparent[3] = 0;
	pPowerFactor[3] = 0;

	//ver0.0.1.1.c - add - fix a critical bug
	pP_accumulate[3] = 0;
	pQ_accumulate[3] = 0;
	pS_accumulate[3] = 0;
	//----------------------------

	float cos_angle = 0;			//ver0.0.1.0
	float Phase_sign = 0;			//ver0.0.1.0
	float nP_sum = 0;				//ver0.0.1.0
	float delta_angle[4];			//ver0.0.1.0

	int PhaseCh_V = 0, PhaseCh_I = 0;
	
	int PhaseCh_V2 = 0;		//ver0.0.1.1.a
	int phase_no2 = 0;

	// Power of Phase R,S,T 
	for (int phase_no = 0; phase_no < PWR_THREE_PHASES_COUNT; phase_no++)
	{
		phase_no2 = (phase_no + 2) % 3;

		if (pMappedAi->MappedInfo[phase_no].bPower_on)
		{
			PhaseCh_V = pMappedAi->MappedInfo[phase_no].V_Ch_Index;
			PhaseCh_I = pMappedAi->MappedInfo[phase_no].I_Ch_Index;

			//#ifdef _DEBUG
			//_DebugString("Neo -  Dll : PwCalcPower() - Phase no[%d], VIndex(%d), I Index(%d)", phase_no, PhaseCh_V, PhaseCh_I);
			//#endif

			//#############################################################################################
			//ver0.0.1.1.a - mod - for delta
			/*
			nP_sum = 0;
			for (k = nFirstIndex; k < nLastIndex; k++)
				nP_sum = nP_sum + (pRawBufAI[PhaseCh_V][k] * pRawBufAI[PhaseCh_I][k]);

			Phase_sign = pFund_angle[PhaseCh_V] - pFund_angle[PhaseCh_I];

			//////////////////    phase P,S    ///////////////////////// 
			pP_Active[phase_no] = (nP_sum / nCalcSize);
			pS_Apparent[phase_no] = prms[PhaseCh_V] * prms[PhaseCh_I];
			/**/

			/**/
			nP_sum = 0;
			if (nGridType == 1)
			{
				PhaseCh_V2 = pMappedAi->MappedInfo[phase_no2].V_Ch_Index;


				//ver0.0.1.2 - add - to fix a critical bug !!!
				if ((PhaseCh_V == -1) || (PhaseCh_V2 == -1) || (PhaseCh_I == -1))
					continue;


				for (k = nFirstIndex; k < nLastIndex; k++)
				{
					nP_sum = nP_sum + ((pRawBufAI[PhaseCh_V][k] - pRawBufAI[PhaseCh_V2][k]) / 3 *
						pRawBufAI[PhaseCh_I][k] * pFirstHarmonicCalVal[phase_no]);
				}

				Phase_sign = v_ang[phase_no] - pFund_angle[PhaseCh_I];

				//////////////////    phase P,S    ///////////////////////// 
				pP_Active[phase_no] = (nP_sum / nCalcSize);
				pS_Apparent[phase_no] = v_mag[phase_no] * prms[PhaseCh_I];
			}
			else
			{
				//ver0.0.1.2 - add - to fix a critical bug !!!
				if ((PhaseCh_V == -1) || (PhaseCh_I == -1))
					continue;


				for (k = nFirstIndex; k < nLastIndex; k++)
				{
					nP_sum = nP_sum + (pRawBufAI[PhaseCh_V][k] * 
										pRawBufAI[PhaseCh_I][k] * pFirstHarmonicCalVal[phase_no]);
				}

				Phase_sign = pFund_angle[PhaseCh_V] - pFund_angle[PhaseCh_I];

				//////////////////    phase P,S    ///////////////////////// 
				pP_Active[phase_no] = (nP_sum / nCalcSize);
				pS_Apparent[phase_no] = prms[PhaseCh_V] * prms[PhaseCh_I];
			}
			/**/
			//ver0.0.1.1.a - mod - for delta - end
			//#############################################################################################


			//////////////////     phase pf    /////////////////////////
			if (pS_Apparent[phase_no] != 0) {
				pPowerFactor[phase_no] = pP_Active[phase_no] / pS_Apparent[phase_no];
			}
			else {
				pPowerFactor[phase_no] = 0;
			}


			///////////////////    phase Q    //////////////////////////
			cos_angle = acos(pPowerFactor[phase_no]) * (180 / M_PI);

			if ((0 <= Phase_sign && Phase_sign < 180) || (-360 <= Phase_sign && Phase_sign < -180)) {  // 1st & 2nd Quadrant
				delta_angle[phase_no] = cos_angle;
			}
			else if ((180 <= Phase_sign && Phase_sign < 360) || (-180 <= Phase_sign && Phase_sign < 0)) {  // 3rd & 4th Quadrant
				delta_angle[phase_no] = cos_angle * -1;
			}

			pQ_Reactive[phase_no] = pS_Apparent[phase_no] * sin(delta_angle[phase_no] * M_PI / 180);
			/////////////////////////////////////////////////////////////////////

			//////////////    Total P,Q,S    /////////////// 
			pP_Active[3] += pP_Active[phase_no];
			pQ_Reactive[3] += pQ_Reactive[phase_no];
			pS_Apparent[3] += pS_Apparent[phase_no];

			//////////////    phase Energy    ////////////// kwh
			pP_accumulate[phase_no] += (pP_Active[phase_no] * 0.2 / 3600);
			pS_accumulate[phase_no] += (pS_Apparent[phase_no] * 0.2 / 3600);
			pQ_accumulate[phase_no] += (pQ_Reactive[phase_no] * 0.2 / 3600);

			//////////////    Total Energy    ////////////// total kwh
			pP_accumulate[3] += pP_accumulate[phase_no];
			pS_accumulate[3] += pS_accumulate[phase_no];
			pQ_accumulate[3] += pQ_accumulate[phase_no];
		}
	}

	////////////////////////////////////////////////////// total pf
	if (pS_Apparent[3] != 0) {
		pPowerFactor[3] = pP_Active[3] / pS_Apparent[3];
	}
	else if (pS_Apparent[3] == 0) {
		pPowerFactor[3] = 0;
	}

	/**
	////////////////////////////////////////////////////////////// measurement time
	if ((pMappedAi->MappedInfo[0].bPower_on) ||
	(pMappedAi->MappedInfo[1].bPower_on) ||
	(pMappedAi->MappedInfo[2].bPower_on))
	{
	nEnergy_count++;

	if (nEnergy_count == 5) {
	nEnergy_sec++;
	nEnergy_count = 0;
	}

	if (nEnergy_sec == 60) {
	nEnergy_min++;
	nEnergy_sec = 0;
	}

	if (nEnergy_min == 60) {
	nEnergy_hour++;
	nEnergy_min = 0;
	}
	}
	*pEnergy_count = nEnergy_count;
	*pEnergy_sec = nEnergy_sec;
	*pEnergy_min = nEnergy_min;
	*pEnergy_hour = nEnergy_hour;
	/**/
	//#ifdef _DEBUG
	//	_DebugString("Neo -  Dll : PwCalcPower() - Total Power Factor (%10.6f)", pPowerFactor[3]);
	//#endif
}
/**/


void PWR_FUNC_EXPORT PwCalcUnbalance(
	void	*pMappedChannelAI,

	float	*pFund_rms,			//ver0.0.1.0		// for unbalance
	float	*pFund_angle,		//ver0.0.1.0		// for P/Q/S & unbalance
	float	*prms,				//ver0.0.1.0

	float	*punb_Vzero,		//ver0.0.1.0
	float	*punb_Vneg,			//ver0.0.1.0
	float	*punb_Izero,		//ver0.0.1.0
	float	*punb_Ineg			//ver0.0.1.0
	)
{
	//PSYS_POWER_AI_CH_SEL_INFO pSelAiInfo = (PSYS_POWER_AI_CH_SEL_INFO)pSelChAi;
	PSYS_POWER_MAPPED_AI pMappedAi = (PSYS_POWER_MAPPED_AI)pMappedChannelAI;

#ifdef _DEBUG
	//for (int k = 0; k < 3; k++)
	//	_DebugString("Neo -  Dll : MappedAI[%d] - On?[%s], V[%d], I[%d]", 
	//		k, pMappedAi->MappedInfo[k].bPower_on ? "ON" : "OFF", 
	//		pMappedAi->MappedInfo[k].VoltageChNo, pMappedAi->MappedInfo[k].CurrentChNo);
	//_DebugString("Neo -  Dll :--------------");
#endif
	int PWR_THREE_PHASES_COUNT = 3;

	// unbalance
	int chno_0 = 0, chno_1 = 0, chno_2 = 0;
	float tmp_re = 0, tmp_im = 0;											//ver0.0.1.0
	float Vzero = 0, Vpos = 0, Vneg = 0, Izero = 0, Ipos = 0, Ineg = 0;		//ver0.0.1.0

	////////////////////////////////////////////////////////////// voltage unbalance
	if ((pMappedAi->MappedInfo[0].VoltageChNo != -1) &&
		(pMappedAi->MappedInfo[1].VoltageChNo != -1) &&
		(pMappedAi->MappedInfo[2].VoltageChNo != -1))
	{
	chno_0 = pMappedAi->MappedInfo[0].V_Ch_Index;
	chno_1 = pMappedAi->MappedInfo[1].V_Ch_Index;
	chno_2 = pMappedAi->MappedInfo[2].V_Ch_Index;


	tmp_re = pFund_rms[chno_0] + (pFund_rms[chno_1] * cos((pFund_angle[chno_1] + 120) * M_PI / 180)) + (pFund_rms[chno_2] * cos((pFund_angle[chno_2] - 120) * M_PI / 180));
	tmp_im = (pFund_rms[chno_1] * sin((pFund_angle[chno_1] + 120) * M_PI / 180)) + (pFund_rms[chno_2] * sin((pFund_angle[chno_2] - 120) * M_PI / 180));
	Vpos = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im)) / 3;

	tmp_re = (pFund_rms[chno_0] + (pFund_rms[chno_1] * cos((pFund_angle[chno_1] - 120) * M_PI / 180)) + (pFund_rms[chno_2] * cos((pFund_angle[chno_2] + 120) * M_PI / 180))) / 3;
	tmp_im = ((pFund_rms[chno_1] * sin((pFund_angle[chno_1] - 120) * M_PI / 180)) + (pFund_rms[chno_2] * sin((pFund_angle[chno_2] + 120) * M_PI / 180))) / 3;
	Vneg = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im));

	tmp_re = (pFund_rms[chno_0] + (pFund_rms[chno_1] * cos(pFund_angle[chno_1] * M_PI / 180)) + (pFund_rms[chno_2] * cos(pFund_angle[chno_2] * M_PI / 180))) / 3;
	tmp_im = ((pFund_rms[chno_1] * sin(pFund_angle[chno_1] * M_PI / 180)) + (pFund_rms[chno_2] * sin(pFund_angle[chno_2] * M_PI / 180))) / 3;
	Vzero = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im));

	*punb_Vzero = (Vzero / Vpos) * 100;
	*punb_Vneg = (Vneg / Vpos) * 100;
	}



	//////////////////////////////////////////////////////////////////////////// current unbalance
	if ((pMappedAi->MappedInfo[0].CurrentChNo != -1) &&
		(pMappedAi->MappedInfo[1].CurrentChNo != -1) &&
		(pMappedAi->MappedInfo[2].CurrentChNo != -1))
	{
	chno_0 = pMappedAi->MappedInfo[0].I_Ch_Index;
	chno_1 = pMappedAi->MappedInfo[1].I_Ch_Index;
	chno_2 = pMappedAi->MappedInfo[2].I_Ch_Index;

	tmp_re = pFund_rms[chno_0] + (pFund_rms[chno_1] * cos((pFund_angle[chno_1] + 120) * M_PI / 180)) + (pFund_rms[chno_2] * cos((pFund_angle[chno_2] - 120) * M_PI / 180));
	tmp_im = (pFund_rms[chno_1] * sin((pFund_angle[chno_1] + 120) * M_PI / 180)) + (pFund_rms[chno_2] * sin((pFund_angle[chno_2] - 120) * M_PI / 180));
	Ipos = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im)) / 3;

	tmp_re = (pFund_rms[chno_0] + (pFund_rms[chno_1] * cos((pFund_angle[chno_1] - 120) * M_PI / 180)) + (pFund_rms[chno_2] * cos((pFund_angle[chno_2] + 120) * M_PI / 180))) / 3;
	tmp_im = ((pFund_rms[chno_1] * sin((pFund_angle[chno_1] - 120) * M_PI / 180)) + (pFund_rms[chno_2] * sin((pFund_angle[chno_2] + 120) * M_PI / 180))) / 3;
	Ineg = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im));

	tmp_re = (pFund_rms[chno_0] + (pFund_rms[chno_1] * cos(pFund_angle[chno_1] * M_PI / 180)) + (pFund_rms[chno_2] * cos(pFund_angle[chno_2] * M_PI / 180))) / 3;
	tmp_im = ((pFund_rms[chno_1] * sin(pFund_angle[chno_1] * M_PI / 180)) + (pFund_rms[chno_2] * sin(pFund_angle[chno_2] * M_PI / 180))) / 3;
	Izero = sqrt((tmp_re * tmp_re) + (tmp_im * tmp_im));

	*punb_Izero = (Izero / Ipos) * 100;
	*punb_Ineg = (Ineg / Ipos) * 100;
	}
}

void PWR_FUNC_EXPORT PwCalcFFTMagnitude(
	int		 nScanSize,
	float   *pInFFT,
	float  *pOutFFT,
	int		nValueType
	)
{
	int i;
	float fTemp;

	for (i = 0; i < (nScanSize / 2); i++)
	{
		fTemp = sqrtf(((pInFFT[2 * i] * pInFFT[2 * i]) + (pInFFT[2 * i + 1] * pInFFT[2 * i + 1])))
			* 2 / nScanSize;

		if (nValueType == 0)	// dB
			pOutFFT[i] = 20 * log10f(fTemp);
		else					// Linear
			pOutFFT[i] = fTemp;
	}
}


//ver0.0.1.0 - no use
// Linear / dB
void PWR_FUNC_EXPORT PwCalcFFTMagnitudeDouble(
	int		 nScanSize,
	float	*pInFFT,
	float	*pOutFFT			//ver0.0.1.0
	)
{
	int i;
	float dTemp;

	for (i = 0; i < (nScanSize / 2); i++)
	{
		dTemp = sqrt(((pInFFT[2 * i] * pInFFT[2 * i]) + (pInFFT[2 * i + 1] * pInFFT[2 * i + 1])))
			* 2 / nScanSize;

		// dB
		pOutFFT[i] = 20 * log10(dTemp);
	}
}

#define CONST_FFT_WINDOW_TYPE_NONE				 0
#define CONST_FFT_WINDOW_TYPE_BARTLETT			 1
#define CONST_FFT_WINDOW_TYPE_HAMMING			 2
#define CONST_FFT_WINDOW_TYPE_HANNING			 3
#define CONST_FFT_WINDOW_TYPE_BLACKMAN			 4
#define CONST_FFT_WINDOW_TYPE_BLACKMAN_HARRIS	 5
#define CONST_FFT_WINDOW_TYPE_WELCH				 6
#define CONST_FFT_WINDOW_TYPE_GAUSSIAN2_5		 7
#define CONST_FFT_WINDOW_TYPE_GAUSSIAN3_5		 8
#define CONST_FFT_WINDOW_TYPE_GAUSSIAN4_5		 9
#define CONST_FFT_WINDOW_TYPE_FLATTOP			 10

float PWR_FUNC_EXPORT PwCalcWindowDensityFactor(
	int nFFTScanSize, 
	int nWinType
	)
{
	int i;
	float fAmpScaleFactor = 1.0;
	float fGausianPriScale = 0;

	float fScanSize = (float)nFFTScanSize;
	float fHalfScanSize = (float)nFFTScanSize / 2;
	float fScanSizeMinusOne = (float)nFFTScanSize - 1;

	switch (nWinType)
	{
		case CONST_FFT_WINDOW_TYPE_NONE:
			break;

		case CONST_FFT_WINDOW_TYPE_BARTLETT:		// Bartlett (triangular) window
			fAmpScaleFactor = 0;

			for (i = 0; i < (nFFTScanSize / 2); i++)
			{
				fAmpScaleFactor = fAmpScaleFactor + (i / fHalfScanSize);
				fAmpScaleFactor = fAmpScaleFactor + (1.0f - (i / fHalfScanSize));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_HAMMING:			// Hamming
			fAmpScaleFactor = 0;
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor + (0.54f - 0.46f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_HANNING:			// Hanning
			fAmpScaleFactor = 0;
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor + (0.50f - 0.50f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_BLACKMAN:		// Blackman
			fAmpScaleFactor = 0;
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor +
					(0.42f - 0.5f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
							0.08f * cosf(4 * M_PI * i / fScanSizeMinusOne));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_BLACKMAN_HARRIS:	// Blackman-Harris
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor +
					(0.35875f - 0.48829f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
								0.14128f * cosf(4 * M_PI * i / fScanSizeMinusOne) -
								0.01168f * cosf(6 * M_PI * i / fScanSizeMinusOne));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_WELCH:			// Welch
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor + (4.0f * i / fScanSize * (1.0f - (i / fScanSize)));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN2_5:		// Gaussian (a=2.5)
			// Precalculate some values, and simplify the fmla to try and reduce overhead
			fGausianPriScale = -2.0f * 2.5f * 2.5f;

			for (i = 0; i <nFFTScanSize; i++)
			{
				// full
				fAmpScaleFactor = fAmpScaleFactor +
					expf(-0.5f * (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)) *
								 (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN3_5:		// Gaussian (a=3.5)
			fGausianPriScale = -2.0f * 3.5f * 3.5f;

			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor +
					expf(fGausianPriScale * (0.25f + ((i / fScanSize) * (i / fScanSize)) - (i / fScanSize)));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN4_5:		// Gaussian (a=4.5)
			fGausianPriScale = -2.0f * 4.5f *4.5f;

			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor +
					expf(fGausianPriScale * (0.25f + ((i / fScanSize) * (i / fScanSize)) - (i / fScanSize)));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		case CONST_FFT_WINDOW_TYPE_FLATTOP:			// Flat Top
			fAmpScaleFactor = 0;
			for (i = 0; i < nFFTScanSize; i++)
			{
				fAmpScaleFactor = fAmpScaleFactor +
					(1.0f - 1.93f  * cosf(2 * M_PI * i / fScanSizeMinusOne) +
							1.29f  * cosf(4 * M_PI * i / fScanSizeMinusOne) +
							0.388f * cosf(6 * M_PI * i / fScanSizeMinusOne) +
							0.028f * cosf(8 * M_PI * i / fScanSizeMinusOne));
			}
			fAmpScaleFactor = nFFTScanSize / fAmpScaleFactor;
			break;

		default:
			break;
	}

	return fAmpScaleFactor;
}


/**
void PWR_FUNC_EXPORT PwWindowing(
	float*	in_, 
	int nFFTScanSize, 
	int nWinType, 
	float fScaleFactor
	)
{
	int i;
	float fGausianPriScale = 0;

	float fScanSize = (float)nFFTScanSize;
	float fHalfScanSize = (float)nFFTScanSize / 2;
	float fScanSizeMinusOne = (float)nFFTScanSize - 1;

	switch (nWinType)
	{
		case CONST_FFT_WINDOW_TYPE_NONE:
			break;

		case CONST_FFT_WINDOW_TYPE_BARTLETT:		// Bartlett (triangular) window
			for (i = 0; i < (nFFTScanSize / 2); i++)
			{
				in_[i*2]   = fScaleFactor * in_[i*2]   * (i / fHalfScanSize);
				in_[i*2+1] = fScaleFactor * in_[i*2+1] * (i / fHalfScanSize);

				in_[(i+(nFFTScanSize/2))*2]	  = fScaleFactor * in_[(i+(nFFTScanSize/2))*2]   * (1.0f - (i / fHalfScanSize));
				in_[(i+(nFFTScanSize/2))*2+1] = fScaleFactor * in_[(i+(nFFTScanSize/2))*2+1] * (1.0f - (i / fHalfScanSize));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_HAMMING:			// Hamming
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2]   = fScaleFactor * in_[i*2]   * (0.54f - 0.46f * cosf(2 * M_PI * i / fScanSizeMinusOne));
				in_[i*2+1] = fScaleFactor * in_[i*2+1] * (0.54f - 0.46f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_HANNING:			// Hanning
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2]   = fScaleFactor * in_[i* 2]  * (0.50f - 0.50f * cosf(2 * M_PI * i / fScanSizeMinusOne));
				in_[i*2+1] = fScaleFactor * in_[i*2+1] * (0.50f - 0.50f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_BLACKMAN:		// Blackman
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2] = fScaleFactor * in_[i*2] *
					(0.42f - 0.5f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 0.08f * cosf(4 * M_PI * i / fScanSizeMinusOne));

				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					(0.42f - 0.5f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 0.08f * cosf(4 * M_PI * i / fScanSizeMinusOne));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_BLACKMAN_HARRIS:	// Blackman-Harris
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2] = fScaleFactor * in_[i*2] *
					(0.35875f - 0.48829f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
					0.14128f * cosf(4 * M_PI * i / fScanSizeMinusOne) -
					0.01168f * cosf(6 * M_PI * i / fScanSizeMinusOne));

				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					(0.35875f - 0.48829f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
					0.14128f * cosf(4 * M_PI * i / fScanSizeMinusOne) -
					0.01168f * cosf(6 * M_PI * i / fScanSizeMinusOne));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_WELCH:			// Welch
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2]   = fScaleFactor * in_[i*2]   * 4 * i / fScanSize *  (1 - (i / fScanSize));
				in_[i*2+1] = fScaleFactor * in_[i*2+1] * 4 * i / fScanSize *  (1 - (i / fScanSize));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN2_5:		// Gaussian (a=2.5)
			fGausianPriScale = -2.0f * 2.5f * 2.5f;

			for (i = 0; i < nFFTScanSize; i++)
			{
				// full
				in_[i*2] = fScaleFactor * in_[i*2] *
					expf(-0.5f * (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)) *
								 (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)));
				
				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					expf(-0.5f * (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)) *
								 (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN3_5:		// Gaussian (a=3.5)
			fGausianPriScale = -2.0f * 3.5f * 3.5f;

			for (i = 0; i < nFFTScanSize; i++)
			{
				// reduced
				in_[i*2] = fScaleFactor * in_[i*2] *
					expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));

				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_GAUSSIAN4_5:		// Gaussian (a=4.5)
			fGausianPriScale = -2.0f * 4.5f * 4.5f;

			for (i = 0; i < nFFTScanSize; i++)
			{
				// reduced
				in_[i*2] = fScaleFactor * in_[i*2] *
					expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));

				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));
			}
			break;

		case CONST_FFT_WINDOW_TYPE_FLATTOP:			// Flat Top
			for (i = 0; i < nFFTScanSize; i++)
			{
				in_[i*2] = fScaleFactor * in_[i*2] *
					(1 - 1.93f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 1.29f  * cosf(4 * M_PI * i / fScanSizeMinusOne) +
						0.388f * cosf(6 * M_PI * i / fScanSizeMinusOne) + 0.028f * cosf(8 * M_PI * i / fScanSizeMinusOne));

				in_[i*2+1] = fScaleFactor * in_[i*2+1] *
					(1 - 1.93f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 1.29f  * cosf(4 * M_PI * i / fScanSizeMinusOne) +
						0.388f * cosf(6 * M_PI * i / fScanSizeMinusOne) + 0.028f * cos(8 * M_PI * i / fScanSizeMinusOne));
			}
			break;

		default:
			break;
	}
}
/**/

/**/
void PWR_FUNC_EXPORT PwWindowing(
	float*	in_,
	int nFFTScanSize,
	int nWinType,
	float fScaleFactor
	)
{
	int i;
	float fGausianPriScale = 0;

	float fScanSize = (float)nFFTScanSize;
	float fHalfScanSize = (float)nFFTScanSize / 2;
	float fScanSizeMinusOne = (float)nFFTScanSize - 1;

	switch (nWinType)
	{
	case CONST_FFT_WINDOW_TYPE_NONE:
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor;
			in_[i * 2 + 1] = fScaleFactor;
		}
		break;

	case CONST_FFT_WINDOW_TYPE_BARTLETT:		// Bartlett (triangular) window
		for (i = 0; i < (nFFTScanSize / 2); i++)
		{
			in_[i * 2] = fScaleFactor * (i / fHalfScanSize);
			in_[i * 2 + 1] = fScaleFactor * (i / fHalfScanSize);

			in_[(i + (nFFTScanSize / 2)) * 2] = fScaleFactor * (1.0f - (i / fHalfScanSize));
			in_[(i + (nFFTScanSize / 2)) * 2 + 1] = fScaleFactor * (1.0f - (i / fHalfScanSize));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_HAMMING:			// Hamming
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor * (0.54f - 0.46f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			in_[i * 2 + 1] = fScaleFactor * (0.54f - 0.46f * cosf(2 * M_PI * i / fScanSizeMinusOne));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_HANNING:			// Hanning
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor * (0.50f - 0.50f * cosf(2 * M_PI * i / fScanSizeMinusOne));
			in_[i * 2 + 1] = fScaleFactor * (0.50f - 0.50f * cosf(2 * M_PI * i / fScanSizeMinusOne));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_BLACKMAN:		// Blackman
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor *
				(0.42f - 0.5f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 0.08f * cosf(4 * M_PI * i / fScanSizeMinusOne));

			in_[i * 2 + 1] = fScaleFactor *
				(0.42f - 0.5f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 0.08f * cosf(4 * M_PI * i / fScanSizeMinusOne));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_BLACKMAN_HARRIS:	// Blackman-Harris
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor * 
				(0.35875f - 0.48829f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
				0.14128f * cosf(4 * M_PI * i / fScanSizeMinusOne) -
				0.01168f * cosf(6 * M_PI * i / fScanSizeMinusOne));

			in_[i * 2 + 1] = fScaleFactor * 
				(0.35875f - 0.48829f * cosf(2 * M_PI * i / fScanSizeMinusOne) +
				0.14128f * cosf(4 * M_PI * i / fScanSizeMinusOne) -
				0.01168f * cosf(6 * M_PI * i / fScanSizeMinusOne));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_WELCH:			// Welch
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor * 4 * i / fScanSize *  (1 - (i / fScanSize));
			in_[i * 2 + 1] = fScaleFactor * 4 * i / fScanSize *  (1 - (i / fScanSize));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_GAUSSIAN2_5:		// Gaussian (a=2.5)
		fGausianPriScale = -2.0f * 2.5f * 2.5f;

		for (i = 0; i < nFFTScanSize; i++)
		{
			// full
			in_[i * 2] = fScaleFactor *
				expf(-0.5f * (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)) *
				(fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)));

			in_[i * 2 + 1] = fScaleFactor *
				expf(-0.5f * (fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)) *
				(fGausianPriScale * ((i - fHalfScanSize) / fHalfScanSize)));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_GAUSSIAN3_5:		// Gaussian (a=3.5)
		fGausianPriScale = -2.0f * 3.5f * 3.5f;

		for (i = 0; i < nFFTScanSize; i++)
		{
			// reduced
			in_[i * 2] = fScaleFactor *
				expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));

			in_[i * 2 + 1] = fScaleFactor *
				expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_GAUSSIAN4_5:		// Gaussian (a=4.5)
		fGausianPriScale = -2.0f * 4.5f * 4.5f;

		for (i = 0; i < nFFTScanSize; i++)
		{
			// reduced
			in_[i * 2] = fScaleFactor *
				expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));

			in_[i * 2 + 1] = fScaleFactor *
				expf(fGausianPriScale * (0.25f + ((i / fScanSize)*(i / fScanSize)) - (i / fScanSize)));
		}
		break;

	case CONST_FFT_WINDOW_TYPE_FLATTOP:			// Flat Top
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor *
				(1 - 1.93f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 1.29f  * cosf(4 * M_PI * i / fScanSizeMinusOne) +
				0.388f * cosf(6 * M_PI * i / fScanSizeMinusOne) + 0.028f * cosf(8 * M_PI * i / fScanSizeMinusOne));

			in_[i * 2 + 1] = fScaleFactor *
				(1 - 1.93f * cosf(2 * M_PI * i / fScanSizeMinusOne) + 1.29f  * cosf(4 * M_PI * i / fScanSizeMinusOne) +
				0.388f * cosf(6 * M_PI * i / fScanSizeMinusOne) + 0.028f * cos(8 * M_PI * i / fScanSizeMinusOne));
		}
		break;

	default:
		for (i = 0; i < nFFTScanSize; i++)
		{
			in_[i * 2] = fScaleFactor;
			in_[i * 2 + 1] = fScaleFactor;
		}
		break;
	}
}
/**/

void PWR_FUNC_EXPORT PwSetInFFTMem_Window(
	int		 nChIndex,

	float  **pRawBufAI,			//ver0.0.1.0
	int		 nDmaScanSize,
	int		 nDmaBlockCount,

	bool	 bIsSeparateMem,
	int		 nRPos,
	int		 nFirstLength,
	int		 nPwrScanSize,

	int		 nFirstIndex,
	int		 nIntIndex,

	float   *pWinFactor,
	float   *pInFFTMem
	)
{
	int k;
	int nPos;

	//LARGE_INTEGER		g_nFreq1, g_nStartTime1, g_nEndTime1;
	//double				g_dPassTimeU;

	//QueryPerformanceFrequency(&g_nFreq1);
	//QueryPerformanceCounter(&g_nStartTime1);

#ifdef _DEBUG
	//_DebugString("Neo -  Dll :DMA scan size(%d), Blocks(%d)", nDmaScanSize, nDmaBlockCount);
	//_DebugString("Neo -  Dll :Separate?(%s), Pos(%d), FirstLen(%d), SCAN(%d), FirstIndex(%d), IntIndex(%d)",
	//	bIsSeparateMem ? "YES" : "NO", nRPos, nFirstLength, nPwrScanSize, nFirstIndex, nIntIndex);

	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll :(%d) (%d) Buf [%d][%d] = %f", nRPos, nFirstIndex, nChIndex, k, pRawBufAI[nChIndex][k]);
	//_DebugString("Neo -  Dll :--------------");
#endif

	int windex = 0;
	for (k = nFirstIndex; k < (nRPos + nPwrScanSize); k++)
	{
		nPos = k;

		if (bIsSeparateMem && (k >= (nDmaBlockCount*nDmaScanSize)))
			nPos = (nRPos + k) - (nDmaBlockCount*nDmaScanSize);

		if (k < nIntIndex) {
			//if (nPos < nIntIndex) {
			pInFFTMem[(k - nFirstIndex) * 2] = (float)pRawBufAI[nChIndex][nPos] * pWinFactor[windex * 2];
			pInFFTMem[(k - nFirstIndex) * 2 + 1] = 0;
		}
		else {
			pInFFTMem[(k - nFirstIndex) * 2] = 0;
			pInFFTMem[(k - nFirstIndex) * 2 + 1] = 0;
		}
		
		windex++;
	}
#ifdef _DEBUG
	//for (k = 0; k < 10; k++)
	//	_DebugString("Neo -  Dll : After In  [%d][%d] = %f", nChIndex, k, pInFFTMem[k]);
	//_DebugString("Neo -  Dll :--------------");
#endif


	//QueryPerformanceCounter(&g_nEndTime1);
	//g_dPassTimeU = (double)(g_nEndTime1.QuadPart - g_nStartTime1.QuadPart) / ((double)(g_nFreq1.QuadPart));
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : ChNo[%d], Set In FFT(%10.7f ms)", nChIndex, g_dPassTimeU * 1000);
#endif

}


//ver0.40
//#####################################################################
//
// Flicker
//
//#####################################################################
void PWR_FUNC_EXPORT PwDownSamples(
	float	*pOrgBuf,			//ver0.0.1.0

	int		 nOrgSize,
	int		 nDownRatio,

	float	*pDownBuf			//ver0.0.1.0
	)
{
	int k;
	int nPos;
#ifdef _DEBUG
	//_DebugString("Neo -  Dll : PwDownSamples() - Org Size(%d), DownRatio(%d), Down Size(%d)",
	//	nOrgSize, nDownRatio, nOrgSize/nDownRatio);

	//for (k = 0; k < 4; k++)
	//	_DebugString("Neo -  Dll : PwDownSamples() - Org Buf(%f)", pOrgBuf[k]);
#endif

	// down samples
	int DownSize = nOrgSize / nDownRatio;

	for (k = 0; k < DownSize; k++)
	{
		pDownBuf[k] = pOrgBuf[k * nDownRatio];
	}

#ifdef _DEBUG
	//for (k = 0; k < 4; k++)
	//	_DebugString("Neo -  Dll : PwDownSamples() - Down Buf(%f)", pDownBuf[k]);
	//_DebugString("Neo -  Dll :--------------");
#endif
}

//##################################
float *flk_que_1[3];			//ver0.0.1.0
float *flk_que_2[3];			//ver0.0.1.0
void PWR_FUNC_EXPORT PwAllocFlickerFilterMem(int nCount1, int nCount2)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (flk_que_1[i] != 0)
			delete flk_que_1[i];

		if (flk_que_2[i] != 0)
			delete flk_que_2[i];
	}

	for (i = 0; i < 3; i++)
	{
		flk_que_1[i] = new float[nCount1];			//ver0.0.1.0

		flk_que_2[i] = new float[nCount2];			//ver0.0.1.0
	}
}

void PWR_FUNC_EXPORT PwFreeFlickerFilterMem()
{
	int i;
	for (i = 0; i < 3; i++)
	{
		delete flk_que_1[i];
		delete flk_que_2[i];
	}
}

void PWR_FUNC_EXPORT PwFlicker_Filter(
	int		 nStep,

	int		 nOrder,
	float 	 fSampleRate,		//ver0.0.1.0
	float 	 fCutOff,			//ver0.0.1.0
	int		 nSampleCount,

	float	*pBuf				//ver0.0.1.0
	)
{
	int i;
	if ((nStep == 1) || (nStep == 2))
	{
		for (i = 0; i < 3; i++)
			memcpy(flk_que_1[i], &pBuf[nSampleCount * i], nSampleCount * sizeof(float));			//ver0.0.1.0

		if (nStep == 1)
		{
			Dsp::SimpleFilter <Dsp::Butterworth::HighPass <1>, 3> v;
			v.setup(nOrder, fSampleRate, fCutOff);							//ver0.0.1.0
			v.process(nSampleCount, flk_que_1);
		}
		else
		{
			Dsp::SimpleFilter <Dsp::Butterworth::LowPass <6>, 3> w;
			w.setup(nOrder, fSampleRate, fCutOff);							//ver0.0.1.0
			w.process(nSampleCount, flk_que_1);
		}

		for (i = 0; i < 3; i++)
			memcpy(&pBuf[nSampleCount * i], flk_que_1[i], nSampleCount * sizeof(float));			//ver0.0.1.0
	}
	else
	{
		for (i = 0; i < 3; i++)
			memcpy(flk_que_2[i], &pBuf[nSampleCount * i], nSampleCount * sizeof(float));			//ver0.0.1.0

		Dsp::SimpleFilter <Dsp::Butterworth::LowPass <1>, 3> h;
		h.setup(nOrder, fSampleRate, fCutOff);								//ver0.0.1.0
		h.process(nSampleCount, flk_que_2);

		for (i = 0; i < 3; i++)
			memcpy(&pBuf[nSampleCount * i], flk_que_2[i], nSampleCount * sizeof(float));			//ver0.0.1.0
	}
}

/**/
void PWR_FUNC_EXPORT PwWFilter(
	float In_Arr[],			//ver0.0.1.0
	int nomV,
	float num[],
	float den[],
	float ScanSize,
	int buffer_size,
	float* Out_Arr			//ver0.0.1.0
	)
{
	rt_nonfinite(In_Arr, buffer_size, num, den, Out_Arr);
	//wFilter(In_Arr, nomV, num, den, ScanSize, buffer_size, Out_Arr);

	return;
}

float PWR_FUNC_EXPORT PwCalc_pst(float* In_Arr, int Pst_time)			//ver0.0.1.0
{
	float Pst = 0;		//ver0.0.1.0

	Pst = rtGetInf(In_Arr, Pst_time);
	//Pst = calc_pst(In_Arr, Pst_time);

	return Pst;
}
/**/