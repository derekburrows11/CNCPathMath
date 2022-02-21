// PolySegFit.h: interface for the CPolySegFit class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYSEGFIT_H__3D50A280_889C_11D8_86C3_0008A15E291C__INCLUDED_)
#define AFX_POLYSEGFIT_H__3D50A280_889C_11D8_86C3_0008A15E291C__INCLUDED_


#include "PointSumInfo.h"



#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000





class CPolySegFit
{
public:
	enum FIT		// Fit type
	{
		FIT_FREE = 0,
		FIT_SET,
		FIT_ATPOINT,
	};
	enum		// Segment Fit Flags
	{
		SFF_POSI = 0x01,
		SFF_POSF = 0x02,
		SFF_DIRI = 0x04,
		SFF_DIRF = 0x08,
	};



	CPolySegFit();
	virtual ~CPolySegFit();

	void Reset();
	void SetPointBuffer(CVector* pBuffer, int sizeBuffer);
	void SetNumPointsInBuffer(int numPts);
	void SetSBuffer(double* pBuffer, int sizeBuffer);
	
	bool AddPoint(const CVector* pPt) { return AddPoint(*pPt); }
	bool AddPoint(const CVector& pt);
	void RemoveFirstPoints(int iNum);
	void SetToInitialPoint();
	void SetToFinalPoint();
	void SetInitialPoint(const CVector& vtPosI, double SInit);
	void SetFinalPoint(const CVector& vtPosF, double SFinal);
	void SetInitialDir(const CVector& vtDirI);
	void SetFinalDir(const CVector& vtDirF);

	void SetSInit(double sInit) { m_SInit = sInit; }
	void SetSFinal(double sFinal) { m_SFinal = sFinal; }
	double GetSInit() { return m_SInit; }
	double GetSFinal() { return m_SFinal; }

	void SetNumPointsToUse(int numPtsToUse);
	void UseAllPoints();
	int GetNumPoints() { return m_numPts; }
	int GetNumPointsStored() { return m_numPtsStored; }

	double GetMinControlSpan();
	double GetMaxResidual() { return m_MaxResidual; }

	void SolveCubic();




// data
protected:
	// point data
	CVector* m_arPts;
	double* m_arS;
	int m_iPtBufferSize;
	int m_iSBufferSize;
	bool m_barPtsFromCArray;
	bool m_barSFromCArray;
	CArray<CVector, CVector&> m_PtArray;
	CArray<double, double> m_SArray;
	int m_numPts;
	int m_numPtsStored;
	int m_numPtsForFit;
	bool m_bUseAllPoints;
	bool m_bPreFit;

	enum FIT m_nFitPosInit;
	enum FIT m_nFitDirInit;
	enum FIT m_nFitPosFinal;
	enum FIT m_nFitDirFinal;

	CPointSPowerData m_PointSums;

	double m_SInit;
	double m_SFinal;
	double m_SSpanOrig;

	// possible boundary conditions
protected:
	CVector m_vtPosInit;
	CVector m_vtPosFinal;
	CVector m_vtDirInitUnit;
	CVector m_vtDirFinalUnit;

	// fitted results
protected:
	double m_dsAvgNorm;
	CVector m_vtSumResidualSq;
	double m_SumResidualSq;		// sum of all axis ResidualSq
	double m_SumResidual;		// sqrt(m_SumResidualSq)
	double m_SumResidualAvg;	// sqrt(m_SumResidualSq / num-?)
	double m_MaxResidual;

public:
	CVector m_Poly[4];
	CVector m_Bez[4];
	CVector m_BezDeriv[4];

	double m_magEnds;
	double m_magControl[2];


// implementation functions
protected:
	void SetEndConditions();
	void SetSValues();
	void NormaliseSValues();
	void AdjustSValues();
	void SetSums();

	void GetResiduals();

	void GetBezierFromPoly();
	void GetBezDerivFromPoly();
	void GetBezDerivFromBezier();
	void GetBezierFromBezDeriv();


	// fitting functions
	bool CheckNumPoints(int numPtsReq, int iLine = 0);
	void FitCubic();
	void FitCubicIPos();
	void FitCubicIPosDir();
	void FitCubicIPosDeriv(double Kv0);
	void FitCubicTwoPos();
	void FitCubicTwoPos0to1();
	void FitCubicTwoPosTwoDeriv0to1(double Kv0, double Kv1);
	void FitCubicTwoPosTwoDir0to1();
	void FitCubicTwoPosIDeriv0to1(double Kv0);
	void FitCubicTwoPosFDeriv0to1(double Kv1);

	friend class CPolySegDblFit;
};

#endif // !defined(AFX_POLYSEGFIT_H__3D50A280_889C_11D8_86C3_0008A15E291C__INCLUDED_)
