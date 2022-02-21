// PolySegInfo.h: interface for the CPolySegFit class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYSEGINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
#define AFX_POLYSEGINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_



#include "PointSumInfo.h"

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


////////////////////////////////




////////////////////////////////

class CPolySegInfo
{
public:
	CPolySegInfo();
//	CPolySegInfo( CPolySegInfo& orig);
	CPolySegInfo& operator=(const CPolySegInfo& orig);

// Functions
	void Reset();
	void AddPoint(CPointData& pt);
	void AddPointUnconfirmed(CPointData& pt);
	void RemoveLatestPoint();
	CPointData& GetFinalPoint();
	CPointData& GetFirstPoint();
	void AppendPoly(CPolySegInfo& psi);

	void ClearAllButFirst();
	double GetMaxXPoint() { return m_PointSums.m_MaxListX; }
	double GetAvgStepSize() { return m_PointSums.GetAvgStepSize(); }
	int GetNumPoints() { return m_PointSums.GetNumPoints(); }
	void SetPolyFrom(CPolySegInfo& src);
	void CalcPointError(CPointData& pt);

	void SetInitValsFromFinal();		// sets initial x, y, dy etc from final
	void CalcInitPV();
	void CalcFinalPV();

	// fit attributes
	double GetSumResidual();
	double GetPolyErrors();
	double GetVariance() { return GetSumResidual() / (GetNumPoints() - m_iFitDOF); }
	double GetVariance(double Sr) { return Sr / (GetNumPoints() - m_iFitDOF); }
	double GetStandardError(double Sr);

	// fitting functions
	void FitPolyToPoints();
	void FitPolyToPointsContPV();
	void FitDblPolyToPointsContPV(CPolySegInfo& nextPoly, double Wa=1, double Wb=1);

// Functions for testing
	int SetPointArray(CVect2* arPts);
	int SetBezierPoints(CVect2* arPts);	// sets 4 points

	void PlotPoints();
	void PlotBezier();
	void PlotPointsAndBezier() { PlotPoints(); PlotBezier(); }

// Data
public:
	double m_coef[4];		// coefficents
	double m_xPolyInit;
	double m_xPolyFinal;
	double m_yValInit[3];	// value[0] and derivatives[1..2] at initial x
	double m_yValFinal[3];	// value[0] and derivatives[1..2] at final x

	bool m_bPolySet;
//	int m_iNumPoints;
	double m_xPointInit, m_xPointFinal;
	double m_Sr, m_Sabs;	// Sum(err^2), Sum(abs(err)) - depends on poly
	char m_iFitDOF;		// DOF used by fitting algorithm
//	int m_nSumRefInit, m_nSumRefFinal;

	CPointSumInfo m_PointSums;	// should have a list of these for different required X intervals
	CPointSumInfoState m_PrevPointSumState;	// to add a point and enable remove
	double m_SrPrev, m_SabsPrev, m_xPointFinalPrev;

};




#endif // !defined(AFX_POLYSEGINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
