// PointSumInfo.h: interface for the CPointSumInfo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POINTSUMINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
#define AFX_POINTSUMINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_


#include <afxtempl.h>
#include "Vector.h"


#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


////////////////////////////////



struct CPointData : CVect2
{
//	stores data for each point when curve fitting
	int n;
	double yPoly, yError, yErrorSq;

	void Zero() { n = 0; yPoly = yError = yErrorSq = 0; }
};

struct CPointPowerData : CPointData		// base class for testing
{
//	stores data sums for each point when curve fitting
// Data
	double Xp[7];
	double YXp[4];
	double Y2;		// used to calc residual
	int num;			// same as Xp[0] but stored as int!

// Functions
	void Zero();
	void SumUpPoint(const CPointPowerData& prevSum, const CVect2& pt);
};

struct CPointSPowerData
{
//	stores data sums for each point when curve fitting
	double Sp[7];
	double SpV[4][3];
	double V2[3];					// used to calc residual
	bool m_bConstValue[3];		// used to check for constant axis value
	double m_Value[3];
	int num;			// same as Sp[0] but stored as int!

// Functions
	void Zero();
	void SumUpPoint(const CPointSPowerData& prevSum, double s, const CVector& pt);
	void SumUpPoint(double s, const CVector& pt);
	void ScaleS(double scale);
};






class CPointSumInfoState
{
public:
	double m_MinListXStep, m_MaxListXStep;
	double m_MinListX, m_MaxListX;
	int m_iNumPoints;		// num points when state saved
};

typedef CList<CPointPowerData, CPointPowerData&> CPointPowerDataList;

class CPointSumInfo
{
public:
	enum {
		OK = 0,
		STEP_TOO_SMALL,
		STEP_TOO_BIG,
		POINT_WITHIN_LIST,
	};
public:
	CPointSumInfo();
	bool AddPoint(CVect2& pt);
	void ClearAll();
//	int ClearAllButFirst();		// not used now
	CPointPowerData& GetFirstSums();
	CPointPowerData& GetLastSums();
	POSITION GetFirstSumsPos();
	CPointPowerData& GetNextSums(POSITION& pos);

	int GetSumsFor(CPointPowerData& sumPoints, int iInit, int iFinal);
	double GetAvgStepSize();
	int GetNumPoints() { return m_SumPointsList.GetCount(); }
	
	void SaveState(CPointSumInfoState& psis);
	bool RetrieveState(CPointSumInfoState& psis);


//Data members
public:	// for test only, normaly   protected:
	CPointPowerDataList m_SumPointsList;

	double m_StepVariation;			// ratio of largest to smallest allowable step size
	double m_MinListXStep, m_MaxListXStep;
	double m_MinListX, m_MaxListX;

};


#endif // !defined(AFX_POINTSUMINFO_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
