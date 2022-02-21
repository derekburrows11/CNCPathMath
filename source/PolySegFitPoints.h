// PolySegFitPoints.h: interface for the CPolySegFit class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYSEGFITPOINTS_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
#define AFX_POLYSEGFITPOINTS_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_



#include "PolySegInfo.h"

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000








/////////////////////////////////////////


class CPolySegFitPoints
{
public:
	enum {
		PSF_OK = 0,
		PSF_REDUCESTEPSIZE,
		PSF_GOTNEXTPOLYSEGMENT,
	};
protected:
/*	enum {
		FITSTAT_NOPOINTS = 0,
		FITSTAT_NOPOLY,
		FITSTAT_CURRSET,			// a current poly is set
		FITSTAT_REDUCESPAN,
		FITSTAT_POLYCHANGE,
	};
*/
public:
	CPolySegFitPoints();
	virtual ~CPolySegFitPoints();

public:
	void Init();
	void SetSmallestXInc(double xSmallInc) { m_SmallestInc = xSmallInc; }
	int NextPoint(double x, double y);
	int NextPoint(double t, CVect2& pt);
	int NextPoint(double t, CVector& pt);
	void BreakPolyAtNextPoint();
	double GetSuggestedXInc() { return m_SuggestedInc; }
	double GetLastGoodXLoc();


protected:
// Fit property settings
	int m_nFitPtsForFirstPoly;
	int m_nFitPtsForPoly;
	int m_nFitPtsForNextAtChange;


// Data
	CPointData m_CurrPt[3];		// to handle 3 axis or 3D vectors

	CPolySegInfo* m_pPrevPoly;
	CPolySegInfo* m_pCurrPoly;
	CPolySegInfo* m_pTempPoly;

	CPolySegInfo m_Poly1, m_Poly2;;

//	CPolySegInfo m_PrevPoly;
//	CPolySegInfo m_CurrPoly;
//	CPolySegInfo m_NextPoly;

	CList<CPolySegInfo, CPolySegInfo&> m_PolySegList;

	double m_TolY;
	double m_TolVariance;
	double m_SmallestInc;
	double m_SuggestedInc;
	int m_nReductions;
	int m_nFitResultFlags;

	bool m_bForcePolyBreak;

protected:
// Functions
	int NextPointMethod1();
	int NextPointMethod2();
	int NextPointMethod3();
	void SetFirstPoint();
	bool CheckPolyFitOK();
	int TryAddPoint1();
	int TryAddPoint3();
	int ForceBreak();
	void StartNextPoly();
	void FindPolyChangeLoc();




};

#endif // !defined(AFX_POLYSEGFITPOINTS_H__6651B2A0_DA1B_11D5_86C3_A499BC9BC324__INCLUDED_)
