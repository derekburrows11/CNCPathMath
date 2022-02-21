// PolySegDblFit.h: interface for the CPolySegDblFit class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYSEGDBLFIT_H__1B3BC120_5E04_11D7_86C3_D7B9E0CF9C26__INCLUDED_)
#define AFX_POLYSEGDBLFIT_H__1B3BC120_5E04_11D7_86C3_D7B9E0CF9C26__INCLUDED_


#include "PolySegFit.h"


#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000




// CPolySegDblFit fits double cubic polys to given points
class CPolySegDblFit : public CPolySegFit
{
public:
	enum		// Segment Fit Flags
	{
		SFF_POSMID  = 0x100,
		SFF_DIRMID  = 0x200,
	};

	CPolySegDblFit();
	virtual ~CPolySegDblFit();

	void Reset();
	void SetInitialStateFromJoin();
	void SetPolyJoinPointRef(double pointRef);
	void RemoveFitAPoints();

	void FitSingleSegment() { CPolySegFit::SolveCubic(); }
	void FitDoubleSegment() { SolveDoubleCubic(); }

	double GetSumResidual() { return m_SumResidual; }

	CVector& GetBaseBezierNode(int iNodeRef) { return m_FitA.m_Bez[iNodeRef]; }
	CVector& GetLeaderBezierNode(int iNodeRef) { return m_FitB.m_Bez[iNodeRef]; }

	double GetBaseMinControlSpan() { return m_FitA.GetMinControlSpan(); }


protected:
	void SolveDoubleCubic();

	// functions for double cubic fitting
	void SetDblCubicPoints();
	void SetSPolyJoin();
	void AdjustSPolyJoin();

	// overidden functions
	void SetSValues();
	void AdjustSValues();
	void SetSums();
	void GetBeziers();
	void GetResiduals();
	void SetAxisFitToConstValue(int ax);

	// using bezier equations
	void FitDblCubicBez();
	void FitDblCubicIPosBez();
	void FitDblCubicIPosDerivBez(double Kv0);
	void FitDblCubicTwoPosBez();
	void FitDblCubicTwoPosIDerivBez(double Kv0);
	void FitDblCubicTwoPosFDerivBez(double Kv1);
	void FitDblCubicTwoPosTwoDerivBez(double Kv0, double Kv1);

	void SetCommonMatricies();			// common matricies for Bezier type fits
	void Solve8by8();

	// using poly equations
	void FitDblCubicIPos();
	void FitDblCubicIPosDeriv(double Kv0);
	void FitDblCubicTwoPos();
	void FitDblCubicTwoPosIDeriv(double Kv0);						// not complete
	void FitDblCubicTwoPosTwoDeriv(double Kv0, double Kv1);	// not complete


// Data
protected:

	// poly fit data for double cubic fitting
	CPolySegFit m_FitA;
	CPolySegFit m_FitB;		// for double cubic fitting

	// common matricies
	CMatrix m_mxSab;
	CMatrix m_mxab2Bez;
	CMatrix m_mxBez2ab;
	CMatrix m_mxSBez;
	CMatrix m_vtSaby;

	// common matricies dependent on BC's
	int m_numDOFs;
	int m_numBCs;
	CMatrix m_mxBCs;
	CMatrix m_mxBCsBez2Bez;
	CMatrix m_mxBCsBez2ab;
	CMatrix m_mxdab;
	CMatrix m_mxdabT;

	// common matricies for 8x8 solve
	CMatrix m_mxBCVals;
	CMatrix m_vtBez;
	CMatrix m_mxAllEqus;
	CMatrix m_mxAllEqus2;
	CMatrix m_vtAllVals;
	CMatrix m_vtAllVals2;
	char* m_arRowArrange;

	// common matricies for reduced DOF solve
	// variable size
	CMatrix m_mxReduEqus;
	CMatrix m_vtReduVals;
	CMatrix m_vtValDiff;
	CMatrix m_vtBCsBez2HR;
	CMatrix m_mxdabS;
	// const size
	CMatrix m_vtab;			// vector(8)


	// common matricies for checking
	CMatrix m_bez88, m_bez44;
	CMatrix m_bez44b;			// just for poly solves


	// data for second poly m_FitB
	double m_pointRefPolyJoin;		// normalised node pos of poly B start
};

#endif // !defined(AFX_POLYSEGDBLFIT_H__1B3BC120_5E04_11D7_86C3_D7B9E0CF9C26__INCLUDED_)
