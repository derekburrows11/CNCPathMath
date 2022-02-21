// PolySegFit.cpp: implementation of the CPolySegFit class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "Matrix.h"
#include "PolyFunc.h"

#include "PolySegFit.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPolySegFit::CPolySegFit()
{
	m_arPts = NULL;
	m_arS = NULL;
	m_barPtsFromCArray = true;
	m_barSFromCArray = true;
	m_iPtBufferSize = 0;;
	m_iSBufferSize = 0;
	Reset();
}

CPolySegFit::~CPolySegFit()
{
}



//////////////////////////////////////////////////////////////////
// Setup functions
//////////////////////////////////////////////////////////////////


void CPolySegFit::SetPointBuffer(CVector* pBuffer, int sizeBuffer)
{
	m_arPts = pBuffer;
	m_iPtBufferSize = sizeBuffer;
	m_barPtsFromCArray = false;
}

void CPolySegFit::SetNumPointsInBuffer(int numPts)
{
	ASSERT(m_barPtsFromCArray == false);
	ASSERT(numPts <= m_iPtBufferSize);
	m_numPtsStored = numPts;
	if (m_bUseAllPoints)
		m_numPts = m_numPtsStored;
}

void CPolySegFit::SetSBuffer(double* pBuffer, int sizeBuffer)
{
	m_arS = pBuffer;
	m_iSBufferSize = sizeBuffer;
	m_barSFromCArray = false;
}

void CPolySegFit::Reset()
{
	m_numPts = 0;
	m_numPtsStored = 0;
	m_PtArray.SetSize(0, 100);
	m_bUseAllPoints = true;
	m_bPreFit = false;
	m_SInit = 0;
	m_SFinal = 0;
	m_nFitPosInit = FIT_FREE;
	m_nFitDirInit = FIT_FREE;
	m_nFitPosFinal = FIT_FREE;
	m_nFitDirFinal = FIT_FREE;

}

void CPolySegFit::SetNumPointsToUse(int numPtsToUse)
{
	m_bUseAllPoints = false;
	m_numPts = numPtsToUse;
}

void CPolySegFit::UseAllPoints()
{
	m_bUseAllPoints = true;
	m_numPts = m_numPtsStored;
}


bool CPolySegFit::AddPoint(const CVector& pt)
{
	if (m_barPtsFromCArray)
	{
		m_PtArray.Add(const_cast<CVector&>(pt));
		ASSERT(++m_numPtsStored == m_PtArray.GetSize());
		m_arPts = m_PtArray.GetData();
		m_numPtsStored = m_PtArray.GetSize();
	}
	else
	{
		ASSERT(m_arPts != NULL);
		ASSERT(m_numPtsStored < m_iPtBufferSize);
		if (m_numPtsStored >= m_iPtBufferSize)
			return false;
		m_arPts[m_numPtsStored++] = pt;
	}
	if (m_bUseAllPoints)
		m_numPts = m_numPtsStored;
	return true;
}

void CPolySegFit::RemoveFirstPoints(int iNum)
{
	if (m_barPtsFromCArray)
	{
		m_PtArray.RemoveAt(0, iNum);
		m_arPts = m_PtArray.GetData();		// shouldn't change!
		m_numPtsStored = m_PtArray.GetSize();
	}
	else
		ASSERT(0);
	if (m_bUseAllPoints)
		m_numPts = m_numPtsStored;
}

void CPolySegFit::SetToInitialPoint()
{
	m_nFitPosInit = FIT_ATPOINT;
}
void CPolySegFit::SetToFinalPoint()
{
	m_nFitPosFinal = FIT_ATPOINT;
}
void CPolySegFit::SetInitialPoint(const CVector& vtPosI, double SInit)
{
	ASSERT(SInit <= 0);		// SInit is SInit value if S first point is 0
	m_vtPosInit = vtPosI;
	m_SInit = SInit;
	m_nFitPosInit = FIT_SET;
}
void CPolySegFit::SetFinalPoint(const CVector& vtPosF, double SFinal)
{
	ASSERT(SFinal >= 0);
	m_vtPosFinal = vtPosF;
	m_SFinal = SFinal;
	m_nFitPosFinal = FIT_SET;
}
void CPolySegFit::SetInitialDir(const CVector& vtDirI)
{
	vtDirI.Unit(m_vtDirInitUnit);
	m_nFitDirInit = FIT_SET;
}
void CPolySegFit::SetFinalDir(const CVector& vtDirF)
{
	vtDirF.Unit(m_vtDirFinalUnit);
	m_nFitDirFinal = FIT_SET;
}


//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

void CPolySegFit::SetEndConditions()
{
	// set init/final points now if applicable
	if (m_nFitPosInit == FIT_ATPOINT)
		m_vtPosInit = m_arPts[0];
	if (m_nFitPosFinal == FIT_ATPOINT)
		m_vtPosFinal = m_arPts[m_numPts-1];
}

void CPolySegFit::SolveCubic()
{

	SetEndConditions();
	SetSValues();
	NormaliseSValues();
	double Kv0 = m_SSpanOrig;		// if normalising S Values
	double Kv1 = m_SSpanOrig;

	int nSegFitFlags = 0;
	if (m_nFitPosInit  != FIT_FREE) nSegFitFlags |= SFF_POSI;
	if (m_nFitPosFinal != FIT_FREE) nSegFitFlags |= SFF_POSF;
	if (m_nFitDirInit  != FIT_FREE) nSegFitFlags |= SFF_DIRI;
	if (m_nFitDirFinal != FIT_FREE) nSegFitFlags |= SFF_DIRF;

	// check dof's
	m_numPtsForFit = m_numPts;
	if (m_nFitPosInit == FIT_SET && m_SInit < m_arS[0])			// if initial point is not at S of first fit point
		m_numPtsForFit++;
	if (m_nFitPosFinal == FIT_SET && m_SFinal > m_arS[m_numPts-1])	// if final point is not at S of last fit point
		m_numPtsForFit++;


	int count = 0;
	for (;;)
	{
		SetSums();						// set SumSp, SumSpV

		m_bPreFit = true;			// don't report >=4 segs
		switch (nSegFitFlags)
		{
		case SFF_POSI | SFF_POSF:
			FitCubicTwoPos0to1();
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRI:
			FitCubicTwoPos0to1();		// just used to get approx end derivatives
			Kv0 = m_Poly[1].Mag();
			FitCubicTwoPosIDeriv0to1(Kv0);
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRF:
			FitCubicTwoPos0to1();		// just used to get approx end derivatives
			Kv1 = (m_Poly[3]*3 + m_Poly[2]*2 + m_Poly[1]).Mag();	// if s=1 at end
			FitCubicTwoPosFDeriv0to1(Kv1);
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRI | SFF_DIRF:
//			FitCubicTwoPosTwoDir0to1();
			FitCubicTwoPos0to1();		// just used to get approx end derivatives
			Kv0 = m_Poly[1].Mag();
			Kv1 = (m_Poly[3]*3 + m_Poly[2]*2 + m_Poly[1]).Mag();	// if s=1 at end
			FitCubicTwoPosTwoDeriv0to1(Kv0, Kv1);
			break;
		default:
			ASSERT(0);			
		}
			

		AdjustSValues();			// if both ends fitted s will still range 0->1

		count++;
		if (m_dsAvgNorm <= 1e-4)
			break;
		if (count >= 2)
			if (m_dsAvgNorm <= 1e-3)
				break;
			else if (count >= 8)
				if (m_dsAvgNorm <= 1e-2)
					break;
				else if (count > 10)
				{
					ASSERT(0);
					break;
				}
	}

	if (m_MaxResidual >= 0.3)
		TRACE1("Large Maximum Residual in CPolySegFit::SolveCubic() of: %g\n", m_MaxResidual);


	GetBezierFromPoly();

}


//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////


void CPolySegFit::SetSValues()
{
	ASSERT(m_numPts <= m_numPtsStored);		// check if not using all points
	if (m_barSFromCArray)
	{
		m_SArray.SetSize(m_numPts);
		m_arS = m_SArray.GetData();
	}
	else
	{
		ASSERT(m_numPts <= m_iSBufferSize);
	}

//	Increase in s value between points is equal to distance between points
// S of first point is always 0 -> arS[0]==0
// if m_nFitPosInit == FIT_SET then SInit must have been set <=0 such that arS[0]==0
// if m_nFitPosInit == FIT_FREE or FIT_ATPOINT SInit will be 0
	double s = 0;
	if (m_nFitPosInit != FIT_SET)
		m_SInit = s;
	else
		ASSERT(m_SInit <= s);		// should have been set to value <= S of first point
	if (m_numPts > 0)
	{
		m_arS[0] = s;
		for (int i = 1; i < m_numPts; i++)
		{
			s += (m_arPts[i] - m_arPts[i-1]).Mag();
			m_arS[i] = s;
		}
	}
	if (m_nFitPosFinal != FIT_SET)
		m_SFinal = s;
	else
		ASSERT(m_SFinal >= s);		// should have been set to value >= S of last point
	m_SSpanOrig = m_SFinal - m_SInit;
}

void CPolySegFit::NormaliseSValues()
{
	// scale so s varies from 0 to 1
	if (m_SInit == 0)
	{
		double fact = 1 / m_SFinal;
		for (int i = 0; i < m_numPts; i++)
			m_arS[i] *= fact;
	}
	else
	{
		double base = m_SInit;
		double fact = 1 / (m_SFinal - base);
		for (int i = 0; i < m_numPts; i++)
			m_arS[i] =  (m_arS[i] - base) * fact;
	}
	m_SInit = 0;
	m_SFinal = 1;
}

void CPolySegFit::AdjustSValues()
{
	CVector vtPos, vtVel, vtAcc;
	CVector vtToPoint;
	double sumResSq = 0;
	double dsSum = 0;
	double magSqToPointMax = 0;
	for (int i = 0; i < m_numPts; i++)
	{
		double s = m_arS[i];
		double sInit = s;
		double ds;
		double magSqToPoint;
		CVector vtPt = m_arPts[i];
		for (;;)
		{
			vtPos = CubicAt(m_Poly, s);
			vtToPoint = vtPt - vtPos;
			magSqToPoint = vtToPoint.MagSq();
			if (magSqToPoint <= 1e-14)		// no change needed
				break;
			//	find ds
			vtVel = CubicD1At(m_Poly, s);
			vtAcc = CubicD2At(m_Poly, s);
			double dsChord, dsCirc, dsTang1, dsTang2;
			double velMagSq = vtVel.MagSq();
			double velMag = sqrt(velMagSq);
			CVector vtCurve = cross(vtVel, vtAcc) / (velMagSq * velMag);
			double curveMagSq = vtCurve.MagSq();

			if (curveMagSq < 1e-6)		// 1000mm radius
			{
				dsTang1 = dot(vtToPoint, vtVel) / velMagSq;			// use straight tangent method
				ds = dsTang1;
			}
			else
			{
				double radius = 1 / sqrt(curveMagSq);

				CVector vtAnormDir = cross(vtCurve, vtVel);
				CVector vtAnormUnit = vtAnormDir / vtAnormDir.Mag();

				CVector vtRadius = vtAnormUnit * -radius;	// from centre to pos
				CVector vtCentre = vtPos - vtRadius;
				CVector vtCentre2Point = vtPt - vtCentre;

			// better tangent method
				CVector vtCentre2Tang = vtCentre2Point * (radius*radius / dot(vtCentre2Point, vtRadius));
				CVector vtTang = vtCentre2Tang - vtRadius;
				dsTang2 = vtTang.Mag() / velMag;
				if (dot(vtVel, vtTang) < 0)
					dsTang2 = -dsTang2;

			// chord method
				CVector vtCentre2Chord = vtCentre2Point * (radius / vtCentre2Point.Mag());
				CVector vtChord = vtCentre2Chord - vtRadius;
				dsChord = vtChord.Mag() / velMag;	// using chord distance
				if (dot(vtVel, vtChord) < 0)
					dsChord = -dsChord;

			// circumference method - doesn't give +/- direction
			//	double sinAngCross = (cross(vtRadius, vtCentre2Point)).Mag() / (radius * vtCentre2Point.Mag());	// doesn't give +/- direction
			// circumference method using chord
				double sinAng = dot(vtCentre2Chord, vtVel) / (velMag * radius);	// gives +/- direction
				double ang = asin(sinAng);
				if (fabs(ang) > 0.1)		// reduce large moves on small radius
					ang *= 0.25;
				double circDist = radius * ang;
				dsCirc = circDist / velMag;

				ds = dsCirc;
			}
			s += ds;
			if (fabs(ds) < 1e-8)
				break;
		}	// for (;;)
		sumResSq += magSqToPoint;
//		double magToPoint = sqrt(magSqToPoint);	// can remove - just for checking
		if (magSqToPoint > magSqToPointMax)
			magSqToPointMax = magSqToPoint;
		ds = s - sInit;
		m_arS[i] = s;
		dsSum += fabs(ds);
	}	// for (int i = 0; i < m_numPts; i++)
	double sSpan = m_arS[m_numPts-1] - m_arS[0];
	ASSERT(sSpan > 0);
	m_dsAvgNorm = dsSum / (m_numPts * sSpan);		// ds averaged over numPts and normalised for total s length
	m_MaxResidual = sqrt(magSqToPointMax);
	m_SumResidualSq = sumResSq;						// sum residual after s points moved
	m_SumResidual = sqrt(sumResSq);
	m_SumResidualAvg = sqrt(sumResSq / (m_numPts-1));

//	m_SInit =					would change only if FIT_FREE but then not used!
//	m_SFinal =					would change only if FIT_FREE but then not used!
}

void CPolySegFit::SetSums()
{
	// set m_PointSums
	CPointSPowerData sums;
	sums.Zero();
	for (int i = 0; i < m_numPts; i++)
		sums.SumUpPoint(m_arS[i], m_arPts[i]);
	m_PointSums = sums;
}

void CPolySegFit::GetBezierFromPoly()
{
	// get Bezier poly (poly span not 0-1 so can't use Cubic2Bezier(m_Bez, m_Poly)
	double sInit = m_SInit;
	double sFinal = m_SFinal;
	double sSpanOn3 = (sFinal - sInit) / 3;
	m_Bez[0] = CubicAt(m_Poly, sInit);
	m_Bez[3] = CubicAt(m_Poly, sFinal);
	m_Bez[1] = m_Bez[0] + CubicD1At(m_Poly, sInit) * sSpanOn3;
	m_Bez[2] = m_Bez[3] - CubicD1At(m_Poly, sFinal) * sSpanOn3;
}

void CPolySegFit::GetBezDerivFromPoly()		// control nodes are relative derivatives
{
	// get Bezier poly (poly span not 0-1 so can't use Cubic2Bezier(m_Bez, m_Poly)
	double sInit = m_SInit;
	double sFinal = m_SFinal;
	m_BezDeriv[0] = CubicAt(m_Poly, sInit);
	m_BezDeriv[3] = CubicAt(m_Poly, sFinal);
	m_BezDeriv[1] = CubicD1At(m_Poly, sInit);
	m_BezDeriv[2] = CubicD1At(m_Poly, sFinal);
	GetBezierFromBezDeriv();
}

void CPolySegFit::GetBezDerivFromBezier()		// control nodes are relative derivatives
{
	double sSpanOn3 = (m_SFinal - m_SInit) / 3;
	m_BezDeriv[0] = m_Bez[0];
	m_BezDeriv[3] = m_Bez[3];
	m_BezDeriv[1] = (m_Bez[1] - m_Bez[0]) / sSpanOn3;
	m_BezDeriv[2] = (m_Bez[3] - m_Bez[2]) / sSpanOn3;
}

void CPolySegFit::GetBezierFromBezDeriv()		// control nodes are relative derivatives
{
	double sSpanOn3 = (m_SFinal - m_SInit) / 3;
	m_Bez[0] = m_BezDeriv[0];
	m_Bez[3] = m_BezDeriv[3];
	m_Bez[1] = m_BezDeriv[0] + m_BezDeriv[1] * sSpanOn3;
	m_Bez[2] = m_BezDeriv[3] - m_BezDeriv[2] * sSpanOn3;
}

double CPolySegFit::GetMinControlSpan()
{
	m_magEnds = (m_Bez[3] - m_Bez[0]).Mag();
	m_magControl[0] = (m_Bez[1] - m_Bez[0]).Mag();
	m_magControl[1] = (m_Bez[3] - m_Bez[2]).Mag();
	double minContRatio = m_magControl[0] < m_magControl[1] ? m_magControl[0] : m_magControl[1];
	minContRatio /= m_magEnds;
	return minContRatio;
}



void CPolySegFit::GetResiduals()
{
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

//	Sr = S(v^2) + (vta' * mxS * vta) - 2*(vta' * vtSv)
//	if (mxS * vta == vtSv)		(only if no other boundary conditions)
//	->  Sr = S(v^2) - (vta' * vtSv)

	CMatrix mxS(4,4);
	CMatrix vtA(4), vtSpV(4);
	CMatrix vtRes1(1,4), res2(1);
	CVector vtSumResidualSq2;

	// set matrix 'mxS' - creates both symmetrical diagonal halves
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
			mxS.elem(r,c) = SumSp[r+c];

	for (int ax = 0; ax < 3; ax++)
	{
		// set vector 'vtA'
		// set vector 'vtSpV'
		for (int i = 0; i < 4; i++)
		{
			vtA[i] = m_Poly[i][ax];
			vtSpV[i] = SumSpV[i][ax];
		}
		vtRes1.Prod(mxS, vtA);
		vtRes1 -= vtSpV * 2;
		res2.ProdTl(vtA, vtRes1);
		m_vtSumResidualSq[ax] = m_PointSums.V2[ax] + res2[0];

		// or Sr = S(v^2) - (vta' * vtSv) method
		// res2.ProdTl(vtA, vtSpV);
		// vtSumResidualSq2[ax] = m_PointSums.V2[ax] - res2[0];
	}
	m_SumResidualSq = m_vtSumResidualSq.Sum();
	m_SumResidual = sqrt(m_SumResidualSq);
}




//////////////////////////////////////////////////////////
// Fitting Functions
//////////////////////////////////////////////////////////


bool CPolySegFit::CheckNumPoints(int numPtsReq, int iLine)
{
	int numPtsExcess = m_numPtsForFit - numPtsReq;
	if (numPtsExcess <= (m_bPreFit ? 0 : 1))
		TRACE2("Low number of excess points(%i) for fit in CPolySegFit line %i\n", numPtsExcess, iLine);
	m_bPreFit = false;
	ASSERT(numPtsExcess >= 0);
	return numPtsExcess >= 0;
}


void CPolySegFit::FitCubic()
{
// Fit parametric cubic to points, no boundary conditions
// Each axis is solved independently

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

/*
		[ S(s^0)  S(s^1)  S(s^2)  S(s^3) ]   a0   S(s^0.v)
		[ S(s^1)  S(s^2)  S(s^3)  S(s^4) ] * a1 = S(s^1.v)
		[ S(s^2)  S(s^3)  S(s^4)  S(s^5) ]   a2   S(s^2.v)
		[ S(s^3)  S(s^4)  S(s^5)  S(s^6) ]   a3   S(s^3.v)
*/

//---- Full solve method ----
	double arA[4];
	double arVal[4];

	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		for (int i = 0; i < 4; i++)
			arVal[i] = SumSpV[i][ax];
		VERIFY(LUFullSymSolve(4, SumSp, arA, arVal));	// solves matrix equ, result in arA
		for (i = 0; i < 4; i++)
			m_Poly[i][ax] = arA[i];
	}

//---- Inverse method ----
	CMatrix mxS(4,4);
	CMatrix vtA(4), vtVal(4);
	// set matrix 'mxS'
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
			mxS.elem(r,c) = SumSp[r+c];
	VERIFY(mxS.Invert() > 1e-8);

	for (ax  = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		for (int i = 0; i < 4; i++)
			vtVal[i] = SumSpV[i][ax];
		vtA.Prod(mxS, vtVal);
		for (i = 0; i < 4; i++)
			m_Poly[i][ax] = vtA[i];
	}
}

void CPolySegFit::FitCubicIPos()
{
// Fit parametric cubic with initial pos at s = 0
// Each axis is solved independently

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	CVector vtPosInit = m_vtPosInit;
/*
		If initial pos and vel are at x=0 then equ's are simplified -> a0 & a1 known
		[ S(s^2)  S(s^3)  S(s^4) ]   a1     S(s^1.y)   [ S(s^1) ]
		[ S(s^3)  S(s^4)  S(s^5) ] * a2  =  S(s^2.y) - [ S(s^2) ] * a0
		[ S(s^4)  S(s^5)  S(s^6) ]   a3     S(s^3.y)   [ S(s^3) ]
*/

//---- Full solve method ----
	double arA[4];
	double arVal[4];
	double a0;

	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		a0 = vtPosInit[ax];
		for (int i = 1; i < 4; i++)
			arVal[i] = SumSpV[i][ax] - SumSp[i] * a0;
		VERIFY(LUFullSymSolve(3, SumSp+2, arA+1, arVal+1));		// solves matrix equ, result in arA

		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = arA[1];
		m_Poly[2][ax] = arA[2];
		m_Poly[3][ax] = arA[3];
	}

//---- Inverse method ----
	CMatrix mxS(3,3);
	CMatrix vtA(3), vtVal(3);
	// set matrix 'mxS'
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			mxS.elem(r,c) = SumSp[2+r+c];
	VERIFY(mxS.Invert() > 1e-8);

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		a0 = vtPosInit[ax];
		for (int i = 1; i < 4; i++)
			vtVal[i-1] = SumSpV[i][ax] - SumSp[i] * a0;
		vtA.Prod(mxS, vtVal);
		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = vtA[0];
		m_Poly[2][ax] = vtA[1];
		m_Poly[3][ax] = vtA[2];
	}
}

void CPolySegFit::FitCubicIPosDir()
{
// Fit parametric cubic with initial pos & direction at s = 0
// Axes are solved together (larger matrix solution)

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	CVector vtPosInit = m_vtPosInit;
	CVector vtDirInit = m_vtDirInitUnit;



/*
	vta = [ k ax2 ax3 ay2 ay3 az2 az3 ]'

	      [ (Vx^2+Vy^2+Vz^2).S(s^2)  Vx.S(s^3)  Vx.S(s^4)  Vy.S(s^3)  Vy.S(s^4)  Vz.S(s^3)  Vz.S(s^4) ]
	      [  Vx.S(s^3)                  S(s^4)     S(s^5)     0          0          0          0      ]
	      [  Vx.S(s^4)                  S(s^5)     S(s^6)     0          0          0          0      ]
	mxS = [  Vy.S(s^3)                  0          0          S(s^4)     S(s^5)     0          0      ]
	      [  Vy.S(s^4)                  0          0          S(s^5)     S(s^6)     0          0      ]
	      [  Vz.S(s^3)                  0          0          0          0          S(s^4)     S(s^5) ]
	      [  Vz.S(s^4)                  0          0          0          0          S(s^5)     S(s^6) ]


	      [ Vx.S(s^1.x)+Vy.S(s^1.y)+Vz.S(s^1.z) - (Vx.x0+Vy.y0+Vz.z0).S(s^1)  ]
	      [ S(s^2.x)       x0.S(s^2)  ]
	      [ S(s^3.x)       x0.S(s^3)  ]
	vtS = [ S(s^2.y)   -   y0.S(s^2)  ]
	      [ S(s^3.y)       y0.S(s^3)  ]
	      [ S(s^2.z)       z0.S(s^2)  ]
	      [ S(s^3.z)       z0.S(s^3)  ]


----------


*/

	// set matrix 'mxS' - creates both symmetrical diagonal halves
	CMatrix mxS(7,7);
	mxS = 0;
	mxS.elem(0,0) = vtDirInit.MagSq() * SumSp[2];
	for (int ax = 0; ax < 3; ax++)
	{
		int ax2 = 2 * ax;
		mxS.elem(0, ax2+1) = mxS.elem(ax2+1,0) = vtDirInit[ax] * SumSp[3];
		mxS.elem(0, ax2+2) = mxS.elem(ax2+2,0) = vtDirInit[ax] * SumSp[4];

		mxS.elem(ax2+1,ax2+1) = SumSp[4];
		mxS.elem(ax2+2,ax2+2) = SumSp[6];
		mxS.elem(ax2+1,ax2+2) = mxS.elem(ax2+2,ax2+1) = SumSp[5];
	}


	// set vector 'vtVals'
	CMatrix vtVals(7);
	double vtVal0 = 0;
	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax])		// constant value, don't need to solve!
		{
			vtVals[2*ax+1] = vtVals[2*ax+2] = 0;
			continue;
		}

		vtVals[2*ax+1] = SumSpV[2][ax] - vtPosInit[ax] * SumSp[2];
		vtVals[2*ax+2] = SumSpV[3][ax] - vtPosInit[ax] * SumSp[3];
		vtVal0 += vtDirInit[ax] * (SumSpV[1][ax] - vtPosInit[ax] * SumSp[1]);
	}
	vtVals[0] = vtVal0;

	// solve
	CMatrix vta(7);
	mxS.LUSolve(vta.GetArray(), vtVals.GetArray());

	for (ax = 0; ax < 3; ax++)
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtDirInit[ax] == 0)		// constant value, don't need to solve!
			vta[2*ax+1] = vta[2*ax+2] = 0;		// eliminate arithmetic error

	// set poly coeffs
	m_Poly[0] = vtPosInit;
	m_Poly[1] = vtDirInit * vta[0];		// vtDirInit * k
	m_Poly[2].Set(vta[1], vta[3], vta[5]);
	m_Poly[3].Set(vta[2], vta[4], vta[6]);

}



void CPolySegFit::FitCubicIPosDeriv(double Kv0)
{
// Fit parametric cubic with initial pos & derivative = Kv * m_vtDirInitUnit at s = 0
// Each axis is solved independently

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(3, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	CVector vtPosInit = m_vtPosInit;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;

/*
		If initial pos and vel are at x=0 then equ's are simplified -> a0 & a1 known
		[ S(s^4)  S(s^5) ] * a2  =  S(s^2.v) - [ S(s^2)  S(s^3) ] * a0
		[ S(s^5)  S(s^6) ]   a3     S(s^3.v)   [ S(s^3)  S(s^4) ]   a1
*/
	// get determinent of 2x2 matrix
	double a0, a1;
	double det, b0, b1;
	det = SumSp[4]*SumSp[6] - SumSp[5]*SumSp[5];
	if (det == 0)
		ASSERT(0);
	det = 1 / det;

	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtDerivInit[ax] == 0)		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		a0 = vtPosInit[ax];
		a1 = vtDerivInit[ax];
		b0 = SumSpV[2][ax] - SumSp[2]*a0 - SumSp[3]*a1;
		b1 = SumSpV[3][ax] - SumSp[3]*a0 - SumSp[4]*a1;
		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = a1;
		m_Poly[2][ax] = det *  ( SumSp[6]*b0 - SumSp[5]*b1);
		m_Poly[3][ax] = det *  (-SumSp[5]*b0 + SumSp[4]*b1);
	}
}


void CPolySegFit::FitCubicTwoPos()
{
// Fit parametric cubic with initial and final pos at sa & sb
// Each axis is solved independently

/*
	[ 1        sa       sa^2     sa^3 ]   a0   Pos(sa)
	[ 1        sb       sb^2     sb^3 ] * a1 = Pos(sb)
	[ 0        1        0        0    ]   a2   A1
	[ 0        0        1        0    ]   a3   A2
	mxa2A * vta = vtA

Invert:
	[ 0                   sa-sb      sa^2-sb^2          sa^3-sb^3 ] = Pos(sa)-Pos(sb)					// r1 - r2
	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)        0    ] = sb^3*Pos(sa)-sa^3*Pos(sb)		// sb^3*r1 - sa^3*r2

	[ sb^3  -sa^3  sb*sa(sa^2-sb^2)  sb^2*sa^2(sa-sb) ]   Pos(sa)   a0*(sb^3-sa^3)
	[ 0        0        1                 0           ] * Pos(sb) = a1
	[ 0        0        0                 1           ]   A1        a2
	[-1        1     sa-sb           sa^2-sb^2        ]   A2        a3*(sb^3-sa^3)
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[sb*sa(sa^2-sb^2)/(sb^3-sa^3)   1   0   (sa-sb)/(sb^3-sa^3)     ]
	[sb^2*sa^2(sa-sb)/(sb^3-sa^3)   0   1   (sa^2-sb^2)/(sb^3-sa^3) ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

*/
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);		// ensure enough data points for a solution, numPtsReq would give exact fit!


	double sa = m_SInit;
	double sb = m_SFinal;
	double sa2 = sa*sa,  sb2 = sb*sb;
	double sa3 = sa*sa2, sb3 = sb*sb2;
	CMatrix mxda(2,4), mxA2aKC(4,2), mxS(4,4);

	// set matrix 'mxda'
	double sb3msa3Inv = 1 / (sb3-sa3);
	double e03 = (sa-sb) * sb3msa3Inv;
	double e13 = (sa+sb) * e03;
	mxda.elem(0,0) = sa*sb*e13;
	mxda.elem(1,0) = sa2*sb2*e03;
	mxda.elem(0,1) = 1;
	mxda.elem(1,1) = 0;
	mxda.elem(0,2) = 0;
	mxda.elem(1,2) = 1;
	mxda.elem(0,3) = e03;
	mxda.elem(1,3) = e13;

	// set matrix 'mxA2aKC'
	mxA2aKC.elem(0,0) = sb3 * sb3msa3Inv;
	mxA2aKC.elem(1,0) = 0;
	mxA2aKC.elem(2,0) = 0;
	mxA2aKC.elem(3,0) = -sb3msa3Inv;
	mxA2aKC.elem(0,1) = -sa3 * sb3msa3Inv;
	mxA2aKC.elem(1,1) = 0;
	mxA2aKC.elem(2,1) = 0;
	mxA2aKC.elem(3,1) = sb3msa3Inv;

	// set matrix 'mxS'
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
			mxS.elem(r,c) = SumSp[r+c];


	CMatrix mxdaS(2,4), mxEqu(2,2);
	mxdaS.Prod(mxda, mxS);
	mxEqu.ProdTr(mxdaS, mxda);
	VERIFY(mxEqu.Invert() > 1e-8);

//		(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
//	or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')

	CMatrix mxSA2aKC(4,2);
	mxSA2aKC.Prod(mxS, mxA2aKC);

	CMatrix vtSv(4), vtAKR(2), vtVal(2), vtAUR(2), vta(4);
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == m_vtPosInit[ax] && m_vtPosInit[ax] == m_vtPosFinal[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		// set vector 'vtSv'
		for (r = 0; r < 4; r++)
			vtSv[r] = SumSpV[r][ax];

		// set vector 'vtAKR'
		vtAKR[0] = m_vtPosInit[ax];
		vtAKR[1] = m_vtPosFinal[ax];

		vtVal.Prod(mxda, (vtSv - mxSA2aKC * vtAKR));
		vtAUR.Prod(mxEqu, vtVal);

		vta.ProdTl(mxda, vtAUR);
		vta += mxA2aKC * vtAKR;

		m_Poly[0][ax] = vta[0];
		m_Poly[1][ax] = vtAUR[0];
		m_Poly[2][ax] = vtAUR[1];
		m_Poly[3][ax] = vta[3];
	}
}

void CPolySegFit::FitCubicTwoPos0to1()
{
// Fit parametric cubic with initial and final pos at s = 0 & 1
// Each axis is solved independently

/*
	[ 1        sa       sa^2     sa^3 ]   a0   Pos(sa)
	[ 1        sb       sb^2     sb^3 ] * a1 = Pos(sb)
	[ 0        1        0        0    ]   a2   A1
	[ 0        0        1        0    ]   a3   A2
	mxa2A * vta = vtA

Invert:
	[ 0                   sa-sb      sa^2-sb^2          sa^3-sb^3 ] = Pos(sa)-Pos(sb)					// r1 - r2
	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)        0    ] = sb^3*Pos(sa)-sa^3*Pos(sb)		// sb^3*r1 - sa^3*r2

	[ sb^3  -sa^3  sb*sa(sa^2-sb^2)  sb^2*sa^2(sa-sb) ]   Pos(sa)   a0*(sb^3-sa^3)
	[ 0        0        1                 0           ] * Pos(sb) = a1
	[ 0        0        0                 1           ]   A1        a2
	[-1        1     sa-sb           sa^2-sb^2        ]   A2        a3*(sb^3-sa^3)
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[sb*sa(sa^2-sb^2)/(sb^3-sa^3)   1   0   (sa-sb)/(sb^3-sa^3)     ]
	[sb^2*sa^2(sa-sb)/(sb^3-sa^3)   0   1   (sa^2-sb^2)/(sb^3-sa^3) ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

	---------------------------
	if sa = 0 & sb = 1
	mxA2a =
	[ 1   0   0   0 ]
	[ 0   0   1   0 ]
	[ 0   0   0   1 ]
	[-1   1  -1  -1 ]
	mxda = mxA2a(uc)' =
	[ 0  1  0 -1 ]
	[ 0  0  1 -1 ]

	mxda*mxS =
	[ S1-S3  S2-S4  S3-S5  S4-S6 ]
	[ S2-S3  S3-S4  S4-S5  S5-S6 ]  

  mxda*mxS*mxda' =
	[ S2-2*S4+S6   S3-S4-S5+S6 ]
	[ S3-S4-S5+S6   S4-2*S5+S6 ]

	mxda * vtSv =
	[ S1v-S3v ]
	[ S2v-S3v ]

	(mxda*mxS) * mxA2a(kc) =
	[ S1-S3-S4+S6    S4-S6 ]
	[ S2-S3-S5+S6    S5-S6 ]

	mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]' =
	[ S1v-S3v - (S1-S3-S4+S6)Pos(sa) - (S4-S6)Pos(sb) ]
	[ S2v-S3v - (S2-S3-S5+S6)Pos(sa) - (S5-S6)Pos(sb) ]
*/
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!


	CMatrix mxEqu(2,2), mxdaSA2aKC(2,2);
	// set matrix 'mxEqu'
	double e01 = SumSp[3] - SumSp[4] - SumSp[5] + SumSp[6];
	mxEqu.elem(0,0) = SumSp[2] - 2*SumSp[4] + SumSp[6];
	mxEqu.elem(0,1) = e01;
	mxEqu.elem(1,0) = e01;
	mxEqu.elem(1,1) = SumSp[4] - 2*SumSp[5] + SumSp[6];
	VERIFY(mxEqu.Invert() > 1e-8);

	// set matrix 'mxdaSA2aKC'
	mxdaSA2aKC.elem(0,0) = SumSp[1] - SumSp[3] - SumSp[4] + SumSp[6];
	mxdaSA2aKC.elem(0,1) = SumSp[4] - SumSp[6];
	mxdaSA2aKC.elem(1,0) = SumSp[2] - SumSp[3] - SumSp[5] + SumSp[6];
	mxdaSA2aKC.elem(1,1) = SumSp[5] - SumSp[6];

	CMatrix vtAKR(2), vtVal(2), vtAUR(2);
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == m_vtPosInit[ax] && m_vtPosInit[ax] == m_vtPosFinal[ax])		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		// set vector 'vtAKR'
		vtAKR[0] = m_vtPosInit[ax];
		vtAKR[1] = m_vtPosFinal[ax];

		// set vector 'vtVal'
		vtVal.Prod(mxdaSA2aKC, vtAKR);
		vtVal[0] = SumSpV[1][ax] - SumSpV[3][ax] - vtVal[0];
		vtVal[1] = SumSpV[2][ax] - SumSpV[3][ax] - vtVal[1];
		vtAUR.Prod(mxEqu, vtVal);

		m_Poly[0][ax] = vtAKR[0];
		m_Poly[1][ax] = vtAUR[0];
		m_Poly[2][ax] = vtAUR[1];
		m_Poly[3][ax] = vtAKR[1] - vtAKR[0] - vtAUR[0] - vtAUR[1];
	}
}



void CPolySegFit::FitCubicTwoPosTwoDeriv0to1(double Kv0, double Kv1)
{
// Fit parametric cubic with initial and final pos & derivative at s = 0 & 1
//	deriv(0) = Kv0 * m_vtDirInitUnit, deriv(1) = Kv1 * m_vtDirFinalUnit
// Each axis is solved independently

	CheckNumPoints(2, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	// cubic is exactly defined
/*
	a2 = (pf-df/3 - pi-2di/3) * 3
	a2 = 3pf-df-3pi-2di
	a3 = bez[3] - bez[0] + (bez[1] - bez[2]) * 3;
	a3 = pf - pi + (pi+di/3-pf+df/3) * 3;
	a3 = 2pi - 2pf + di + df;
	a3 = pf - a2 - a1 - a0
*/
	m_Poly[0] = m_vtPosInit;
	m_Poly[1] = m_vtDirInitUnit * Kv0;
	m_Poly[3] = (m_vtPosInit - m_vtPosFinal) * 2 + m_Poly[1] + m_vtDirFinalUnit * Kv1;
	m_Poly[2] = m_vtPosFinal - m_Poly[3] - m_Poly[1] - m_Poly[0];
//	m_Bez[0] = m_vtPosInit;
//	m_Bez[1] = m_vtPosInit + m_vtDirInitUnit * (Kv0 / 3);
//	m_Bez[2] = m_vtPosFinal - m_vtDirFinalUnit * (Kv1 / 3);
//	m_Bez[3] = m_vtPosFinal;
//	Bezier2Cubic(m_Poly, m_Bez);
}

void CPolySegFit::FitCubicTwoPosTwoDir0to1()
{
// Fit parametric cubic with initial and final pos & direction at s = 0 & 1
//	m_vtDirInitUnit, m_vtDirFinalUnit
// Axes are solved together; Kv0 & Kv1 are the only dof's so a 2x2 inverse

/*
	3 axis:
	vta = [ax0 ax1 ax2 ax3 ay0 ay1 ay2 ay3 az0 az1 az2 az3]'
	vtA = [x0 x1 y0 y1 z0 z1 k0 k1]'

	mxA2a =
	[ 1    0    0    0    0    0    0    0   ]
	[ 0    0    0    0    0    0   Vx0   0   ]
	[-3    3    0    0    0    0 -2Vx0 -Vx1  ]
	[ 2   -2    0    0    0    0   Vx0  Vx1  ]
	[ 0    0    1    0    0    0    0    0   ]
	[ 0    0    0    0    0    0   Vy0   0   ]
	[ 0    0   -3    3    0    0 -2Vy0 -Vy1  ]
	[ 0    0    2   -2    0    0   Vy0  Vy1  ]
	[ 0    0    0    0    1    0    0    0   ]
	[ 0    0    0    0    0    0   Vz0   0   ]
	[ 0    0    0    0   -3    3 -2Vz0 -Vz1  ]
	[ 0    0    0    0    2   -2   Vz0  Vz1  ]

	mxda = transpose of unknown columns of mxA2a
	[ 0  Vx0 -2Vx0  Vx0   0  Vy0 -2Vy0  Vy0   0  Vz0 -2Vz0  Vz0 ]
	[ 0    0  -Vx1  Vx1   0    0  -Vy1  Vy1   0    0  -Vz1  Vz1 ]

	mxda * mxS * vta  =  mxda * vtSv
	(mxda*mxS) * mxda' * vtA(ur)  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * vtA(kr)

	mxda * mxS =
	[ Vx0(S1-2S2+S3)  Vx0(S2-2S3+S4)  Vx0(S3-2S4+S5)  Vx0(S4-2S5+S6)      Vy0(S1-2S2+S3)  ..... ]
	[ Vx1(S3-S2)      Vx1(S4-S3)      Vx1(S5-S4)      Vx1(S6-S5)          Vy1(S3-S2)      ..... ]

	mxda*mxS * mxda' =
	[ (S2-4S3+6S4-4S5+S6)*(Vx0^2+Vy0^2+Vz0^2)    (S6-3S5+3S4-S3)*(Vx0*Vx1+Vy0*Vy1+Vz0*Vz1) ]
	[ (S6-3S5+3S4-S3)*(Vx0*Vx1+Vy0*Vy1+Vz0*Vz1)            (S4-2S5+S6)*(Vx1^2+Vy1^2+Vz1^2) ]

	mxda * vtSv =
	[ Vx0(S1x-2S2x+S3x) + Vy0(S1y-2S2y+S3y) + Vz0(S1z-2S2z+S3z) ]
	[ Vx1(S3x-S2x) + Vy1(S3y-S2y) + Vz1(S3z-S2z)                ]

	(mxda*mxS) * mxA2a(kc) =
	[ Vx0(S1-2S2-2S3+8S4-7S5+2S6)  Vx0(3S3-8S4+7S5-2S6)  Vy0(S1-2S2-2S3+8S4-7S5+2S6)  Vy0(3S3-8S4+7S5-2S6)  Vz0(S1-2S2-2S3+8S4-7S5+2S6)  Vz0(3S3-8S4+7S5-2S6) ]
	[ Vx1(2S6-5S5+3S4+S3-S2)       Vx1(-2S6+5S5-3S4)     Vy1(2S6-5S5+3S4+S3-S2)       Vy1(-2S6+5S5-3S4)     Vz1(2S6-5S5+3S4+S3-S2)       Vz1(-2S6+5S5-3S4)    ]

*/
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(4, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDirInit = m_vtDirInitUnit;
	CVector vtDirFinal = m_vtDirFinalUnit;

	CMatrix mxEqu(2,2), mxdaSA2aKC(2,2);
	// set matrix 'mxEqu'
	double e01 = (SumSp[6] - 3*SumSp[5] + 3*SumSp[4] - SumSp[3]) * dot(vtDirInit, vtDirFinal);
	mxEqu.elem(0,0) = (SumSp[2] - 4*SumSp[3] + 6*SumSp[4] - 4*SumSp[5] + SumSp[6]) * vtDirInit.MagSq();
	mxEqu.elem(0,1) = e01;
	mxEqu.elem(1,0) = e01;
	mxEqu.elem(1,1) = (SumSp[4] - 2*SumSp[5] + SumSp[6]) * vtDirFinal.MagSq();
	VERIFY(mxEqu.Invert() > 1e-8);

	// set matrix 'mxdaSA2aKC'
	e01 = 3*SumSp[3] - 8*SumSp[4] + 7*SumSp[5] - 2*SumSp[6];
	double e11 = -3*SumSp[4] + 5*SumSp[5] - 2*SumSp[6];
	mxdaSA2aKC.elem(0,0) = SumSp[1] - 2*SumSp[2] + SumSp[3] - e01;	// * Vxyz0
	mxdaSA2aKC.elem(0,1) = e01;												// * Vxyz0
	mxdaSA2aKC.elem(1,0) = SumSp[3] - SumSp[2] - e11;					// * Vxyz1
	mxdaSA2aKC.elem(1,1) = e11;												// * Vxyz1

	CMatrix vtAKR(2), vtValKR(2), vtVal(2), vtAUR(2);
	///////////////////////////
	// vtVal is calculated for each axis
	vtVal = 0;
	for (int ax = 0; ax < 3; ax++)
	{
		// set vector 'vtAKR'
		vtAKR[0] = vtPosInit[ax];
		vtAKR[1] = vtPosFinal[ax];

		// set vector 'vtVal'
		vtValKR.Prod(mxdaSA2aKC, vtAKR);
		vtVal[0] += (SumSpV[3][ax] - 2*SumSpV[2][ax] + SumSpV[1][ax] - vtValKR[0]) * vtDirInit[ax];
		vtVal[1] += (SumSpV[3][ax] - SumSpV[2][ax] - vtValKR[1]) * vtDirFinal[ax];
	}
	vtAUR.Prod(mxEqu, vtVal);		// vtAUR = [k0 k1]'
	if (vtAUR[0] <= 0 || vtAUR[1] <= 0)
		TRACE("Dir coeff <= 0 in CPolySegFit::FitCubicTwoPosTwoDir0to1()");

	m_Poly[0] = vtPosInit;
	m_Poly[1] = vtDirInit * vtAUR[0];
	m_Poly[3] = (vtPosInit - vtPosFinal)*2 + m_Poly[1] + vtDirFinal * vtAUR[1];
	m_Poly[2] = vtPosFinal - m_Poly[3] - m_Poly[1] - m_Poly[0];
}

void CPolySegFit::FitCubicTwoPosIDeriv0to1(double Kv0)
{
// Fit parametric cubic with initial and final pos & initial derivative at s = 0 & 1
//	deriv(0) = Kv0 * m_vtDirInitUnit
// Each axis is solved independently

/*
	[ 1        sa       sa^2     sa^3 ]   a0   Pos(sa)
	[ 0        1       2sa      3sa^2 ] * a1 = Der(sa)
	[ 1        sb       sb^2     sb^3 ]   a2   Pos(sb)
	[ 0        0        1        0    ]   a3   A2
	mxa2A * vta = vtA

Invert: not done yet
	[ 0                   sa-sb      sa^2-sb^2          sa^3-sb^3 ] = Pos(sa)-Pos(sb)					// r1 - r2
	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)        0    ] = sb^3*Pos(sa)-sa^3*Pos(sb)		// sb^3*r1 - sa^3*r2

	[ sb^3  -sa^3  sb*sa(sa^2-sb^2)  sb^2*sa^2(sa-sb) ]   Pos(sa)   a0*(sb^3-sa^3) = a0(sb-sa)(sb^2+sa.sb+sb^2)
	[ 0        0        1                 0           ] * Pos(sb) = a1
	[ 0        0        0                 1           ]   A1        a2
	[-1        1     sa-sb           sa^2-sb^2        ]   A2        a3*(sb^3-sa^3)
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[sb^2*sa^2(sa-sb)/(sb^3-sa^3)   0   1   (sa^2-sb^2)/(sb^3-sa^3)   ]
  =[-sb^2*sa^2/(sa^2+sa.sb+sb^2)   0   1  -(sa+sb)/(sa^2+sa.sb+sb^2) ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

	---------------------------
	if sa = 0 & sb = 1
	mxa2A =
	[ 1   0   0   0 ]
	[ 0   1   0   0 ]
	[ 1   1   1   1 ]
	[ 0   0   1   0 ]
	mxA2a =
	[ 1   0   0   0 ]
	[ 0   1   0   0 ]
	[ 0   0   0   1 ]
	[-1  -1   1  -1 ]
	mxda = mxA2a(uc)' =
	[ 0   0   1  -1 ]

	mxda*mxS =
	[ S2-S3  S3-S4  S4-S5  S5-S6 ]

  mxda*mxS*mxda' =
	[ S4-2*S5+S6 ]

	mxda * vtSv =
	[ S2v-S3v ]

	(mxda*mxS) * mxA2a(kc) =
	[ S2-S3-S5+S6   S3-S4-S5+S6   S5-S6 ]

	mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
*/

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(3, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	double equInv = SumSp[4] - 2*SumSp[5] + SumSp[6];
	equInv = 1 / equInv;

	// set matrix 'mxdaSA2aKC'
	double mxdaSA2aKC0 = SumSp[2] - SumSp[3] - SumSp[5] + SumSp[6];
	double mxdaSA2aKC1 = SumSp[3] - SumSp[4] - SumSp[5] + SumSp[6];
	double mxdaSA2aKC2 = SumSp[5] - SumSp[6];

	double vtAKR0, vtAKR1, vtAKR2, val, a2;
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == m_vtPosInit[ax] && m_vtPosInit[ax] == m_vtPosFinal[ax]
			&& m_vtDirInitUnit[ax] == 0)		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		// set vector 'vtAKR'
		vtAKR0 = m_vtPosInit[ax];
		vtAKR1 = m_vtDirInitUnit[ax] * Kv0;
		vtAKR2 = m_vtPosFinal[ax];

		// set vector 'vtVal'
		val = mxdaSA2aKC0 * vtAKR0 + mxdaSA2aKC1 * vtAKR1 + mxdaSA2aKC2 * vtAKR2;
		val = SumSpV[2][ax] - SumSpV[3][ax] - val;
		a2 = equInv * val;

		m_Poly[0][ax] = vtAKR0;
		m_Poly[1][ax] = vtAKR1;
		m_Poly[2][ax] = a2;
		m_Poly[3][ax] = vtAKR2 - vtAKR0 - vtAKR1 - a2;
	}
}




void CPolySegFit::FitCubicTwoPosFDeriv0to1(double Kv1)
{
// Fit parametric cubic with initial and final pos & final derivative at s = 0 & 1
//	deriv(1) = Kv1 * m_vtDirFinalUnit
// Each axis is solved independently

/*
	[ 1        sa       sa^2     sa^3 ]   a0   Pos(sa)
	[ 0        1       2sb      3sb^2 ] * a1 = Der(sb)
	[ 1        sb       sb^2     sb^3 ]   a2   Pos(sb)
	[ 0        1        0        0    ]   a3   A1
	mxa2A * vta = vtA

Invert:
	[ 0                   sa-sb      sa^2-sb^2          sa^3-sb^3 ] = r1 - r3
*	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)        0    ] = sb^3*r1 - sa^3*r3

	[ 0       2sb(sa-sb)-sa^2+sb^2      0      2sb(sa^3-sb^3)-3sb^2(sa^2-sb^2)                 ] = 2sb(r1-r3) - (sa^2-sb^2)r2
	[ 0       -sb^2+2sa.sb-sa^2         0      2sb(sa-sb)(sa^2+sa.sb+sb^2)-3sb^2(sa-sb)(sa+sb) ] = 2sb(r1-r3) - (sa^2-sb^2)r2
	[ 0           -(sb-sa)^2            0      sb(sa-sb)(2sa^2-sa.sb-sb^2)                     ] = 2sb(r1-r3) - (sa^2-sb^2)r2
	[ 0           -(sb-sa)^2            0      sb(2sa+sb)(sb-sa)^2                             ] = 2sb(r1-r3) - (sa^2-sb^2)r2

	[ 0    3sb^2(sa-sb)-sa^3+sb^3                 3sb^2(sa^2-sb^2)-2sb(sa^3-sb^3)   0 ] = 3sb^2(r1-r3) - (sa^3-sb^3)r2
	[ 0   -3sb^2(sb-sa)+(sb-sa)(sa^2+sa.sb+sb^2)  -sb(2sa+sb)(sb-sa)^2              0 ] = 3sb^2(r1-r3) - (sa^3-sb^3)r2
	[ 0    (sb-sa)(sa^2+sa.sb-2sb^2)              -sb(2sa+sb)(sb-sa)^2              0 ] = 3sb^2(r1-r3) - (sa^3-sb^3)r2
*	[ 0   -(sa+2sb)(sb-sa)^2                      -sb(2sa+sb)(sb-sa)^2              0 ] = 3sb^2(r1-r3) - (sa^3-sb^3)r2

//	[ 2sb   2sb.sa-sa^2    0    2sb.sa^3-3sa^2.sb^2  ] = 2sb*r1 - sa^2*r2
//	[ 2sb   sa(2sb-sa)     0    sa^2.sb(2sa-3sb)     ] = 2sb*r1 - sa^2*r2

	[ (sb-sa)(2sa+sb)(sb^3-sa^3)   sb.sa(sb-sa)(2sa+sb)(sb^2-sa^2) - sb.sa^2(sa+2sb)(sb-sa)^2   0     0  ] = (sb-sa)(2sa+sb).(sb^3*r1 - sa^3*r3) + (sb*sa^2).(3sb^2(r1-r3) - (sa^3-sb^3)r2)
	[ (sb-sa)^2.(2sa+sb)(sb^2+sa.sb+sa^2)     sb.sa(sb-sa)^2( (2sa+sb)(sa+sb) - sa(sa+2sb))     0     0  ] = (sb-sa)(2sa+sb).(sb^3*r1 - sa^3*r3) + sb.sa^2.3sb^2(r1-r3) - sb.sa^2.(sa^3-sb^3)r2
	[ (2sa+sb)(sb-sa)^2.(sa^2+sa.sb+sb^2)     sb.sa(sb-sa)^2.(sa^2+sa.sb+sb^2)                  0     0  ] = (sb-sa)(2sa+sb)sb^3*r1 - (sb-sa)(2sa+sb)sa^3*r3 + sb.sa^2.3sb^2.r1-sb.sa^2.3sb^2.r3 - sb.sa^2.(sa^3-sb^3)r2
	[ (2sa+sb)(sb-sa)^2.(sa^2+sa.sb+sb^2)     sb.sa(sb-sa)^2.(sa^2+sa.sb+sb^2)                  0     0  ] = ((sb-sa)(2sa+sb) + 3sa^2)sb^3.r1 - ((sb-sa)(2sa+sb)sa + 3sb^3)sa^2.r3 - sb.sa^2.(sa^3-sb^3)r2
	[ (2sa+sb)(sb-sa)^2.(sa^2+sa.sb+sb^2)     sb.sa(sb-sa)^2.(sa^2+sa.sb+sb^2)                  0     0  ] = (sa^2+sa.sb+sb^2)sb^3.r1 - (-2sa^3+sa^2.sb+sa.sb^2 + 3sb^3)sa^2.r3 - sb.sa^2.(sa^3-sb^3)r2
	[ (2sa+sb)(sb-sa)^2.(sa^2+sa.sb+sb^2)     sb.sa(sb-sa)^2.(sa^2+sa.sb+sb^2)                  0     0  ] = (sa^2+sa.sb+sb^2)sb^3.r1 - (sa^2+sa.sb+sb^2)(3sb-2sa)sa^2.r3 - sb.sa^2.(sa^2+sa.sb+sb^2)(sa-sb)r2
	[ (2sa+sb)(sb-sa)^2                       sb.sa(sb-sa)^2                                    0     0  ] = sb^3.r1 - (3sb-2sa)sa^2.r3 - sb.sa^2.(sa-sb)r2


	[ sb^3    sb.sa^2.(sb-sa)   (2sa-3sb)sa^2      -sa.sb(sb-sa)^2  ]   Pos(sa)   a0*(2sa+sb)(sb-sa)^2
	[ 0            0              0                         1       ] * Der(sb) = a1
	[-3sb^2   sa^3-sb^3          3sb^2          -(sa+2sb)(sb-sa)^2  ]   Pos(sb)   a2*sb(2sa+sb)(sb-sa)^2
	[ 2sb     sb^2-sa^2         -2sb                     (sb-sa)^2  ]   A1        a3*sb(2sa+sb)(sb-sa)^2
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[-sa.sb/(2sa+sb)   1   -(sa+2sb)/((2sa+sb)sb)   1/((2sa+sb)sb)  ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

	---------------------------
	if sa = 0 & sb = 1
	mxa2A =
	[ 1   0   0   0 ]
	[ 0   1   2   3 ]
	[ 1   1   1   1 ]
	[ 0   1   0   0 ]

	[-2  -1   0   1 ]		r2-2*r3
	[ 3   2   1   0 ]		3r3-r2
	 
	mxA2a =
	[ 1   0   0   0 ]
	[ 0   0   0   1 ]
	[-3  -1   3  -2 ]
	[ 2   1  -2   1 ]
	mxda = mxA2a(uc)' =
	[ 0   1  -2   1 ]

	mxda*mxS =
	[ S1-2*S2+S3  S2-2*S3+S4  S3-2*S4+S5  S4-2*S5+S6 ]

  mxda*mxS*mxda' =
	[ S2-4*S3+6*S4-4*S5+S6 ]

	mxda * vtSv =
	[ S1v-2*S2v+S3v ]

	(mxda*mxS) * mxA2a(kc) =
	[ S1-2*S2-2*S3+8*S4-7*S5+2*S6  -S3+3*S4-3*S5+S6   3*S3-8*S4+7*S5-2*S6 ]

	mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
*/

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CheckNumPoints(3, __LINE__);			// ensure enough data points for a solution, numPtsReq would give exact fit!

	double equInv = SumSp[2] - 4*SumSp[3] + 6*SumSp[4] - 4*SumSp[5] + SumSp[6];
	equInv = 1 / equInv;

	// set matrix 'mxdaSA2aKC'
	double mxdaSA2aKC1 =  -SumSp[3] + 3*SumSp[4] - 3*SumSp[5] +   SumSp[6];
	double mxdaSA2aKC2 = 3*SumSp[3] - 8*SumSp[4] + 7*SumSp[5] - 2*SumSp[6];
	double mxdaSA2aKC0 = SumSp[1] - 2*SumSp[2] + SumSp[3] - mxdaSA2aKC2;

	double vtAKR0, vtAKR1, vtAKR2, val, a1;
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == m_vtPosInit[ax] && m_vtPosInit[ax] == m_vtPosFinal[ax]
			&& m_vtDirFinalUnit[ax] == 0)		// constant value, don't need to solve!
		{
			m_Poly[0][ax] = m_PointSums.m_Value[ax];
			m_Poly[1][ax] = m_Poly[2][ax] = m_Poly[3][ax] = 0;
			continue;
		}

		// set vector 'vtAKR'
		vtAKR0 = m_vtPosInit[ax];
		vtAKR1 = m_vtDirFinalUnit[ax] * Kv1;
		vtAKR2 = m_vtPosFinal[ax];

		// set vector 'vtVal'
		val = mxdaSA2aKC0 * vtAKR0 + mxdaSA2aKC1 * vtAKR1 + mxdaSA2aKC2 * vtAKR2;
		val = SumSpV[1][ax] - 2*SumSpV[2][ax] + SumSpV[3][ax] - val;
		a1 = equInv * val;

		m_Poly[0][ax] = vtAKR0;
		m_Poly[1][ax] = a1;
		m_Poly[2][ax] =-3*vtAKR0 - vtAKR1 + 3*vtAKR2 - 2*a1;
		m_Poly[3][ax] = 2*vtAKR0 + vtAKR1 - 2*vtAKR2 + a1;
	}
}



