// PolySegInfo.cpp: implementation of the CPolySegInfo class.
//
//////////////////////////////////////////////////////////////////////


#include "stdafx.h"


#include "Matrix.h"
#include "PolyFunc.h"
#include "..\Graph\Plot.h"		// only for CVect2 for testing ???

#include "PolySegInfo.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif









//////////////////////////////////
// CPolySegInfo functions
//////////////////////////////////

CPolySegInfo::CPolySegInfo()
{
	Reset();
}

CPolySegInfo& CPolySegInfo::operator=(const CPolySegInfo& orig)
{
	for (int i = 0; i < 3; i++)
	{
		m_coef[i] = orig.m_coef[i];
		m_yValInit[i] = orig.m_yValInit[i];
		m_yValFinal[i] = orig.m_yValFinal[i];
	}
	m_coef[3] = orig.m_coef[3];		// one more of coef!
	m_xPolyInit = orig.m_xPolyInit;
	m_xPolyFinal = orig.m_xPolyFinal;
	m_xPointInit = orig.m_xPointInit;
	m_xPointFinal = orig.m_xPointFinal;
//	GetNumPoints() = orig.GetNumPoints();
//	m_PointSums = orig.m_PointSums;		// a class with a list
	return *this;
}

void CPolySegInfo::Reset()
{
	m_PointSums.ClearAll();
	m_PrevPointSumState.m_iNumPoints = -1;

	m_Sr = m_Sabs = 0;
	m_bPolySet = false;
//	m_nSumRefInit = m_nSumRefFinal = 0;
}

void CPolySegInfo::ClearAllButFirst()
{
	// Clears all but the first point
	CPointData pd = GetFirstPoint();
	Reset();
	AddPoint(pd);
}

void CPolySegInfo::SetInitValsFromFinal()
{
	m_xPolyInit = m_xPolyFinal;
	m_yValInit[0] = CubicAt(m_coef, m_xPolyFinal);
	m_yValInit[1] = CubicD1At(m_coef, m_xPolyFinal);
}

void CPolySegInfo::CalcInitPV()
{
	// xPolyInit & coef's must be set
	m_yValInit[0] = CubicAt(m_coef, m_xPolyInit);
	m_yValInit[1] = CubicD1At(m_coef, m_xPolyInit);
}
void CPolySegInfo::CalcFinalPV()
{
	// xPolyFinal & coef's must be set
	m_yValFinal[0] = CubicAt(m_coef, m_xPolyFinal);
	m_yValFinal[1] = CubicD1At(m_coef, m_xPolyFinal);
}

void CPolySegInfo::AddPoint(CPointData& pt)
{
	// This function adds the current point to the CPointSumData member
	if (!m_PointSums.AddPoint(pt))
	{
//		ASSERT(0);		// step not right!
		return;
	}

	if (GetNumPoints() == 1)
		m_xPointInit = m_xPointFinal = pt.x;
	m_Sr += pt.yErrorSq;
	m_Sabs += fabs(pt.yError);
	if (pt.x > m_xPointFinal)
		m_xPointFinal = pt.x;
}

void CPolySegInfo::AddPointUnconfirmed(CPointData& pt)
{
	// This function adds the current point to the CPointSumData member
	// backup of current state is done so removal is possible

	m_PointSums.SaveState(m_PrevPointSumState);
	if (!m_PointSums.AddPoint(pt))
		ASSERT(0);

	if (GetNumPoints() == 0)
		m_xPointInit = m_xPointFinal = pt.x;
	//	store in case of removing point
	m_SrPrev = m_Sr;
	m_SabsPrev = m_Sabs;
	m_xPointFinalPrev = m_xPointFinal;
	//
	m_Sr += pt.yErrorSq;
	m_Sabs += fabs(pt.yError);
	if (pt.x > m_xPointFinal)
		m_xPointFinal = pt.x;
}

void CPolySegInfo::AppendPoly(CPolySegInfo& psi)
{
	// This function appends psi
	POSITION pos = psi.m_PointSums.GetFirstSumsPos();
	while (pos)
		AddPoint(psi.m_PointSums.GetNextSums(pos));


}

void CPolySegInfo::RemoveLatestPoint()
{
	if ((m_PrevPointSumState.m_iNumPoints == GetNumPoints() - 1)
		&& (m_PointSums.RetrieveState(m_PrevPointSumState)))
		{
			m_Sr = m_SrPrev;
			m_Sabs = m_SabsPrev;
			m_xPointFinal = m_xPointFinalPrev;
		}
	else
		ASSERT(0);
}

CPointData& CPolySegInfo::GetFirstPoint()
{
	return m_PointSums.GetFirstSums();
}

CPointData& CPolySegInfo::GetFinalPoint()
{
	return m_PointSums.GetLastSums();
}

void CPolySegInfo::SetPolyFrom(CPolySegInfo& src)
{
	for (int i = 0; i < sizeof(m_coef)/sizeof(*m_coef); i++)
		m_coef[i] = src.m_coef[i];
	m_bPolySet = src.m_bPolySet;
}

void CPolySegInfo::CalcPointError(CPointData& pt)
{
	pt.yPoly = CubicAt(m_coef, pt.x);
	pt.yError = pt.yPoly - pt.y;
	pt.yErrorSq = pt.yError * pt.yError;
}


double CPolySegInfo::GetSumResidual()
{
/*	Sr = Sum((y - Y)^2)		where y is actual value and Y is cubic value
	Sr = S(y^2) + a0^2.S(x^0) + a1^2.S(x^2) + a2^2.S(x^4) + a3^2.S(x^6)
		+ 2{-a0.S(x^0.y) - a1.S(x^1.y) - a2.S(x^2.y) - a3.S(x^3.y) }
		+ 2{ a0.a1.S(x^1) + a0.a2.S(x^2) + (a0.a3+a1.a2).S(x^3) + a1.a3.S(x^4) + a2.a3.S(x^5) }
*/
	const CPointPowerData& sumPoints = m_PointSums.GetLastSums();		// based on poly range
	const double* SumXp = sumPoints.Xp;
	const double* SumYXp = sumPoints.YXp;
	double a0 = m_coef[0], a1 = m_coef[1], a2 = m_coef[2], a3 = m_coef[3];		// faster??
	double Sr;
	Sr = a0*a1*SumXp[1] + a0*a2*SumXp[2] + (a0*a3+a1*a2)*SumXp[3] + a1*a3*SumXp[4] + a2*a3*SumXp[5];
	Sr -= a0*SumYXp[0] + a1*SumYXp[1] + a2*SumYXp[2] + a3*SumYXp[3];
	Sr *= 2;
	Sr += a0*a0*SumXp[0] + a1*a1*SumXp[2] + a2*a2*SumXp[4] + a3*a3*SumXp[6];
	Sr += sumPoints.Y2;
	return Sr;
	// Standard Error is sqrt(Variance) = sqrt(Sr / (n-4))   (4, or num d.o.f. in fit)
}

double CPolySegInfo::GetStandardError(double Sr)	// Standard Error is sqrt(Variance) = sqrt(Sr / (n-4))   (4, or num d.o.f. in fit)
{
	return sqrt(Sr / (GetNumPoints() - m_iFitDOF));
}

double CPolySegInfo::GetPolyErrors()
{
	CPointPowerDataList& spList = m_PointSums.m_SumPointsList;
	double SumRes = 0;
	POSITION pos = spList.GetTailPosition();
	while (pos)
	{
		CPointPowerData& pd = spList.GetPrev(pos);
		pd.yPoly = CubicAt(m_coef, pd.x);
		pd.yError = pd.yPoly - pd.y;
		pd.yErrorSq = pd.yError * pd.yError;
		SumRes += pd.yErrorSq;
	}
//	arPD[i].yErrorSq = SumRes;		// store it in next element!
	return SumRes;
}

int CPolySegInfo::SetPointArray(CVect2* arPts)
{
	CPointPowerDataList& spList = m_PointSums.m_SumPointsList;
	POSITION pos = spList.GetTailPosition();
	int numPt = 0;
	while (pos)
	{
		CPointPowerData& pd = spList.GetPrev(pos);
		arPts[numPt].x = pd.x;
		arPts[numPt].y = pd.y;
		numPt++;
	}
	return numPt;
}

int CPolySegInfo::SetBezierPoints(CVect2* arPts)		// sets 4 points
{
	double spanOn3 = (m_xPolyFinal - m_xPolyInit) / 3;
	arPts[0].x = m_xPolyInit;
	arPts[0].y = CubicAt(m_coef, m_xPolyInit);
	arPts[1].x = m_xPolyInit + spanOn3;
	arPts[1].y = arPts[0].y + CubicD1At(m_coef, m_xPolyInit) * spanOn3;

	arPts[3].x = m_xPolyFinal;
	arPts[3].y = CubicAt(m_coef, m_xPolyFinal);
	arPts[2].x = m_xPolyFinal - spanOn3;
	arPts[2].y = arPts[3].y - CubicD1At(m_coef, m_xPolyFinal) * spanOn3;
	return 4;		// 4 points set
}

void CPolySegInfo::PlotPoints()
{
	// plot points
	int numPts = GetNumPoints();
	CVect2* arPts = new CVect2[numPts];

	SetPointArray(arPts);
	::Plot(arPts, numPts, GPS_POINT);
	delete[] arPts;
}

void CPolySegInfo::PlotBezier()
{
	// plot bezier poly
	CVect2 arPts[4];
	SetBezierPoints(arPts);
	::Plot(arPts, 4, GPS_BEZIER);
}


// Fit poly to current points - not continuous from previous poly
void CPolySegInfo::FitPolyToPoints()
{
	const CPointPowerData& sumPoints = m_PointSums.GetLastSums();
	// fill matrix for poly soltion
/*	int mxIdx = 0;
	for (int j = 0; j < 4; j++)			// row
	{
		int xpIdx = j;
		for (int k = 0; k < 4; k++)		// column
			matrix[mxIdx++] = SumXp[xpIdx++];
	}
*/
	// LUFullSymSolve uses a vector instead of matrix for a fully sym matrix
	LUFullSymSolve(4, sumPoints.Xp, m_coef, sumPoints.YXp);		// solves matrix equ, result in m_arFitPoly
	m_iFitDOF = 4;
}

// Fit poly to current points with Pos and Vel continuous from previous poly
void CPolySegInfo::FitPolyToPointsContPV()
{
	const CPointPowerData& SumPts = m_PointSums.GetLastSums();
	const double (&SumXp)[7] = SumPts.Xp;
	const double (&SumYXp)[4] = SumPts.YXp;

	double det, b[2];
	double Xa = m_xPolyInit;		// initial x

	if (Xa == 0)		// Xa is x location of initial pos & vel
	{
/*
		If initial pos and vel are at x=0 then equ's are simplified -> a0 & a1 known
		[ S(x^4)  S(x^5) ] * a2  =  S(x^2.y) - [ S(x^2)  S(x^3) ] * a0
		[ S(x^5)  S(x^6) ]   a3     S(x^3.y)   [ S(x^3)  S(x^4) ]   a1
*/
		// get determinent of 2x2 matrix
		det = SumXp[4]*SumXp[6] - SumXp[5]*SumXp[5];
		if (det == 0)
			ASSERT(0);
		det = 1 / det;

		m_coef[0] = m_yValInit[0];	// derivative[0] is just value!
		m_coef[1] = m_yValInit[1];
		// m_arFitPoly[0,1] are already set
		b[0] = SumYXp[2] - SumXp[2]*m_coef[0] - SumXp[3]*m_coef[1];
		b[1] = SumYXp[3] - SumXp[3]*m_coef[0] - SumXp[4]*m_coef[1];
		m_coef[2] = det * ( SumXp[6]*b[0] - SumXp[5]*b[1]);
		m_coef[3] = det * (-SumXp[5]*b[0] + SumXp[4]*b[1]);
//		return;
	}
double c0[4], c22[4], c44[4];
for (int i = 0; i<4; i++) c0[i] = m_coef[i];

/*
	Using a two equation solution with know pos and vel

	Let the 4x4 defined below be called mxSBCs
	mxSBC = [ mxda * mxS * vta ]
	        [      mxBCs * vta ]

	mxSBCs * vta = vtSyBCs			from below

	coef's for shifted equ are vtb (instead of vta) where x is 0 at start

	[ 1   -x0    x0^2  -x0^3 ]   b0   a0
	[ 0    1   -2x0    3x0^2 ] * b1 = a1
	[ 0    0     1    -3x0   ]   b2   a2
	[ 0    0     0      1    ]   b3   a3
		mxBtoA * vtb = vta

	[ 1    x0    x0^2   x0^3 ]   a0   b0 = pos init
	[ 0    1    2x0    3x0^2 ] * a1 = b1 = vel init
	[ 0    0     1     3x0   ]   a2   b2
	[ 0    0     0      1    ]   a3   b3
		mxAtoB * vta = vtb

	->  mxSBCs * mxBtoA * vtb = vtSyBCs
	As b0 & b1 are know (pos & vel at x0) b2 & b3 can be solved by rearranging as a 2x2
	Then use:
		vta = mxBtoA * vtb

OR:
	         [ 1    x0    x0^2   x0^3 ]
	mxAtoB = [ 0    1    2x0    3x0^2 ]
	         [ 0    0     1      0    ]
	         [ 0    0     0      1    ]
		and
	         [ 1   -x0    x0^2  2x0^3 ]
	mxBtoA = [ 0    1   -2x0   -3x0^2 ]
	         [ 0    0     1      0    ]
	         [ 0    0     0      1    ]

*/

/*
	If initial x = xa
	[ 1        xa       xa^2     xa^3 ]   a0    Pos(xa)
	[ 0        1       2xa      3xa^2 ] * a1 =  Vel(xa)
	[ 0        0        1        0    ]   a2    A2
	[ 0        0        0        1    ]   a3    A3
	mxa2A * vta  =  vtA

	eliminate a0, a1 (solve for these!)
	step 1:
	[ 1        0       -xa^2   -2xa^3 ]   a0    Pos(xa) - xa*Vel(xa)		//  -xa * row2
	[ 0        1       2xa      3xa^2 ] * a1 =  Vel(xa)
	[ 0        0        1        0    ]   a2    A2
	[ 0        0        0        1    ]   a3    A3

	[ 1       -xa       xa^2    2xa^3 ]   Pos(xa)   a0
	[ 0        1      -2xa     -3xa^2 ] * Vel(xa) = a1
	[ 0        0        1        0    ]   A2        a2
	[ 0        0        0        1    ]   A3        a3
	mxA2a * vtA = vta



	dSr/da2 = dSra/da2 + dSra/da0 * da0/da2 + dSra/da1 * da1/da2 = 0		// from all dependent params
	-->
	dSr/da2   [ da0/da2  da1/da2    1    0 ]   dSra/da0     0
	dSr/da3 = [ da0/da3  da1/da3    0    1 ] * dSra/da1  =  0
	                                           dSra/da2
	                                           dSra/da3

	-> vtdSrRed = mxda * vtdSra = 0

	-> mxda = [  xa^2   -2xa     1    0 ]
	          [ 2xa^3   -3xa^2   0    1 ]
	     
	  mxda * vtdSr = 0
	& vtdSra = mxS * vta  -  vtSy

	mxda * mxS * vta  =  mxda * vtSy		// 2 equ's		(note: mxda doesn't have an inverse - not square!)
	     mxBCs * vta  =  vtBCs				// 2 equ's
	  ( mxSBCs * vta  =  vtSyBCs )

	(mxda * mxS * mxA2a) * vtA  =  mxda * vtSy		// 2 equ's for the 2 unknowns in vtA
	use columns of (mxda * mxS * mxA2a) associated with unknowns for inverse
	use columns of (mxda * mxS * mxA2a) associated with knowns to multipy and subtract from rhs

	(mxda*mxS) * mxA2a(uc) * vtA(ur)  =  mxda * vtSy - (mxda*mxS) * mxA2a(kc) * vtA(kr)
	uc = unknown columns
	kc = known columns

	note:  mxda = mxA2a(uc)'
	(mxda*mxS) * mxda' * vtA(ur)  =  mxda * vtSy - (mxda*mxS) * mxA2a(kc) * vtA(kr)
	mxda * mxS * mxda' * vtA(ur)  =  mxda * (vtSy - mxS * mxA2a(kc) * vtA(kr))

*/


	CMatrix mxda(2,4), mxS(4,4), mxBCs(2,4);
	CMatrix vtSy(4), vtBCs(2);

	double XaP;
	// set matrix 'mxda'
	mxda.elem(0, 1) = -2*Xa;
	XaP = Xa*Xa;
	mxda.elem(0, 0) = XaP;
	mxda.elem(1, 1) = -3*XaP;
	mxda.elem(1, 0) = 2*Xa*XaP;
	mxda.elem(0, 2) = 1;
	mxda.elem(1, 3) = 1;
	mxda.elem(0, 3) = 0;
	mxda.elem(1, 2) = 0;

	// set matrix 'mxS'
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
			mxS.elem(r, c) = SumXp[r+c];

	// set vector 'vtSy'
	for (r = 0; r < 4; r++)
		vtSy[r] = SumYXp[r];

	// set matrix 'mxBCs'
	XaP = 1;
	for (int c = 0; ; )
	{
		mxBCs.elem(0, c) = XaP;
		if (++c == 4) break;
		mxBCs.elem(1, c) = c * XaP;
		XaP *= Xa;
	}
	mxBCs.elem(1, 0) = 0;

	// set vector 'vtBCs'
	vtBCs[0] = m_yValInit[0];		// Pos(Xa)
	vtBCs[1] = m_yValInit[1];		// Vel(Xa)

/*
	vta = [ a0  a1  a2  a3 ]'

	       [  xa^2   -2xa     1    0 ]
	mxda = [ 2xa^3   -3xa^2   0    1 ]

	      [ Sa(x^0)  Sa(x^1)  Sa(x^2)  Sa(x^3) ]
	mxS = [ Sa(x^1)  Sa(x^2)  Sa(x^3)  Sa(x^4) ]
	      [ Sa(x^2)  Sa(x^3)  Sa(x^4)  Sa(x^5) ]
	      [ Sa(x^3)  Sa(x^4)  Sa(x^5)  Sa(x^6) ]

	       [ S(x^0.y) ]
	vtSy = [ S(x^1.y) ]
	       [ S(x^2.y) ]
	       [ S(x^3.y) ]

	        [ 1        xa       xa^2     xa^3 ]
	mxBCs = [ 0        1       2xa      3xa^2 ]

	vtBCs = [ Pos(xa) ]
	        [ Vel(xa) ]

	------------------------------
	mxda * mxS * vta  =  mxda * vtSy		// 2 equ's		(note: mxda doesn't have an inverse - not square!)
	     mxBCs * vta  =  vtBCs				// 2 equ's


	         [ 1   -xa    xa^2  -xa^3 ]
	mxBtoA = [ 0    1   -2xa    3xa^2 ]
	         [ 0    0     1    -3xa   ]
	         [ 0    0     0      1    ]

	mxda * mxS * vta  =  mxda * vtSy		// 2 equ's		(note: mxda doesn't have an inverse - not square!)
	     mxBCs * vta  =  vtBCs				// 2 equ's

*/

//#define TwoByTwo
//#ifndef TwoByTwo
{	
	CMatrix mxSBCs(4,4);
	CMatrix vtSyBCs(4);

	mxSBCs.ProdPart(mxda, mxS);
	vtSyBCs.ProdPart(mxda, vtSy);
	mxBCs.CopyToDestLoc(mxSBCs, 2,0);
	vtBCs.CopyToDestLoc(vtSyBCs, 2,0);

	mxSBCs.LUSolve(m_coef, vtSyBCs.GetArray());

for (i = 0; i<4; i++)
	c44[i] = m_coef[i];
}
//#endif
//#if defined TwoByTwo

	CMatrix mxSBCs(2,4);			// only two top rows required!
	CMatrix vtSyBCs(2);			// only two top rows required!

	mxSBCs.Prod(mxda, mxS);
	vtSyBCs.Prod(mxda, vtSy);

afxDump << "mxda: " << mxda;
afxDump << "mxS: " << mxS;
afxDump << "vtSy: " << vtSy;

afxDump << "mxSBCs: " << mxSBCs;
afxDump << "vtSyBCs: " << vtSyBCs;

	CMatrix mxBtoA(4,4);
	mxBtoA.Identity();
	mxBtoA.elem(0,1) = -Xa;
	mxBtoA.elem(1,2) = -2*Xa;
	mxBtoA.elem(2,3) = -3*Xa;
	XaP = Xa*Xa;
	mxBtoA.elem(0,2) = XaP;
	mxBtoA.elem(1,3) = 3*XaP;
	mxBtoA.elem(0,3) = -Xa*XaP;

	CMatrix mxbSBCs(2,4);
	mxbSBCs.Prod(mxSBCs, mxBtoA);		// only two top rows required!
//	vtbSyBCs = vtSyBCs;
	CMatrix mxb22(2,2);
	CMatrix mxTemp22(2,2);
	CMatrix vtb22(2);

	mxb22.CopyFromSrcLoc(mxbSBCs, 0,2);
	mxb22.Invert();
	mxTemp22.CopyFromSrcLoc(mxbSBCs, 0,0);

	vtb22 = vtSyBCs - mxTemp22 * vtBCs;

afxDump << "mxBtoA: " << mxBtoA;
afxDump << "mxbSBCs: " << mxbSBCs;
afxDump << "mxb22: " << mxb22;
afxDump << "mxTemp22: " << mxTemp22;
afxDump << "vtb22: " << vtb22;

	CMatrix vta(4), vtb(4);
	vtb.ProdPart(mxb22, vtb22);
	vtb[2] = vtb[0];			// move result to elements 2,3
	vtb[3] = vtb[1];
	vtb[0] = vtBCs[0];
	vtb[1] = vtBCs[1];

	vta.Prod(mxBtoA, vtb);
	vta.CopyToArray(m_coef);

afxDump << "vta: " << vta;
afxDump << "vtb: " << vtb;

for (i = 0; i<4; i++)
	c22[i] = m_coef[i];

//#endif

	m_iFitDOF = 2;
}



void CPolySegInfo::FitDblPolyToPointsContPV(CPolySegInfo& nextPoly, double Wa, double Wb)
{
// Fit two polys to current and extrapolated points with continuous Pos and Vel at initial and middle knot

	const CPointPowerData& SumPtsA = m_PointSums.GetLastSums();
	const CPointPowerData& SumPtsB = nextPoly.m_PointSums.GetLastSums();
	const double (&SumAXp)[7] = SumPtsA.Xp;			// better array display for debuging!
	const double (&SumAYXp)[4] = SumPtsA.YXp;
	const double (&SumBXp)[7] = SumPtsB.Xp;
	const double (&SumBYXp)[4] = SumPtsB.YXp;

//	double* fitPoly = poly.coef;


	CMatrix mxdab(4,8), mxSab(8,8), mxBCs(4,8);
	CMatrix vtSaby(8), vtBCs(4);

	double Xa = m_xPolyInit;				// initial x of first poly
	double Xb = nextPoly.m_xPolyInit;	// initial x of next poly
	double XaP, XbP;

	// set matrix 'mxdab'
	mxdab = 0;
	mxdab.elem(0  ,1  ) = -2*Xa;
	mxdab.elem(0+2,1+4) = -2*Xb;
	mxdab.elem(0  ,1+4) = -2*(Xa-Xb);
	XaP = Xa*Xa;
	XbP = Xb*Xb;
	mxdab.elem(0  ,0  ) = XaP;
	mxdab.elem(0+2,0+4) = XbP;
	mxdab.elem(0  ,0+4) = (XaP-XbP);
	mxdab.elem(1  ,1  ) = -3*XaP;
	mxdab.elem(1+2,1+4) = -3*XbP;
	mxdab.elem(1  ,1+4) = -3*(XaP-XbP);
	XaP *= Xa;
	XbP *= Xb;
	mxdab.elem(1  ,0  ) = 2*XaP;
	mxdab.elem(1+2,0+4) = 2*XbP;
	mxdab.elem(1  ,0+4) = 2*(XaP-XbP);
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(0+2,2+4) = 1;
	mxdab.elem(1+2,3+4) = 1;

	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumPtsA.Xp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumPtsB.Xp[r+c];
		}

	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumPtsA.YXp[r];
		vtSaby[r+4] = Wb * SumPtsB.YXp[r];
	}

	// set matrix 'mxBCs'
	XaP = XbP = 1;
	int c = 0;
	while (1)
	{
		mxBCs.elem(0, c+0) = XaP;
		mxBCs.elem(2, c+0) = XbP;
		mxBCs.elem(2, c+4) = -XbP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+4) = 0;
		if (++c == 4) break;
		mxBCs.elem(1, c+0) = c * XaP;
		mxBCs.elem(3, c+0) = c * XbP;
		mxBCs.elem(3, c+4) = -c * XbP;
		XaP *= Xa;
		XbP *= Xb;
	}
	mxBCs.elem(1, 0) = 0;
	mxBCs.elem(3, 0) = 0;
	mxBCs.elem(3, 4) = 0;

	// set vector 'vtBCs'
	vtBCs[0] = m_yValInit[0];		// Pos(Xa)
	vtBCs[1] = m_yValInit[1];		// Vel(Xa)
	vtBCs[2] = 0;
	vtBCs[3] = 0;


/*
	vtab = [ a0  a1  a2  a3  b0  b1  b2  b3 ]'

	        [  xa^2   -2xa     1    0     xa^2-xb^2     2xb-2xa     0    0 ]
	mxdab = [ 2xa^3   -3xa^2   0    1   2xa^3-2xb^3   3xb^2-3xa^2   0    0 ]
	        [  0        0      0    0          xb^2        -2xb     1    0 ]
	        [  0        0      0    0         2xb^3        -3xb^2   0    1 ]

	        [ Sa(x^0)  Sa(x^1)  Sa(x^2)  Sa(x^3)  0        0        0        0       ]
	mxSab = [ Sa(x^1)  Sa(x^2)  Sa(x^3)  Sa(x^4)  0        0        0        0       ]
	        [ Sa(x^2)  Sa(x^3)  Sa(x^4)  Sa(x^5)  0        0        0        0       ]
	        [ Sa(x^3)  Sa(x^4)  Sa(x^5)  Sa(x^6)  0        0        0        0       ]
	        [ 0        0        0        0        Sb(x^0)  Sb(x^1)  Sb(x^2)  Sb(x^3) ]
	        [ 0        0        0        0        Sb(x^1)  Sb(x^2)  Sb(x^3)  Sb(x^4) ]
	        [ 0        0        0        0        Sb(x^2)  Sb(x^3)  Sb(x^4)  Sb(x^5) ]
	        [ 0        0        0        0        Sb(x^3)  Sb(x^4)  Sb(x^5)  Sb(x^6) ]

	         [ Sa(x^0.y) ]
	vtSaby = [ Sa(x^1.y) ]
	         [ Sa(x^2.y) ]
	         [ Sa(x^3.y) ]
	         [ Sb(x^0.y) ]
	         [ Sb(x^1.y) ]
	         [ Sb(x^2.y) ]
	         [ Sb(x^3.y) ]

	        [ 1        xa       xa^2     xa^3     0        0        0        0    ]
	mxBCs = [ 0        1       2xa      3xa^2     0        0        0        0    ]
	        [ 1        xb       xb^2     xb^3    -1       -xb      -xb^2    -xb^3 ]
	        [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]

	        [ Pos(xa) ]
	vtBCs = [ Vel(xa) ]
	        [     0   ]
	        [     0   ]


	------------------------------
	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's



	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

	vtAB = [ Pos(xa)  Vel(xa)  0  0  A2  A3  B2  B3 ]'

	          [ 1        xa       xa^2     xa^3     0        0        0        0    ]
	mxab2AB = [ 0        1       2xa      3xa^2     0        0        0        0    ]
	          [ 1        xb       xb^2     xb^3    -1       -xb      -xb^2    -xb^3 ]
	          [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]

	          [ 1       -xa       0        0        xa^2    2xa^3       0        0    ]
	mxAB2ab = [ 0        1        0        0      -2xa     -3xa^2       0        0    ]
	          [ 0        0        0        0        1        0          0        0    ]
	          [ 0        0        0        0        0        1          0        0    ]
	          [ 1       -xa      -1        xb  xa^2-xb^2  2xa^3-2xb^3   xb^2    2xb^3 ]
	          [ 0        1        0       -1   2xb-2xa    3xb^2-3xa^2 -2xb     -3xb^2 ]
	          [ 0        0        0        0        0        0          1        0    ]
	          [ 0        0        0        0        0        0          0        1    ]

  if
	vtAB = [ Pos(xa)  Vel(xa)  A2  A3  0  0  B2  B3 ]'

	          [ 1       -xa       xa^2     2xa^3       0      0      0        0    ]
	mxAB2ab = [ 0        1      -2xa      -3xa^2       0      0      0        0    ]
	          [ 0        0        1         0          0      0      0        0    ]
	          [ 0        0        0         1          0      0      0        0    ]
	          [ 1       -xa  xa^2-xb^2   2xa^3-2xb^3  -1      xb     xb^2    2xb^3 ]
	          [ 0        1   2xb-2xa     3xb^2-3xa^2   0     -1    -2xb     -3xb^2 ]
	          [ 0        0        0         0          0      0      1        0    ]
	          [ 0        0        0         0          0      0      0        1    ]

*/

double c88[8], c44[8];

//#define FourByFour
//#ifndef FourByFour				// solves an 8x8 matrix
{	
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	vtAllVals.ProdPart(mxdab, vtSaby);
	mxBCs.CopyToDestLoc(mxAllEqus, 4,0);
	vtBCs.CopyToDestLoc(vtAllVals, 4,0);

afxDump << "mxdab: " << mxdab;
afxDump << "mxSab: " << mxSab;
afxDump << "vtSaby: " << vtSaby;
afxDump << "mxBCs: " << mxBCs;
afxDump << "vtBCs: " << vtBCs;

afxDump << "mxAllEqus: " << mxAllEqus;
afxDump << "vtAllVals: " << vtAllVals;

// rearrange rows
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);
//	int arRowArrange[] = {4, 5, 0, 1, 6, 7, 2, 3}; 
//                       0, 1, 5, 4, 2, 6, 3, 7   matlab p
	int arRowArrange[] = {4, 5, 7, 6, 0, 2, 1, 3};
	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}
afxDump << "mxAllEqus2: " << mxAllEqus2;
afxDump << "vtAllVals2: " << vtAllVals2;

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

afxDump << "vtab: " << vtab;

	CMatrix vtValSolu;
	vtValSolu = mxAllEqus2 * vtab;
	CMatrix vtDiff;
	vtDiff = vtValSolu - vtAllVals2;
	afxDump << "vtValSolu: " << vtValSolu;
	afxDump << "vtDiff: " << vtDiff;

	for (int i = 0; i < 4; i++)
	{
		m_coef[i] = vtab[i];
		nextPoly.m_coef[i] = vtab[i+4];
	}

for (i = 0; i<8; i++) c88[i] = vtab[i];

}
//#endif
//#if defined FourByFour			// solves a 4x4 matrix


/*
	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's

	mxdab * mxSab * mxAB2ab * vtAB  =  mxdab * vtSaby		// 4 equ's

	mxdab * mxSab * mxAB2ab[2C] * vtAB[2R]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1C] * vtAB[1R]
	[1C] = 1st half of Columns etc.
	
*/
// Doesn't use mxBCs or vtBCs

	// set mxAB2ab2HC
	CMatrix mxAB2ab2HC(8,4);
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

	CMatrix mxdabS(4,8);
	CMatrix vtdabSy(4);
	CMatrix mxAllEqus(4,4);
	CMatrix vtAllVals(4);

	mxdabS.Prod(mxdab, mxSab);
	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

	CMatrix vtab(8);
//	vtab = mxAB2ab[1C] * vtAB[1R]
	vtab = 0;														// partly setup vtab
	vtab[0] = vtab[4] = m_yValInit[0] - Xa*m_yValInit[1];		// Pos(Xa) - Xa*Vel(Xa)
	vtab[1] = vtab[5] = m_yValInit[1];							// Vel(Xa)

	CMatrix vtValDiff(4);
	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

afxDump << "mxdab: " << mxdab;
afxDump << "mxSab: " << mxSab;
afxDump << "vtSaby: " << vtSaby;
afxDump << "mxdabS: " << mxdabS;
afxDump << "vtdabSy: " << vtdabSy;

afxDump << "mxAllEqus: " << mxAllEqus;
afxDump << "vtAllVals: " << vtAllVals;
afxDump << "partly set vtab: " << vtab;

	CMatrix vtAB2HR(4);
	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
//	vtab = mxAB2ab[1C] * vtAB[1R] + mxAB2ab[2C] * vtAB[2R]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

afxDump << "vtab: " << vtab;

	for (int i = 0; i < 4; i++)
	{
		m_coef[i] = vtab[i];
		nextPoly.m_coef[i] = vtab[i+4];
	}

for (i = 0; i<8; i++) c44[i] = vtab[i];

//#endif

	m_iFitDOF = 2;
	nextPoly.m_iFitDOF = 2;
}




