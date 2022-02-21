// PolySegParaFit.cpp: implementation of the CPolySegParaFit class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "Matrix.h"
#include "PolyFunc.h"


#include "PolySegParaFit.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

enum		// Segment Fit Flags
{
	SFF_POSINIT  = 0x01,
	SFF_POSFINAL = 0x02,
	SFF_DIRINIT  = 0x04,
	SFF_DIRFINAL = 0x08,
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPolySegParaFit::CPolySegParaFit()
{
	Reset();
}

CPolySegParaFit::~CPolySegParaFit()
{
}


//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

void CPolySegParaFit::Reset()
{
	m_numPts = 0;
	m_nFlags = 0;
}
void CPolySegParaFit::AddPoint(CVector* pPt)
{
	m_arPts[m_numPts++] = *pPt;
	ASSERT(m_numPts <= sizeof(m_arPts)/sizeof(*m_arPts));
}
void CPolySegParaFit::AddPoint(CVector& Pt)
{
	m_arPts[m_numPts++] = Pt;
	ASSERT(m_numPts <= sizeof(m_arPts)/sizeof(*m_arPts));
}
void CPolySegParaFit::SetToInitPoint()
{
	m_vtPosInit = m_arPts[0];
	m_nFlags |= SFF_POSINIT;
}
void CPolySegParaFit::SetToFinalPoint()
{
	m_vtPosFinal = m_arPts[m_numPts-1];
	m_nFlags |= SFF_POSFINAL;
}
void CPolySegParaFit::SetInitPos(CVector& vtPosI)
{
	m_vtPosInit = vtPosI;
	m_nFlags |= SFF_POSINIT;
}
void CPolySegParaFit::SetFinalPos(CVector& vtPosF)
{
	m_vtPosFinal = vtPosF;
	m_nFlags |= SFF_POSFINAL;
}
void CPolySegParaFit::SetInitDir(CVector& vtDirI)
{
	m_vtDirInit = vtDirI;
	m_nFlags |= SFF_DIRINIT;
}
void CPolySegParaFit::SetFinalDir(CVector& vtDirF)
{
	m_vtDirFinal = vtDirF;
	m_nFlags |= SFF_DIRFINAL;
}

void CPolySegParaFit::Fit()
{
	double Kv0, Kv1;
	int count = 0;

	SetSValues();
	NormaliseSValues();
	while (1)
	{
		SetSums();					// set SumSp, SumSpV
//		m_PointSums.ScaleS(1/m_arS[m_numPts-1]);		// set as if max s was 1, m_arS still need scaling

		FitCubicPos();
		GetBezier();
		
		FitCubicTwoPos();		// just used to get initial derivative
		GetBezier();

		FitCubicTwoPos0to1();
		GetBezier();
		
		
		Kv0 = m_Poly[1].Mag() / m_vtDirInit.Mag();
		Kv1 = (m_Poly[3]*3 + m_Poly[2]*2 + m_Poly[1]).Mag() / m_vtDirFinal.Mag();
//		FitCubicPosDeriv(Kv);

		GetBezier();
		AdjustSValues();
		m_SumResidual = m_SumResidual;
		count++;
		if (m_dsAvgNorm <= 1e-4)
			break;
	}

	GetBezier();


}

void CPolySegParaFit::Main()
{
	// set some points!
	CVector pts[7];
	pts[0].Set(0,0,0);
	pts[1].Set(0,1,0);
	pts[2].Set(6,10,0);
	pts[3].Set(8,7,0);
	pts[4].Set(10,3,0);
	pts[5].Set(9,2,0);
	pts[6].Set(6,-2,0);
	m_vtPosInit.Set(0,0,0);
	m_vtDirInit.Set(-0.5,3,0);
	m_numPts = sizeof(pts) / sizeof(*pts);
	for (int i = 0; i < m_numPts; i++)
		m_arPts[i] = pts[i];

	LARGE_INTEGER liFreq;
	LARGE_INTEGER liStart;
	LARGE_INTEGER liEnd;
	
	VERIFY(QueryPerformanceFrequency(&liFreq));
	double countPeriod = 1e3 / liFreq.LowPart;		// in ms
	QueryPerformanceCounter(&liStart);
	
	int count1 = 0;
	SetSValues();
	while (1)
	{
		SetSums();				// set SumSp, SumSpV
		FitCubicPosDir();
		GetBezier();
		AdjustSValues();
		m_SumResidual = m_SumResidual;
		count1++;
		if (m_dsAvgNorm <= 1e-4)
			break;
	}
	double sr1 = m_SumResidual;

	QueryPerformanceCounter(&liEnd);
	double t1 = (liEnd.LowPart - liStart.LowPart) * countPeriod;
	QueryPerformanceCounter(&liStart);


	int count2 = 0;
	double Kv = 1;
	SetSValues();
	while (1)
	{
		SetSums();					// set SumSp, SumSpV
		FitCubicPos();		// just used to get initial derivative
		Kv = m_Poly[1].Mag() / m_vtDirInit.Mag();
		FitCubicPosDeriv(Kv);

		GetBezier();
		AdjustSValues();
		m_SumResidual = m_SumResidual;
		count2++;
		if (m_dsAvgNorm <= 1e-4)
			break;
	}
	double sr2 = m_SumResidual;

	QueryPerformanceCounter(&liEnd);
	double t2 = (liEnd.LowPart - liStart.LowPart) * countPeriod;

// need to set m_PointSumsA & m_PointSumsB!
//	FitDblCubicPosDeriv(1);
//	FitDblCubicPos();

}

void CPolySegParaFit::SetSValues()
{
//	Increase in s value between points is equal to distance between points
	double s = 0;
	m_arS[0] = s;
	for (int i = 1; i < m_numPts; i++)
	{
		s += (m_arPts[i] - m_arPts[i-1]).Mag();
		m_arS[i] = s;
	}
	m_SInit = 0;
	m_SFinal = s;
}

void CPolySegParaFit::NormaliseSValues()
{
	// scale so s varies from 0 to 1
	double fact = 1 / m_arS[m_numPts-1];
	for (int i = 0; i < m_numPts; i++)
		m_arS[i] *= fact;
	m_SInit = m_arS[0];					// should be 0
	m_SFinal = m_arS[m_numPts-1];		// should be 1
}

void CPolySegParaFit::AdjustSValues()
{
	CVector vtPos, vtVel, vtAcc;
	CVector vtToPoint;
	double sumResSq = 0;
	double dsSum = 0;
	for (int i = 0; i < m_numPts; i++)
	{
		double s = m_arS[i];
		double sInit = s;
		double ds;
		double magSqToPoint;
		CVector vtPt = m_arPts[i];
		while (1)
		{
			vtPos = CubicAt(m_Poly, s);
			vtToPoint = vtPt - vtPos;
			magSqToPoint = vtToPoint.MagSq();
			if (magSqToPoint == 0)		// no change needed
				break;
			//	find ds
			vtVel = CubicD1At(m_Poly, s);
			vtAcc = CubicD2At(m_Poly, s);
			double dsTang, dsChord, dsCirc, dsTang2;
			double velMagSq = vtVel.MagSq();
			double velMag = sqrt(velMagSq);
			CVector vtCurve = cross(vtVel, vtAcc) / (velMagSq * velMag);
			double curveMagSq = vtCurve.MagSq();
				dsTang = dot(vtToPoint, vtVel) / velMagSq;	// use tangent method
			if (dsTang == 0)			// at best location
				break;
			if (curveMagSq < 1e-8)
				// tangent method
				ds = dot(vtToPoint, vtVel) / velMagSq;	// use tangent method
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
				double sinAngCross = (cross(vtRadius, vtCentre2Point)).Mag() / (radius * vtCentre2Point.Mag());	// doesn't give +/- direction
			// circumference method using chord
				double sinAng = dot(vtCentre2Chord, vtVel) / (velMag * radius);	// gives +/- direction
				double ang = asin(sinAng);
				double circDist = radius * ang;
					dsCirc = circDist / velMag;

				ds = dsCirc;
			}
			s += ds;
			if (fabs(ds) < 1e-8)
				break;
		}	// while (1)
		sumResSq += magSqToPoint;
		ds = s - sInit;
		m_arS[i] = s;
		dsSum += fabs(ds);
	}	// for (int i = 0; i < m_numPts; i++)
	double sSpan = m_arS[m_numPts-1] - m_arS[0];
	ASSERT(sSpan > 0);
	m_dsAvgNorm = dsSum / (m_numPts * sSpan);		// ds averaged over numPts and normalised for total s length
//	m_SumResidual = sumResSq;		// or sqrt(sumResSq / (m_numPts - ?))
	m_SumResidual = sqrt(sumResSq / (m_numPts-1));

	m_SInit = m_arS[0];
	m_SFinal = m_arS[m_numPts-1];
}


void CPolySegParaFit::SetSums()
{
	// set m_PointSums
	CPointSPowerData sums;
	sums.Zero();
	for (int i = 0; i < m_numPts; i++)
		sums.SumUpPoint(m_arS[i], m_arPts[i]);
	m_PointSums = sums;
}

void CPolySegParaFit::GetBezier()
{
	// get Bezier poly (poly span not 0-1 so can't use Cubic2Bezier(m_Bez, m_Poly)
	double sInit = m_arS[0];
	double sFinal = m_arS[m_numPts - 1];
	double sSpanOn3 = (sFinal - sInit) / 3;
	m_Bez[0] = CubicAt(m_Poly, sInit);
	m_Bez[3] = CubicAt(m_Poly, sFinal);
	m_Bez[1] = m_Bez[0] + CubicD1At(m_Poly, sInit) * sSpanOn3;
	m_Bez[2] = m_Bez[3] - CubicD1At(m_Poly, sFinal) * sSpanOn3;
}



void CPolySegParaFit::FitCubicPosDir()
{
// Fit parametric cubic with initial pos & direction at s = 0
// Axes are solved together (larger matrix solution)

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtDirInit = m_vtDirInit;



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
		vtVals[2*ax+1] = SumSpV[2][ax] - vtPosInit[ax] * SumSp[2];
		vtVals[2*ax+2] = SumSpV[3][ax] - vtPosInit[ax] * SumSp[3];
		vtVal0 += vtDirInit[ax] * (SumSpV[1][ax] - vtPosInit[ax] * SumSp[1]);
	}
	vtVals[0] = vtVal0;

	// solve
	CMatrix vta(7);
	mxS.LUSolve(vta.GetArray(), vtVals.GetArray());

	// set poly coeffs
	m_Poly[0] = vtPosInit;
	m_Poly[1] = vtDirInit * vta[0];		// vtDirInit * k
	m_Poly[2].Set(vta[1], vta[3], vta[5]);
	m_Poly[3].Set(vta[2], vta[4], vta[6]);

}

void CPolySegParaFit::FitCubicPosDeriv(double Kv)
{
// Fit parametric cubic with initial pos & derivative = Kv * m_vtDirInit at s = 0
// Each axis is solved independently

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtDirInit = m_vtDirInit;

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
		a0 = vtPosInit[ax];
		a1 = vtDirInit[ax] * Kv;
		b0 = SumSpV[2][ax] - SumSp[2]*a0 - SumSp[3]*a1;
		b1 = SumSpV[3][ax] - SumSp[3]*a0 - SumSp[4]*a1;
		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = a1;
		m_Poly[2][ax] = det *  ( SumSp[6]*b0 - SumSp[5]*b1);
		m_Poly[3][ax] = det *  (-SumSp[5]*b0 + SumSp[4]*b1);
	}
}

void CPolySegParaFit::FitCubicPos()
{
// Fit parametric cubic with initial pos at s = 0
// Each axis is solved independently

	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
/*
		If initial pos and vel are at x=0 then equ's are simplified -> a0 & a1 known
		[ S(s^2)  S(s^3)  S(s^4) ]   a1     S(s^1.y)   [ S(s^1) ]
		[ S(s^3)  S(s^4)  S(s^5) ] * a2  =  S(s^2.y) - [ S(s^2) ] * a0
		[ S(s^4)  S(s^5)  S(s^6) ]   a3     S(s^3.y)   [ S(s^3) ]
*/
	double arA[4];
	double arVal[4];
	double a0;

	for (int ax  = 0; ax < 3; ax++)
	{
		a0 = vtPosInit[ax];
		for (int i  = 1; i < 4; i++)
			arVal[i] = SumSpV[i][ax] - SumSp[i] * a0;
		VERIFY(LUFullSymSolve(3, SumSp+2, arA+1, arVal+1));		// solves matrix equ, result in m_arFitPoly

		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = arA[1];
		m_Poly[2][ax] = arA[2];
		m_Poly[3][ax] = arA[3];
	}

	CMatrix mxS(3,3);
	CMatrix vtA(3), vtVal(3);
	// set matrix 'mxS'
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			mxS.elem(r,c) = SumSp[2+r+c];
	mxS.Invert();

	for (ax  = 0; ax < 3; ax++)
	{
		a0 = vtPosInit[ax];
		for (int i  = 1; i < 4; i++)
			vtVal[i-1] = SumSpV[i][ax] - SumSp[i] * a0;
		vtA.Prod(mxS, vtVal);
		m_Poly[0][ax] = a0;
		m_Poly[1][ax] = vtA[0];
		m_Poly[2][ax] = vtA[1];
		m_Poly[3][ax] = vtA[2];
	}
}

void CPolySegParaFit::FitCubicTwoPos()
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
	[ 0     sa-sb  sa^2-sb^2  sa^3-sb^3 ]    = Pos(sa)-Pos(sb)					// r1 - r2
	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)   0 ] = sb^3*Pos(sa)-sa^3*Pos(sb)		// sb^3*r1 - sa^3*r2

	[ sb^3  -sa^3  sb*sa(sa^2-sb^2)  sb^2*sa^2(sa-sb) ]   Pos(sa)   a0*(sb^3-sa^3)
	[ 0        1        0                 0           ] * Pos(sb) = a1
	[ 0        0        1                 0           ]   A1        a2
	[-1        1     sa-sb           sa^2-sb^2        ]   A2        a3*(sb^3-sa^3)
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[sb*sa(sa^2-sb^2)/(sb^3-sa^3)   0   1   (sa-sb)/(sb^3-sa^3)     ]
	[sb^2*sa^2(sa-sb)/(sb^3-sa^3)   0   0   (sa^2-sb^2)/(sb^3-sa^3) ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

*/
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;
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
	mxda.elem(0,1) = 0;
	mxda.elem(1,1) = 0;
	mxda.elem(0,2) = 1;
	mxda.elem(1,2) = 0;
	mxda.elem(0,3) = e03;
	mxda.elem(1,3) = e13;

	// set matrix 'mxA2aKC'
	mxA2aKC.elem(0,0) = sb3 * sb3msa3Inv;
	mxA2aKC.elem(1,0) = 0;
	mxA2aKC.elem(2,0) = 0;
	mxA2aKC.elem(3,0) = -sb3msa3Inv;
	mxA2aKC.elem(0,1) = -sa3 * sb3msa3Inv;
	mxA2aKC.elem(1,1) = 1;
	mxA2aKC.elem(2,1) = 0;
	mxA2aKC.elem(3,1) = sb3msa3Inv;

	// set matrix 'mxS'
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
			mxS.elem(r,c) = SumSp[r+c];


	CMatrix mxdaS(2,4), mxEqu(2,2);
	mxdaS.Prod(mxda, mxS);
	mxEqu.ProdTr(mxdaS, mxda);
	mxEqu.Invert();

//		(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
//	or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')

	CMatrix mxSA2aKC(4,2);
	mxSA2aKC.Prod(mxS, mxA2aKC);

	CMatrix vtSv(4), vtAKR(2), vtVal(2), vtAUR(2), vta(4);
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
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

void CPolySegParaFit::FitCubicTwoPos0to1()
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
	[ 0     sa-sb  sa^2-sb^2  sa^3-sb^3 ]    = Pos(sa)-Pos(sb)					// r1 - r2
	[ sb^3-sa^3   sb*sa(sb^2-sa^2)   sb^2*sa^2(sb-sa)   0 ] = sb^3*Pos(sa)-sa^3*Pos(sb)		// sb^3*r1 - sa^3*r2

	[ sb^3  -sa^3  sb*sa(sa^2-sb^2)  sb^2*sa^2(sa-sb) ]   Pos(sa)   a0*(sb^3-sa^3)
	[ 0        1        0                 0           ] * Pos(sb) = a1
	[ 0        0        1                 0           ]   A1        a2
	[-1        1     sa-sb           sa^2-sb^2        ]   A2        a3*(sb^3-sa^3)
	mxA2a * vtA = vta

	mxda =    mxA2a(uc)' =
	[sb*sa(sa^2-sb^2)/(sb^3-sa^3)   0   1   (sa-sb)/(sb^3-sa^3)     ]
	[sb^2*sa^2(sa-sb)/(sb^3-sa^3)   0   0   (sa^2-sb^2)/(sb^3-sa^3) ]

	(mxda*mxS) * mxda' * [A1 A2]'  =  mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]'
or	mxda * mxS * mxda' * [A1 A2]'  =  mxda * (vtSv - mxS * mxA2a(kc) * [Pos(sa) Pos(sb)]')
		uc = unknown columns
		kc = known columns

	---------------------------
	if sa = 0 & sb = 1
	mxA2a =
	[ 1   0   0   0 ]
	[ 0   1   0   0 ]
	[ 0   0   1   0 ]
	[-1   1  -1  -1 ]
	mxda =
	[ 0  0  1 -1 ]
	[ 0  0  0 -1 ]

	mxda*mxS =
	[ S2-S3  S3-S4  S4-S5  S5-S6 ]  
	[   -S3    -S4    -S5    -S6 ]
	mxda*mxS*mxda' =
	[ S4-2*S5+S6   S6-S5 ]
	[    S6-S5        S6 ]

	mxda * vtSv =
	[ S2v-S3v ]
	[    -S3v ]

	(mxda*mxS) * mxA2a(kc) =
	[ S2-S3-S5+S6   S3-S4+S5-S6 ]
	[    S6-S3        -S4-S6    ]

	mxda * vtSv - (mxda*mxS) * mxA2a(kc) * [Pos(sa) Pos(sb)]' =
	[ S2v-S3v - (S2-S3-S5+S6)Pos(sa) - (S3-S4+S5-S6)Pos(sb) ]
	[    -S3v - (S6-S3)Pos(sa) - (-S4-S6)Pos(sb)             ]
*/
	const double (&SumSp)[7] = m_PointSums.Sp;
	double (&SumSpV)[4][3] = m_PointSums.SpV;

	CMatrix mxEqu(2,2), mxdaSA2aKC(2,2);
	// set matrix 'mxEqu'
	mxEqu.elem(0,0) = SumSp[4] - 2*SumSp[5] + SumSp[6];
	mxEqu.elem(0,1) = SumSp[6] - SumSp[5];
	mxEqu.elem(1,0) = SumSp[6] - SumSp[5];
	mxEqu.elem(1,1) = SumSp[6];
	mxEqu.Invert();

	// set matrix 'mxdaSA2aKC'
	mxdaSA2aKC.elem(0,0) = SumSp[2] - SumSp[3] - SumSp[5] + SumSp[6];
	mxdaSA2aKC.elem(0,1) = SumSp[3] - SumSp[4] + SumSp[5] - SumSp[6];
	mxdaSA2aKC.elem(1,0) = SumSp[6] - SumSp[3];
	mxdaSA2aKC.elem(1,1) =-SumSp[4] - SumSp[6];

	CMatrix vtAKR(2), vtVal(2), vtAUR(2);
	///////////////////////////
	// vtVal is calculated for each axis
	for (int ax = 0; ax < 3; ax++)
	{
		// set vector 'vtAKR'
		vtAKR[0] = m_vtPosInit[ax];
		vtAKR[1] = m_vtPosFinal[ax];

		// set vector 'vtVal'
		vtVal.Prod(mxdaSA2aKC, vtAKR);
		vtVal[0] = SumSpV[2][ax] - SumSpV[3][ax] - vtVal[0];
		vtVal[1] =-SumSpV[3][ax] - vtVal[1];
		vtAUR.Prod(mxEqu, vtVal);

		m_Poly[0][ax] = vtAKR[0];
		m_Poly[1][ax] = vtAUR[0];
		m_Poly[2][ax] = vtAUR[1];
		m_Poly[3][ax] = vtAKR[1] - vtAKR[0] - vtAUR[0] - vtAUR[1];
	}
}

void CPolySegParaFit::FitCubicTwoPosDeriv0to1(double Kv0, double Kv1)
{
// Fit parametric cubic with initial and final pos & derivative at s = 0 & 1
//	deriv(0) = Kv0 * m_vtDirInit, deriv(1) = Kv1 * m_vtDirFinal
// Each axis is solved independently

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
	m_Poly[1] = m_vtDirInit * Kv0;
	m_Poly[3] = (m_vtPosInit - m_vtPosFinal) * 2 - m_Poly[1] - m_vtDirFinal * Kv1;
	m_Poly[2] = m_vtPosFinal - m_Poly[3] - m_Poly[1] - m_Poly[0];
//	m_Bez[0] = m_vtPosInit;
//	m_Bez[1] = m_vtPosInit + m_vtDirInit * (Kv0 / 3);
//	m_Bez[2] = m_vtPosFinal - m_vtDirFinal * (Kv1 / 3);
//	m_Bez[3] = m_vtPosFinal;
//	Bezier2Cubic(m_Poly, m_Bez);
}

void CPolySegParaFit::FitDblCubicPosDeriv(double Kv)
{
// Fit two parametric cubics with initial pos & derivative = Kv * m_vtDirInit at Sa
// Pos and deriv are continuous at knot at Sb
// Each axis is solved independently
	
	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_PointSumsA.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_PointSumsA.SpV;
	const double (&SumBSp)[7] = m_PointSumsB.Sp;
	double (&SumBSpV)[4][3] = m_PointSumsB.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtDerivInit = m_vtDirInit * Kv;


	CMatrix mxdab(4,8), mxSab(8,8), mxBCs(4,8);
	CMatrix vtSaby(8), vtBCs(4);

	double Sa = m_SPolyAInit;				// initial s of first poly
	double Sb = m_SPolyBInit;				// initial s of next poly
	double SaP, SbP;

	// set matrix 'mxdab'
	mxdab = 0;
	mxdab.elem(0  ,1  ) = -2*Sa;
	mxdab.elem(0+2,1+4) = -2*Sb;
	mxdab.elem(0  ,1+4) = -2*(Sa-Sb);
	SaP = Sa*Sa;
	SbP = Sb*Sb;
	mxdab.elem(0  ,0  ) = SaP;
	mxdab.elem(0+2,0+4) = SbP;
	mxdab.elem(0  ,0+4) = (SaP-SbP);
	mxdab.elem(1  ,1  ) = -3*SaP;
	mxdab.elem(1+2,1+4) = -3*SbP;
	mxdab.elem(1  ,1+4) = -3*(SaP-SbP);
	SaP *= Sa;
	SbP *= Sb;
	mxdab.elem(1  ,0  ) = 2*SaP;
	mxdab.elem(1+2,0+4) = 2*SbP;
	mxdab.elem(1  ,0+4) = 2*(SaP-SbP);
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(0+2,2+4) = 1;
	mxdab.elem(1+2,3+4) = 1;

	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumASp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumBSp[r+c];
		}

	// set matrix 'mxBCs'
	SaP = SbP = 1;
	int c = 0;
	while (1)
	{
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(2, c+0) = SbP;
		mxBCs.elem(2, c+4) = -SbP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+4) = 0;
		if (++c == 4) break;
		mxBCs.elem(1, c+0) = c * SaP;
		mxBCs.elem(3, c+0) = c * SbP;
		mxBCs.elem(3, c+4) = -c * SbP;
		SaP *= Sa;
		SbP *= Sb;
	}
	mxBCs.elem(1, 0) = 0;
	mxBCs.elem(3, 0) = 0;
	mxBCs.elem(3, 4) = 0;




/*
	vtab = [ a0  a1  a2  a3  b0  b1  b2  b3 ]'

	        [  sa^2   -2xa     1    0     sa^2-sb^2     2xb-2xa     0    0 ]
	mxdab = [ 2xa^3   -3xa^2   0    1   2xa^3-2xb^3   3xb^2-3xa^2   0    0 ]
	        [  0        0      0    0          sb^2        -2xb     1    0 ]
	        [  0        0      0    0         2xb^3        -3xb^2   0    1 ]

	        [ Sa(s^0)  Sa(s^1)  Sa(s^2)  Sa(s^3)  0        0        0        0       ]
	mxSab = [ Sa(s^1)  Sa(s^2)  Sa(s^3)  Sa(s^4)  0        0        0        0       ]
	        [ Sa(s^2)  Sa(s^3)  Sa(s^4)  Sa(s^5)  0        0        0        0       ]
	        [ Sa(s^3)  Sa(s^4)  Sa(s^5)  Sa(s^6)  0        0        0        0       ]
	        [ 0        0        0        0        Sb(s^0)  Sb(s^1)  Sb(s^2)  Sb(s^3) ]
	        [ 0        0        0        0        Sb(s^1)  Sb(s^2)  Sb(s^3)  Sb(s^4) ]
	        [ 0        0        0        0        Sb(s^2)  Sb(s^3)  Sb(s^4)  Sb(s^5) ]
	        [ 0        0        0        0        Sb(s^3)  Sb(s^4)  Sb(s^5)  Sb(s^6) ]

	         [ Sa(s^0.y) ]
	vtSaby = [ Sa(s^1.y) ]
	         [ Sa(s^2.y) ]
	         [ Sa(s^3.y) ]
	         [ Sb(s^0.y) ]
	         [ Sb(s^1.y) ]
	         [ Sb(s^2.y) ]
	         [ Sb(s^3.y) ]

	        [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxBCs = [ 0        1       2xa      3xa^2     0        0        0        0    ]
	        [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	        [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]

	        [ Pos(sa) ]
	vtBCs = [ Vel(sa) ]
	        [     0   ]
	        [     0   ]


	------------------------------
	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's



	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

	vtAB = [ Pos(sa)  Vel(sa)  0  0  A2  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxab2AB = [ 0        1       2xa      3xa^2     0        0        0        0    ]
	          [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	          [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]

	          [ 1       -sa       0        0        sa^2    2xa^3       0        0    ]
	mxAB2ab = [ 0        1        0        0      -2xa     -3xa^2       0        0    ]
	          [ 0        0        0        0        1        0          0        0    ]
	          [ 0        0        0        0        0        1          0        0    ]
	          [ 1       -sa      -1        sb  sa^2-sb^2  2xa^3-2xb^3   sb^2    2xb^3 ]
	          [ 0        1        0       -1   2xb-2xa    3xb^2-3xa^2 -2xb     -3xb^2 ]
	          [ 0        0        0        0        0        0          1        0    ]
	          [ 0        0        0        0        0        0          0        1    ]


	if vtAB = [ Pos(sa)  Vel(sa)  A2  A3  0  0  B2  B3 ]'
	then:
	          [ 1       -sa       sa^2     2xa^3       0      0      0        0    ]
	mxAB2ab = [ 0        1      -2xa      -3xa^2       0      0      0        0    ]
	          [ 0        0        1         0          0      0      0        0    ]
	          [ 0        0        0         1          0      0      0        0    ]
	          [ 1       -sa  sa^2-sb^2   2xa^3-2xb^3  -1      sb     sb^2    2xb^3 ]
	          [ 0        1   2xb-2xa     3xb^2-3xa^2   0     -1    -2xb     -3xb^2 ]
	          [ 0        0        0         0          0      0      1        0    ]
	          [ 0        0        0         0          0      0      0        1    ]

*/

double c88[8], c44[8];
/////////////////////////////////
// 8x8 matrix
//
{	
for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	// set vector 'vtBCs'
	vtBCs[0] = vtPosInit[ax];			// Pos(Sa)
	vtBCs[1] = vtDerivInit[ax];		// Vel(Sa)
	vtBCs[2] = 0;
	vtBCs[3] = 0;
	
	
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	vtAllVals.ProdPart(mxdab, vtSaby);
	mxBCs.CopyToDestLoc(mxAllEqus, 4,0);
	vtBCs.CopyToDestLoc(vtAllVals, 4,0);
/*
afxDump << "mxdab: " << mxdab;
afxDump << "mxSab: " << mxSab;
afxDump << "vtSaby: " << vtSaby;
afxDump << "mxBCs: " << mxBCs;
afxDump << "vtBCs: " << vtBCs;

afxDump << "mxAllEqus: " << mxAllEqus;
afxDump << "vtAllVals: " << vtAllVals;
*/
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

//afxDump << "vtab: " << vtab;

	CMatrix vtValSolu;
	vtValSolu = mxAllEqus2 * vtab;
	CMatrix vtDiff;
	vtDiff = vtValSolu - vtAllVals2;
	afxDump << "vtValSolu: " << vtValSolu;
	afxDump << "vtDiff: " << vtDiff;

	for (int i = 0; i < 4; i++)
	{
		m_PolyA[i][ax] = vtab[i];
		m_PolyB[i][ax] = vtab[i+4];
	}

for (i = 0; i<8; i++)
	c88[i] = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}








///////////////////////////////////////////////
// 4x4 matrix
// Calculate by changing vtab to vtAB with only 4 unknowns and rearranging to solve
// a 4x4 matrix equation
/*
	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's

	mxdab * mxSab * mxAB2ab * vtAB  =  mxdab * vtSaby		// 4 equ's

	where mxAB2ab is inv(mxab2AB) and is obtained easily algebraicly!
	Note:  mxab2AB * vtab = vtAB
	is derived from: mxBCs * vtab = vtBcs with ones added for extra vtab elements to make a square matrix

	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	[1HC] = 1st half of Columns etc.
	vtAB[2HR] is 2nd half of rows which are the only 4 unknowns

*/
{
for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}

	
	// set mxAB2ab2HC - '2HC' => 2nd half of columns
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
//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax] - Sa*vtDerivInit[ax];		// Pos(Sa) - Sa*Vel(Sa)
	vtab[1] = vtab[5] = vtDerivInit[ax];								// Vel(Sa)

	CMatrix vtValDiff(4);
	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;
/*
afxDump << "mxdab: " << mxdab;
afxDump << "mxSab: " << mxSab;
afxDump << "vtSaby: " << vtSaby;
afxDump << "mxdabS: " << mxdabS;
afxDump << "vtdabSy: " << vtdabSy;

afxDump << "mxAllEqus: " << mxAllEqus;
afxDump << "vtAllVals: " << vtAllVals;
afxDump << "partly set vtab: " << vtab;
*/
	CMatrix vtAB2HR(4);
	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;

	for (int i = 0; i < 4; i++)
	{
		m_PolyA[i][ax] = vtab[i];
		m_PolyB[i][ax] = vtab[i+4];
	}

for (i = 0; i<8; i++)
	c44[i] = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}




}

/*
	        [  sa^2   -2xa     1    0     sa^2-sb^2     2xb-2xa     0    0 ]
	mxdab = [ 2xa^3   -3xa^2   0    1   2xa^3-2xb^3   3xb^2-3xa^2   0    0 ]
	        [  0        0      0    0          sb^2        -2xb     1    0 ]
	        [  0        0      0    0         2xb^3        -3xb^2   0    1 ]

	        [ Sa(s^0)  Sa(s^1)  Sa(s^2)  Sa(s^3)  0        0        0        0       ]
	mxSab = [ Sa(s^1)  Sa(s^2)  Sa(s^3)  Sa(s^4)  0        0        0        0       ]
	        [ Sa(s^2)  Sa(s^3)  Sa(s^4)  Sa(s^5)  0        0        0        0       ]
	        [ Sa(s^3)  Sa(s^4)  Sa(s^5)  Sa(s^6)  0        0        0        0       ]
	        [ 0        0        0        0        Sb(s^0)  Sb(s^1)  Sb(s^2)  Sb(s^3) ]
	        [ 0        0        0        0        Sb(s^1)  Sb(s^2)  Sb(s^3)  Sb(s^4) ]
	        [ 0        0        0        0        Sb(s^2)  Sb(s^3)  Sb(s^4)  Sb(s^5) ]
	        [ 0        0        0        0        Sb(s^3)  Sb(s^4)  Sb(s^5)  Sb(s^6) ]

mxdAB[1HR](4x8)
4x8 * 8x8 * 8x4
Find mxdab * msSab * mxdab':
mxdab * msSab = 

sa^2*Sa0-2xa*Sa1+Sa2       sa^2*Sa1-2xa*Sa2+Sa3       sa^2*Sa2-2xa*Sa3+Sa4       sa^2*Sa3-2xa*Sa4+Sa5       (sa^2-sb^2)Sb0+2(sb-sa)Sb1         (sa^2-sb^2)Sb1+2(sb-sa)Sb2         (sa^2-sb^2)Sb2+2(sb-sa)Sb3         (sa^2-sb^2)Sb3+2(sb-sa)Sb4
2xa^3*Sa0-3xa^2*Sa1+Sa3    2xa^3*Sa1-3xa^2*Sa2+Sa4    2xa^3*Sa2-3xa^2*Sa3+Sa5    2xa^3*Sa3-3xa^2*Sa4+Sa6    2(sa^3-sb^3)Sb0+3(sb^2-sa^2)Sb1    2(sa^3-sb^3)Sb1+3(sb^2-sa^2)Sb2    2(sa^3-sb^3)Sb2+3(sb^2-sa^2)Sb3    2(sa^3-sb^3)Sb3+3(sb^2-sa^2)Sb4
  0                          0                          0                          0                        sb^2*Sb0-2xb*Sb1+Sb2               sb^2*Sb1-2xb*Sb2+Sb3               sb^2*Sb2-2xb*Sb3+Sb4               sb^2*Sb3-2xb*Sb4+Sb5
  0                          0                          0                          0                        2xb^3*Sb0-3xb^2*Sb1+Sb3            2xb^3*Sb1-3xb^2*Sb2+Sb4            2xb^3*Sb2-3xb^2*Sb3+Sb5            2xb^3*Sb3-3xb^2*Sb4+Sb6

(mxdab * msSab) * mxdab'
elem(0,0) =
sa^4*Sa0-2xa^3*Sa1+sa^2*Sa2 - (2xa^3*Sa1-4xa^2*Sa2+2xa*Sa3) + sa^2*Sa2-2xa*Sa3+Sa4 + 0 + (sa^2-sb^2)^2*Sb0+2(sb-sa)(sa^2-sb^2)Sb1 + 2(sb-sa)(sa^2-sb^2)Sb1+4(sb-sa)^2*Sb2 + 0 + 0
=> sa^4*Sa0 - 4xa^3*Sa1 + 6xa^2*Sa2 - 4xa*Sa3 + Sa4 + (sa^2-sb^2)^2*Sb0 + 4(sb-sa)(sa^2-sb^2)Sb1 + 4(sb-sa)^2*Sb2




*/


void CPolySegParaFit::FitDblCubicPos()
{
// Fit two parametric cubics with initial pos at Sa
// Pos and deriv are continuous at knot at Sb
// Each axis is solved independently

	
/*
	vtab = [ a0  a1  a2  a3  b0  b1  b2  b3 ]'

	        [ Sa(s^0)  Sa(s^1)  Sa(s^2)  Sa(s^3)  0        0        0        0       ]
	mxSab = [ Sa(s^1)  Sa(s^2)  Sa(s^3)  Sa(s^4)  0        0        0        0       ]
	        [ Sa(s^2)  Sa(s^3)  Sa(s^4)  Sa(s^5)  0        0        0        0       ]
	        [ Sa(s^3)  Sa(s^4)  Sa(s^5)  Sa(s^6)  0        0        0        0       ]
	        [ 0        0        0        0        Sb(s^0)  Sb(s^1)  Sb(s^2)  Sb(s^3) ]
	        [ 0        0        0        0        Sb(s^1)  Sb(s^2)  Sb(s^3)  Sb(s^4) ]
	        [ 0        0        0        0        Sb(s^2)  Sb(s^3)  Sb(s^4)  Sb(s^5) ]
	        [ 0        0        0        0        Sb(s^3)  Sb(s^4)  Sb(s^5)  Sb(s^6) ]

	         [ Sa(s^0.y) ]
	vtSaby = [ Sa(s^1.y) ]
	         [ Sa(s^2.y) ]
	         [ Sa(s^3.y) ]
	         [ Sb(s^0.y) ]
	         [ Sb(s^1.y) ]
	         [ Sb(s^2.y) ]
	         [ Sb(s^3.y) ]

	        [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxBCs = [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]	       
	        [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]

	        [ Pos(sa) ]
	vtBCs = [     0   ]
	        [     0   ]


------------------------------

	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

-----------------
	if vtAB = [ Pos(sa)  A1  A2  A3  0  0  B2  B3 ]'
	then:
	mxab2AB = [ 1       sa        sa^2         sa^3     0       0        0        0    ]
	          [ 0       1         0            0        0       0        0        0    ]
	          [ 0       0         1            0        0       0        0        0    ]
	          [ 0       0         0            1        0       0        0        0    ]
	          [ 1       sb        sb^2         sb^3    -1      -sb      -sb^2    -sb^3 ]
	          [ 0       1        2xb          3xb^2     0      -1      -2xb     -3xb^2 ]
	          [ 0       0         0            0        0       0        1        0    ]
	          [ 0       0         0            0        0       0        0        1    ]

to invert:
	for b1    [ 0      -1       -2xb         -3xb^2     0       1       2xb      3xb^2 ]	-vtAB(6)

	for b0	 [-1      -sb       -sb^2        -sb^3     1       sb       sb^2     sb^3 ]	-vtAB(5)
				 [ 0    sa-sb   sa^2-sb^2    sa^3-sb^3     1       sb       sb^2     sb^3 ]	+vtAB(1)
				 [ 0       sa   sa^2+sb^2   sa^3+2xb^3     1       0       -sb^2   -2xb^3 ]	+sb*vtAB(6)

	mxAB2ab = [ 1      -sa       -sa^2        -sa^3     0       0        0        0    ]
	          [ 0       1         0            0        0       0        0        0    ]
	          [ 0       0         1            0        0       0        0        0    ]
	          [ 0       0         0            1        0       0        0        0    ]
	          [ 1      -sa  -sa^2-sb^2  -sa^3-2xb^3    -1       sb       sb^2    2xb^3 ]
	          [ 0       1        2xb          3xb^2     0      -1      -2xb     -3xb^2 ]
	          [ 0       0         0            0        0       0        1        0    ]
	          [ 0       0         0            0        0       0        0        1    ]

-----------------

	Rearrange vtAB to group unknows:

	vtAB = [ Pos(sa)  0  0  A1  A2  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxab2AB = [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	          [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]
	          [ 0        1        0        0        0        0        0        0    ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]

	          [ 1        0        0       -sa      -sa^2    -sa^3       0        0    ]
	mxAB2ab = [ 0        0        0        1        0        0          0        0    ]
	          [ 0        0        0        0        1        0          0        0    ]
	          [ 0        0        0        0        0        1          0        0    ]
	          [ 1       -1        sb      -sa -sa^2-sb^2  -sa^3-2xb^3   sb^2    2xb^3 ]
	          [ 0        0       -1        1       2xb      3xb^2     -2xb     -3xb^2 ]
	          [ 0        0        0        0        0        0          1        0    ]
	          [ 0        0        0        0        0        0          0        1    ]

  mxdab is transpose of mxAB2ab(:,4:8)

	mxdab = [ -sa      1     0     0        -sa        1       0    0 ]
	        [ -sa^2    0     1     0   -sa^2-sb^2     2xb      0    0 ]
	        [ -sa^3    0     0     1  -sa^3-2xb^3     3xb^2    0    0 ]
	        [  0       0     0     0         sb^2    -2xb      1    0 ]
	        [  0       0     0     0        2xb^3    -3xb^2    0    1 ]


						mxdab
	dSr/dA2   [ da0/dA2  da1/dA2  da2/dA2  da3/dA2 ]   dSra/da0     0
	dSr/dA3 = [ da0/dA3  da1/dA3  da2/dA3  da3/dA3 ] * dSra/da1  =  0
	                                                   dSra/da2
	                                                   dSra/da3
	elements of mxAB2ab are da(row)/dA(col) so transpose gives mxdab
	mxdab gives dSumRes/d(vtAB) from dSumRes/d(vtab)

-----------------
	mxdab * mxSab * vtab  =  mxdab * vtSaby		- 5 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					- 3 equ's
	OR
	mxSab * mxAB2ab * vtAB  =  vtSaby				Incorporates BC's
	vtab = mxAB2ab * vtAB
	
-----------------



*/

	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_PointSumsA.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_PointSumsA.SpV;
	const double (&SumBSp)[7] = m_PointSumsB.Sp;
	double (&SumBSpV)[4][3] = m_PointSumsB.SpV;

	CVector vtPosInit = m_vtPosInit;

	CMatrix mxdab(5,8), mxBCs(3,8);
	CMatrix mxSab(8,8), vtSaby(8);

	double Sa = m_SPolyAInit;				// initial s of first poly
	double Sb = m_SPolyBInit;				// initial s of next poly
	double SaP, SbP;

	// set matrix 'mxdab'
	mxdab = 0;
	mxdab.elem(0  ,1  ) = 1;
	mxdab.elem(1  ,2  ) = 1;
	mxdab.elem(2  ,3  ) = 1;
	mxdab.elem(0  ,5  ) = 1;
	mxdab.elem(3  ,6  ) = 1;
	mxdab.elem(4  ,7  ) = 1;
	mxdab.elem(0  ,0  ) = -Sa;
	mxdab.elem(0  ,0+4) = -Sa;
	mxdab.elem(1  ,1+4) = 2*Sb;
	mxdab.elem(1+2,1+4) = -2*Sb;
	SaP = Sa*Sa;
	SbP = Sb*Sb;
	mxdab.elem(1  ,0  ) = -SaP;
	mxdab.elem(1  ,0+4) = -SaP-SbP;
	mxdab.elem(1+2,0+4) = SbP;
	mxdab.elem(2  ,1+4) = 3*SbP;
	mxdab.elem(2+2,1+4) = -3*SbP;
	SaP *= Sa;
	SbP *= Sb;
	mxdab.elem(2  ,0  ) = -SaP;
	mxdab.elem(2  ,0+4) = -SaP-2*SbP;
	mxdab.elem(2+2,0+4) = 2*SbP;

	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumASp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumBSp[r+c];
		}

	// set matrix 'mxBCs'
	SaP = SbP = 1;
	int c = 0;
	while (1)
	{
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+0) = SbP;
		mxBCs.elem(1, c+4) = -SbP;
		if (++c == 4) break;
		mxBCs.elem(2, c+0) = c * SbP;
		mxBCs.elem(2, c+4) = -c * SbP;
		SaP *= Sa;
		SbP *= Sb;
	}
	mxBCs.elem(2, 0) = 0;
	mxBCs.elem(2, 4) = 0;





double c88[8][3], c44[8][3], c44b[8][3];

/////////////////////////////////
// 8x8 matrix
//
{	
	CMatrix vtBCs(3);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	// set vector 'vtBCs'
	vtBCs[0] = vtPosInit[ax];			// Pos(Sa)
	vtBCs[1] = 0;
	vtBCs[2] = 0;
	
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	vtAllVals.ProdPart(mxdab, vtSaby);
	mxBCs.CopyToDestLoc(mxAllEqus, 5,0);
	vtBCs.CopyToDestLoc(vtAllVals, 5,0);
/*
afxDump << "mxdab: " << mxdab;
afxDump << "mxSab: " << mxSab;
afxDump << "vtSaby: " << vtSaby;
afxDump << "mxBCs: " << mxBCs;
afxDump << "vtBCs: " << vtBCs;

afxDump << "mxAllEqus: " << mxAllEqus;
afxDump << "vtAllVals: " << vtAllVals;
*/
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
/*
afxDump << "mxAllEqus2: " << mxAllEqus2;
afxDump << "vtAllVals2: " << vtAllVals2;
*/
	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

//afxDump << "vtab: " << vtab;

	CMatrix vtValSolu;
	vtValSolu = mxAllEqus2 * vtab;
	CMatrix vtDiff;
	vtDiff = vtValSolu - vtAllVals2;
	afxDump << "vtValSolu: " << vtValSolu;
	afxDump << "vtDiff: " << vtDiff;

	for (int i = 0; i < 4; i++)
	{
		m_PolyA[i][ax] = vtab[i];
		m_PolyB[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c88[i][ax] = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}








///////////////////////////////////////////////
// 4x4 matrix
//

// Calculate by changing vtab to vtAB with only 4 unknowns and rearranging to solve
// a 4x4 matrix equation
/*
	mxdab * mxSab * vtab  =  mxdab * vtSaby		- 5 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					- 3 equ's

	mxSab * mxAB2ab * vtAB  =  vtSaby				Incorporates BC's
	vtab = mxAB2ab * vtAB

	where mxAB2ab is inv(mxab2AB) and is obtained algebraicly from mxab2AB!
	Note:  mxab2AB * vtab = vtAB
	mxab2AB is derived from: mxBCs * vtab = vtBCs with ones added for extra vtab elements to make a square matrix

	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	[1HC] = 1st half of Columns etc.
	vtAB[2HR] is 2nd half of rows which are the only 4 unknowns

*/
{
for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	
	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	CMatrix mxAB2ab2HC(8,5);
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

	CMatrix mxdabS(5,8);
	CMatrix vtdabSy(5);
	CMatrix mxAllEqus(5,5);
	CMatrix vtAllVals(5);

	// mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	mxdabS.Prod(mxdab, mxSab);
	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

	CMatrix vtab(8);
	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax];					// Pos(Sa)

	CMatrix vtValDiff(5);
	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

	CMatrix vtAB2HR(5);
	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;

	for (int i = 0; i < 4; i++)
	{
		m_PolyA[i][ax] = vtab[i];
		m_PolyB[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c44[i][ax] = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}



// without mxdab *
{
for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	
	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	CMatrix mxAB2ab2HC(8,5);
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

//	CMatrix mxdabS(5,8);
//	CMatrix vtdabSy(5);
	CMatrix mxAllEqus(5,5);
	CMatrix vtAllVals(5);

	// without mxdab *
	// mxSab * mxAB2ab * vtAB  =  vtSaby
	// mxSab[2HR] * mxAB2ab[2HC] * vtAB[2HR]  =  vtSaby - mxSab[2HR] * mxAB2ab[1HC] * vtAB[1HR]
	CMatrix mxSab2HR(5,8);
	mxSab2HR.CopyFromSrcLoc(mxSab, 3,0);
//	mxdabS.Prod(mxdab, mxSab);
//	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxSab2HR, mxAB2ab2HC);

	CMatrix vtab(8);
	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax];					// Pos(Sa)

	CMatrix vtValDiff(5);
	vtValDiff.Prod(mxSab2HR, vtab);
	vtAllVals = vtSaby - vtValDiff;

	CMatrix vtAB2HR(5);
	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;


for (int i = 0; i < 8; i++)
	c44b[i][ax] = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}








}




