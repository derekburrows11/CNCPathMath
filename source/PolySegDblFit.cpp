// PolySegDblFit.cpp: implementation of the CPolySegDblFit class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "Matrix.h"
#include "PolyFunc.h"


#include "PolySegDblFit.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPolySegDblFit::CPolySegDblFit()
{
	Reset();


}

CPolySegDblFit::~CPolySegDblFit()
{
}


//////////////////////////////////////////////////////////////////
// Setup functions
//////////////////////////////////////////////////////////////////

void CPolySegDblFit::Reset()
{
	CPolySegFit::Reset();	// call base
//	m_numPtsStored = 0;
	m_FitA.Reset();
	m_FitB.Reset();
	m_pointRefPolyJoin = -1;


	// Set constant size matricies
	m_mxSab.SetSize(8,8);
	m_vtSaby.SetSize(8);
	m_mxab2Bez.SetSize(8,8);
	m_mxBez2ab.SetSize(8,8);
	m_mxSBez.SetSize(8,8);
	m_mxBCsBez2Bez.SetSize(8,8);

	m_mxBCVals.SetSize(3,8);
	m_vtBez.SetSize(8);
	m_vtab.SetSize(8);
	m_mxAllEqus.SetSize(8,8);
	m_mxAllEqus2.SetSize(8,8);
	m_vtAllVals.SetSize(8);
	m_vtAllVals2.SetSize(8);

	m_bez88.SetSize(3,8);		// columns of vectors
	m_bez44.SetSize(3,8);
	m_bez44b.SetSize(3,8);

}

//////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////

void CPolySegDblFit::SetCommonMatricies()			// for bezier solutions
{
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


////////////////////////////////////////

	With derivative beziers
	vtab  = [  a0   a1   a2   a3   b0   b1   b2   b3 ]'
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxab2Bez * vtab = vtBez
	mxBez2ab * vtBez = vtab

	mxab2Bez =
	[ 1        sa       sa^2           sa^3          0        0        0              0      ] ...[1]
	[ 0        1       2sa            3sa^2          0        0        0              0      ] ...[2]
	[ 0        1       2sb            3sb^2          0        0        0              0      ] ...[3]
	[ 1        sb       sb^2           sb^3          0        0        0              0      ] ...[4]
	[ 0        0        0              0             1        sb       sb^2           sb^3   ]
	[ 0        0        0              0             0        1       2sb            3sb^2   ]
	[ 0        0        0              0             0        1       2sc            3sc^2   ]
	[ 0        0        0              0             1        sc       sc^2           sc^3   ]

	Invert:
	[ 1        0      -sa^2          -2sa^3     ] ... [1] - sa[2]
	[ 0      sb-sa   sb^2-sa^2      sb^3-sa^3   ] ... [4] - [1]
	[ 0        0     2(sb-sa)     3(sb^2-sa^2)  ] ... [3] - [2]
	[ 1        0      -sb^2          -2sb^3     ] ... [4] - sb[3]
	[ 0        2     2(sb+sa)     3(sb^2+sa^2)  ] ... [2] + [3]

	[ 0   sb^2-sa^2   2sa.sb(sb-sa)      0      ] ... sb^2[2] - sa^2[3]
	[ 0      sb-sa      0        -3sa.sb(sb-sa) ] ... sb[2] - sa[3]



	[ 0        0     sb^2-sa^2              2(sb^3-sa^3)        ] ... [1] - sa[2] - ([4] - sb[3])
	[ 0        0   (sb-sa)(sb+sa)    2(sb-sa)(sb^2+sa.sb+sa^2)  ] ... [1] - [4] - sa[2] + sb[3]
a3:[ 0        0        0        4(sb-sa)(sb^2+sa.sb+sa^2) - 3(sb+sa)^2(sb-sa)  ] ... 2([1] - [4] - sa[2] + sb[3]) - (sb+sa)([3] - [2])
	[ 0        0        0        (sb-sa)(4(sb^2+sa.sb+sa^2) - 3(sb+sa)^2)  ] ... 2([1] - [4] - sa[2] + sb[3]) - (sb+sa)([3] - [2])
	[ 0        0        0        (sb-sa)(4sb^2+4sa.sb+4sa^2 -3sb^2 -3sa^2 -6sa.sb)  ] ... 2([1] - [4] - sa[2] + sb[3]) - (sb+sa)([3] - [2])
	[ 0        0        0        (sb-sa)(sb^2 -2sa.sb +sa^2)  ] ... 2([1] - [4] - sa[2] + sb[3]) - (sb+sa)([3] - [2])
	[ 0        0        0        (sb-sa)^3  ] ... 2([1] - [4]) + (sb-sa)([2] + [3])

  	[ 0        0     sb^2-sa^2                       2(sb^3-sa^3)       ] ... [1] - sa[2] - ([4] - sb[3])
	[ 0        0   (sb-sa)(sb+sa)             2(sb-sa)(sb^2+sa.sb+sa^2) ] ... [1] - [4] - sa[2] + sb[3]
a2:[ 0        0  3(sb-sa)(sb+sa)^2 - 4(sb^2+sa.sb+sa^2)(sb-sa)    0    ] ... 3(sb+sa)([1] - [4] - sa[2] + sb[3]) - 2(sb^2+sa.sb+sa^2)([3] - [2])
	[ 0        0   (sb-sa)(3(sb+sa)^2 - 4(sb^2+sa.sb+sa^2))    0    ] ... 3(sb+sa)([1] - [4]) + 3(sb+sa)( - sa[2] + sb[3]) - 2(sb^2+sa.sb+sa^2)([3] - [2])
	[ 0        0   (sb-sa)(3sb^2+3sa^2+6sa.sb - 4sb^2 -4sa.sb -4sa^2)    0    ] ... 3(sb+sa)([1] - [4]) + 3(sb^2+sa.sb)[3] - 3(sa^2+sa.sb)[2] - 2(sb^2+sa.sb+sa^2)([3] - [2])
	[ 0        0   (sb-sa)(-sb^2 +2sa.sb -sa^2)    0    ] ... 3(sb+sa)([1] - [4]) + (sb^2+sa.sb-2sa^2)[3] - (sa^2+sa.sb-2sb^2)[2]
	[ 0        0  -(sb-sa)^3    0    ] ... 3(sb+sa)([1] - [4]) + (sb^2+sa.sb+sa^2-3sa^2)[3] - (sa^2-sb^2+sa.sb-sb^2)[2]
	[ 0        0  -(sb-sa)^3    0    ] ... 3(sb+sa)([1] - [4]) + ((sb-sa)(sb+sa)+sa(sb-sa))[3] + ((sb-sa)(sb+sa)+sb(sb-sa))[2]
	[ 0        0  -(sb-sa)^3    0    ] ... 3(sb+sa)([1] - [4]) + (sb-sa)((2sb+sa)[2] + (sb+2sa)[3])

a1:[ 0    (sb-sa)^3(sb+sa)    0      0    ] ... 2sa.sb{3(sb+sa)([1] - [4]) + (sb-sa)((2sb+sa)[2] + (sb+2sa)[3])} + (sb-sa)^2(sb^2[2] - sa^2[3])
	[ 0    (sb-sa)^3(sb+sa)    0      0    ] ... 6sa.sb(sb+sa)([1] - [4]) + (sb-sa){(2sa.sb(2sb+sa)+(sb-sa)sb^2)[2] + (2sa.sb(sb+2sa)-(sb-sa)sa^2)[3]}
	[ 0    (sb-sa)^3(sb+sa)    0      0    ] ... 6sa.sb(sb+sa)([1] - [4]) + (sb-sa){sb(sb^2+3sa.sb+2sa^2)[2] + sa(sa^2+3sa.sb+2sb^2)[3]}
	[ 0    (sb-sa)^3(sb+sa)    0      0    ] ... 6sa.sb(sb+sa)([1] - [4]) + (sb-sa){sb(sb+sa)(sb+2sa)[2] + sa(sb+sa)(2sb+sa)[3]}
	[ 0    (sb-sa)^3           0      0    ] ... 6sa.sb([1] - [4]) + (sb-sa){sb(sb+2sa)[2] + sa(2sb+sa)[3]}

a0:[ (sb-sa)^3        0      -sa^2(sb-sa)^3   -2sa^3(sb-sa)^3     ] ... (sb-sa)^3([1] - sa[2])
	[ (sb-sa)^3        0      -sa^2(sb-sa)^3     0     ] ... (sb-sa)^3([1] - sa[2]) + 2sa^3{ 2([1] - [4]) + (sb-sa)([2] + [3]) }
	[ (sb-sa)^3        0            0            0     ] ... (sb-sa)^3([1] - sa[2]) + 2sa^3{ 2([1] - [4]) + (sb-sa)([2] + [3]) } - sa^2{ 3(sb+sa)([1] - [4]) + (sb-sa)((2sb+sa)[2] + (sb+2sa)[3]) }
	[ (sb-sa)^3        0            0            0     ] ... (sb-sa)^3([1] - sa[2]) + 4sa^3[1] - 4sa^3[4] + 2sa^3(sb-sa)[2] + 2sa^3(sb-sa)[3]  -  3sa^2(sb+sa)[1] + 3sa^2(sb+sa)[4]) - sa^2(sb-sa)(2sb+sa)[2] - sa^2(sb-sa)(sb+2sa)[3]
	[ (sb-sa)^3        0            0            0     ] ... (sb-sa)^3([1] - sa[2]) + (sa^3-3sb.sa^2)[1] + (3sb.sa^2 - sa^3)[4] + 2sa^3(sb-sa)[2] + 2sa^3(sb-sa)[3] - sa^2{ (sb-sa)((2sb+sa)[2] + (sb+2sa)[3]) }
	[ (sb-sa)^3        0            0            0     ] ... (sb-sa)^3([1] - sa[2]) + sa^2(sa-3sb)[1] + sa^2(3sb-sa)[4] + (sb-sa)(2sa^3-sa^2(2sb+sa))[2] + (sb-sa)(2sa^3-sa^2(sb+2sa))[3]
	[ (sb-sa)^3        0            0            0     ] ... (sb-sa)^3([1] - sa[2]) + sa^2(sa-3sb)[1] + sa^2(3sb-sa)[4] + (sb-sa)(sa^3-2sb.sa^2)[2] - sb.sa^2(sb-sa)[3]
	[ (sb-sa)^3        0            0            0     ] ... sb^2(sb-3sa)[1] + sa^2(3sb-sa)[4] - sa.sb^2(sb-sa)[2] - sb.sa^2(sb-sa)[3]



	mxBez2ab =
	[ sb^2(sb-3sa)   -(sb-sa)sa.sb^2    -(sb-sa)sb.sa^2      sa^2(3sb-sa)  0                   0                  0               0            ] / (sb-sa)^3
	[ 6sa.sb        sb(sb-sa)(sb+2sa)  sa(sb-sa)(2sb+sa)   -6sa.sb         0                   0                  0               0            ] / (sb-sa)^3
	[-3(sb+sa)       -(sb-sa)(2sb+sa)   -(sb-sa)(sb+2sa)    3(sb+sa)       0                   0                  0               0            ] / (sb-sa)^3
	[ 2               (sb-sa)            (sb-sa)           -2              0                   0                  0               0            ] / (sb-sa)^3
	[ 0                   0                  0              0              sc^2(sc-3sb)   -(sc-sb)sb.sc^2    -(sc-sb)sc.sb^2      sb^2(sb+3sc) ] / (sc-sb)^3
	[ 0                   0                  0              0              6sb.sc        sc(sc-sb)(sc+2sb)  sb(sc-sb)(2sc+sb)   -6sb.sc        ] / (sc-sb)^3
	[ 0                   0                  0              0             -3(sc+sb)       -(sc-sb)(2sc+sb)   -(sc-sb)(sc+2sb)    3(sc+sb)      ] / (sc-sb)^3
	[ 0                   0                  0              0              2               (sc-sb)            (sc-sb)           -2             ] / (sc-sb)^3


----------------------

						mxdab
	dSr/dA2   [ da0/dA2  da1/dA2  da2/dA2  da3/dA2 ]   dSra/da0     0
	dSr/dA3 = [ da0/dA3  da1/dA3  da2/dA3  da3/dA3 ] * dSra/da1  =  0
	                                                   dSra/da2
	                                                   dSra/da3
	elements of mxAB2ab are da(row)/dA(col) so transpose gives mxdab
	mxdab gives dSumRes/d(vtAB) from dSumRes/d(vtab)


*/
	m_arRowArrange = NULL;

	// All 8x8 matricies
	CMatrix& mxSab = m_mxSab;
	CMatrix& mxab2Bez = m_mxab2Bez;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& mxSBez = m_mxSBez;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double sc = m_FitB.m_SFinal;				// final s of next poly
	double SaP, SbP, ScP;

	// set matrix 'mxSab'
	// double Wa = 1;
	// double Wb = 1;
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = SumASp[r+c];			// * Wa ?? do to 'vtSaby' as well!
			mxSab.elem(r+4, c+4) = SumBSp[r+c];		// * Wb ?? do to 'vtSaby' as well!
		}

	// set matrix 'mxab2Bez'
	mxab2Bez = 0;
	SaP = SbP = ScP = 1;
	for (int c = 0; c < 4; c++)
	{
		mxab2Bez.elem(1, c+0) = c * SaP;
		mxab2Bez.elem(2, c+0) = c * SbP;
		mxab2Bez.elem(5, c+4) = c * SbP;
		mxab2Bez.elem(6, c+4) = c * ScP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
			ScP *= sc;
		}
		mxab2Bez.elem(0, c+0) = SaP;
		mxab2Bez.elem(3, c+0) = SbP;
		mxab2Bez.elem(4, c+4) = SbP;
		mxab2Bez.elem(7, c+4) = ScP;
	}

	// set matrix mxBez2ab
	mxBez2ab = 0;
	for (int group = 0; group < 2; group++)
	{
		double sA, sB;			// local selections
		if (group == 0)
		{
			sA = sa;
			sB = sb;
		}
		else
		{
			sA = sb;
			sB = sc;
		}

		double sBmsA = sB - sA;
		double OnsBmsA = 1 / sBmsA;
		double OnsBmsAp2 = OnsBmsA * OnsBmsA;
		double OnsBmsAp3 = OnsBmsAp2 * OnsBmsA;
		double sAp2 = sA*sA;
		double sBp2 = sB*sB;

		int br = group * 4;
		int bc = group * 4;
		mxBez2ab.elem(br+0,bc+0) =  OnsBmsAp3 * sBp2*(sB-3*sA);
		mxBez2ab.elem(br+0,bc+1) = -OnsBmsAp2 * sA*sBp2;
		mxBez2ab.elem(br+0,bc+2) = -OnsBmsAp2 * sB*sAp2;
		mxBez2ab.elem(br+0,bc+3) =  OnsBmsAp3 * sAp2*(3*sB-sA);

		mxBez2ab.elem(br+2,bc+0) = -OnsBmsAp3 * 3*(sB+sA);
		mxBez2ab.elem(br+2,bc+1) = -OnsBmsAp2 * (2*sB+sA);
		mxBez2ab.elem(br+2,bc+2) = -OnsBmsAp2 * (sB+2*sA);
		mxBez2ab.elem(br+2,bc+3) = -mxBez2ab.elem(br+2,bc+0);

		mxBez2ab.elem(br+1,bc+0) =  OnsBmsAp3 * 6*sA*sB;
		mxBez2ab.elem(br+1,bc+1) = -mxBez2ab.elem(br+2,bc+2) * sB;
		mxBez2ab.elem(br+1,bc+2) = -mxBez2ab.elem(br+2,bc+1) * sA;
		mxBez2ab.elem(br+1,bc+3) = -mxBez2ab.elem(br+1,bc+0);

		mxBez2ab.elem(br+3,bc+0) =  OnsBmsAp3 * 2;
		mxBez2ab.elem(br+3,bc+1) =  OnsBmsAp2;
		mxBez2ab.elem(br+3,bc+2) =  OnsBmsAp2;
		mxBez2ab.elem(br+3,bc+3) = -OnsBmsAp3 * 2;
	}

	// set matrix mxSBez
	mxSBez.Prod(mxSab, mxBez2ab);
}

/////////////////////

void CPolySegDblFit::Solve8by8()
{
// solve as 8x8 matrix
	return;

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBCVals = m_mxBCVals;
	CMatrix& mxAllEqus = m_mxAllEqus;
	CMatrix& mxAllEqus2 = m_mxAllEqus2;
	CMatrix& vtAllVals = m_vtAllVals;
	CMatrix& vtAllVals2 = m_vtAllVals2;
	CMatrix& vtBez = m_vtBez;

	mxAllEqus.ProdPart(m_mxdab, m_mxSBez);
	m_mxBCs.CopyToDestLoc(m_mxAllEqus, m_numDOFs,0);

	for (int ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax])		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];
			vtSaby[r+4] = SumBSpV[r][ax];
		}

		// set vector 'vtAllVals'
		vtAllVals.ProdPart(m_mxdab, vtSaby);
		int idxVal = m_numDOFs;
		for (int idxBC = 0; idxBC < m_numBCs; idxBC++)
			vtAllVals[idxVal++] = mxBCVals.elem(ax, idxBC);

//if (ax == 0)
//	afxDump << "mxAllEqus: " << mxAllEqus;
	// rearrange rows to work OK with LUSolve() !!!
		if (m_arRowArrange == NULL)
			mxAllEqus.LUSolve(vtBez.GetArray(), vtAllVals.GetArray());
		else
		{
			for (int rowD = 0; rowD < 8; rowD++)
			{
				int rowS = m_arRowArrange[rowD];
				for (int col = 0; col < 8; col++)
					mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
				vtAllVals2[rowD] = vtAllVals[rowS];
			}
			mxAllEqus2.LUSolve(vtBez.GetArray(), vtAllVals2.GetArray());
		}

		for (int i = 0; i < 4; i++)
		{
			m_FitA.m_BezDeriv[i][ax] = vtBez[i];
			m_FitB.m_BezDeriv[i][ax] = vtBez[i+4];
		}

		// Get Polys if required
		CMatrix& vtab = m_vtab;
		vtab.Prod(m_mxBez2ab, vtBez);
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez88.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez88.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez88: " << m_bez88;
}


////////////////////////////////////////////


void CPolySegDblFit::SolveDoubleCubic()
{
	LARGE_INTEGER liFreq;
	LARGE_INTEGER liStart;
	LARGE_INTEGER liEnd;
	VERIFY(QueryPerformanceFrequency(&liFreq));
	double countPeriod = 1e3 / liFreq.LowPart;		// in ms
	QueryPerformanceCounter(&liStart);




	int count = 0;
	double Kv0 = 1;
	double Kv1 = 1;

	SetDblCubicPoints();
	SetEndConditions();
	SetSValues();

	int nSegFitFlags = 0;
	if (m_nFitPosInit  != FIT_FREE) nSegFitFlags |= SFF_POSI;
	if (m_nFitPosFinal != FIT_FREE) nSegFitFlags |= SFF_POSF;
	if (m_nFitDirInit  != FIT_FREE) nSegFitFlags |= SFF_DIRI;
	if (m_nFitDirFinal != FIT_FREE) nSegFitFlags |= SFF_DIRF;

	for (;;)
	{
		SetSums();					// set m_PointSumsA, m_PointSumsB

		switch (nSegFitFlags)
		{
		case 0:
			FitDblCubicBez();
			break;
		case SFF_POSI:
			FitDblCubicIPosBez();
			break;
		case SFF_POSI | SFF_DIRI:
			FitDblCubicIPosBez();
			Kv0 = m_FitA.m_BezDeriv[1].Mag();			// keep mag of deriv same, just rotate to set direction
			FitDblCubicIPosDerivBez(Kv0);
			break;
		case SFF_POSI | SFF_POSF:
			FitDblCubicTwoPosBez();
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRI:
			FitDblCubicTwoPosBez();
			Kv0 = m_FitA.m_BezDeriv[1].Mag();			// keep mag of deriv same, just rotate to set direction
			FitDblCubicTwoPosIDerivBez(Kv0);
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRF:
			FitDblCubicTwoPosBez();
			Kv1 = m_FitB.m_BezDeriv[2].Mag();			// keep mag of deriv same, just rotate to set direction
			FitDblCubicTwoPosFDerivBez(Kv1);
			break;
		case SFF_POSI | SFF_POSF | SFF_DIRI | SFF_DIRF:
			FitDblCubicTwoPosBez();
			Kv0 = m_FitA.m_BezDeriv[1].Mag();			// keep mag of deriv same, just rotate to set direction
			Kv1 = m_FitB.m_BezDeriv[2].Mag();			// keep mag of deriv same, just rotate to set direction
			FitDblCubicTwoPosTwoDerivBez(Kv0, Kv1);
			break;
		default:
			ASSERT(0);			
		}

	
		//GetResiduals();
		AdjustSValues();

		count++;
		if (m_dsAvgNorm <= 1e-4)
			break;
		if (count >= 2)
			if (m_dsAvgNorm <= 1e-3)
				break;
			else if (count >= 8)
				if (m_dsAvgNorm <= 1e-2)
					break;
				else if (count >= 20)
					if (m_dsAvgNorm <= 5e-2)
						break;
					else
					{
						ASSERT(0);
						break;
					}
	}

	double maxA = m_FitA.m_MaxResidual;
	double maxB = m_FitB.m_MaxResidual;
	m_MaxResidual = max(maxA, maxB);

	if (m_MaxResidual >= 0.3)
		TRACE1("Large Maximum Residual in CPolySegDblFit::SolveDoubleCubic() of: %g\n", m_MaxResidual);

	QueryPerformanceCounter(&liEnd);
	double t1 = (liEnd.LowPart - liStart.LowPart) * countPeriod;
}

void CPolySegDblFit::SetInitialStateFromJoin()
{
	double SInit = m_FitB.GetSInit() - m_FitB.m_arS[0];			// get SInit relative to first point in FitB
	SetInitialPoint(m_FitA.m_BezDeriv[3], SInit);		// must set SInit too!
	SetInitialDir(m_FitA.m_BezDeriv[2]);
}

void CPolySegDblFit::SetPolyJoinPointRef(double pointRef)
{
	m_pointRefPolyJoin = pointRef;
}

void CPolySegDblFit::RemoveFitAPoints()
{
	RemoveFirstPoints(m_FitA.GetNumPoints());
}

void CPolySegDblFit::SetDblCubicPoints()
{
	int numPtsA = (int)ceil(m_pointRefPolyJoin);			// if join is on point, point is included in m_FitB
	int numPtsB = m_numPts - numPtsA;
	m_FitA.SetPointBuffer(m_arPts, numPtsA);
	m_FitB.SetPointBuffer(m_arPts + numPtsA, numPtsB);
	m_FitA.SetNumPointsInBuffer(numPtsA);
	m_FitB.SetNumPointsInBuffer(numPtsB);
}

void CPolySegDblFit::SetSValues()		// set s values for double cubic fit
{
	CPolySegFit::SetSValues();
	int numPtsA = m_FitA.GetNumPoints();
	int numPtsB = m_FitB.GetNumPoints();
	m_FitA.SetSBuffer(m_arS, numPtsA);
	m_FitB.SetSBuffer(m_arS + numPtsA, numPtsB);

	m_FitA.SetSInit(GetSInit());
	m_FitB.SetSFinal(GetSFinal());
	SetSPolyJoin();
}

void CPolySegDblFit::SetSPolyJoin()
{
	// Find S of poly join from normalised node pos
	int pointRefInt = (int)m_pointRefPolyJoin;
	double pointRefFract = m_pointRefPolyJoin - pointRefInt;
	double sPolyJoin = m_arS[pointRefInt] + pointRefFract * (m_arS[pointRefInt+1] - m_arS[pointRefInt]);
	m_FitA.SetSFinal(sPolyJoin);
	m_FitB.SetSInit(sPolyJoin);
}

void CPolySegDblFit::AdjustSValues()
{
	m_FitA.AdjustSValues();
	m_FitB.AdjustSValues();
	AdjustSPolyJoin();

	m_dsAvgNorm = m_FitA.m_dsAvgNorm + m_FitB.m_dsAvgNorm;
	m_MaxResidual = max(m_FitA.m_MaxResidual, m_FitB.m_MaxResidual);
	m_vtSumResidualSq = m_FitA.m_vtSumResidualSq + m_FitB.m_vtSumResidualSq;
	m_SumResidualSq = m_FitA.m_SumResidualSq + m_FitB.m_SumResidualSq;
	m_SumResidual = sqrt(m_SumResidualSq);
}

void CPolySegDblFit::AdjustSPolyJoin()
{
	SetSPolyJoin();		// just set for now
}

void CPolySegDblFit::GetResiduals()
{
	m_FitA.GetResiduals();
	m_FitB.GetResiduals();
	m_vtSumResidualSq = m_FitA.m_vtSumResidualSq + m_FitB.m_vtSumResidualSq;
	m_SumResidualSq = m_FitA.m_SumResidualSq + m_FitB.m_SumResidualSq;
	m_SumResidual = sqrt(m_SumResidualSq);
}

void CPolySegDblFit::SetSums()		// set num points and sums for double cubic fit
{
	// SetDblCubics() must be called first
	m_FitA.SetSums();
	m_FitB.SetSums();

	for (int ax = 0; ax < 3; ax++)
	{
		m_PointSums.m_bConstValue[ax] =
			(m_FitA.m_PointSums.m_bConstValue[ax] && m_FitB.m_PointSums.m_bConstValue[ax]		// constant value, don't need to solve!
			&& m_FitA.m_PointSums.m_Value[ax] == m_FitB.m_PointSums.m_Value[ax]);				// same value
		m_PointSums.m_Value[ax] = m_FitA.m_PointSums.m_Value[ax];
	}

}

void CPolySegDblFit::GetBeziers()
{
	m_FitA.GetBezierFromPoly();
	m_FitB.GetBezierFromPoly();
}

void CPolySegDblFit::SetAxisFitToConstValue(int ax)
{
	double val = m_FitA.m_PointSums.m_Value[ax];
	for (int i = 1; i < 4; i++)
	{
		m_FitA.m_BezDeriv[i][ax] = 0;
		m_FitB.m_BezDeriv[i][ax] = 0;
		m_FitA.m_Poly[i][ax] = 0;
		m_FitB.m_Poly[i][ax] = 0;
	}
	m_FitA.m_BezDeriv[0][ax] = val;
	m_FitB.m_BezDeriv[0][ax] = val;
	m_FitA.m_BezDeriv[3][ax] = val;
	m_FitB.m_BezDeriv[3][ax] = val;
	m_FitA.m_Poly[0][ax] = val;
	m_FitB.m_Poly[0][ax] = val;
}











void CPolySegDblFit::FitDblCubicBez()
{
/* Fit two parametric cubics with:
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]

	vtBCs =
	[     0   ]
	[     0   ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ 0  0  na0  da0  da1  na1  db1  nb1 ]'

	mxBez2BCsBez =
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	Invert:
	mxBCsBez2Bez =
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[-1        0        0        0       0        1        0        0    ]
	[ 0       -1        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 6 equ's
	                   mxBCs * vtBez  =  vtBCs					- 2 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 6 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 2 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 6 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 6 equ's

-------------------------------

	mxBCsBez2Bez[1HC] * vtBCsBez[1HR] =
	[ 0  0  0  0  0  0  0  0 ]'
*/

	SetCommonMatricies();

	int numBCs = 2;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,3) = 1;
	mxBCs.elem(0,4) = -1;
	mxBCs.elem(1,2) = 1;
	mxBCs.elem(1,5) = -1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,2) = 1;
	mxBCsBez2Bez.elem(1,3) = 1;
	mxBCsBez2Bez.elem(2,4) = 1;
	mxBCsBez2Bez.elem(3,5) = 1;
	mxBCsBez2Bez.elem(4,0) = -1;
	mxBCsBez2Bez.elem(4,5) = 1;
	mxBCsBez2Bez.elem(5,1) = -1;
	mxBCsBez2Bez.elem(5,4) = 1;
	mxBCsBez2Bez.elem(6,6) = 1;
	mxBCsBez2Bez.elem(7,7) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

//	char arRowArrange[] = {0, 1, 2, 3, 4, 5, 6, 7};
//	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = 0;
		m_mxBCVals.elem(ax,1) = 0;
	}
	Solve8by8();


/////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax])		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ 0  0  0  0  0  0  0  0 ]'
		CMatrix& vtab = m_vtab;

		m_mxReduEqus.LUSolve(vtBCsBez2HR.GetArray(), m_vtReduVals.GetArray());

		// vtBCsBez = [ 0  0  na0  da0  da1  na1  db1  nb1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[1][ax] = vtBCsBez2HR[1];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[2];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[3];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[3];
		m_FitB.m_BezDeriv[2][ax] = vtBCsBez2HR[4];
		m_FitB.m_BezDeriv[3][ax] = vtBCsBez2HR[5];;


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab = mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (int i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;
}


void CPolySegDblFit::FitDblCubicIPosBez()
{
/* Fit two parametric cubics with:
	Initial pos at Sa
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]

	vtBCs =
	[ Pos(sa) ]
	[     0   ]
	[     0   ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  0  0  da0  da1  na1  db1  nb1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[ 0       -1        0        0       0        1        0        0    ]
	[ 0        0       -1        0       1        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 5 equ's
	                   mxBCs * vtBez  =  vtBCs					- 3 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 5 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 3 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 5 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 5 equ's

-------------------------------

	mxBCsBez2Bez[1HC] * vtBCsBez[1HR] =
	[ Pos(sa)  0  0  0  0  0  0  0 ]'
*/

	SetCommonMatricies();

	int numBCs = 3;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,3) = 1;
	mxBCs.elem(1,4) = -1;
	mxBCs.elem(2,2) = 1;
	mxBCs.elem(2,5) = -1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,3) = 1;
	mxBCsBez2Bez.elem(2,4) = 1;
	mxBCsBez2Bez.elem(3,5) = 1;
	mxBCsBez2Bez.elem(4,1) = -1;
	mxBCsBez2Bez.elem(4,5) = 1;
	mxBCsBez2Bez.elem(5,2) = -1;
	mxBCsBez2Bez.elem(5,4) = 1;
	mxBCsBez2Bez.elem(6,6) = 1;
	mxBCsBez2Bez.elem(7,7) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {5, 0, 1, 2, 3, 4, 6, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = 0;
		m_mxBCVals.elem(ax,2) = 0;
	}
	Solve8by8();


/////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax])		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  0  0  0  0  0  0  0 ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax];			// * Pos(sa)
			vtab[4+i] = 0;
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		m_mxReduEqus.LUSolve(vtBCsBez2HR.GetArray(), m_vtReduVals.GetArray());

		// vtBCsBez = [ Pos(sa)  0  0  da0  da1  na1  db1  nb1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[1];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[2][ax] = vtBCsBez2HR[3];
		m_FitB.m_BezDeriv[3][ax] = vtBCsBez2HR[4];;


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;
}


void CPolySegDblFit::FitDblCubicIPosDerivBez(double Kv0)
{
/* Fit two parametric cubics with:
	Initial pos & derivative = Kv0 * m_vtDirInitUnit at Sa
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]

	vtBCs =
	[ Pos(sa) ]
	[ Vel(sa) ]
	[     0   ]
	[     0   ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  da1  na1  db1  nb1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[ 0        0       -1        0       0        1        0        0    ]
	[ 0        0        0       -1       1        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 4 equ's
	                   mxBCs * vtBez  =  vtBCs					- 4 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 4 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 4 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 4 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 4 equ's

-------------------------------

	mxBCsBez2Bez[1HC] * vtBCsBez[1HR] =
	[ Pos(sa)  Vel(sa)  0  0  0  0  0  0 ]'
*/

	SetCommonMatricies();

	int numBCs = 4;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,1) = 1;
	mxBCs.elem(2,3) = 1;
	mxBCs.elem(2,4) = -1;
	mxBCs.elem(3,2) = 1;
	mxBCs.elem(3,5) = -1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,1) = 1;
	mxBCsBez2Bez.elem(2,4) = 1;
	mxBCsBez2Bez.elem(3,5) = 1;
	mxBCsBez2Bez.elem(4,2) = -1;
	mxBCsBez2Bez.elem(4,5) = 1;
	mxBCsBez2Bez.elem(5,3) = -1;
	mxBCsBez2Bez.elem(5,4) = 1;
	mxBCsBez2Bez.elem(6,6) = 1;
	mxBCsBez2Bez.elem(7,7) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {4, 5, 0, 1, 2, 3, 6, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = vtDerivInit[ax];		// Vel(sa)
		m_mxBCVals.elem(ax,2) = 0;
		m_mxBCVals.elem(ax,3) = 0;
	}
	Solve8by8();


/////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax]
			&& vtDerivInit[ax] == 0)		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  Vel(sa)  0  0  0  0  0  0 ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax] + mxBez2ab.elem(0+i,1) * vtDerivInit[ax];	// * Pos(sa), * Vel(sa)
			vtab[4+i] = 0;
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		m_mxReduEqus.LUSolve(vtBCsBez2HR.GetArray(), m_vtReduVals.GetArray());

		// vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  da1  na1  db1  nb1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtDerivInit[ax];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[0];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[2][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[3][ax] = vtBCsBez2HR[3];;


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;
}


void CPolySegDblFit::FitDblCubicTwoPosBez()
{
/* Fit two parametric cubics with:
	Initial pos at Sa
	Final pos at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	vtBCs =
	[ Pos(sa) ]
	[     0   ]
	[     0   ]
	[ Pos(sc) ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  0  0  Pos(sc)  da0  da1  na1  db1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0       -1        0        0       0        0        1        0    ]
	[ 0        0       -1        0       0        1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        0        0        1       0        0        0        0    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 4 equ's
	                   mxBCs * vtBez  =  vtBCs					- 4 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 4 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 4 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 4 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 4 equ's

-------------------------------

	mxBCsBez2Bez[1HC] * vtBCsBez[1HR] =
	[ Pos(sa)  0  0  0  0  0  0  Pos(sc) ]'
*/

	SetCommonMatricies();

	int numBCs = 4;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,3) = 1;
	mxBCs.elem(1,4) = -1;
	mxBCs.elem(2,2) = 1;
	mxBCs.elem(2,5) = -1;
	mxBCs.elem(3,7) = 1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,4) = 1;
	mxBCsBez2Bez.elem(2,5) = 1;
	mxBCsBez2Bez.elem(3,6) = 1;
	mxBCsBez2Bez.elem(4,1) = -1;
	mxBCsBez2Bez.elem(4,6) = 1;
	mxBCsBez2Bez.elem(5,2) = -1;
	mxBCsBez2Bez.elem(5,5) = 1;
	mxBCsBez2Bez.elem(6,7) = 1;
	mxBCsBez2Bez.elem(7,3) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {4, 0, 1, 2, 3, 5, 6, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = 0;
		m_mxBCVals.elem(ax,2) = 0;
		m_mxBCVals.elem(ax,3) = vtPosFinal[ax];		// Pos(sc)
	}
	Solve8by8();


/////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtPosInit[ax] == vtPosFinal[ax])		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  0  0  0  0  0  0  Pos(sc) ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax];			// * Pos(sa)
			vtab[4+i] = mxBez2ab.elem(4+i,7) * vtPosFinal[ax];			// * Pos(sc)
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		m_mxReduEqus.LUSolve(vtBCsBez2HR.GetArray(), m_vtReduVals.GetArray());

		// vtBCsBez = [ Pos(sa)  0  0  Pos(sc)  da0  da1  na1  db1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[1];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[2][ax] = vtBCsBez2HR[3];
		m_FitB.m_BezDeriv[3][ax] = vtPosFinal[ax];


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;
}





void CPolySegDblFit::FitDblCubicTwoPosIDerivBez(double Kv0)
{
/* Fit two parametric cubics with:
	Initial pos & derivative = Kv0 * m_vtDirInitUnit at Sa
	Final pos at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	vtBCs =
	[ Pos(sa) ]
	[ Vel(sa) ]
	[     0   ]
	[     0   ]
	[ Pos(sc) ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  Pos(sc)  da1  na1  db1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0       -1        0       0        0        1        0    ]
	[ 0        0        0       -1       0        1        0        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        0        0        0       1        0        0        0    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 3 equ's
	                   mxBCs * vtBez  =  vtBCs					- 5 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 3 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 5 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 3 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 3 equ's

-------------------------------

	mxBCsBez2Bez[1HC]) * vtBCsBez[1HR] =
	[ Pos(sa)  Vel(sa)  0  0  0  0  0  Pos(sc) ]'
*/

	SetCommonMatricies();

	int numBCs = 5;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,1) = 1;
	mxBCs.elem(2,3) = 1;
	mxBCs.elem(2,4) = -1;
	mxBCs.elem(3,2) = 1;
	mxBCs.elem(3,5) = -1;
	mxBCs.elem(4,7) = 1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,1) = 1;
	mxBCsBez2Bez.elem(2,5) = 1;
	mxBCsBez2Bez.elem(3,6) = 1;
	mxBCsBez2Bez.elem(4,2) = -1;
	mxBCsBez2Bez.elem(4,6) = 1;
	mxBCsBez2Bez.elem(5,3) = -1;
	mxBCsBez2Bez.elem(5,5) = 1;
	mxBCsBez2Bez.elem(6,7) = 1;
	mxBCsBez2Bez.elem(7,4) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {3, 4, 6, 5, 0, 1, 2, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = vtDerivInit[ax];		// Vel(sa)
		m_mxBCVals.elem(ax,2) = 0;
		m_mxBCVals.elem(ax,3) = 0;
		m_mxBCVals.elem(ax,4) = vtPosFinal[ax];		// Pos(sc)
	}

	Solve8by8();


///////////////////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);
	VERIFY(m_mxReduEqus.Invert() > 1e-8);			// is a 3x3 matrix

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtPosInit[ax] == vtPosFinal[ax]
			&& vtDerivInit[ax] == 0)		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC] * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  Vel(sa)  0  0  0  0  0  Pos(sc) ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax] + mxBez2ab.elem(0+i,1) * vtDerivInit[ax];	// * Pos(sa), * Vel(sa)
			vtab[4+i] = mxBez2ab.elem(4+i,7) * vtPosFinal[ax];			// * Pos(sc)
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		vtBCsBez2HR.Prod(m_mxReduEqus, m_vtReduVals);

		// vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  Pos(sc)  da1  na1  db1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtDerivInit[ax];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[0];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[2][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[3][ax] = vtPosFinal[ax];


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;

}	// end CPolySegDblFit::FitDblCubicTwoPosIDerivBez(double Kv0)



void CPolySegDblFit::FitDblCubicTwoPosFDerivBez(double Kv1)
{
/* Fit two parametric cubics with:
	Initial pos at Sa
	Final pos & derivative = Kv0 * m_vtDirInitUnit at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	vtBCs =
	[ Pos(sa) ]
	[     0   ]
	[     0   ]
	[ Vel(sc) ]
	[ Pos(sc) ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  0  0  Vel(sc)  Pos(sc)  da0  da1  na1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0       -1        0        0       0        0        0        1    ]
	[ 0        0       -1        0       0        0        1        0    ]
	[ 0        0        0        1       0        0        0        0    ]
	[ 0        0        0        0       1        0        0        0    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 3 equ's
	                   mxBCs * vtBez  =  vtBCs					- 5 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 3 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 5 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 3 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 3 equ's

-------------------------------

	mxBCsBez2Bez[1HC]) * vtBCsBez[1HR] =
	[ Pos(sa)  0  0  0  0  0  Vel(sc)  Pos(sc) ]'
*/

	SetCommonMatricies();

	int numBCs = 5;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDerivFinal = m_vtDirFinalUnit * Kv1;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,3) = 1;
	mxBCs.elem(1,4) = -1;
	mxBCs.elem(2,2) = 1;
	mxBCs.elem(2,5) = -1;
	mxBCs.elem(3,6) = 1;
	mxBCs.elem(4,7) = 1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,5) = 1;
	mxBCsBez2Bez.elem(2,6) = 1;
	mxBCsBez2Bez.elem(3,7) = 1;
	mxBCsBez2Bez.elem(4,1) = -1;
	mxBCsBez2Bez.elem(4,7) = 1;
	mxBCsBez2Bez.elem(5,2) = -1;
	mxBCsBez2Bez.elem(5,6) = 1;
	mxBCsBez2Bez.elem(6,3) = 1;
	mxBCsBez2Bez.elem(7,4) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {3, 0, 5, 4, 1, 2, 6, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = 0;
		m_mxBCVals.elem(ax,2) = 0;
		m_mxBCVals.elem(ax,3) = vtDerivFinal[ax];		// Vel(sc)
		m_mxBCVals.elem(ax,4) = vtPosFinal[ax];		// Pos(sc)
	}

	Solve8by8();


///////////////////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);
	VERIFY(m_mxReduEqus.Invert() > 1e-8);			// is a 3x3 matrix

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtPosInit[ax] == vtPosFinal[ax]
			&& vtDerivFinal[ax] == 0)		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC] * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  0  0  0  0  0  Vel(sc)  Pos(sc) ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax];																// * Pos(sa)
			vtab[4+i] = mxBez2ab.elem(4+i,7) * vtPosFinal[ax] + mxBez2ab.elem(4+i,6) * vtDerivFinal[ax];		// * Pos(sc), Vel(sc)
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		vtBCsBez2HR.Prod(m_mxReduEqus, m_vtReduVals);

		// vtBCsBez = [ Pos(sa)  0  0  Vel(sc)  Pos(sc)  da0  da1  na1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[1];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[2];
		m_FitB.m_BezDeriv[2][ax] = vtDerivFinal[ax];
		m_FitB.m_BezDeriv[3][ax] = vtPosFinal[ax];


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;

}	// end CPolySegDblFit::FitDblCubicTwoPosFDerivBez(double Kv1)



void CPolySegDblFit::FitDblCubicTwoPosTwoDerivBez(double Kv0, double Kv1)
{
/* Fit two parametric cubics with:
	Initial pos & derivative = Kv0 * m_vtDirInitUnit at Sa
	Final pos at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
/*
	vtBez = [ na0  da0  da1  na1  nb0  db0  db1  nb1 ]'

	mxBCs * vtBez = vtBCs

	mxBCs =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]

	vtBCs =
	[ Pos(sa) ]
	[ Vel(sa) ]
	[     0   ]
	[     0   ]
	[ Vel(sc) ]
	[ Pos(sc) ]

------------------------------

	Arrange vtBCsBez from vtBCs to group unknows:

	vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  Vel(sc)  Pos(sc)  da1  na1 ]'

	mxBez2BCsBez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        1      -1        0        0        0    ]
	[ 0        0        1        0       0       -1        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        0        1        0       0        0        0        0    ]
	[ 0        0        0        1       0        0        0        0    ]

	Invert:
	mxBCsBez2Bez =
	[ 1        0        0        0       0        0        0        0    ]
	[ 0        1        0        0       0        0        0        0    ]
	[ 0        0        0        0       0        0        1        0    ]
	[ 0        0        0        0       0        0        0        1    ]
	[ 0        0       -1        0       0        0        0        1    ]
	[ 0        0        0       -1       0        0        1        0    ]
	[ 0        0        0        0       1        0        0        0    ]
	[ 0        0        0        0       0        1        0        0    ]

	mxBez2BCsBez * vtBez = vtBCsBez
	mxBCsBez2Bez * vtBCsBez = vtBez

	mxBez2BCsBez * mxab2Bez * vtab = vtBCsBez
	mxBez2ab * mxBCsBez2Bez * vtBCsBez = vtab

	mxBez2ab * mxBCsBez2Bez = mxBCsBez2ab
	mxdab is transpose of mxBCsBez2ab(columns of Bez DOF's - [2HC])

-----------------

	mxdab * mxSab * mxBez2ab * vtBez  =  mxdab * vtSaby		- 2 equ's
	                   mxBCs * vtBez  =  vtBCs					- 6 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez * vtBCsBez  =  mxdab * vtSaby		- 2 equ's
	                   mxBCs * mxBCsBez2Bez * vtBCsBez  =  vtBCs					- 6 equ's - redundent

	mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBez2ab * mxBCsBez2Bez[1HC]) * vtBCsBez[1HR]		- 2 equ's
	mxdab * mxSab * mxBCsBez2ab[2HC] * vtBCsBez[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxBCsBez2ab[1HC]) * vtBCsBez[1HR]		- 2 equ's

-------------------------------

	mxBCsBez2Bez[1HC]) * vtBCsBez[1HR] =
	[ Pos(sa)  Vel(sa)  0  0  0  0  Vel(sc)  Pos(sc) ]'
*/

	SetCommonMatricies();

	int numBCs = 6;
	int numDOFs = 8 - numBCs;
	m_numBCs = numBCs;
	m_numDOFs = numDOFs;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;
	CVector vtDerivFinal = m_vtDirFinalUnit * Kv1;

	CMatrix& mxBCs = m_mxBCs;
	CMatrix& mxBCsBez2Bez = m_mxBCsBez2Bez;
	CMatrix& mxBCsBez2ab = m_mxBCsBez2ab;
	CMatrix& mxdab = m_mxdab;
	CMatrix& mxdabT = m_mxdabT;

	// set matrix mxBCs - a part of Bez2BCsBez
	mxBCs.SetSize(m_numBCs,8);
	mxBCs = 0;
	mxBCs.elem(0,0) = 1;
	mxBCs.elem(1,1) = 1;
	mxBCs.elem(2,3) = 1;
	mxBCs.elem(2,4) = -1;
	mxBCs.elem(3,2) = 1;
	mxBCs.elem(3,5) = -1;
	mxBCs.elem(4,6) = 1;
	mxBCs.elem(5,7) = 1;

	// set matrix mxBCsBez2Bez
	mxBCsBez2Bez = 0;
	mxBCsBez2Bez.elem(0,0) = 1;
	mxBCsBez2Bez.elem(1,1) = 1;
	mxBCsBez2Bez.elem(2,6) = 1;
	mxBCsBez2Bez.elem(3,7) = 1;
	mxBCsBez2Bez.elem(4,2) = -1;
	mxBCsBez2Bez.elem(4,7) = 1;
	mxBCsBez2Bez.elem(5,3) = -1;
	mxBCsBez2Bez.elem(5,6) = 1;
	mxBCsBez2Bez.elem(6,4) = 1;
	mxBCsBez2Bez.elem(7,5) = 1;


	// set matrix mxBCsBez2Bez
	mxBCsBez2ab.Prod(m_mxBez2ab, mxBCsBez2Bez);


	// set matrix mxdab
	mxdabT.SetSize(8,numDOFs);
	mxdabT.CopyFromSrcLoc(mxBCsBez2ab, 0, numBCs);
	mxdabT.Transpose(mxdab);


/////////////////////////////////
// 8x8 matrix

	char arRowArrange[] = {2, 3, 5, 4, 0, 1, 6, 7};
	m_arRowArrange = arRowArrange;

	// set matrix 'm_mxBCVals'
	for (int ax = 0; ax < 3; ax++)
	{
		m_mxBCVals.elem(ax,0) = vtPosInit[ax];			// Pos(sa)
		m_mxBCVals.elem(ax,1) = vtDerivInit[ax];		// Vel(sa)
		m_mxBCVals.elem(ax,2) = 0;
		m_mxBCVals.elem(ax,3) = 0;
		m_mxBCVals.elem(ax,4) = vtDerivFinal[ax];		// Vel(sc)
		m_mxBCVals.elem(ax,5) = vtPosFinal[ax];		// Pos(sc)
	}

	Solve8by8();


///////////////////////////////////////////////
// Reduced DOF matrix

	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CMatrix& vtSaby = m_vtSaby;
	CMatrix& mxBez2ab = m_mxBez2ab;
	CMatrix& vtBCsBez2HR = m_vtBCsBez2HR;
	vtBCsBez2HR.SetSize(numDOFs);

	CMatrix& mxBCsBez2ab2HC = mxdabT;		// mxBCsBez2ab2HC is just transpose of mxdab!!
	m_mxdabS.Prod(mxdab, m_mxSab);
	m_mxReduEqus.Prod(m_mxdabS, mxBCsBez2ab2HC);
	VERIFY(m_mxReduEqus.Invert() > 1e-8);			// is a 2x2 matrix

	for (ax = 0; ax < 3; ax++)
	{
		if (m_PointSums.m_bConstValue[ax] && m_PointSums.m_Value[ax] == vtPosInit[ax] && vtPosInit[ax] == vtPosFinal[ax]
			&& vtDerivInit[ax] == 0 && vtDerivFinal[ax] == 0)		// constant value, don't need to solve!
		{
			SetAxisFitToConstValue(ax);
			continue;
		}

		// set vector 'vtSaby'
		for (int r = 0; r < 4; r++)
		{
			vtSaby[r] = SumASpV[r][ax];				// * Wa ??
			vtSaby[r+4] = SumBSpV[r][ax];				// * Wb ??
		}
		m_vtReduVals.Prod(mxdab, vtSaby);

		//	vtab(part) = mxBez2ab * mxBCsBez2Bez[1HC] * vtBCsBez[1HR]			- partly setup vtab
		// mxBCsBez2Bez[1HC] * vtBCsBez[1HR] = [ Pos(sa)  Vel(sa)  0  0  0  0  Vel(sc)  Pos(sc) ]'
		CMatrix& vtab = m_vtab;
		for (int i = 0; i < 4; i++)
		{
			vtab[0+i] = mxBez2ab.elem(0+i,0) * vtPosInit[ax]  + mxBez2ab.elem(0+i,1) * vtDerivInit[ax];		// * Pos(sa), * Vel(sa)
			vtab[4+i] = mxBez2ab.elem(4+i,7) * vtPosFinal[ax] + mxBez2ab.elem(4+i,6) * vtDerivFinal[ax];		// * Pos(sc), * Vel(sc)
		}
		m_vtValDiff.Prod(m_mxdabS, m_vtab);
		m_vtReduVals -= m_vtValDiff;

		vtBCsBez2HR.Prod(m_mxReduEqus, m_vtReduVals);

		// vtBCsBez = [ Pos(sa)  Vel(sa)  0  0  Vel(sc)  Pos(sc)  da1  na1 ]'
		m_FitA.m_BezDeriv[0][ax] = vtPosInit[ax];
		m_FitA.m_BezDeriv[1][ax] = vtDerivInit[ax];
		m_FitA.m_BezDeriv[2][ax] = vtBCsBez2HR[0];
		m_FitB.m_BezDeriv[1][ax] = vtBCsBez2HR[0];
		m_FitA.m_BezDeriv[3][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[0][ax] = vtBCsBez2HR[1];
		m_FitB.m_BezDeriv[2][ax] = vtDerivFinal[ax];
		m_FitB.m_BezDeriv[3][ax] = vtPosFinal[ax];


		// Get Polys if required
		//	vtab = mxBCsBez2ab[1HC] * vtBCsBez[1HR] + mxBCsBez2ab[2HC] * vtBCsBez[2HR]
		vtab += mxBCsBez2ab2HC * vtBCsBez2HR;		// vtab is now complete!
		for (i = 0; i < 4; i++)
		{
			m_FitA.m_Poly[i][ax] = vtab[i];
			m_FitB.m_Poly[i][ax] = vtab[i+4];
		}

		for (i = 0; i < 4; i++)
		{
			m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
			m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
		}
	}	// for (int ax = 0; ax < 3; ax++)
	m_FitA.GetBezierFromBezDeriv();
	m_FitB.GetBezierFromBezDeriv();
//	afxDump << "bez44: " << m_bez44;

}	// end CPolySegDblFit::FitDblCubicTwoPosTwoDerivBez(double Kv0, double Kv1)

















void CPolySegDblFit::FitDblCubicIPos()
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
	        [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]

	        [ Pos(sa) ]
	vtBCs = [     0   ]
	        [     0   ]


------------------------------

	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

-----------------

	Arrange vtAB from vtBCs to group unknows:

	vtAB = [ Pos(sa)  0  0  A1  A2  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxab2AB = [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	          [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	          [ 0        1        0        0        0        0        0        0    ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]

to invert:
	for a0    [ 1       sa        sa^2         sa^3     0        0       0        0    ]   [1]

	for b1    [ 0      -1       -2sb         -3sb^2     0       1       2sb      3sb^2 ]	-[3]

	for b0	 [-1      -sb       -sb^2        -sb^3     1       sb       sb^2     sb^3 ]	-[2]
				 [ 0    sa-sb   sa^2-sb^2    sa^3-sb^3     1       sb       sb^2     sb^3 ]	[1] - [2]
				 [ 0       sa   sa^2+sb^2   sa^3+2sb^3     1       0       -sb^2   -2sb^3 ]	[1] - [2] + sb*[3]


	          [ 1        0        0       -sa      -sa^2    -sa^3       0        0    ]
	mxAB2ab = [ 0        0        0        1        0        0          0        0    ]
	          [ 0        0        0        0        1        0          0        0    ]
	          [ 0        0        0        0        0        1          0        0    ]
	          [ 1       -1        sb      -sa -sa^2-sb^2  -sa^3-2sb^3   sb^2    2sb^3 ]
	          [ 0        0       -1        1       2sb      3sb^2     -2sb     -3sb^2 ]
	          [ 0        0        0        0        0        0          1        0    ]
	          [ 0        0        0        0        0        0          0        1    ]



  mxdab is transpose of mxAB2ab(:,4:8)

	mxdab = [ -sa      1     0     0        -sa        1       0    0 ]
	        [ -sa^2    0     1     0   -sa^2-sb^2     2sb      0    0 ]
	        [ -sa^3    0     0     1  -sa^3-2sb^3     3sb^2    0    0 ]
	        [  0       0     0     0         sb^2    -2sb      1    0 ]
	        [  0       0     0     0        2sb^3    -3sb^2    0    1 ]

part of mxdab * mxSab:
   0		0		0		0	 sb^2*Sb(s^0)-2sb  *Sb(s^1)+Sb(s^2)		 sb^2*Sb(s^1)-2sb  *Sb(s^2)+Sb(s^3)		 sb^2*Sb(s^2)-2sb  *Sb(s^3)+Sb(s^4)		 sb^2*Sb(s^3)-2sb  *Sb(s^4)+Sb(s^5)
	0		0		0		0	2sb^3*Sb(s^0)-3sb^2*Sb(s^1)+Sb(s^3)		2sb^3*Sb(s^1)-3sb^2*Sb(s^2)+Sb(s^4)		2sb^3*Sb(s^2)-3sb^2*Sb(s^3)+Sb(s^5)		2sb^3*Sb(s^3)-3sb^2*Sb(s^4)+Sb(s^6)

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

	OR reduce to 5x5 matrix

	mxdab * mxSab * mxAB2ab * vtAB  =  mxdab * vtSaby		- 5 equ's
	        mxBCs * mxAB2ab * vtAB  =  vtBCs					- 3 equ's - redundent

-----------------



*/

	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;

	CMatrix mxdab(5,8), mxBCs(3,8);
	CMatrix mxSab(8,8), vtSaby(8);

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double SaP, SbP;

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
	for (int c = 0; c < 4; c++)
	{
		mxBCs.elem(2, c+0) = c * SbP;
		mxBCs.elem(2, c+4) = -c * SbP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
		}
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+0) = SbP;
		mxBCs.elem(1, c+4) = -SbP;
	}

	
	// set matrix 'mxdab'
	mxdab = 0;
	mxdab.elem(0  ,1  ) = 1;
	mxdab.elem(1  ,2  ) = 1;
	mxdab.elem(2  ,3  ) = 1;
	mxdab.elem(0  ,5  ) = 1;
	mxdab.elem(3  ,6  ) = 1;
	mxdab.elem(4  ,7  ) = 1;
	mxdab.elem(0  ,0  ) = -sa;
	mxdab.elem(0  ,0+4) = -sa;
	mxdab.elem(1  ,1+4) = 2*sb;
	mxdab.elem(1+2,1+4) = -2*sb;
	SaP = sa*sa;
	SbP = sb*sb;
	mxdab.elem(1  ,0  ) = -SaP;
	mxdab.elem(1  ,0+4) = -SaP-SbP;
	mxdab.elem(1+2,0+4) = SbP;
	mxdab.elem(2  ,1+4) = 3*SbP;
	mxdab.elem(2+2,1+4) = -3*SbP;
	SaP *= sa;
	SbP *= sb;
	mxdab.elem(2  ,0  ) = -SaP;
	mxdab.elem(2  ,0+4) = -SaP-2*SbP;
	mxdab.elem(2+2,0+4) = 2*SbP;







/////////////////////////////////
// 8x8 matrix
//
{	
	CMatrix vtBCs(3);
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	mxBCs.CopyToDestLoc(mxAllEqus, 5,0);

	// set vector 'vtBCs'
	vtBCs[1] = 0;
	vtBCs[2] = 0;

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

	vtAllVals.ProdPart(mxdab, vtSaby);
	vtBCs.CopyToDestLoc(vtAllVals, 5,0);

// rearrange rows to work OK with LUSolve() !!!
	int arRowArrange[] = {5, 7, 6, 0, 1, 2, 3, 4};
	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez88.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez88.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

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
	CMatrix vtdabSy(5);
	CMatrix mxAllEqus(5,5);
	CMatrix vtAllVals(5);
	CMatrix mxdabS(5,8);
	CMatrix mxAB2ab2HC(8,5);
	CMatrix vtab(8);
	CMatrix vtValDiff(5);
	CMatrix vtAB2HR(5);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!
	// mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	mxdabS.Prod(mxdab, mxSab);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	vtdabSy.Prod(mxdab, vtSaby);

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax];					// Pos(Sa)

	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!


	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

}	// for (int ax = 0; ax < 3; ax++)
}



// without mxdab *
{
	CMatrix mxAllEqus(5,5);
	CMatrix vtAllVals(5);
	CMatrix mxAB2ab2HC(8,5);
	CMatrix mxSab2HR(5,8);
	CMatrix vtab(8);
	CMatrix vtValDiff(5);
	CMatrix vtSaby2HR(5);
	CMatrix vtAB2HR(5);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

	// without mxdab *
	// mxSab * mxAB2ab * vtAB  =  vtSaby
	// mxSab[2HR] * mxAB2ab[2HC] * vtAB[2HR]  =  vtSaby[2HR] - mxSab[2HR] * mxAB2ab[1HC] * vtAB[1HR]
	mxSab2HR.CopyFromSrcLoc(mxSab, 3,0);
//	mxdabS.Prod(mxdab, mxSab);
//	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxSab2HR, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax];					// Pos(Sa)

	vtValDiff.Prod(mxSab2HR, vtab);
	vtSaby2HR.CopyFromSrcLoc(vtSaby, 3,0);
	vtAllVals = vtSaby2HR - vtValDiff;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez44b.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez44b.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

}	// for (int ax = 0; ax < 3; ax++)
}

afxDump << "bez88: " << m_bez88;
afxDump << "bez44: " << m_bez44;
afxDump << "bez44b without mxdab*: " << m_bez44b;

}




void CPolySegDblFit::FitDblCubicIPosDeriv(double Kv0)
{
// Fit two parametric cubics with initial pos & derivative = Kv * m_vtDirInitUnit at Sa
// Pos and deriv are continuous at knot at Sb
// Each axis is solved independently
	

/*
	vtab = [ a0  a1  a2  a3  b0  b1  b2  b3 ]'

	        [  sa^2   -2sa     1    0     sa^2-sb^2     2sb-2sa     0    0 ]
	mxdab = [ 2sa^3   -3sa^2   0    1   2sa^3-2sb^3   3sb^2-3sa^2   0    0 ]
	        [  0        0      0    0          sb^2        -2sb     1    0 ]
	        [  0        0      0    0         2sb^3        -3sb^2   0    1 ]

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
	mxBCs = [ 0        1       2sa      3sa^2     0        0        0        0    ]
	        [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	        [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]

	        [ Pos(sa) ]
	vtBCs = [ Vel(sa) ]
	        [     0   ]
	        [     0   ]


	------------------------------
	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's



	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

	Arrange vtAB from vtBCs to group unknows:

	vtAB = [ Pos(sa)  Vel(sa)  0  0  A2  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxab2AB = [ 0        1       2sa      3sa^2     0        0        0        0    ]
	          [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	          [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]


to invert:
	for a0    [ 1        0       -sa^2       -2sa^3     0        0       0        0    ]   [1] - sa[2]

	for a1    [ 0        1       2sa          3sa^2     0        0       0        0    ]   [2]

	for b1    [ 0        0   2sa-2sb    3sa^2-3sb^2     0       1       2sb      3sb^2 ]	[2] - [4]

	for b0    [ 0    sa-sb   sa^2-sb^2    sa^3-sb^3     1       sb       sb^2     sb^3 ]	[1] - [3]
				 [ 0       sa   sa^2+sb^2   sa^3+2sb^3     1       0       -sb^2   -2sb^3 ]	[1] - [3] + sb[4]
				 [ 0        0   sb^2-sa^2  2sb^3-2sa^3     1       0       -sb^2   -2sb^3 ]	[1] - [3] + sb[4] - sa[2]




	          [ 1       -sa       0        0        sa^2    2sa^3       0        0    ]
	mxAB2ab = [ 0        1        0        0      -2sa     -3sa^2       0        0    ]
	          [ 0        0        0        0        1        0          0        0    ]
	          [ 0        0        0        0        0        1          0        0    ]
	          [ 1       -sa      -1        sb  sa^2-sb^2  2sa^3-2sb^3   sb^2    2sb^3 ]
	          [ 0        1        0       -1   2sb-2sa    3sb^2-3sa^2 -2sb     -3sb^2 ]
	          [ 0        0        0        0        0        0          1        0    ]
	          [ 0        0        0        0        0        0          0        1    ]

  mxdab is transpose of mxAB2ab(:,5:8)

	mxdab = [  sa^2  -2sa    1     0    sa^2-sb^2     2sb-2sa      0    0 ]
	        [ 2sa^3  -3sa^2  0     1  2sa^3-2sb^3     3sb^2-3sa^2  0    0 ]
	        [  0       0     0     0         sb^2    -2sb          1    0 ]
	        [  0       0     0     0        2sb^3    -3sb^2        0    1 ]


*/

	double Wa = 1;				// ??
	double Wb = 1;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;


	CMatrix mxdab(4,8), mxSab(8,8), mxBCs(4,8);
	CMatrix vtSaby(8), vtBCs(4);

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double SaP, SbP;


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
	for (int c = 0; c < 4; c++)
	{
		mxBCs.elem(1, c+0) = c * SaP;
		mxBCs.elem(3, c+0) = c * SbP;
		mxBCs.elem(3, c+4) = -c * SbP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
		}
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+4) = 0;
		mxBCs.elem(2, c+0) = SbP;
		mxBCs.elem(2, c+4) = -SbP;
	}


	// set matrix 'mxdab'
	mxdab = 0;
	mxdab.elem(0  ,1  ) = -2*sa;
	mxdab.elem(0+2,1+4) = -2*sb;
	mxdab.elem(0  ,1+4) = -2*(sa-sb);
	SaP = sa*sa;
	SbP = sb*sb;
	mxdab.elem(0  ,0  ) = SaP;
	mxdab.elem(0+2,0+4) = SbP;
	mxdab.elem(0  ,0+4) = (SaP-SbP);
	mxdab.elem(1  ,1  ) = -3*SaP;
	mxdab.elem(1+2,1+4) = -3*SbP;
	mxdab.elem(1  ,1+4) = -3*(SaP-SbP);
	SaP *= sa;
	SbP *= sb;
	mxdab.elem(1  ,0  ) = 2*SaP;
	mxdab.elem(1+2,0+4) = 2*SbP;
	mxdab.elem(1  ,0+4) = 2*(SaP-SbP);
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(0+2,2+4) = 1;
	mxdab.elem(1+2,3+4) = 1;






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

	// rearrange rows
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);
	int arRowArrange[] = {4, 5, 7, 6, 0, 2, 1, 3};
	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez88.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez88.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

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

	CMatrix mxAllEqus(4,4);
	CMatrix vtAllVals(4);
	CMatrix mxAB2ab2HC(8,4);
	
	CMatrix vtab(8);
	CMatrix vtValDiff(4);
	CMatrix vtAB2HR(4);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

{
for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}

	CMatrix mxdabS(4,8);
	CMatrix vtdabSy(4);

	mxdabS.Prod(mxdab, mxSab);
	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = vtPosInit[ax] - sa*vtDerivInit[ax];		// Pos(Sa) - sa*Vel(Sa)
	vtab[1] = vtab[5] = vtDerivInit[ax];								// Vel(Sa)

	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

}	// for (int ax = 0; ax < 3; ax++)
}

afxDump << "bez88: " << m_bez88;
afxDump << "bez44: " << m_bez44;

}









void CPolySegDblFit::FitDblCubicTwoPos()
{
/* Fit two parametric cubics with:
	Initial pos at Sa
	Final pos at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
	
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
	        [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	        [ 0        0        0        0        1        sc       sc^2     sc^3 ]

	        [ Pos(sa) ]
	vtBCs = [     0   ]
	        [     0   ]
	        [ Pos(sc) ]


------------------------------

	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

-----------------

	Arrange vtAB from vtBCs to group unknows:

	vtAB = [ Pos(sa)  0  0  Pos(sc)  A2  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ]
	mxab2AB = [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]
	          [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	          [ 0        0        0        0        1        sc       sc^2     sc^3 ]
	          [ 0        0        1        0        0        0        0        0    ]
	          [ 0        0        0        1        0        0        0        0    ]
	          [ 0        0        0        0        0        0        1        0    ]
	          [ 0        0        0        0        0        0        0        1    ]
	to invert:
	          [ 0      sa-sb  sa^2-sb^2  sa^3-sb^3  1        sb        sb^2      sb^3 ]  ...[1]-[2]  ...([x] => from vtAB[x])
	          [ 1      0  sa(sa-2sb) sa(sa^2-3sb^2) 0        sa       2sa.sb  3sa.sb^2]  ...[1]-sa[3]

	          [ 1      sb       sb^2     sb^3       0     sc-sb   sc^2-sb^2  sc^3-sb^3]  ...[2]+[4]
	          [sa-sb  0  sa.sb(sb-sa) sa.sb(sb^2-sa^2) 0  sa(sc-sb) sa(sc^2-sb^2) sa(sc^3-sb^3)]  ...sa([2]+[4])-sb[1]
	          [sa-sc  0  sa(sb(sb-sa)-(sc-sb)(sa-2sb))  sa(sb(sb^2-sa^2)-(sc-sb)(sa^2-3sb^2))    0     0   sa((sc^2-sb^2)-2sb(sc-sb))  sa((sc^3-sb^3)-3sb^2(sc-sb))]  ...sa([2]+[4])-sb[1] - (sc-sb)([1]-sa[3])
	          [sa-sc  0  sa(-sa.sc+2sb.sc-sb^2)  sa(-sc.sa^2+3sc.sb^2-2sb^3)    0     0   sa(sc^2-2sb.sc+sb^2)  sa(sc^3-3sb^2.sc+2sb^3))]  ...sa([2]+[4])-sb[1] - (sc-sb)([1]-sa[3])
	          [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc(sc^2-sb^2)+2sb^2(sb-sc))]  ...sa([2]+[4])-sb[1] - (sc-sb)([1]-sa[3])
	for a0    [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([2]+[4])-sb[1] - (sc-sb)([1]-sa[3]) ==> sa([2]+[4]) - sc[1] + sa(sc-sb)[3]

	          [ 0  sa(sc-sa)  sa(sc(sb-sa)+sb(sc-sb)+sa(sc-sa))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([2]+[4]) - sc[1] + sa(sc-sb)[3] - (sa-sc)[1] ==> sa([2]+[4]) - sa[1] + sa(sc-sb)[3]
	          [ 0  sc-sa  sc(sb-sa)+sb(sc-sb)+sa(sc-sa)  sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa)    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[2] + [4] - [1] + (sc-sb)[3]
	for a1    [ 0  sc-sa  sb(2sc-sb)-sa^2  sb^2(3sc-2sb)-sa^3    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[2] + [4] - [1] + (sc-sb)[3]

	          [ 0     0  sb(2sc-sb)-sa^2-2sb(sc-sa)  sb^2(3sc-2sb)-3sb^2(sc-sa)-sa^3    0     sc-sa   (sc-sb)^2+2sb(sc-sa)  (sc+2sb)(sc-sb)^2+3sb^2(sc-sa)]  ...[2] + [4] - [1] + (sc-sb)[3] - (sc-sa)[3] ==> [2] + [4] - [1] + (sa-sb)[3]
	for b1    [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3                 0     sc-sa    sb(sb-2sa)+sc^2       sb^2(2sb-3sa)+sc^3]  ...[2] + [4] - [1] + (sa-sb)[3]

	          [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sc^3-sc^2(sc-sa)       sc.sb^2(2sb-3sa)+sc^4-sc^3(sc-sa)]  ...sc([2] + [4] - [1] + (sa-sb)[3]) - (sc-sa)[4] ==> sc[2] + sa[4] - sc[1] + sc(sa-sb)[3])
	for b0    [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sa.sc^2       sc.sb^2(2sb-3sa)+sa.sc^3]  ...sc[2] + sa[4] - sc[1] + sc(sa-sb)[3])


	          [ sc      -sa  -sa(sc-sb)   -sa  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))     sa(sc-sb)^2        sa(sc+2sb)(sc-sb)^2   ] / (sc-sa)
	mxAB2ab = [-1        1       sc-sb     1     sa^2-sb(2sc-sb)        sa^3-sb^2(3sc-2sb)                  -(sc-sb)^2         -(sc+2sb)(sc-sb)^2   ] / (sc-sa)
	          [ 0        0        0        0        1                      0                                     0                       0          ]
	          [ 0        0        0        0        0                      1                                     0                       0          ]
	          [ sc      -sc   sc(sb-sa)   -sa    -sc(sa-sb)^2           sc(sb^2(3sa-2sb)-sa^3)      sb.sc(sb-2sa)+sa.sc^2   sc.sb^2(2sb-3sa)+sa.sc^3] / (sc-sa)
	          [-1        1       sa-sb     1      (sa-sb)^2             sa^3-sb^2(3sa-2sb)             sb(2sa-sb)-sc^2         sb^2(3sa-2sb)-sc^3   ] / (sc-sa)
	          [ 0        0        0        0        0                      0                                     1                       0          ]
	          [ 0        0        0        0        0                      0                                     0                       1          ]



  mxdab is transpose of mxAB2ab(:,5:8)

	        [  sa(sc(sb-sa)+sb(sc-sb))          sa^2-sb(2sc-sb)        1     0  -sc(sb-sa)^2                           (sb-sa)^2   0    0 ]
	mxdab = [  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))   sa^3-sb^2(3sc-2sb)     0     1  -sc(sb^2(2sb-3sa)+sa^3)       sb^2(2sb-3sa)+sa^3   0    0 ]
	        [  sa(sc-sb)^2                              -(sc-sb)^2     0     0   sc(sb  ( sb-2sa)+sa.sc)        -sb(sb-2sa)-sc^2   1    0 ]
	        [  sa(sc+2sb)(sc-sb)^2              -(sc+2sb)(sc-sb)^2     0     0   sc(sb^2(2sb-3sa)+sa.sc^2)   -sb^2(2sb-3sa)-sc^3   0    1 ]
	               /(sc-sa)                             /(sc-sa)                       /(sc-sa)                     /(sc-sa)

Regular expersion method:
	        [  sa(2sb.sc-sa.sc-sb^2)          sa^2-sb.sc-sb(sc-sb)        1     0  -sc(sb^1(1(sb-sa)-sa)+sa^2)        sb^1(1(sb-sa)-sa)+sa^2    0    0 ]
	mxdab = [  sa(3sc.sb^2-sc.sa^2-2sb^3)     sa^3-sc.sb^2-2sb^2(sc-sb)   0     1  -sc(sb^2(2(sb-sa)-sa)+sa^3)        sb^2(2(sb-sa)-sa)+sa^3    0    0 ]
	        [  sa(sc^2-2sb.sc+sb^2)             -(sc-sb)^2                0     0   sc(sb^1(1(sb-sa)-sa)+sa.sc^1)    -sb^1(1(sb-sa)-sa)-sc^2    1    0 ]
	        [  sa(sc^3-3sc.sb^2+2sb^3)           -(sc+2sb)(sc-sb)^2       0     0   sc(sb^2(2(sb-sa)-sa)+sa.sc^2)    -sb^2(2(sb-sa)-sa)-sc^3    0    1 ]
	               /(sc-sa)                             /(sc-sa)                       /(sc-sa)                           /(sc-sa)

Best Calculating method:
	        [  sa(sc(sb-sa)+sb(sc-sb))            -sc(sb-sa)-sb(sc-sb)-sa(sc-sa)              1     0  -sc(sb-sa)^2                       (sb-sa)^2                                0    0 ]
	mxdab = [  sa(sc(sb+sa)(sb-sa)+2sb^2(sc-sb))  -sc(sb+sa)(sb-sa)-2sb^2(sc-sb)-sa^2(sc-sa)  0     1  -sc(2sb^2(sb-sa)-sa(sb^2-sa^2))    2sb^2(sb-sa)-sa(sb^2-sa^2)               0    0 ]
	        [  sa(sc-sb)^2                                    -(sc-sb)^2                      0     0   sc(sb(sb-sa)+sa(sc-sb))          -sb(sb-sa)-sa(sc-sb)-sc(sc-sa)            1    0 ]
	        [  sa(sc+2sb)(sc-sb)^2                    -(sc+2sb)(sc-sb)^2                      0     0   sc(2sb^2(sb-sa)+sa(sc^2-sb^2))   -2sb^2(sb-sa)-sa(sc^2-sb^2)-sc^2(sc-sa)   0    1 ]
	               /(sc-sa)                             /(sc-sa)                                                 /(sc-sa)                           /(sc-sa)


						mxdab
	dSr/dA2   [ da0/dA2  da1/dA2  da2/dA2  da3/dA2 ]   dSra/da0     0
	dSr/dA3 = [ da0/dA3  da1/dA3  da2/dA3  da3/dA3 ] * dSra/da1  =  0
	                                                   dSra/da2
	                                                   dSra/da3
	elements of mxAB2ab are da(row)/dA(col) so transpose gives mxdab
	mxdab gives dSumRes/d(vtAB) from dSumRes/d(vtab)

-----------------
	mxdab * mxSab * vtab  =  mxdab * vtSaby		- 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					- 4 equ's

	OR reduce to 4x4 matrix

	mxdab * mxSab * mxAB2ab * vtAB  =  mxdab * vtSaby		- 4 equ's
	        mxBCs * mxAB2ab * vtAB  =  vtBCs					- 4 equ's - redundent

	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxAB2ab[1HC]) * vtAB[1HR]			- 4 equ's

-----------------



*/

	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;

	int numBCs = 4;
	int numDOFs = 8 - numBCs;

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double sc = m_FitB.m_SFinal;				// final s of next poly
	double SaP, SbP, ScP;

	CMatrix mxSab(8,8), vtSaby(8);
	CMatrix mxdab(numDOFs,8);
	CMatrix mxBCs(numBCs,8);

	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumASp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumBSp[r+c];
		}


	// set matrix 'mxBCs'
	SaP = SbP = ScP = 1;
	for (int c = 0; c < 4; c++)
	{
		mxBCs.elem(2, c+0) = c * SbP;
		mxBCs.elem(2, c+4) = -c * SbP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
			ScP *= sc;
		}
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(1, c+0) = SbP;
		mxBCs.elem(1, c+4) = -SbP;
		mxBCs.elem(3, c+0) = 0;
		mxBCs.elem(3, c+4) = ScP;
	}


	// set matrix 'mxdab'
	double scmsa = sc - sa;
	double scmsb = sc - sb;
	double sbmsa = sb - sa;
	double scpsa = sc + sa;
	double scpsb = sc + sb;
	double sbpsa = sb + sa;

	double Onscmsa = 1 / scmsa;
	double saOnscmsa = sa / scmsa;
	double scOnscmsa = sc / scmsa;

	mxdab = 0;
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(2  ,6  ) = 1;
	mxdab.elem(3  ,7  ) = 1;

	double base = sc*sbmsa + sb*scmsb;
	mxdab.elem(0  ,0  ) = saOnscmsa*base;
	mxdab.elem(0  ,1  ) = -Onscmsa*(base + sa*scmsa);

	base = sc*sbpsa*sbmsa + 2*sb*sb*scmsb;
	mxdab.elem(1  ,0  ) = saOnscmsa*base;
	mxdab.elem(1  ,1  ) = -Onscmsa*(base + sa*sa*scmsa);

	mxdab.elem(2  ,1  ) = -Onscmsa*scmsb*scmsb;
	mxdab.elem(2  ,0  ) = -sa*mxdab.elem(2,1);

	mxdab.elem(3  ,1  ) = (sc+2*sb) * mxdab.elem(2,1);
	mxdab.elem(3  ,0  ) = -sa*mxdab.elem(3,1);
	
	mxdab.elem(0  ,1+4) = Onscmsa*sbmsa*sbmsa;
	mxdab.elem(0  ,0+4) = -sc*mxdab.elem(0,1+4);

	double base1 = 2*sb*sb*sbmsa;
	mxdab.elem(1  ,1+4) = Onscmsa*(base1 - sa*sbpsa*sbmsa);
	mxdab.elem(1  ,0+4) = -sc*mxdab.elem(1,1+4);

	double base2 = sb*sbmsa + sa*scmsb;
	mxdab.elem(2  ,0+4) =  scOnscmsa*base2;
	mxdab.elem(2  ,1+4) = -Onscmsa*(base2 + sc*scmsa);

	base1 += sa*scpsb*scmsb;
	mxdab.elem(3  ,0+4) =  scOnscmsa*base1;
	mxdab.elem(3  ,1+4) = -Onscmsa*(base1 + sc*sc*scmsa);
	








/////////////////////////////////
// 8x8 matrix
//
{	
	CMatrix vtBCs(numBCs);
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	mxBCs.CopyToDestLoc(mxAllEqus, numDOFs,0);

	// set vector 'vtBCs'
	vtBCs[1] = 0;
	vtBCs[2] = 0;

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
	vtBCs[3] = vtPosFinal[ax];			// Pos(Sc)

	vtAllVals.ProdPart(mxdab, vtSaby);
	vtBCs.CopyToDestLoc(vtAllVals, numDOFs,0);

// rearrange rows to work OK with LUSolve() !!!
	int arRowArrange[] = {4, 5, 6, 0, 1, 2, 3, 7};

	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez88.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez88.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

}	// for (int ax = 0; ax < 3; ax++)
}








///////////////////////////////////////////////
// 4x4 matrix
//

// Calculate by changing vtab to vtAB with only 4 unknowns and rearranging to solve
// a 4x4 matrix equation

//	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - (mxdab * mxSab * mxAB2ab[1HC]) * vtAB[1HR]			- 4 equ's


	CMatrix mxAllEqus(numDOFs,numDOFs);
	CMatrix vtAllVals(numDOFs);
	CMatrix mxAB2ab2HC(8,numDOFs);
	
	CMatrix vtab(8);
	CMatrix vtValDiff(numDOFs);
	CMatrix vtAB2HR(numDOFs);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

{
	CMatrix vtdabSy(numDOFs);
	CMatrix mxdabS(numDOFs,8);

	// mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	mxdabS.Prod(mxdab, mxSab);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	vtdabSy.Prod(mxdab, vtSaby);

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = Onscmsa * (sc*vtPosInit[ax] - sa*vtPosFinal[ax]);			// Pos(sa) & Pos(sc)
	vtab[1] = vtab[5] = Onscmsa * (  -vtPosInit[ax] + vtPosFinal[ax]);				// Pos(sa) & Pos(sc)

	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

	m_FitA.GetBezDerivFromPoly();
	m_FitB.GetBezDerivFromPoly();
	for (i = 0; i < 4; i++)
	{
		m_bez44.elem(ax,0+i) = m_FitA.m_BezDeriv[i][ax];
		m_bez44.elem(ax,4+i) = m_FitB.m_BezDeriv[i][ax];
	}

}	// for (int ax = 0; ax < 3; ax++)
}

afxDump << "bez88: " << m_bez88;
afxDump << "bez44: " << m_bez44;

}	// end CPolySegDblFit::FitDblCubicTwoPos()




void CPolySegDblFit::FitDblCubicTwoPosIDeriv(double Kv0)
{
/* Fit two parametric cubics with:
	Initial pos & derivative = Kv0 * m_vtDirInitUnit at Sa
	Final pos at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/
	
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
	mxBCs = [ 0        1       2sa      3sa^2     0        0        0        0    ]
	        [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]	       
	        [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	        [ 0        0        0        0        1        sc       sc^2     sc^3 ]

	        [ Pos(sa) ]
	vtBCs = [ Vel(sa) ]
	        [     0   ]
	        [     0   ]
	        [ Pos(sc) ]


------------------------------

	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

-----------------

	Arrange vtAB from vtBCs to group unknows:

	vtAB = [ Pos(sa)  Vel(sa)  0  0  Pos(sc)  A3  B2  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ] ...1  prev 1
	mxab2AB = [ 0        1       2sa      3sa^2     0        0        0        0    ] ...2
	          [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ] ...3  prev 2
	          [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ] ...4  prev 3
	          [ 0        0        0        0        1        sc       sc^2     sc^3 ] ...5  prev 4
	          [ 0        0        0        1        0        0        0        0    ] ...6
	          [ 0        0        0        0        0        0        1        0    ] ...7
	          [ 0        0        0        0        0        0        0        1    ] ...8

	to invert:
	          [ 0      sa-sb  sa^2-sb^2  sa^3-sb^3  1        sb        sb^2      sb^3 ]  ...[1]-[3]  ...([x] => from vtAB[x])
	          [ 1      0  sa(sa-2sb) sa(sa^2-3sb^2) 0        sa       2sa.sb  3sa.sb^2]  ...[1]-sa[4]


	          [ 1      sb       sb^2     sb^3       0     sc-sb   sc^2-sb^2  sc^3-sb^3]  ...[3]+[5]
	          [sa-sb  0  sa.sb(sb-sa)                   sa.sb(sb^2-sa^2)                         0   sa(sc-sb)  sa(sc^2-sb^2)                sa(sc^3-sb^3)]  ...sa([3]+[5])-sb[1]
	          [sa-sc  0  sa(sb(sb-sa)-(sc-sb)(sa-2sb))  sa(sb(sb^2-sa^2)-(sc-sb)(sa^2-3sb^2))    0     0        sa((sc^2-sb^2)-2sb(sc-sb))   sa((sc^3-sb^3)-3sb^2(sc-sb))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	          [sa-sc  0  sa(-sa.sc+2sb.sc-sb^2)  sa(-sc.sa^2+3sc.sb^2-2sb^3)    0     0   sa(sc^2-2sb.sc+sb^2)  sa(sc^3-3sb^2.sc+2sb^3))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	          [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc(sc^2-sb^2)+2sb^2(sb-sc))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	for a0    [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4]) ==> sa([3]+[5]) - sc[1] + sa(sc-sb)[4]

	          [ 0  sa(sc-sa)  sa(sc(sb-sa)+sb(sc-sb)+sa(sc-sa))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([3]+[5]) - sc[1] + sa(sc-sb)[4] - (sa-sc)[1] ==> sa([3]+[5]) - sa[1] + sa(sc-sb)[4]
	          [ 0  sc-sa  sc(sb-sa)+sb(sc-sb)+sa(sc-sa)  sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa)    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]
	for a1    [ 0  sc-sa  sb(2sc-sb)-sa^2  sb^2(3sc-2sb)-sa^3    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]

	          [ 0     0  sb(2sc-sb)-sa^2-2sb(sc-sa)  sb^2(3sc-2sb)-3sb^2(sc-sa)-sa^3    0     sc-sa   (sc-sb)^2+2sb(sc-sa)  (sc+2sb)(sc-sb)^2+3sb^2(sc-sa)]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[4] ==> [3] + [5] - [1] + (sa-sb)[4]
	for b1    [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3                 0     sc-sa    sb(sb-2sa)+sc^2       sb^2(2sb-3sa)+sc^3]  ...[3] + [5] - [1] + (sa-sb)[4]

	          [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sc^3-sc^2(sc-sa)       sc.sb^2(2sb-3sa)+sc^4-sc^3(sc-sa)]  ...sc([3] + [5] - [1] + (sa-sb)[4]) - (sc-sa)[5] ==> sc[3] + sa[5] - sc[1] + sc(sa-sb)[4])
	for b0    [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sa.sc^2       sc.sb^2(2sb-3sa)+sa.sc^3]  ...sc[3] + sa[5] - sc[1] + sc(sa-sb)[4])
------------------ above from TwoPos

	from a1   [ 0    sc-sa   sb(2sc-sb)-sa^2             sb^2(3sc-2sb)-sa^3                          0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]
	for a2    [ 0     0      sb(2sc-sb)-sa^2-2sa(sc-sa)  sb^2(3sc-2sb)-sa^3-3sa^2(sc-sa)             0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      sb(2sc-sb)-sa(2sc-sa)       sb^2(3sc-2sb)-sa^2(3sc-2sa)                 0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      sa^2-sb^2+2sc(sb-sa)        2(sa^3-sb^3)+3sc(sb^2-sa^2)                 0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(3sc.sb+3sc.sa-2sa^2-2sa.sb-2sb^2)   0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(sb(3sc-sa-2sb)+sa(3sc-sb-2sa))      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(sb(2(sc-sb)+sc-sa)+sa(2(sc-sa)+sc-sb))  0 0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]

	for a1    [ 0  -(sb-sa)((sc-sa)+(sc-sb))  0    sa(sb-sa)(4sb(sc-sb)+2sb(sc-sa))+sa(sb-sa))      0     0    2sa(sc-sb)^2    2sa(sc+2sb)(sc-sb)^2]  ...2sa[3] + 2sa[5] - 2sa[1] + 2sa(sc-sb)[4] - 2sa(sc-sa)[2] - (sb-sa)((sc-sa)+(sc-sb)) [2]
	          [ 0  -(sb-sa)((sc-sa)+(sc-sb))  0    sa(sb-sa)(6sc.sb - 4sb^2 - sa.sb - sa^2))        0     0    2sa(sc-sb)^2    2sa(sc+2sb)(sc-sb)^2]  ...2sa[3] + 2sa[5] - 2sa[1] + 2sa(sc-sb)[4] - ((sb+sa)(sc-sa) + (sb-sa)(sc-sb)) [2] ==> (-2bc+2ab+aa-2ab+bb)[2] ==> +((b-a)^2-2b(c-a))[2]
	          [ 0  -(sb-sa)((sc-sa)+(sc-sb))  0    sa(sb-sa)(6sc.sb - 4sb^2 - sa.sb - sa^2))        0     0    2sa(sc-sb)^2    2sa(sc+2sb)(sc-sb)^2]  ...2sa[3] + 2sa[5] - 2sa[1] + 2sa(sc-sb)[4] + ((sb-sa)^2-2sb(sc-sa))[2]


	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(sb(3sc-sa-2sb)+sa(3sc-sb-2sa))      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...([3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]) / (sb-sa)((sc-sa)+(sc-sb))

	from a0   [sa-sc  0       sa(sc(sb-sa)+sb(sc-sb))     sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...(sa([3]+[5]) - sc[1] + sa(sc-sb)[4]) / sa(sc(sb-sa)+sb(sc-sb))

(sb-sa)((sc-sa)+(sc-sb)) = sa^2-sb^2+2sb.sc-2sa.sc
sc(sb-sa)+sb(sc-sb) = 2sb.sc - sa.sc - sb^2 = (sa^2-sb^2+2sb.sc-2sa.sc) + sa(sc - sa)

	          [ sc      -sa                -sa(sc-sb)    -sa       sa(sc(sb-sa)+sb(sc-sb))           sa(sc(sb^2-sa^2)+2sb^2(sc-sb))     sa(sc-sb)^2        sa(sc+2sb)(sc-sb)^2   ] / (sc-sa)
	mxAB2ab = [2sa    2sb(sc-sa)-(sb-sa)^2    -2sa      -2sa(sc-sb)      -2sa            sa(sb-sa)(4sb(sc-sb)+2sb(sc-sa))+sa(sb-sa))   2sa(sc-sb)^2       2sa(sc+2sb)(sc-sb)^2   ] / (sb-sa)((sc-sa)+(sc-sb))
	          [-1       -(sc-sa)                1            sc-sb         1                 -(sb-sa)(sb(3sc-sa-2sb)+sa(3sc-sb-2sa))     -(sc-sb)^2          (sc+2sb)(sc-sb)^2   ] / (sb-sa)((sc-sa)+(sc-sb))
	          [ 0        0        0        0        0                      1                                     0                       0          ]
	          [ sc      -sc   sc(sb-sa)   -sa    -sc(sa-sb)^2           sc(sb^2(3sa-2sb)-sa^3)      sb.sc(sb-2sa)+sa.sc^2   sc.sb^2(2sb-3sa)+sa.sc^3] / (sc-sa)
	          [-1        1       sa-sb     1      (sa-sb)^2             sa^3-sb^2(3sa-2sb)             sb(2sa-sb)-sc^2         sb^2(3sa-2sb)-sc^3   ] / (sc-sa)
	          [ 0        0        0        0        0                      0                                     1                       0          ]
	          [ 0        0        0        0        0                      0                                     0                       1          ]

NOT FINISHED YET

  mxdab is transpose of mxAB2ab(:,6:8)

	mxdab =




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
	ASSERT(0);	// NOT FINISHED YET

	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;

	CMatrix mxSab(8,8), vtSaby(8);
	CMatrix mxdab(4,8), mxBCs(4,8);

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double sc = m_FitB.m_SFinal;				// final s of next poly
	double SaP, SbP, ScP;


	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumASp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumBSp[r+c];
		}


	// set matrix 'mxBCs'
	SaP = SbP = ScP = 1;
	for (int c = 0; c < 4; c++)
	{
		mxBCs.elem(1, c+0) = c * SaP;
		mxBCs.elem(1, c+4) = 0;
		mxBCs.elem(3, c+0) = c * SbP;
		mxBCs.elem(3, c+4) = -c * SbP;
		mxBCs.elem(5, c+0) = 0;
		mxBCs.elem(5, c+4) = c * ScP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
			ScP *= sc;
		}
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(2, c+0) = SbP;
		mxBCs.elem(2, c+4) = -SbP;
		mxBCs.elem(4, c+0) = 0;
		mxBCs.elem(4, c+4) = ScP;
	}


	// set matrix 'mxdab'
	double scmsa = sc - sa;
	double scmsb = sc - sb;
	double sbmsa = sb - sa;
	double scpsa = sc + sa;
	double scpsb = sc + sb;
	double sbpsa = sb + sa;

	double Onscmsa = 1 / scmsa;
	double saOnscmsa = sa / scmsa;
	double scOnscmsa = sc / scmsa;

	mxdab = 0;
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(2  ,6  ) = 1;
	mxdab.elem(3  ,7  ) = 1;

	double base = sc*sbmsa + sb*scmsb;
	mxdab.elem(0  ,0  ) = saOnscmsa*base;
	mxdab.elem(0  ,1  ) = -Onscmsa*(base + sa*scmsa);

	base = sc*sbpsa*sbmsa + 2*sb*sb*scmsb;
	mxdab.elem(1  ,0  ) = saOnscmsa*base;
	mxdab.elem(1  ,1  ) = -Onscmsa*(base + sa*sa*scmsa);

	mxdab.elem(2  ,1  ) = -Onscmsa*scmsb*scmsb;
	mxdab.elem(2  ,0  ) = -sa*mxdab.elem(2,1);

	mxdab.elem(3  ,1  ) = (sc+2*sb) * mxdab.elem(2,1);
	mxdab.elem(3  ,0  ) = -sa*mxdab.elem(3,1);
	
	mxdab.elem(0  ,1+4) = Onscmsa*sbmsa*sbmsa;
	mxdab.elem(0  ,0+4) = -sc*mxdab.elem(0,1+4);

	double base1 = 2*sb*sb*sbmsa;
	mxdab.elem(1  ,1+4) = Onscmsa*(base1 - sa*sbpsa*sbmsa);
	mxdab.elem(1  ,0+4) = -sc*mxdab.elem(1,1+4);

	double base2 = sb*sbmsa + sa*scmsb;
	mxdab.elem(2  ,0+4) =  scOnscmsa*base2;
	mxdab.elem(2  ,1+4) = -Onscmsa*(base2 + sc*scmsa);

	base1 += sa*scpsb*scmsb;
	mxdab.elem(3  ,0+4) =  scOnscmsa*base1;
	mxdab.elem(3  ,1+4) = -Onscmsa*(base1 + sc*sc*scmsa);
	








//double c88[8][3], c44[8][3], c44b[8][3];
CMatrix c88(3,8), c44(3,8), c44b(3,8);

/////////////////////////////////
// 8x8 matrix
//
{	
	CMatrix vtBCs(4);
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	mxBCs.CopyToDestLoc(mxAllEqus, 4,0);

	// set vector 'vtBCs'
	vtBCs[1] = 0;
	vtBCs[2] = 0;

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
	vtBCs[3] = vtPosFinal[ax];			// Pos(Sc)

	vtAllVals.ProdPart(mxdab, vtSaby);
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
// rearrange rows to work OK with LUSolve() !!!
//	int arRowArrange[] = {4, 5, 0, 1, 6, 7, 2, 3}; 
//                       0, 1, 5, 4, 2, 6, 3, 7   matlab p
//	int arRowArrange[] = {4, 5, 7, 6, 0, 2, 1, 3};
	int arRowArrange[] = {5, 7, 6, 0, 1, 2, 3, 4};
	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}

//afxDump << "8x8 mxAllEqus2: " << mxAllEqus2;
//afxDump << "vtAllVals2: " << vtAllVals2;

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

	CMatrix vtValSolu;
	vtValSolu = mxAllEqus2 * vtab;
	CMatrix vtDiff;
	vtDiff = vtValSolu - vtAllVals2;
	afxDump << "vtValSolu: " << vtValSolu;
	afxDump << "vtDiff: " << vtDiff;

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c88.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}








///////////////////////////////////////////////
// 4x4 matrix
//

// Calculate by changing vtab to vtAB with only 4 unknowns and rearranging to solve
// a 4x4 matrix equation
/*
	mxdab * mxSab * vtab  =  mxdab * vtSaby		- 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					- 4 equ's

	mxSab * mxAB2ab * vtAB  =  vtSaby				Incorporates BC's
	vtab = mxAB2ab * vtAB

	where mxAB2ab is inv(mxab2AB) and is obtained algebraicly from mxab2AB!
	Note:  mxab2AB * vtab = vtAB
	mxab2AB is derived from: mxBCs * vtab = vtBCs with ones added for extra vtab elements to make a square matrix

	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	[1HC] = 1st half of Columns etc.
	vtAB[2HR] is 2nd half of rows which are the only 4 unknowns

*/

	CMatrix mxAllEqus(4,4);
	CMatrix vtAllVals(4);
	CMatrix mxAB2ab2HC(8,4);
	
	CMatrix vtab(8);
	CMatrix vtValDiff(4);
	CMatrix vtAB2HR(4);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

{
	CMatrix vtdabSy(4);
	CMatrix mxdabS(4,8);

	// mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	mxdabS.Prod(mxdab, mxSab);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	vtdabSy.Prod(mxdab, vtSaby);

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = sc*vtPosInit[ax] - sa*vtPosFinal[ax];			// Pos(sa) & Pos(sc)
	vtab[1] = vtab[5] =   -vtPosInit[ax] + vtPosFinal[ax];					// Pos(sa) & Pos(sc)

	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

afxDump << "4x4 mxAllEqus: " << mxAllEqus;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c44.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}



// without mxdab *
{
	CMatrix mxSab2HR(4,8);
	CMatrix vtSaby2HR(4);

	// without mxdab *
	// mxSab * mxAB2ab * vtAB  =  vtSaby
	// mxSab[2HR] * mxAB2ab[2HC] * vtAB[2HR]  =  vtSaby[2HR] - mxSab[2HR] * mxAB2ab[1HC] * vtAB[1HR]
	mxSab2HR.CopyFromSrcLoc(mxSab, 4,0);
//	mxdabS.Prod(mxdab, mxSab);
//	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxSab2HR, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = sc*vtPosInit[ax] - sa*vtPosFinal[ax];			// Pos(sa) & Pos(sc)
	vtab[1] = vtab[5] =   -vtPosInit[ax] + vtPosFinal[ax];					// Pos(sa) & Pos(sc)

	vtValDiff.Prod(mxSab2HR, vtab);
	vtSaby2HR.CopyFromSrcLoc(vtSaby, 4,0);
	vtAllVals = vtSaby2HR - vtValDiff;

afxDump << "4x4 without mxdab* mxAllEqus: " << mxAllEqus;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;


for (int i = 0; i < 8; i++)
	c44b.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}

afxDump << "c88: " << c88;
afxDump << "c44: " << c44;
afxDump << "c44b without mxdab*: " << c44b;

}	// end CPolySegDblFit::FitDblCubicTwoPosIDeriv(double Kv0)





void CPolySegDblFit::FitDblCubicTwoPosTwoDeriv(double Kv0, double Kv1)
{
/* Fit two parametric cubics with:
	Initial pos & derivative = Kv0 * m_vtDirInitUnit at Sa
	Final pos & derivative = Kv1 * m_vtDirFinalUnit at Sc
	Pos and deriv are continuous at knot at Sb
	Each axis is solved independently
*/

	
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
	mxBCs = [ 0        1       2sa      3sa^2     0        0        0        0    ]
	        [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ]	       
	        [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ]
	        [ 0        0        0        0        1        sc       sc^2     sc^3 ]
	        [ 0        0        0        0        0        1       2sc      3sc^2 ]

	        [ Pos(sa) ]
	vtBCs = [ Vel(sa) ]
	        [     0   ]
	        [     0   ]
	        [ Pos(sc) ]
	        [ Vel(sc) ]


------------------------------

	mxab2AB * vtab = vtAB
	mxAB2ab * vtAB = vtab

-----------------

	Arrange vtAB from vtBCs to group unknows:

	vtAB = [ Pos(sa)  Vel(sa)  0  0  Pos(sc)  Vel(sc)  A3  B3 ]'

	          [ 1        sa       sa^2     sa^3     0        0        0        0    ] ...1  prev 1
	mxab2AB = [ 0        1       2sa      3sa^2     0        0        0        0    ] ...2
	          [ 1        sb       sb^2     sb^3    -1       -sb      -sb^2    -sb^3 ] ...3  prev 2
	          [ 0        1       2sb      3sb^2     0       -1      -2sb     -3sb^2 ] ...4  prev 3
	          [ 0        0        0        0        1        sc       sc^2     sc^3 ] ...5  prev 4
	          [ 0        0        0        0        0        1       2sc      3sc^2 ] ...6
	          [ 0        0        0        1        0        0        0        0    ] ...7
	          [ 0        0        0        0        0        0        0        1    ] ...8

	to invert:
	          [ 0      sa-sb  sa^2-sb^2  sa^3-sb^3  1        sb        sb^2      sb^3 ]  ...[1]-[3]  ...([x] => from vtAB[x])
	          [ 1      0  sa(sa-2sb) sa(sa^2-3sb^2) 0        sa       2sa.sb  3sa.sb^2]  ...[1]-sa[4]



	          [ 1      sb       sb^2     sb^3       0     sc-sb   sc^2-sb^2  sc^3-sb^3]  ...[3]+[5]
	          [sa-sb  0  sa.sb(sb-sa) sa.sb(sb^2-sa^2) 0  sa(sc-sb) sa(sc^2-sb^2) sa(sc^3-sb^3)]  ...sa([3]+[5])-sb[1]
	          [sa-sc  0  sa(sb(sb-sa)-(sc-sb)(sa-2sb))  sa(sb(sb^2-sa^2)-(sc-sb)(sa^2-3sb^2))    0     0   sa((sc^2-sb^2)-2sb(sc-sb))  sa((sc^3-sb^3)-3sb^2(sc-sb))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	          [sa-sc  0  sa(-sa.sc+2sb.sc-sb^2)  sa(-sc.sa^2+3sc.sb^2-2sb^3)    0     0   sa(sc^2-2sb.sc+sb^2)  sa(sc^3-3sb^2.sc+2sb^3))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	          [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc(sc^2-sb^2)+2sb^2(sb-sc))]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4])
	for a0    [sa-sc  0  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([3]+[5])-sb[1] - (sc-sb)([1]-sa[4]) ==> sa([3]+[5]) - sc[1] + sa(sc-sb)[4]

	          [ 0  sa(sc-sa)  sa(sc(sb-sa)+sb(sc-sb)+sa(sc-sa))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa))    0     0   sa(sc-sb)^2  sa(sc+2sb)(sc-sb)^2]  ...sa([3]+[5]) - sc[1] + sa(sc-sb)[4] - (sa-sc)[1] ==> sa([3]+[5]) - sa[1] + sa(sc-sb)[4]
	          [ 0  sc-sa  sc(sb-sa)+sb(sc-sb)+sa(sc-sa)  sc(sb^2-sa^2)+2sb^2(sc-sb)+sa^2(sc-sa)    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]
	for a1    [ 0  sc-sa  sb(2sc-sb)-sa^2  sb^2(3sc-2sb)-sa^3    0     0   (sc-sb)^2  (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]

	          [ 0     0  sb(2sc-sb)-sa^2-2sb(sc-sa)  sb^2(3sc-2sb)-3sb^2(sc-sa)-sa^3    0     sc-sa   (sc-sb)^2+2sb(sc-sa)  (sc+2sb)(sc-sb)^2+3sb^2(sc-sa)]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[4] ==> [3] + [5] - [1] + (sa-sb)[4]
	for b1    [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3                 0     sc-sa    sb(sb-2sa)+sc^2       sb^2(2sb-3sa)+sc^3]  ...[3] + [5] - [1] + (sa-sb)[4]

	          [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sc^3-sc^2(sc-sa)       sc.sb^2(2sb-3sa)+sc^4-sc^3(sc-sa)]  ...sc([3] + [5] - [1] + (sa-sb)[4]) - (sc-sa)[5] ==> sc[3] + sa[5] - sc[1] + sc(sa-sb)[4])
	for b0    [ 0     0    -sc(sa-sb)^2             sc(sb^2(3sa-2sb)-sa^3)        -(sc-sa)    0    sb.sc(sb-2sa)+sa.sc^2       sc.sb^2(2sb-3sa)+sa.sc^3]  ...sc[3] + sa[5] - sc[1] + sc(sa-sb)[4])
------------------ above from TwoPos
	for b1    [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3       0     sc-sa    sb(sb-2sa)+sc^2               sb^2(2sb-3sa)+sc^3]  ...[3] + [5] - [1] + (sa-sb)[4]
	          [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3       0     0        sb(sb-2sa)+sc^2-2sc(sc-sa)    sb^2(2sb-3sa)+sc^3-3sc^2(sc-sa)]  ...[3] + [5] - [1] + (sa-sb)[4] - (sc-sa)[6]
	          [ 0     0    -(sa-sb)^2                sb^2(3sa-2sb)-sa^3       0     0        sb(sb-2sa)-sc(sc-2sa)         sb^2(2sb-3sa)-sc^2(2sc-3sa)]  ...[3] + [5] - [1] + (sa-sb)[4] - (sc-sa)[6]

	for a1    [ 0    sc-sa   sb(2sc-sb)-sa^2             sb^2(3sc-2sb)-sa^3       0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4]
	          [ 0     0      sb(2sc-sb)-sa^2-2sa(sc-sa)  sb^2(3sc-2sb)-sa^3-3sa^2(sc-sa)      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      sb(2sc-sb)-sa(2sc-sa)       sb^2(3sc-2sb)-sa^2(3sc-2sa)      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      sa^2-sb^2+2sc(sb-sa)        2(sa^3-sb^3)+3sc(sb^2-sa^2)      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(3sc.sb+3sc.sa-2sa^2-2sa.sb-2sb^2)      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(sb(3sc-sa-2sb)+sa(3sc-sb-2sa))      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	          [ 0     0      (sb-sa)((sc-sa)+(sc-sb))    (sb-sa)(sb(2(sc-sb)+sc-sa)+sa(2(sc-sa)+sc-sb))      0     0      (sc-sb)^2      (sc+2sb)(sc-sb)^2]  ...[3] + [5] - [1] + (sc-sb)[4] - (sc-sa)[2]
	              sa^2-sb^2+2sc(sb-sa)
	              (sa-sb)(sa+sb)+2sc(sb-sa)
	              (sb-sa)((sc-sa)+(sc-sb))


NOT FINISHED YET


	          [ sc      -sa  -sa(sc-sb)   -sa  sa(sc(sb-sa)+sb(sc-sb))  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))     sa(sc-sb)^2        sa(sc+2sb)(sc-sb)^2   ] / (sc-sa)
	mxAB2ab = [-1        1       sc-sb     1     sa^2-sb(2sc-sb)        sa^3-sb^2(3sc-2sb)                  -(sc-sb)^2         -(sc+2sb)(sc-sb)^2   ] / (sc-sa)
	          [ 0        0        0        0        1                      0                                     0                       0          ]
	          [ 0        0        0        0        0                      1                                     0                       0          ]
	          [ sc      -sc   sc(sb-sa)   -sa    -sc(sa-sb)^2           sc(sb^2(3sa-2sb)-sa^3)      sb.sc(sb-2sa)+sa.sc^2   sc.sb^2(2sb-3sa)+sa.sc^3] / (sc-sa)
	          [-1        1       sa-sb     1      (sa-sb)^2             sa^3-sb^2(3sa-2sb)             sb(2sa-sb)-sc^2         sb^2(3sa-2sb)-sc^3   ] / (sc-sa)
	          [ 0        0        0        0        0                      0                                     1                       0          ]
	          [ 0        0        0        0        0                      0                                     0                       1          ]



  mxdab is transpose of mxAB2ab(:,3:8)

	        [  sa(sc(sb-sa)+sb(sc-sb))          sa^2-sb(2sc-sb)        1     0  -sc(sb-sa)^2                           (sb-sa)^2   0    0 ]
	mxdab = [  sa(sc(sb^2-sa^2)+2sb^2(sc-sb))   sa^3-sb^2(3sc-2sb)     0     1  -sc(sb^2(2sb-3sa)+sa^3)       sb^2(2sb-3sa)+sa^3   0    0 ]
	        [  sa(sc-sb)^2                              -(sc-sb)^2     0     0   sc(sb  ( sb-2sa)+sa.sc)        -sb(sb-2sa)-sc^2   1    0 ]
	        [  sa(sc+2sb)(sc-sb)^2              -(sc+2sb)(sc-sb)^2     0     0   sc(sb^2(2sb-3sa)+sa.sc^2)   -sb^2(2sb-3sa)-sc^3   0    1 ]
	               /(sc-sa)                             /(sc-sa)                       /(sc-sa)                     /(sc-sa)

Regular expersion method:
	        [  sa(2sb.sc-sa.sc-sb^2)          sa^2-sb.sc-sb(sc-sb)        1     0  -sc(sb^1(1(sb-sa)-sa)+sa^2)        sb^1(1(sb-sa)-sa)+sa^2    0    0 ]
	mxdab = [  sa(3sc.sb^2-sc.sa^2-2sb^3)     sa^3-sc.sb^2-2sb^2(sc-sb)   0     1  -sc(sb^2(2(sb-sa)-sa)+sa^3)        sb^2(2(sb-sa)-sa)+sa^3    0    0 ]
	        [  sa(sc^2-2sb.sc+sb^2)             -(sc-sb)^2                0     0   sc(sb^1(1(sb-sa)-sa)+sa.sc^1)    -sb^1(1(sb-sa)-sa)-sc^2    1    0 ]
	        [  sa(sc^3-3sc.sb^2+2sb^3)           -(sc+2sb)(sc-sb)^2       0     0   sc(sb^2(2(sb-sa)-sa)+sa.sc^2)    -sb^2(2(sb-sa)-sa)-sc^3    0    1 ]
	               /(sc-sa)                             /(sc-sa)                       /(sc-sa)                           /(sc-sa)

Best Calculating method:
	        [  sa(sc(sb-sa)+sb(sc-sb))            -sc(sb-sa)-sb(sc-sb)-sa(sc-sa)              1     0  -sc(sb-sa)^2                       (sb-sa)^2                                0    0 ]
	mxdab = [  sa(sc(sb+sa)(sb-sa)+2sb^2(sc-sb))  -sc(sb+sa)(sb-sa)-2sb^2(sc-sb)-sa^2(sc-sa)  0     1  -sc(2sb^2(sb-sa)-sa(sb^2-sa^2))    2sb^2(sb-sa)-sa(sb^2-sa^2)               0    0 ]
	        [  sa(sc-sb)^2                                    -(sc-sb)^2                      0     0   sc(sb(sb-sa)+sa(sc-sb))          -sb(sb-sa)-sa(sc-sb)-sc(sc-sa)            1    0 ]
	        [  sa(sc+2sb)(sc-sb)^2                    -(sc+2sb)(sc-sb)^2                      0     0   sc(2sb^2(sb-sa)+sa(sc^2-sb^2))   -2sb^2(sb-sa)-sa(sc^2-sb^2)-sc^2(sc-sa)   0    1 ]
	               /(sc-sa)                             /(sc-sa)                                                 /(sc-sa)                           /(sc-sa)


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
	ASSERT(0);		// NOT FINISHED YET

	double Wa = 1;
	double Wb = 1;

	const double (&SumASp)[7] = m_FitA.m_PointSums.Sp;			// better array display for debuging!
	double (&SumASpV)[4][3] = m_FitA.m_PointSums.SpV;
	const double (&SumBSp)[7] = m_FitB.m_PointSums.Sp;
	double (&SumBSpV)[4][3] = m_FitB.m_PointSums.SpV;

	CVector vtPosInit = m_vtPosInit;
	CVector vtPosFinal = m_vtPosFinal;
	CVector vtDerivInit = m_vtDirInitUnit * Kv0;
	CVector vtDerivFinal = m_vtDirFinalUnit * Kv1;

	CMatrix mxSab(8,8), vtSaby(8);
	CMatrix mxdab(4,8), mxBCs(4,8);

	double sa = m_FitA.m_SInit;				// initial s of first poly
	double sb = m_FitB.m_SInit;				// initial s of next poly
	double sc = m_FitB.m_SFinal;				// final s of next poly
	double SaP, SbP, ScP;


	// set matrix 'mxSab'
	mxSab = 0;
	for (int r = 0; r < 4; r++)
		for (int c = 0; c < 4; c++)
		{
			mxSab.elem(r, c) = Wa * SumASp[r+c];
			mxSab.elem(r+4, c+4) = Wb * SumBSp[r+c];
		}


	// set matrix 'mxBCs'
	SaP = SbP = ScP = 1;
	for (int c = 0; c < 4; c++)
	{
		mxBCs.elem(1, c+0) = c * SaP;
		mxBCs.elem(1, c+4) = 0;
		mxBCs.elem(3, c+0) = c * SbP;
		mxBCs.elem(3, c+4) = -c * SbP;
		mxBCs.elem(5, c+0) = 0;
		mxBCs.elem(5, c+4) = c * ScP;
		if (c > 0)
		{
			SaP *= sa;
			SbP *= sb;
			ScP *= sc;
		}
		mxBCs.elem(0, c+0) = SaP;
		mxBCs.elem(0, c+4) = 0;
		mxBCs.elem(2, c+0) = SbP;
		mxBCs.elem(2, c+4) = -SbP;
		mxBCs.elem(4, c+0) = 0;
		mxBCs.elem(4, c+4) = ScP;
	}


	// set matrix 'mxdab'
	double scmsa = sc - sa;
	double scmsb = sc - sb;
	double sbmsa = sb - sa;
	double scpsa = sc + sa;
	double scpsb = sc + sb;
	double sbpsa = sb + sa;

	double Onscmsa = 1 / scmsa;
	double saOnscmsa = sa / scmsa;
	double scOnscmsa = sc / scmsa;

	mxdab = 0;
	mxdab.elem(0  ,2  ) = 1;
	mxdab.elem(1  ,3  ) = 1;
	mxdab.elem(2  ,6  ) = 1;
	mxdab.elem(3  ,7  ) = 1;

	double base = sc*sbmsa + sb*scmsb;
	mxdab.elem(0  ,0  ) = saOnscmsa*base;
	mxdab.elem(0  ,1  ) = -Onscmsa*(base + sa*scmsa);

	base = sc*sbpsa*sbmsa + 2*sb*sb*scmsb;
	mxdab.elem(1  ,0  ) = saOnscmsa*base;
	mxdab.elem(1  ,1  ) = -Onscmsa*(base + sa*sa*scmsa);

	mxdab.elem(2  ,1  ) = -Onscmsa*scmsb*scmsb;
	mxdab.elem(2  ,0  ) = -sa*mxdab.elem(2,1);

	mxdab.elem(3  ,1  ) = (sc+2*sb) * mxdab.elem(2,1);
	mxdab.elem(3  ,0  ) = -sa*mxdab.elem(3,1);
	
	mxdab.elem(0  ,1+4) = Onscmsa*sbmsa*sbmsa;
	mxdab.elem(0  ,0+4) = -sc*mxdab.elem(0,1+4);

	double base1 = 2*sb*sb*sbmsa;
	mxdab.elem(1  ,1+4) = Onscmsa*(base1 - sa*sbpsa*sbmsa);
	mxdab.elem(1  ,0+4) = -sc*mxdab.elem(1,1+4);

	double base2 = sb*sbmsa + sa*scmsb;
	mxdab.elem(2  ,0+4) =  scOnscmsa*base2;
	mxdab.elem(2  ,1+4) = -Onscmsa*(base2 + sc*scmsa);

	base1 += sa*scpsb*scmsb;
	mxdab.elem(3  ,0+4) =  scOnscmsa*base1;
	mxdab.elem(3  ,1+4) = -Onscmsa*(base1 + sc*sc*scmsa);
	








//double c88[8][3], c44[8][3], c44b[8][3];
CMatrix c88(3,8), c44(3,8), c44b(3,8);

/////////////////////////////////
// 8x8 matrix
//
{	
	CMatrix vtBCs(4);
	CMatrix mxAllEqus(8,8);
	CMatrix vtAllVals(8);
	CMatrix vtab(8);
	CMatrix mxAllEqus2(8,8);
	CMatrix vtAllVals2(8);

	mxAllEqus.ProdPart(mxdab, mxSab);
	mxBCs.CopyToDestLoc(mxAllEqus, 4,0);

	// set vector 'vtBCs'
	vtBCs[1] = 0;
	vtBCs[2] = 0;

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
	vtBCs[3] = vtPosFinal[ax];			// Pos(Sc)

	vtAllVals.ProdPart(mxdab, vtSaby);
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
// rearrange rows to work OK with LUSolve() !!!
//	int arRowArrange[] = {4, 5, 0, 1, 6, 7, 2, 3}; 
//                       0, 1, 5, 4, 2, 6, 3, 7   matlab p
//	int arRowArrange[] = {4, 5, 7, 6, 0, 2, 1, 3};
	int arRowArrange[] = {5, 7, 6, 0, 1, 2, 3, 4};
	for (int rowD = 0; rowD < 8; rowD++)
	{
		int rowS = arRowArrange[rowD];
		for (int col = 0; col < 8; col++)
			mxAllEqus2.elem(rowD, col) = mxAllEqus.elem(rowS, col);
		vtAllVals2[rowD] = vtAllVals[rowS];
	}

//afxDump << "8x8 mxAllEqus2: " << mxAllEqus2;
//afxDump << "vtAllVals2: " << vtAllVals2;

	mxAllEqus2.LUSolve(vtab.GetArray(), vtAllVals2.GetArray());

	CMatrix vtValSolu;
	vtValSolu = mxAllEqus2 * vtab;
	CMatrix vtDiff;
	vtDiff = vtValSolu - vtAllVals2;
	afxDump << "vtValSolu: " << vtValSolu;
	afxDump << "vtDiff: " << vtDiff;

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c88.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}








///////////////////////////////////////////////
// 4x4 matrix
//

// Calculate by changing vtab to vtAB with only 4 unknowns and rearranging to solve
// a 4x4 matrix equation
/*
	mxdab * mxSab * vtab  =  mxdab * vtSaby		- 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					- 4 equ's

	mxSab * mxAB2ab * vtAB  =  vtSaby				Incorporates BC's
	vtab = mxAB2ab * vtAB

	where mxAB2ab is inv(mxab2AB) and is obtained algebraicly from mxab2AB!
	Note:  mxab2AB * vtab = vtAB
	mxab2AB is derived from: mxBCs * vtab = vtBCs with ones added for extra vtab elements to make a square matrix

	mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	[1HC] = 1st half of Columns etc.
	vtAB[2HR] is 2nd half of rows which are the only 4 unknowns

*/

	CMatrix mxAllEqus(4,4);
	CMatrix vtAllVals(4);
	CMatrix mxAB2ab2HC(8,4);
	
	CMatrix vtab(8);
	CMatrix vtValDiff(4);
	CMatrix vtAB2HR(4);

	// set mxAB2ab2HC - '2HC' => 2nd half of columns
	mxdab.Transpose(mxAB2ab2HC);		// mxAB2ab2HC is just transpose of mxdab!!

{
	CMatrix vtdabSy(4);
	CMatrix mxdabS(4,8);

	// mxdab * mxSab * mxAB2ab[2HC] * vtAB[2HR]  =  mxdab * vtSaby - mxdab * mxSab * mxAB2ab[1HC] * vtAB[1HR]
	mxdabS.Prod(mxdab, mxSab);
	mxAllEqus.Prod(mxdabS, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}
	vtdabSy.Prod(mxdab, vtSaby);

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = sc*vtPosInit[ax] - sa*vtPosFinal[ax];			// Pos(sa) & Pos(sc)
	vtab[1] = vtab[5] =   -vtPosInit[ax] + vtPosFinal[ax];					// Pos(sa) & Pos(sc)

	vtValDiff.Prod(mxdabS, vtab);
	vtAllVals = vtdabSy - vtValDiff;

afxDump << "4x4 mxAllEqus: " << mxAllEqus;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;

	for (int i = 0; i < 4; i++)
	{
		m_FitA.m_Poly[i][ax] = vtab[i];
		m_FitB.m_Poly[i][ax] = vtab[i+4];
	}

for (i = 0; i < 8; i++)
	c44.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}



// without mxdab *
{
	CMatrix mxSab2HR(4,8);
	CMatrix vtSaby2HR(4);

	// without mxdab *
	// mxSab * mxAB2ab * vtAB  =  vtSaby
	// mxSab[2HR] * mxAB2ab[2HC] * vtAB[2HR]  =  vtSaby[2HR] - mxSab[2HR] * mxAB2ab[1HC] * vtAB[1HR]
	mxSab2HR.CopyFromSrcLoc(mxSab, 4,0);
//	mxdabS.Prod(mxdab, mxSab);
//	vtdabSy.Prod(mxdab, vtSaby);
	mxAllEqus.Prod(mxSab2HR, mxAB2ab2HC);

for (int ax = 0; ax < 3; ax++)
{
	// set vector 'vtSaby'
	for (r = 0; r < 4; r++)
	{
		vtSaby[r] = Wa * SumASpV[r][ax];
		vtSaby[r+4] = Wb * SumBSpV[r][ax];
	}

	//	vtab(part) = mxAB2ab[1HC] * vtAB[1HR]				// partly setup vtab
	vtab = 0;
	vtab[0] = vtab[4] = sc*vtPosInit[ax] - sa*vtPosFinal[ax];			// Pos(sa) & Pos(sc)
	vtab[1] = vtab[5] =   -vtPosInit[ax] + vtPosFinal[ax];					// Pos(sa) & Pos(sc)

	vtValDiff.Prod(mxSab2HR, vtab);
	vtSaby2HR.CopyFromSrcLoc(vtSaby, 4,0);
	vtAllVals = vtSaby2HR - vtValDiff;

afxDump << "4x4 without mxdab* mxAllEqus: " << mxAllEqus;

	mxAllEqus.LUSolve(vtAB2HR.GetArray(), vtAllVals.GetArray());
	//	vtab = mxAB2ab[1HC] * vtAB[1HR] + mxAB2ab[2HC] * vtAB[2HR]
	vtab += mxAB2ab2HC * vtAB2HR;		// vtab is now complete!

//afxDump << "vtab: " << vtab;


for (int i = 0; i < 8; i++)
	c44b.elem(ax,i) = vtab[i];

}	// for (int ax = 0; ax < 3; ax++)
}



afxDump << "c88: " << c88;
afxDump << "c44: " << c44;
afxDump << "c44b without mxdab*: " << c44b;

}	// end CPolySegDblFit::FitDblCubicTwoPosTwoDeriv(double Kv0, double Kv1)






