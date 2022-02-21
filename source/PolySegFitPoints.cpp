// PolySegFitPoints.cpp: implementation of the CPolySegFitPoints class.
//
//////////////////////////////////////////////////////////////////////


//#include "stdafx.h"

#include <math.h>


#include "PolySegFitPoints.h"



#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


/*
	double x, y, xPow;

	Num++;
	SumYXp[0] += y;
	SumXp[1]  += (xPow = x);
	SumYXp[1] += xPow * y;
	SumXp[2]  += (xPow *= x);
	SumYXp[2] += xPow * y;
	SumXp[3]  += (xPow *= x);
	SumYXp[3] += xPow * y;
	SumXp[4]  += (xPow *= x);
	SumXp[5]  += (xPow *= x);
	SumXp[6]  += (xPow * x);
	SumY2 += y * y;		// needed to calc residual error (sum of squares)

------------

	Best cubic fit
	Y = a0 + a1.x + a2.x^2 + a3.x^3

	Sum of residuals squared = Sr
		Sr uses all the sums S(x^0..6) plus S(y^2)
	Sr = Sum((y - Y)^2)		where y is actual value and Y is cubic value

	Sr = S(y^2) + a0^2.S(x^0) + a1^2.S(x^2) + a2^2.S(x^4) + a3^2.S(x^6)
		+ 2{-a0.S(x^0.y) - a1.S(x^1.y) - a2.S(x^2.y) - a3.S(x^3.y) }
		+ 2{ a0.a1.S(x^1) + a0.a2.S(x^2) + (a0.a3+a1.a2).S(x^3) + a1.a3.S(x^4) + a2.a3.S(x^5) }

		S(...) = sum(...)
		note S(x^0) = N	number of points

	dSr/da0     [ S(x^0)  S(x^1)  S(x^2)  S(x^3) ]   a0     S(x^0.y)
	dSr/da1 = 2*[ S(x^1)  S(x^2)  S(x^3)  S(x^4) ] * a1 - 2*S(x^1.y)
	dSr/da2     [ S(x^2)  S(x^3)  S(x^4)  S(x^5) ]   a2     S(x^2.y)
	dSr/da3     [ S(x^3)  S(x^4)  S(x^5)  S(x^6) ]   a3     S(x^3.y)

	For minimum Sr all dSr/da are 0.  Therefor:
	[ S(x^0)  S(x^1)  S(x^2)  S(x^3) ]   a0   S(x^0.y)
	[ S(x^1)  S(x^2)  S(x^3)  S(x^4) ] * a1 = S(x^1.y)
	[ S(x^2)  S(x^3)  S(x^4)  S(x^5) ]   a2   S(x^2.y)
	[ S(x^3)  S(x^4)  S(x^5)  S(x^6) ]   a3   S(x^3.y)

	--> vtdSra = 2 * mxS * vta  -  2 * vtSy
	->  mxS * vta  =  vtSy

-------------------------
	With end boundary conditions


	Initial Pos and Vel fixed and initial x = 0:
	a0 = P0 & a1 = V0

	[ S(x^4)  S(x^5) ] * a2  =  S(x^2.y) - [ S(x^2)  S(x^3) ] * a0
	[ S(x^5)  S(x^6) ]   a3     S(x^3.y)   [ S(x^3)  S(x^4) ]   a1

	If initial x != 0, but x = xa
	[ 1        xa       xa^2     xa^3 ]   a0    Pos(xa)
	[ 0        1       2xa      3xa^2 ] * a1 =  Vel(xa)
	                                      a2
	                                      a3
	mxBCs * vta  =  vtBCs

	eliminate a0, a1 (solve for these!)
	step 1:
	[ 1        0       -xa^2   -2xa^3 ]   a0    Pos(xa) - xa*Vel(xa)		//  -xa * row2
	[ 0        1       2xa      3xa^2 ] * a1 =  Vel(xa)
	                                      a2
	                                      a3
	dSr/da2 = dSra/da2 + dSra/da0 * da0/da2 + dSra/da1 * da1/da2 = 0		// from all dependent params
	-->
	dSr/da2   [ da0/da2  da1/da2    1    0 ]   dSra/da0     0
	dSr/da3 = [ da0/da3  da1/da3    0    1 ] * dSra/da1  =  0
	                                           dSra/da2
	                                           dSra/da3

	-> vtdSrRed = mxda * vtdSra = 0

	          [  xa^2   -2xa     1    0 ]
	-> mxda = [ 2xa^3   -3xa^2   0    1 ]
	     
	  mxda * vtdSra = 0
	& vtdSra = 2 * mxS * vta  -  2 * vtSy

	mxda * mxS * vta  =  mxda * vtSy		// 2 equ's		(note: mxda doesn't have an inverse - not square!)
	     mxBCs * vta  =  vtBCs				// 2 equ's







--------------------------
	Optimise knot of two polys with continuous pos and vel
	  P(a0 a1 a2 a3)    P(b0 b1 b2 b3)
	|-----------------|-----------------|
	xa						xb

	assuming x starts at 0 for poly B !!:
  b0 = a0 + a1.xb + a2.xb^2 + a3.xb^3		continuous pos
  b1 = a1 + 2a2.xb + 3a3.xb^2					continuous vel

  for any other starting x for poly B:
  0 = a0-b0 + (a1-b1).xb + (a2-b2).xb^2 + (a3-b3).xb^3		continuous pos
  0 = (a1-b1) + 2(a2-b2).xb + 3(a3-b3).xb^2						continuous vel


		First set of sums is over poly A
	Sr = S(y^2) + a0^2.S(x^0) + a1^2.S(x^2) + a2^2.S(x^4) + a3^2.S(x^6)
		+ 2{-a0.S(x^0.y) - a1.S(x^1.y) - a2.S(x^2.y) - a3.S(x^3.y) }
		+ 2{ a0.a1.S(x^1) + a0.a2.S(x^2) + (a0.a3+a1.a2).S(x^3) + a1.a3.S(x^4) + a2.a3.S(x^5) }

		Following sums are over poly B
		+ S(y^2) + b0^2.S(x^0) + b1^2.S(x^2) + b2^2.S(x^4) + b3^2.S(x^6)
		+ 2{-b0.S(x^0.y) - b1.S(x^1.y) - b2.S(x^2.y) - b3.S(x^3.y) }
		+ 2{ b0.b1.S(x^1) + b0.b2.S(x^2) + (b0.b3+b1.b2).S(x^3) + b1.b3.S(x^4) + b2.b3.S(x^5) }





	Initial pos & vel of poly a is set
	pos & vel between polys a & b is continous

	Find minimum of Sr = Wa*Sra + Wb*Srb		with weightings for Sra and Srb
			dSr = Wa*dSra + Wb*dSrb = 0
	dSra/da0    [ Sa(x^0)  Sa(x^1)  Sa(x^2)  Sa(x^3) ]   a0    [ Sa(x^0.y) ]
	dSra/da1 = 2[ Sa(x^1)  Sa(x^2)  Sa(x^3)  Sa(x^4) ] * a1 - 2[ Sa(x^1.y) ]
	dSra/da2    [ Sa(x^2)  Sa(x^3)  Sa(x^4)  Sa(x^5) ]   a2    [ Sa(x^2.y) ]
	dSra/da3    [ Sa(x^3)  Sa(x^4)  Sa(x^5)  Sa(x^6) ]   a3    [ Sa(x^3.y) ]

	dSrb/db0    [ Sb(x^0)  Sb(x^1)  Sb(x^2)  Sb(x^3) ]   b0    [ Sb(x^0.y) ]
	dSrb/db1 = 2[ Sb(x^1)  Sb(x^2)  Sb(x^3)  Sb(x^4) ] * b1 - 2[ Sb(x^1.y) ]
	dSrb/db2    [ Sb(x^2)  Sb(x^3)  Sb(x^4)  Sb(x^5) ]   b2    [ Sb(x^2.y) ]
	dSrb/db3    [ Sb(x^3)  Sb(x^4)  Sb(x^5)  Sb(x^6) ]   b3    [ Sb(x^3.y) ]

	-->
	dSra/da0    [ Sa(x^0)  Sa(x^1)  Sa(x^2)  Sa(x^3)  0        0        0        0       ]   a0    [ Sa(x^0.y) ]
	dSra/da1 = 2[ Sa(x^1)  Sa(x^2)  Sa(x^3)  Sa(x^4)  0        0        0        0       ] * a1 - 2[ Sa(x^1.y) ]
	dSra/da2    [ Sa(x^2)  Sa(x^3)  Sa(x^4)  Sa(x^5)  0        0        0        0       ]   a2    [ Sa(x^2.y) ]
	dSra/da3    [ Sa(x^3)  Sa(x^4)  Sa(x^5)  Sa(x^6)  0        0        0        0       ]   a3    [ Sa(x^3.y) ]
	dSrb/db0    [ 0        0        0        0        Sb(x^0)  Sb(x^1)  Sb(x^2)  Sb(x^3) ]   b0    [ Sb(x^0.y) ]
	dSrb/db1    [ 0        0        0        0        Sb(x^1)  Sb(x^2)  Sb(x^3)  Sb(x^4) ]   b1    [ Sb(x^1.y) ]
	dSrb/db2    [ 0        0        0        0        Sb(x^2)  Sb(x^3)  Sb(x^4)  Sb(x^5) ]   b2    [ Sb(x^2.y) ]
	dSrb/db3    [ 0        0        0        0        Sb(x^3)  Sb(x^4)  Sb(x^5)  Sb(x^6) ]   b3    [ Sb(x^3.y) ]
	-> vtdSrab = 2 * mxSab * vtab  -  2 * vtSaby

	At xb:
	Apos(xb)   [ 1  xb  xb^2  xb^3  ]   a0
	Avel(xb) = [ 0   1  2xb   3xb^2 ] * a1
	                                    a2
	                                    a3

	Bpos(xb)   [ 1  xb  xb^2  xb^3  ]   b0
	Bvel(xb) = [ 0   1  2xb   3xb^2 ] * b1
	                                    b2
	                                    b3


boundary condition equations:
	[ 1        xb       xb^2     xb^3    -1       -xb      -xb^2    -xb^3  ]   a..    0				// Pos(xb) cont.
	[ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2  ] * a.. =  0				// Vel(xb) cont.
	[ 1        xa       xa^2     xa^3     0        0        0        0     ]   b..    Pos(xa)		// Pos(xa)
	[ 0        1       2xa      3xa^2     0        0        0        0     ]   b..    Vel(xa)		// Vel(xa)
	-> mxBCs * vtab = vtBCs		// 4 equ's

	eliminate a0, a1, b0, b1 (solve for these!)
	step 1:
	[ 1        0 xb^2-2xa*xb xb^3-3xb*xa^2 -1       -xb      -xb^2    -xb^3    ]   a..    0 - xb*Vel(xa)				// -xb * row4
	[ 0        0   2(xb-xa)  3(xb^2-xa^2)   0       -1      -2xb     -3xb^2    ] * a.. =  0 - Vel(xa)					//  -1 * row4
	[ 1        0       -xa^2   -2xa^3       0        0        0        0       ]   b..    Pos(xa) - xa*Vel(xa)		// -xa * row4
	[ 0        1       2xa      3xa^2       0        0        0        0       ]   b..    Vel(xa)
	step 2:
	[ 0        0   (xb-xa)^2 xb^3-3xb*xa^2+2xa^3 -1       -xb      -xb^2    -xb^3    ]   a..    (xa-xb)*Vel(xa) - Pos(xa)	//  -1 * row3
	[ 0        0   2(xb-xa)   3(xb^2-xa^2)        0       -1      -2xb     -3xb^2    ] * a.. =  -Vel(xa)
	[ 1        0       -xa^2   -2xa^3             0        0        0        0       ]   b..    Pos(xa) - xa*Vel(xa)
	[ 0        1       2xa      3xa^2             0        0        0        0       ]   b..    Vel(xa)
	step 3:
	[ 0        0   xa^2-xb^2  2xa^3-2xb^3    -1        0       xb^2     2xb^3    ]   a..    xa*Vel(xa) - Pos(xa)	// -xb * row2
	[ 0        0   2xb-2xa    3xb^2-3xa^2     0       -1      -2xb     -3xb^2    ] * a.. =  -Vel(xa)
	[ 1        0       -xa^2   -2xa^3         0        0        0        0       ]   b..    Pos(xa) - xa*Vel(xa)
	[ 0        1       2xa      3xa^2         0        0        0        0       ]   b..    Vel(xa)

	solution:
	a0   [      xa^2        2xa^3     0        0         1     -xa ] * [ a2  a3  b2  b3  Pos(xa)  Vel(xa) ]'
	a1 = [    -2xa         -3xa^2     0        0         0      1  ]
	b0   [ xa^2-xb^2  2xa^3-2xb^3     xb^2    2xb^3      1     -xa ]
	b1   [ 2xb-2xa    3xb^2-3xa^2   -2xb     -3xb^2      0      1  ]



	dSr/da2 = dSra/da2 + dSra/da0 * da0/da2 + dSra/da1 * da1/da2 + dSrb/db0 * db0/da2 + dSrb/db1 * db1/da2 = 0		// from all dependent params
	-->
	dSr/da2   [ da0/da2  da1/da2    1    0   db0/da2  db1/da2   0    0 ]   dSra/da0     0
	dSr/da3 = [ da0/da3  da1/da3    0    1   db0/da3  db1/da3   0    0 ] * dSra/da1  =  0
	dSr/db2   [ da0/db2  da1/db2    0    0   db0/db2  db1/db2   1    0 ]   dSra/da2     0
	dSr/db3   [ da0/db3  da1/db3    0    0   db0/db3  db1/db3   0    1 ]   dSra/da3     0
	                                                                       dSrb/db0
	                                                                       dSrb/db1
	                                                                       dSrb/db2
	                                                                       dSrb/db3
	-> vtdSrRed = mxdab * vtdSrab = 0

	           [  xa^2   -2xa     1    0     xa^2-xb^2     2xb-2xa     0    0 ]
	-> mxdab = [ 2xa^3   -3xa^2   0    1   2xa^3-2xb^3   3xb^2-3xa^2   0    0 ]
	           [  0        0      0    0          xb^2        -2xb     1    0 ]
	           [  0        0      0    0         2xb^3        -3xb^2   0    1 ]


	  mxdab * vtdSrab = 0
	& vtdSrab = 2*mxSab * vtab  -  2*vtSaby

	mxdab * mxSab * vtab  =  mxdab * vtSaby		// 4 equ's		(note: mxdab doesn't have an inverse - not square!)
	        mxBCs * vtab  =  vtBCs					// 4 equ's

	------------------------------

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

	        [ 1        xb       xb^2     xb^3    -1       -xb      -xb^2    -xb^3 ]
	mxBCs = [ 0        1       2xb      3xb^2     0       -1      -2xb     -3xb^2 ]
	        [ 1        xa       xa^2     xa^3     0        0        0        0    ]
	        [ 0        1       2xa      3xa^2     0        0        0        0    ]

	        [     0   ]
	vtBCs = [     0   ]
	        [ Pos(xa) ]
	        [ Vel(xa) ]

	----------------------------------



  




**************************************
	If only pos is known at xb:

	a0 + a1.xb + a2.xb^2 + a3.xb^3 = val(xb)
	a0 = val(xb) - a1.xb - a2.xb^2 - a3.xb^3

	da1/da0 = -1/xb
	da2/da0 = -1/xb^2
	da3/da0 = -1/xb^3

	da0/da1 = -xb
	da0/da2 = -xb^2
	da0/da3 = -xb^3

	dSr'/da1 = dSr/da1 + dSr/da0 * da0/da1
	eliminate a0:
	dSr/da0   [ S(x^0)  S(x^1)  S(x^2)  S(x^3) ]   a0   S(x^0.y)
	dSr/da1 = [ S(x^1)  S(x^2)  S(x^3)  S(x^4) ] * a1 - S(x^1.y)
	dSr/da2   [ S(x^2)  S(x^3)  S(x^4)  S(x^5) ]   a2   S(x^2.y)
	dSr/da3   [ S(x^3)  S(x^4)  S(x^5)  S(x^6) ]   a3   S(x^3.y)

	  [         0            0            0            0   ]   a0
	- [ S(x^0).xb    S(x^1).xb    S(x^2).xb    S(x^3).xb   ] * a1
	  [ S(x^0).xb^2  S(x^1).xb^2  S(x^2).xb^2  S(x^3).xb^2 ]   a2
	  [ S(x^0).xb^3  S(x^1).xb^3  S(x^2).xb^3  S(x^3).xb^3 ]   a3

  
	 
	Sr = S(y^2) + (val - a1.xb - a2.xb^2 - a3.xb^3)^2.S(x^0) + a1^2.S(x^2) + a2^2.S(x^4) + a3^2.S(x^6)
		+ 2{-(val - a1.xb - a2.xb^2 - a3.xb^3).S(x^0.y) - a1.S(x^1.y) - a2.S(x^2.y) - a3.S(x^3.y) }
		+ 2{ a0.a1.S(x^1) + a0.a2.S(x^2) + (a0.a3+a1.a2).S(x^3) + a1.a3.S(x^4) + a2.a3.S(x^5) }


	[      0       0       0       0 ]   [         1           xb           xb^2         xb^3 ]   a0   -val(xb)
	[ S(x^1)  S(x^2)  S(x^3)  S(x^4) ] - [ S(x^0).xb    S(x^1).xb    S(x^2).xb    S(x^3).xb   ] * a1 = S(x^1.y)
	[ S(x^2)  S(x^3)  S(x^4)  S(x^5) ]   [ S(x^0).xb^2  S(x^1).xb^2  S(x^2).xb^2  S(x^3).xb^2 ]   a2   S(x^2.y)
	[ S(x^3)  S(x^4)  S(x^5)  S(x^6) ]   [ S(x^0).xb^3  S(x^1).xb^3  S(x^2).xb^3  S(x^3).xb^3 ]   a3   S(x^3.y)


*/









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// CPolySegFitPoints Functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

CPolySegFitPoints::CPolySegFitPoints()
{
	Init();
}

CPolySegFitPoints::~CPolySegFitPoints()
{

}


void CPolySegFitPoints::Init()
{
	// Fit property settings
	m_nFitPtsForFirstPoly = 10;
	m_nFitPtsForPoly = 10;
	m_nFitPtsForNextAtChange = 10;


	m_TolY = 1e-4;
//	m_TolFirstSr = 1e-12;	// squared value!
	m_TolVariance = m_TolY * m_TolY;	// squared value!

//	m_nFitStatus = FITSTAT_NOPOINTS;
	m_SmallestInc = 0;
	m_SuggestedInc = 0;
	m_nReductions = 0;

	m_bForcePolyBreak = false;

	m_pPrevPoly = &m_Poly1;
	m_pCurrPoly = &m_Poly2;

	m_CurrPt[0].Zero();
}

//////////////////////////////////////////////////////////////////////

/*
parameters:
	ErrorTol -> Maximum error of any point = ~0.05mm
	InitialNumPoints -> Initial points to use = ~8
	InitialResidualTol -> maximum on initial fit = ~ ErrorTol^2

Start:

	Add InitialNumPoints points to CurrPoly and fit poly
	if (Residual > InitialResidualTol)
		reduce step size (1/2) and redo --^

	SumResidual = Residual
	SumResidualMax = Residual
Main Loop:
		Get Next Point
		Calc error of next point from CurrPoly
		if (Error <= ErrorTol)
		{
			SumResidualMax += Error^2
			store point to CurrPoly data, goto 'Main Loop'
		}
		else if (Error Not OK)
		{
			Refit CurrPoly with points up to but not including new point
			if (Residual OK)
				get error for new point
				if (Error OK)
					store point to CurrPoly data, goto 'Main Loop'
				else if (Error Not OK)
					don't include new point in CurrPoly
					goto 'Find Poly Change Location'
			else if (Residual Not OK)
				go back a little ??
				goto 'Find Poly Change Location'
		}
		Finish off

Find Poly Change Location:
		Set some NextPoly points after current point
		Find best change point between CurrPoly and NextPoly

		Advance CurrPoly to PrevPoly, NextPoly to CurrPoly
		Calc next point for CurrPoly
		goto 'Main Loop'

Smple method:
Main Loop:
{
	Get Next Point
	Calc error of next point from CurrPoly
	if (Error <= ErrorTol)
	{
		SumResidualMax += Error^2
		store point to CurrPoly data, goto 'Main Loop'
	}
	else if (Error Not OK)
	{
		Refit CurrPoly with points including new point
		Check if fit is OK
		if OK
			goto 'MainLoop'
		else
			remove last point
			Check if OK if not tested with previous point set - should be!
			end CurrPoly and start new from prev point
	}
}

Finish off




*/

// Main entry function for adding points to poly
int CPolySegFitPoints::NextPoint(double x, double y)
{
	m_CurrPt[0].x = x;
	m_CurrPt[1].x = m_CurrPt[2].x = 0;
	m_CurrPt[0].y = y;
	m_CurrPt[1].y = m_CurrPt[2].y = 0;
	return NextPointMethod3();
}

int CPolySegFitPoints::NextPoint(double t, CVector& pt)
{
	m_CurrPt[0].x = m_CurrPt[1].x = m_CurrPt[2].x = t;
	m_CurrPt[0].y = pt.x;
	m_CurrPt[1].y = pt.y;
	m_CurrPt[2].y = pt.z;
	return NextPointMethod3();
}

int CPolySegFitPoints::NextPoint(double t, CVect2& pt)
{
	m_CurrPt[0].x = m_CurrPt[1].x = m_CurrPt[2].x = t;
	m_CurrPt[0].y = pt.x;
	m_CurrPt[1].y = pt.y;
	m_CurrPt[2].y = 0;
	return NextPointMethod3();
}

void CPolySegFitPoints::BreakPolyAtNextPoint()
{
	m_bForcePolyBreak = true;
}


#if 0
int CPolySegFitPoints::NextPointMethod1()
{
	int retVal = PSF_OK;
	if (m_nFitStatus == FITSTAT_CURRSET)
		TryAddPoint1();
	else if (m_nFitStatus == FITSTAT_POLYCHANGE)
	{
//		m_NextPoly.AddPoint(m_CurrPt[0]);
//		if (m_NextPoly.GetNumPoints() == m_nFitPtsForNextAtChange)
//			FindPolyChangeLoc();
	}
	else if (m_nFitStatus <= FITSTAT_NOPOLY)		// no initial poly is established yet
	{
		if (m_nFitStatus == FITSTAT_NOPOINTS)
		{
			SetFirstPoint();
			m_nFitStatus = FITSTAT_NOPOLY;
		}
		m_pCurrPoly->AddPoint(m_CurrPt[0]);
//		AddFitPoint();
		if (m_pCurrPoly->GetNumPoints() == m_nFitPtsForPoly)	// calc first poly attempt once enough points
		{
			m_pCurrPoly->FitPolyToPointsContPV();

//			double SumErrors = m_pCurrPoly->GetPolyErrors();		// for testing
//			double SumResidual = m_pCurrPoly->GetSumResidual();
			double Variance = m_pCurrPoly->GetVariance();
			if (Variance <= m_TolVariance)			// Sr is squared value
			{
				m_nFitStatus = FITSTAT_CURRSET;
			}
			else
			{
				m_nFitStatus = FITSTAT_NOPOLY;
				m_SuggestedInc = m_pCurrPoly->GetAvgStepSize() / 10;
				m_pCurrPoly->ClearAllButFirst();
				m_nReductions++;
				retVal = PSF_REDUCESTEPSIZE;
			}
		}
	}
	else
		ASSERT(0);

	return retVal;
}

int CPolySegFitPoints::NextPointMethod2()
{
	int retVal = PSF_OK;
	if (m_nFitStatus == FITSTAT_CURRSET)
		TryAddPoint1();
	else if (m_nFitStatus <= FITSTAT_NOPOLY)		// no initial poly is established yet
	{
		if (m_nFitStatus == FITSTAT_NOPOINTS)
		{
			SetFirstPoint();
			m_nFitStatus = FITSTAT_NOPOLY;
		}
		m_pCurrPoly->AddPoint(m_CurrPt[0]);
//		AddFitPoint();
		if (m_pCurrPoly->GetNumPoints() == m_nFitPtsForPoly)	// calc first poly attempt once enough points
		{
			m_pCurrPoly->FitPolyToPoints();

//			double SumErrors = m_pCurrPoly->GetPolyErrors();		// for testing
//			double SumResidual = m_pCurrPoly->GetSumResidual();
			double Variance = m_pCurrPoly->GetVariance();
			if (Variance <= m_TolVariance)			// Sr is squared value
			{
				m_nFitStatus = FITSTAT_CURRSET;
			}
			else
			{
				m_SuggestedInc = m_pCurrPoly->GetAvgStepSize();
				if (m_SuggestedInc > 1.2 * m_SmallestInc)		// increment is at least 20% larger than min
				{
					m_SuggestedInc /= 2;
					if (m_SuggestedInc < m_SmallestInc)
						m_SuggestedInc = m_SmallestInc;
					m_pCurrPoly->ClearAllButFirst();
					m_nReductions++;
					m_nFitStatus = FITSTAT_NOPOLY;
					retVal = PSF_REDUCESTEPSIZE;
				}
				else
					m_nFitStatus = FITSTAT_CURRSET;
			}
		}
	}
	else
		ASSERT(0);

	return retVal;
}
#endif

int CPolySegFitPoints::TryAddPoint1()
{
// Current poly is set already
//	int iPointResult = 0;

	m_pCurrPoly->CalcPointError(m_CurrPt[0]);
	if ((fabs(m_CurrPt[0].yError) <= m_TolY)
		&& (m_CurrPt[0].yErrorSq <= 1.5 * m_pCurrPoly->GetVariance()))		// will find a sudden change in standard error
	{
		// add this point to current set   
		m_pCurrPoly->AddPoint(m_CurrPt[0]);
		return 1;
	}
// refit poly to points up until this then check error again
//	m_pCurrPoly->FitPolyToPointsContPV();
	m_pCurrPoly->FitPolyToPoints();

////////
	double SumErrors = m_pCurrPoly->GetPolyErrors();		// for testing
	double SumResidual = m_pCurrPoly->GetSumResidual();
///////
	m_pCurrPoly->CalcPointError(m_CurrPt[0]);
	if (fabs(m_CurrPt[0].yError) > m_TolY)
	{
		StartNextPoly();
		return -1;					// if error too large finish poly
	}
	double Variance = m_pCurrPoly->GetVariance();
	if (Variance > m_TolVariance)			// Variance is a squared value
	{
		StartNextPoly();
		return -2;					// if variance too large finish poly
	}

	if (m_CurrPt[0].yErrorSq > 1.5 * Variance)		// will find a sudden change in standard error
	{
		StartNextPoly();
		return -3;					// if error increases suddenly
	}
	// add this point to current poly segment   
	m_pCurrPoly->AddPoint(m_CurrPt[0]);
	return 2;
}



/*
Simple method:
Start:
Main Loop:
{
	Get Next Point
	if (poly not fitted yet)
	{
		Add point to CurrPoly
		if (enough points)
		{
			fit poly
			if (!CheckPolyFitOK())
				reduce step size (1/2) and remove all but first point
			else
				SumResidual = SumResidualMax = Residual
		}
		goto Main Loop
	}
	else (if poly is fitted)
	 
	Calc error of next point from CurrPoly
	if (Error <= ErrorTol)	// [quick check OK]
	{
		Add point to CurrPoly, sum up other error data
	}
	else if (Error Not OK)	// [quick check not OK]
	{
		Add point to CurrPoly
		Refit CurrPoly with points including new point
		if (fit is not OK)
		{
			remove last point from CurrPoly
			Refit CurrPoly with points not including new point
			Check if OK if not already tested with previous point set - should be!
			Start next Poly from end of CurrPoly
		}
	}
	if (NumPoints in CurrPoly > ~20)
		Reduce points and increase step size

}
Finish off
*/

int CPolySegFitPoints::NextPointMethod3()
{
	int retVal = PSF_OK;
	if (m_bForcePolyBreak)
	{
		bool bBroken = (ForceBreak() != 0);
		m_bForcePolyBreak = false;
		if (bBroken)
			return retVal;
	}
	if (m_pCurrPoly->m_bPolySet)
		TryAddPoint3();
	else							// no initial poly is established yet
	{
		if (m_pCurrPoly->GetNumPoints() == 0)
			SetFirstPoint();			// the very first
		m_pCurrPoly->AddPoint(m_CurrPt[0]);
		if (m_pCurrPoly->GetNumPoints() == m_nFitPtsForPoly)	// calc first poly attempt once enough points
		{
			m_pCurrPoly->FitPolyToPoints();
			if (CheckPolyFitOK())
				m_pCurrPoly->m_bPolySet = true;
			else
			{
				m_SuggestedInc = m_pCurrPoly->GetAvgStepSize();
				if (m_SuggestedInc > 1.2 * m_SmallestInc)		// increment is at least 20% larger than min
				{
					m_SuggestedInc /= 2;
					if (m_SuggestedInc < m_SmallestInc)
						m_SuggestedInc = m_SmallestInc;
					m_pCurrPoly->ClearAllButFirst();
					m_nReductions++;
					retVal = PSF_REDUCESTEPSIZE;
				}
				else
					m_pCurrPoly->m_bPolySet = true;
			}
		}
	}
	return retVal;
}

int CPolySegFitPoints::TryAddPoint3()
{
/* returns:
	-1	- New poly started, point added to new poly
	 1 - Added, point within tolerence of existing fitted poly
	 2 - Added, point within tolerence of refitted poly
*/
	ASSERT(m_pCurrPoly->m_bPolySet);			// Current poly is set already
	m_pCurrPoly->CalcPointError(m_CurrPt[0]);
	if (fabs(m_CurrPt[0].yError) <= m_TolY)
	{
		m_pCurrPoly->AddPoint(m_CurrPt[0]);			// add this point to current set   
		return 1;
	}
// refit poly to points including this point then check error again
	m_pCurrPoly->AddPointUnconfirmed(m_CurrPt[0]);			// add this point such that it can be removed
	m_pCurrPoly->FitPolyToPoints();

	if (!CheckPolyFitOK())
	{
		m_pCurrPoly->RemoveLatestPoint();
		m_pCurrPoly->FitPolyToPoints();
		StartNextPoly();
		m_pCurrPoly->AddPoint(m_CurrPt[0]);		// final point of prev poly segment as added in StartNextPoly()
		return -1;
	}
	else
		return 2;
}

int CPolySegFitPoints::ForceBreak()
{
	if (m_pCurrPoly->m_bPolySet)			// Current poly is set already
	{
		m_pCurrPoly->AddPoint(m_CurrPt[0]);			// add this point to current set   
//		m_pCurrPoly->FitPolyToPoints();		// not needed, will be fit in StartNextPoly()

		StartNextPoly();
		return 1;		// break successful
	}
	else if (m_pPrevPoly->m_bPolySet)		// add currpoly to prevpoly
	{
		m_pPrevPoly->AppendPoly(*m_pCurrPoly);		// add currpoly to prevpoly
		return 1;
	}
	else
		ASSERT(0); // trying to break before very first poly is set!
						// should reduce step size to get enough points to make a poly!!
	return 0;
}


////////////////////////////////

double CPolySegFitPoints::GetLastGoodXLoc()
{
	return m_pCurrPoly->GetMaxXPoint();
}

void CPolySegFitPoints::SetFirstPoint()
{
	// either initial pos and slope are set before
	// or not fixed and determined from curve fit
	m_pCurrPoly->m_yValInit[0] = m_CurrPt[0].y;		// pos
	m_pCurrPoly->m_yValInit[1] = 0;					// vel
	m_pCurrPoly->m_xPolyInit = m_CurrPt[0].x;
}


bool CPolySegFitPoints::CheckPolyFitOK()
{
///////// for testing
	double SumErrors = m_pCurrPoly->GetPolyErrors();		// for testing
	double SumResidual = m_pCurrPoly->GetSumResidual();
//////////

	m_nFitResultFlags = 0;
	m_pCurrPoly->CalcPointError(m_CurrPt[0]);
	if (fabs(m_CurrPt[0].yError) > m_TolY)
		m_nFitResultFlags |= 0x01;					// if error too large
	double Variance = m_pCurrPoly->GetVariance();
	if (Variance > m_TolVariance)			// Variance is a squared value
		m_nFitResultFlags |= 0x02;					// if variance too large

	return (m_nFitResultFlags == 0);
}


//if (AfxMessageBox("Send another plot?", MB_OKCANCEL | MB_ICONQUESTION) == IDCANCEL)
//	ExitThread(0);			// OK! apart from exit code?
//	_endthread();
//	ExitProcess(2);		// dosn't deallocate
//	PostQuitMessage(2);	// will return
//	AfxAbort();				// doesn't deallocate
//	AfxGetMainWnd()->SendMessage(WM_CLOSE);	// closes main wnd but then continues!
//	AfxGetMainWnd()->SendMessage(WM_QUIT);		// doesn't close window!
//	AfxEndThread(2);		// can't be used on CWinApp thread

void CPolySegFitPoints::StartNextPoly()
{
	// start a new poly segment with this point
	// First set coef's of prevPoly (if not first past!) with DblPoly fit to prev/curr polys
	m_pCurrPoly->m_xPolyFinal = m_pCurrPoly->m_xPointFinal;
	if (m_pPrevPoly->GetNumPoints() > 0)		// if not on first poly
	{
		m_pCurrPoly->CalcInitPV();		// for testing
		m_pPrevPoly->FitDblPolyToPointsContPV(*m_pCurrPoly);	// Initial pos and vel of prevPoly must be set!
		m_pPrevPoly->CalcFinalPV();	// for testing
		m_pCurrPoly->CalcInitPV();
		double vi0 = m_pPrevPoly->m_yValInit[0], vi1 = m_pPrevPoly->m_yValInit[1];	// for testing 
		m_pPrevPoly->CalcInitPV();		// for testing
		m_pPrevPoly->m_yValInit[0] = vi0; m_pPrevPoly->m_yValInit[1] = vi1; 

// more tests
		double SumErrC = m_pCurrPoly->GetPolyErrors();
		double SumResC = m_pCurrPoly->GetSumResidual();
		double SumErrP = m_pPrevPoly->GetPolyErrors();
		double SumResP = m_pPrevPoly->GetSumResidual();

		m_PolySegList.AddHead(*m_pPrevPoly);		// store previous poly segment

		static bool bPlotCurrentAlso = false;		// for testing


		m_pPrevPoly->PlotPointsAndBezier();		// for testing
		if (bPlotCurrentAlso)
			m_pCurrPoly->PlotPointsAndBezier();		// for testing
	}

	m_pTempPoly = m_pPrevPoly;			// advance polys
	m_pPrevPoly = m_pCurrPoly;
	m_pCurrPoly = m_pTempPoly;

	m_pCurrPoly->Reset();
	m_pCurrPoly->m_xPolyInit = m_pPrevPoly->m_xPolyFinal;
//	m_pCurrPoly->SetInitValsFromFinal();
	m_pCurrPoly->AddPoint(m_pPrevPoly->GetFinalPoint());		// add final point of prev poly segment

}

void CPolySegFitPoints::FindPolyChangeLoc()
{
	ASSERT(0);
/*		NextPoly has some points now
		Find best change point between CurrPoly and NextPoly

		Advance CurrPoly to PrevPoly, NextPoly to CurrPoly
		Calc next point for CurrPoly
*/
	// find join between m_CurrPoly & m_NextPoly
//	m_pCurrPoly->FitDblPolyToPointsContPV(m_NextPoly);


}









//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Global functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//const double pi = acos(-1.0);

int DoFFT()			// Fast Fourier Transform inplementation
{
	// Ref: Numerical Methods for Engineers, 2nd Edition, p.424
	const int M = 10;
	const int N = 1 << M;
	double x[N];		// real values - initial real data in this array
	double y[N];		// imaginary values - initial imaginary data in this array -> 0!
	double xtmp, ytmp;
	double mag[N];

	// First part is FFT calculation
	int N2 = N;
	for (int k = 1; k <= M; k++)
	{
		int N1 = N2;
		N2 /= 2;
		double angle = 0;
		double arg = 2*pi / N1;
		for (int j = 0; j < N2; j++)
		{
			double c = cos(angle);
			double s = -sin(angle);
			for (int i = j; i < N; i += N1)
			{
				int l = i + N2;
				xtmp = x[i] - x[l];
				x[i] += x[l];
				ytmp = y[i] - y[l];
				y[i] += y[l];
				x[l] = xtmp * c - ytmp * s;
				y[l] = ytmp * c + xtmp * s;
			}
			angle = (j+1) * arg;
		}
	}

	// Second part is bit-reversal routine to unscramble
	// the order of the resulting Fourier coefficients
	int j = 0;
	for (int i = 0; i <= N-2; i++)
	{
		if (i < j)
		{
			xtmp = x[j];
			x[j] = x[i];
			x[i] = xtmp;
			ytmp = y[j];
			y[j] = y[i];
			y[i] = ytmp;
		}
		int k = N / 2;
		while (k <= j)
		{
			j -= k;
			k /= 2;
		}
		j += k;
	}
	for (i = 0; i < N; i++)
	{
		x[i] /= N;
		y[i] /= N;
		mag[i] = sqrt(x[i]*x[i] + y[i]*y[i]);
	}

	return 1;
}


/*

	B-spline with uneven spans

			  /  \		  /	\
			/		 \		/		  \
	  a0=0	a1		a2		a3		a4=0
			t1		t2		t3		t4


	dv = a0.t + j.t^2 / 2

		= (a0 + a1).t / 2

	dp = v0.t + a0.t^2 / 2 + j.t^3 / 6

		= v0.t + (2.a0 + a1).t^2 / 6


	a1, a2, a3 are the d.o.f.
	a0 = a4 = 0
	v0 = v4 = 0
	p0 = p4 = 0
	p1 + p2 + p3 = 1

	2.v1 = a1.t1
	2.v2 = a1.t1 + (a1+a2)t2							= -a3.t4 - (a2+a3)t3
	2.v3 = a1.t1 + (a1+a2)t2 + (a2+a3)t3			= -a3.t4
	2.v4 = a1.t1 + (a1+a2)t2 + (a2+a3)t3 + a3.t4 = 0

	6.p1 = a1.t1^2
	6.p2 = a1.t1^2 + 3.a1.t1.t2 + (2.a1+a2).t2^2
	6.p3 = a1.t1^2 + 3.a1.t1.t2 + (2.a1+a2).t2^2 + 3(a1.t1 + (a1+a2)t2)t3 + (2.a2+a3).t3^2
	6.p4 = a1.t1^2 + 3.a1.t1.t2 + (2.a1+a2).t2^2 + 3(a1.t1 + (a1+a2)t2)t3 + (2.a2+a3).t3^2
												+ 3(a1.t1 + (a1+a2)t2 + (a2+a3)t3)t4 + 2.a3.t4^2 = 0

	6.p4 = a1.t1^2 + 3.a1.t1.t2 + (2.a1+a2).t2^2 + 3(a1.t1 + (a1+a2)t2)t3 + (2.a2+a3).t3^2
												- a3.t4^2 = 0

	set a1 = 1
	gives 2 equ and 2 dof a2, a3:

	t1 + t2 + a2.t2 + (a2+a3)t3 + a3.t4 = 0
	t1^2 + 3.t1.t2 + (2+a2).t2^2 + 3(t1 + (1+a2)t2)t3 + (2.a2+a3).t3^2 - a3.t4^2 = 0

	a2 (t2 + t3) + a3 (t3 + t4)  =  - (t1 + t2)

	a2 (t2^2 + 3.t2.t3 + 2.t3^2) + a3 (t3^2 - t4^2)  =  - t1^2 - 3.t1.t2 - 2.t2^2 - 3(t1+t2)t3
	a2 (t2 + t3)(t2 + 2.t3) + a3 (t3 + t4)(t3 - t4)  =  - (t1 + t2)(t1 + 2.t2 + 3.t3)

	|T| = (ad-bc) = (t2 + t3)(t3 + t4)(t3 - t4) - (t2 + t3)(t2 + 2.t3)(t3 + t4)
	|T| = - (t2 + t3)(t3 + t4)(t2 + t3 + t4)

	[a2] = [ (t3 + t4)(t3 - t4)    -(t3 + t4) ] * [ -(t1 + t2)                   ] / |T|
	[a3]   [-(t2 + t3)(t2 + 2.t3)   (t2 + t3) ]   [ -(t1 + t2)(t1 + 2.t2 + 3.t3) ]


	a2 =  (t1 + t2)(t3 + t4)(t1 + 2.t2 + 2.t3 + t4)  / |T|
	a3 = -(t1 + t2)(t2 + t3)(t1 + t2 + t3)           / |T|


	calculate p1, p2, p3
	6.p1 = a1.t1^2
		  = t1^2

	6.p2 = t1^2 + 3.t1.t2 + (2 + a2).t2^2
	6.p3 = - 3(t1 + (1+a2)t2 + (a2+a3)t3)t4 - 2.a3.t4^2


*/
