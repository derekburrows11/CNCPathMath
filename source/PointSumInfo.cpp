// PointSumInfo.cpp: implementation of the CPointSumInfo class.
//
//////////////////////////////////////////////////////////////////////


//#include "stdafx.h"


#include "PointSumInfo.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif





//////////////////////////////////
// SPointPowerData functions
//////////////////////////////////

void CPointPowerData::Zero()
{
	num = 0;
	for (int i = 0; i < 7; i++)
		Xp[i] = 0;
	for (i = 0; i < 4; i++)
		YXp[i] = 0;
	Y2 = 0;
}

void CPointPowerData::SumUpPoint(const CPointPowerData& prevSum, const CVect2& pt)
{
	double xPow;
	double x = pt.x, y = pt.y;		// create local copies	
	num = prevSum.num + 1;
	Xp[0]  = prevSum.Xp[0] + 1;
	YXp[0] = prevSum.YXp[0] + y;
	Xp[1]  = prevSum.Xp[1] + (xPow = x);
	YXp[1] = prevSum.YXp[1] + (xPow * y);
	Xp[2]  = prevSum.Xp[2] + (xPow *= x);
	YXp[2] = prevSum.YXp[2] + (xPow * y);
	Xp[3]  = prevSum.Xp[3] + (xPow *= x);
	YXp[3] = prevSum.YXp[3] + (xPow * y);
	Xp[4]  = prevSum.Xp[4] + (xPow *= x);
	Xp[5]  = prevSum.Xp[5] + (xPow *= x);
	Xp[6]  = prevSum.Xp[6] + (xPow * x);
	Y2 = prevSum.Y2 + (y * y);		// needed to calc residual error (sum of squares)

	this->x = x;		// store point values for checking
	this->y = y;
}

///////////////////////////////////////////

void CPointSPowerData::Zero()
{
	num = 0;
	for (int i = 0; i < 7; i++)
		Sp[i] = 0;
	for (int ax = 0; ax < 3; ax++)
	{
		for (i = 0; i < 4; i++)
			SpV[i][ax] = 0;
		V2[ax] = 0;
		m_bConstValue[ax] = true;		// starts true
		m_Value[ax] = 0;
	}
}

void CPointSPowerData::SumUpPoint(const CPointSPowerData& prevSum, double s, const CVector& pt)
{
	double sPow, sPow2, sPow3;
	num = prevSum.num + 1;
	Sp[0] = prevSum.Sp[0] + 1;
	Sp[1] = prevSum.Sp[1] + (sPow = s);
	Sp[2] = prevSum.Sp[2] + (sPow *= s);
	sPow2 = sPow;
	Sp[3] = prevSum.Sp[3] + (sPow *= s);
	sPow3 = sPow;
	Sp[4] = prevSum.Sp[4] + (sPow *= s);
	Sp[5] = prevSum.Sp[5] + (sPow *= s);
	Sp[6] = prevSum.Sp[6] + (sPow * s);
	for (int ax = 0; ax < 3; ax++)
	{
		double v = pt[ax];
		SpV[0][ax] = prevSum.SpV[0][ax] + v;
		SpV[1][ax] = prevSum.SpV[1][ax] + v * s;
		SpV[2][ax] = prevSum.SpV[2][ax] + v * sPow2;
		SpV[3][ax] = prevSum.SpV[3][ax] + v * sPow3;
		V2[ax] = prevSum.V2[ax] + v*v;				// needed to calc residual error (sum of squares)
		// check if all values of this axis are the same!
		m_bConstValue[ax] = prevSum.m_bConstValue[ax];
		m_Value[ax] = v;
		if (m_bConstValue[ax])
			if (num != 1 && v != prevSum.m_Value[ax])	// just set first time
				m_bConstValue[ax] = false;
	}
}

void CPointSPowerData::SumUpPoint(double s, const CVector& pt)
{
	double sPow, sPow2, sPow3;
	num++;
	Sp[0] += 1;
	Sp[1] += (sPow = s);
	Sp[2] += (sPow *= s);
	sPow2 = sPow;
	Sp[3] += (sPow *= s);
	sPow3 = sPow;
	Sp[4] += (sPow *= s);
	Sp[5] += (sPow *= s);
	Sp[6] += (sPow * s);
	for (int ax = 0; ax < 3; ax++)
	{
		double v = pt[ax];
		SpV[0][ax] += v;
		SpV[1][ax] += v * s;
		SpV[2][ax] += v * sPow2;
		SpV[3][ax] += v * sPow3;
		V2[ax] += v*v;				// needed to calc residual error (sum of squares)
		// check if all values of this axis are the same!
		if (m_bConstValue[ax])
			if (num == 1)	// just set first time
				m_Value[ax] = v;
			else if (v != m_Value[ax])
				m_bConstValue[ax] = false;
	}
}

void CPointSPowerData::ScaleS(double scale)
{
	double scalePow, scalePow2, scalePow3;
	Sp[1] *= (scalePow = scale);
	Sp[2] *= (scalePow *= scale);
	scalePow2 = scalePow;
	Sp[3] *= (scalePow *= scale);
	scalePow3 = scalePow;
	Sp[4] *= (scalePow *= scale);
	Sp[5] *= (scalePow *= scale);
	Sp[6] *= (scalePow * scale);
	for (int ax = 0; ax < 3; ax++)
	{
		SpV[1][ax] *= scale;
		SpV[2][ax] *= scalePow2;
		SpV[3][ax] *= scalePow3;
	}
}


//////////////////////////////////
// CPointSumInfo functions
//////////////////////////////////

CPointSumInfo::CPointSumInfo() : m_SumPointsList(100)
{
	m_StepVariation = 4;

/*	m_iStart = 0;
	m_iEnd = 0;
	m_iLength = 0;
	m_iBufferSize = sizeof(m_SumPoints) / sizeof(m_SumPoints[0]);
	m_iStartRef = 0;
*/
}


bool CPointSumInfo::AddPoint(CVect2& pt)
{
	CPointPowerData pntData;
	if (GetNumPoints() == 0)			// set first point
	{
		pntData.Zero();
		pntData.SumUpPoint(pntData, pt);
		m_SumPointsList.AddHead(pntData);

		m_MinListX = m_MaxListX = pt.x;
		return true;
	}

	double stepSize = 0;
	int nLocStep = 0;			// step between min & max
	// find new point location relative to list and set New Point flags
	if (pt.x > m_MaxListX)		// usual case
		stepSize = pt.x - m_MaxListX;
	else if (pt.x < m_MinListX)
		stepSize = m_MinListX - pt.x;

	if (GetNumPoints() <= 1)
		m_MinListXStep = m_MaxListXStep = stepSize;
	else if (stepSize != 0)
		if (stepSize > m_MaxListXStep)
		{
			nLocStep = 1;
			if (stepSize > m_StepVariation * m_MinListXStep)
				nLocStep = 2;		// step to large
		}
		else if (stepSize < m_MinListXStep)
		{
			nLocStep = -1;
			if (stepSize * m_StepVariation < m_MaxListXStep)
				nLocStep = -2;	// step too small
		}

	// add point to list - or could test if step not good!
	pntData.SumUpPoint(m_SumPointsList.GetHead(), pt);
	m_SumPointsList.AddHead(pntData);

	if (pt.x > m_MaxListX)		// usual case
		m_MaxListX = pt.x;
	else if (pt.x < m_MinListX)
		m_MinListX = pt.x;
	if (nLocStep > 0)
		m_MaxListXStep = stepSize;
	else if (nLocStep < 0)
		m_MinListXStep = stepSize;

	if (stepSize != 0 && abs(nLocStep) <= 1)
		return true;
	return false;	
}

void CPointSumInfo::SaveState(CPointSumInfoState& psis)
{
	psis.m_MaxListX = m_MaxListX;
	psis.m_MinListX = m_MinListX;
	psis.m_MaxListXStep = m_MaxListXStep;
	psis.m_MinListXStep = m_MinListXStep;
	psis.m_iNumPoints = GetNumPoints();
}

bool CPointSumInfo::RetrieveState(CPointSumInfoState& psis)
{
	while (GetNumPoints() > psis.m_iNumPoints)
		m_SumPointsList.RemoveHead();
	if (GetNumPoints() != psis.m_iNumPoints)
	{
		ASSERT(0);
		return false;
	}
	m_MaxListX = psis.m_MaxListX;
	m_MinListX = psis.m_MinListX;
	m_MaxListXStep = psis.m_MaxListXStep;
	m_MinListXStep = psis.m_MinListXStep;
	return true;
}

void CPointSumInfo::ClearAll()
{
	m_SumPointsList.RemoveAll();
}
/*
int CPointSumInfo::ClearAllButFirst()
{
	if (m_SumPointsList.IsEmpty())
		return 0;

	while (m_SumPointsList.GetCount() > 1)			// avoid removing all and reallocating memory
		m_SumPointsList.RemoveHead();
	m_MinListX = m_MaxListX = m_SumPointsList.GetTail().Xp[1];
	
	return 1;
}
*/
double CPointSumInfo::GetAvgStepSize()
{
	if (GetNumPoints() >= 2)
		return (m_MaxListX - m_MinListX) / (GetNumPoints() - 1);
	return 0;
}

CPointPowerData& CPointSumInfo::GetFirstSums()
{
	return m_SumPointsList.GetTail();
}
POSITION CPointSumInfo::GetFirstSumsPos()
{
	return m_SumPointsList.GetTailPosition();
}
CPointPowerData& CPointSumInfo::GetNextSums(POSITION& pos)
{
	return m_SumPointsList.GetPrev(pos);
}

CPointPowerData& CPointSumInfo::GetLastSums()
{
	return m_SumPointsList.GetHead();
}

int CPointSumInfo::GetSumsFor(CPointPowerData& sumPoints, int iInit, int iFinal)
{
	if (iInit == 0 && iFinal == -1)
		sumPoints = m_SumPointsList.GetHead();
	else
		ASSERT(0);
	return 1;
}

