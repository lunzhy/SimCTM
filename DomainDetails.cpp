/**
* @file DomainDetails.cpp
* @brief This file contains the detailed implementation of the methods declared in the header file.
*
*
*
* @author
* @version 
* @date 2013-7-4   19:41
* @note
* @todo
*/

#include "DomainDetails.h"
#include "SctmMath.h"
#include "SctmUtils.h"

double FDVertex::Distance( FDVertex *vertex1, FDVertex *vertex2 )
{
	return SctmMath::sqrt(SctmMath::square(vertex1->X - vertex2->X) + SctmMath::square(vertex1->Y - vertex2->Y));
}

bool FDVertex::IsAtBoundary()
{
	return this->BoundaryCond.Valid();
}

void FDBoundary::SetBndCond(BndCond bndtype, double bndvalue1, double bndvalue2)
{
	this->valid = true;
	this->bndType = bndtype;
	if (bndtype == BC_Artificial)
	{
		bndvalue1 = 0;
		bndvalue2 = 0;
	}
	this->bndValue1 = bndvalue1;
	this->bndValue2 = bndvalue2;
}

double FDBoundary::BndValue2() const
{
	if (bndType == BC_Dirichlet || bndType == BC_Artificial)
	{
		return 0;
	}
	else
	{
		return bndValue2;
	}
}

double FDBoundary::BndValue1() const
{
	if (bndType == BC_Artificial)
	{
		return 0;
	}
	else
	{
		return bndValue1;
	}
}

double FDBoundary::BndValuePotential() const
{
	double ret = 0;
	if (bndType == BC_Dirichlet)
	{
		ret = bndValue1;
	}
	else
	{
		SCTM_ASSERT(true, 9);
	}
	return ret;
}

double FDBoundary::BndValueElecFieldWestEast() const
{
	double ret = 0;
	switch (bndType)
	{
	case BC_Dirichlet:
		SCTM_ASSERT(true, 9);
		break;
	case BC_Neumann:
		ret = bndValue1;
		break;
	case BC_Artificial:
		ret = 0;
	}
	return ret;
}

double FDBoundary::BndValueElecFieldSouthNorth() const
{
	double ret = 0;
	switch (bndType)
	{
	case BC_Dirichlet:
		SCTM_ASSERT(true, 9);
		break;
	case BC_Neumann:
		ret = bndValue2;
	case BC_Artificial:
		ret = 0;
	}
	return ret;
}
