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

double FDVertex::Distance( FDVertex *vertex1, FDVertex *vertex2 )
{
	return SctmMath::sqrt(SctmMath::square(vertex1->X - vertex2->X) + SctmMath::square(vertex1->Y - vertex2->Y));
}

bool FDVertex::IsAtBoundary()
{
	return this->BoundaryCond.Valid;
}

void FDBoundary::SetBndCond(BndCond bndType, double bndValue)
{
	this->Valid = true;
	this->Type = bndType;
	this->Value = bndValue;
}
