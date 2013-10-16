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

bool FDVertex::IsAtBoundary(FDBoundary::BCName bcName)
{
	return this->BndCond.Valid(bcName);
}

void FDBoundary::SetBndCond(BCName bcName, BCType bcType, double bcValue1, double bcValue2)
{
	//this->valid = true;
	bc_valid.insert(map<BCName, bool>::value_type(bcName, true));
	this->bc_types[bcName] = bcType;

	switch (bcType)
	{
	case BC_Dirichlet:
		this->bc_values[bcName] = bcValue1;
		this->bc_values_second[bcName] = bcValue2; //BC_Dirichlet might have the second value (in case of current)
		break;
	case BC_Neumann:
		this->bc_values[bcName] = bcValue1;
		this->bc_values_second[bcName] = bcValue2;
		break;
	case BC_Artificial:
		this->bc_values[bcName] = 0;
		this->bc_values_second[bcName] = 0;
		break;
	}
}

FDBoundary::BCType FDBoundary::GetBCType(BCName bcName)
{
	map<BCName, BCType>::iterator iter;
	iter = this->bc_types.find(bcName);
	SCTM_ASSERT(iter!=this->bc_types.end(), 10010);
	return iter->second;
}

double FDBoundary::GetBCValue(BCName bcName)
{
	map<BCName, double>::iterator iter;
	iter = this->bc_values.find(bcName);
	SCTM_ASSERT(iter!=this->bc_values.end(), 10010);
	return iter->second;
}

double FDBoundary::GetBCValueWestEast(BCName bcName)
{
	return GetBCValue(bcName);
}

double FDBoundary::GetBCValueSouthNorth(BCName bcName)
{
	map<BCName, double>::iterator iter;
	iter = this->bc_values_second.find(bcName);
	SCTM_ASSERT(iter!=this->bc_values_second.end(), 10010);
	return iter->second;
}

bool FDBoundary::Valid(BCName bcName)
{
	return this->bc_valid[bcName];
	//map<BCName, bool>::iterator iter;
	//iter = this->bc_valid.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_valid.end(), 10010);
	//return iter->second;
}

void FDBoundary::RefreshBndCondValue(BCName bcName, double bcValue1, double bcValue2)
{
	//check if the boundary condition exist.
	map<BCName, BCType>::iterator iter;
	iter = this->bc_types.find(bcName);
	SCTM_ASSERT(iter!=this->bc_types.end(), 10014);

	this->bc_values[bcName] = bcValue1;
	this->bc_values_second[bcName] = bcValue2;
}

FDElement::FDElement(unsigned int _id, FDVertex *_swVertex, FDVertex *_seVertex, FDVertex *_neVertex, FDVertex *_nwVertex)
	:id(_id), SouthwestVertex(_swVertex), SoutheastVertex(_seVertex), NortheastVertex(_neVertex), NorthwestVertex(_nwVertex)
{
	WestLength =FDVertex::Distance(NorthwestVertex, SouthwestVertex);
	SouthLength = FDVertex::Distance(SouthwestVertex, SoutheastVertex);
	EastLength = FDVertex::Distance(NortheastVertex, SoutheastVertex);
	NorthLength =FDVertex::Distance(NorthwestVertex, NortheastVertex);
	Area = WestLength * SouthLength;
	///TODO: judge if the corresponding lengths equal to each other
	double relError = (WestLength - EastLength) / WestLength;
	SCTM_ASSERT(relError < 0.01, 10003);
	relError = (SouthLength - NorthLength) / SouthLength;
	SCTM_ASSERT(relError < 0.01, 10003);
}
