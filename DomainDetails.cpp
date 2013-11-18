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
#include "SctmPhys.h"
using SctmPhys::PhysProperty;

double FDVertex::Distance( FDVertex *vertex1, FDVertex *vertex2 )
{
	return SctmMath::sqrt(SctmMath::square(vertex1->X - vertex2->X) + SctmMath::square(vertex1->Y - vertex2->Y));
}

bool FDVertex::IsAtBoundary(FDBoundary::BCName bcName)
{
	return this->BndCond.Valid(bcName);
}

FDVertex::FDVertex(unsigned _id, double _x, double _y) : X(_x), Y(_y), id(_id)
{
	//These pointers will be set to NULL outside when the setting of adjacency of a vertex 
	//and initializing the contact.
	EastVertex = NULL;
	WestVertex = NULL;
	NorthVertex = NULL;
	SouthVertex = NULL;
	NortheastElem = NULL;
	NorthwestElem = NULL;
	SoutheastElem = NULL;
	SouthwestElem = NULL;
	Contact = NULL;
	Phys = new PhysProperty();
}

void FDBoundary::SetBndCond(BCName bcName, BCType bcType, double bcValue, VectorValue bcNormVec /*= VectorValue(0, 0)*/)
{
	SCTM_ASSERT( this->bc_valid.find(bcName) == this->bc_valid.end(), 10017 ); //make sure bcName not exists
	SCTM_ASSERT( bcType == BC_Dirichlet || bcNormVec.DirectionValid(), 10016 );

	//bc_valid.insert(map<BCName, bool>::value_type(bcName, true));
	this->bc_valid[bcName] = true;
	this->bc_types[bcName] = bcType;
	this->bc_values[bcName] = bcValue;
	if (bcNormVec.DirectionValid())
	{
		bcNormVec.Normalize();
	}
	this->bc_normVector[bcName] = bcNormVec;
}

FDBoundary::BCType FDBoundary::GetBCType(BCName bcName)
{
	//map<BCName, BCType>::iterator iter;
	//iter = this->bc_types.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_types.end(), 10010);
	//return iter->second;

	//the same with above method
	SCTM_ASSERT(bc_types.find(bcName)!=bc_types.end(), 10010);
	return bc_types[bcName];

}

double FDBoundary::GetBCValue(BCName bcName)
{
	//map<BCName, double>::iterator iter;
	//iter = this->bc_values.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_values.end(), 10010);
	//return iter->second;

	SCTM_ASSERT(bc_values.find(bcName)!=bc_values.end(), 10010);
	return bc_values[bcName];
}

bool FDBoundary::Valid(BCName bcName)
{
	return bc_valid[bcName];
	//map<BCName, bool>::iterator iter;
	//iter = this->bc_valid.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_valid.end(), 10010);
	//return iter->second;
}

void FDBoundary::RefreshBndCond(BCName bcName, double newValue)
{
	SCTM_ASSERT(bc_valid.find(bcName)!=bc_valid.end(), 10014); //make sure bcName exists

	bc_values[bcName] = newValue;
}

void FDBoundary::RefreshBndCond(BCName bcName, BCType bcType, double bcVal, VectorValue bcNormVec /*= VectorValue(0, 0)*/)
{
	SCTM_ASSERT(bc_valid.find(bcName)!=bc_valid.end(), 10014); //make sure bcName exists
	SCTM_ASSERT( bcType == BC_Dirichlet || bcNormVec.DirectionValid(), 10016 );

	bc_types[bcName] = bcType;
	bc_values[bcName] = bcVal;
	if (bcNormVec.DirectionValid())
	{
		bcNormVec.Normalize();
	}
	bc_normVector[bcName] = bcNormVec;
}

VectorValue & FDBoundary::GetBCNormVector(BCName bcName)
{
	SCTM_ASSERT(bc_normVector.find(bcName)!=bc_normVector.end(), 10010);
	return bc_normVector[bcName];
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
